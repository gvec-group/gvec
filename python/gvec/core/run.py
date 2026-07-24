# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""run gvec from python"""

import logging
import re
import shutil
import time
import warnings
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Literal

import numpy as np
import xarray as xr
from pandas import read_csv
import scipy

import gvec
from gvec.core.state import State
from gvec.errors import catch_gvec_errors
from gvec.lib import _binding, _run, _globals
from gvec.util import CaseInsensitiveDict as cidict

DEFAULT_MINIMIZE_TOL = 1e-6  # different to default `minimize_tol` in fortran
DEFAULT_TOTALITER = 10**5  # different to default `maxIter` in fortran
AUTO_STAGE_IOTA_INIT_MAXITER = 10
AUTO_STAGE_IOTA_INIT_IOTA_TOL = 1e-3
AUTO_PICARD_CURRENT_IOTA_TOL = 1e-5
AUTO_PICARD_CURRENT_IOTA_TOL_TARGET = 1e-8  # least-squares fit: limited to single precision!
AUTO_MINIMIZE_TOL_MIN = 1e-3
DEFAULT_PICARD_CURRENT_MAXRESTARTS = 50
DIAGNOSTICS_NPOINTS = 256 + 1
IOTA_BSPLINE_DEGREE = 5
IOTA_SAMPLING_FACTOR = 4
IOTA_KNOTS_PER_ELEMENT_MIN = 2
IOTA_KNOTS_PER_ELEMENT_MAX = 256
IOTA_COEFS_MAX = 768  # at 20 characters per coefficient: 15360 characters + overhead < 2^14
IOTA_MIN_RHO = 1e-2  # ignore near-axis region for Δiota & Δfit, TODO: fix accuracy problems


def run(
    parameters: Mapping | Path | str,
    restartstate: Path | State | None = None,
    runpath: Path | str | None = None,
    redirect_gvec_stdout: bool = True,
    quiet: bool = False,
    parameter_format: Literal["toml", "yaml"] = "toml",
    keep_intermediates: Literal["all", "stages"] | None = None,
    loglevel: Literal["WARNING", "INFO", "DEBUG"] | None = None,
):
    """Run GVEC with the provided parameters.

    Parameters
    ----------
    parameters : Mapping | Path | str
        GVEC parameter dictionary or parameter file.
    restartstate : Path | State | None, optional
        Path to restart file or restart State object, by default None
    runpath : Path | str | None, optional
        Path to the directory where GVEC is executed. Overwrites existing directories or cerates a new directory if it does not exist.
        The default is None.
    redirect_gvec_stdout : bool, optional
        Whether to redirect GVEC's stdout. The default is True.
    quiet : bool, optional
        Whether to suppress all output. The default is False.
    parameter_format : Literal["toml", "yaml"], optional
        Format of the parameter file automatically generated if `picard_current="auto"`, by default "toml"
    keep_intermediates : Literal["all", "stages"] | None, optional
        Whether to keep intermediate results of GVEC. With `"all"`, all intermediate results are kept. With `"stages"`,
        only the final restarts from each stage are kept. With `None`, all intermediate results are deleted. The default is None.
    loglevel : Literal["WARNING", "INFO", "DEBUG"] | None, optional
        Set the loglevel for the gvec logger.

    Returns
    -------
    gvec.Run
        Run-object containing history and current state of the GVEC run.
    """

    # logger setup
    gvec.util.logging_setup()
    logger = logging.getLogger("gvec.run")
    if loglevel is not None:
        logger.setLevel(loglevel)

    # rundirectory setup
    if runpath is None:
        runpath = Path.cwd()
    else:
        runpath = Path(runpath)

        if not runpath.exists():
            runpath.mkdir()
            logger.info(f"created run directory {runpath}")
        else:
            logger.info(f"run directory {runpath} is overwritten")

    # transform parameters into hierarchical cidict.
    if isinstance(parameters, (Path, str)):
        parameters = Path(parameters)
        params = gvec.util.read_parameters(parameters)
    else:
        params = gvec.util.stack_parameters(parameters)

    # Run the case
    with gvec.util.chdir(runpath):
        run_instance = Run(
            params=params,
            state=restartstate,
            redirect_gvec_stdout=redirect_gvec_stdout,
            quiet=quiet,
            parameter_format=parameter_format,
        )

        run_instance.run(keep_intermediates=keep_intermediates)

    return run_instance


class Run:
    def __init__(
        self,
        params: Mapping,
        state: Path | State | None = None,
        redirect_gvec_stdout: bool = True,
        quiet: bool = True,
        parameter_format: Literal["toml", "yaml"] = "toml",
    ):
        """
        State of a GVEC run during a stage, e.g. a picard iteration during a current optimization run.

        Parameters
        ----------
        params : Mapping
            GVEC parameter dictionary.
        state : Path | gvec.State | None, optional
            Statefile or State to restart from. The default is None.
        redirect_gvec_stdout : bool, optional
            Whether to redirect GVEC's stdout. The default is True.
        quiet : bool, optional
            Whether to suppress all output. The default is False.

        Raises
        ------
        ValueError
            If neither params nor a state is provided.
        """

        if params:
            self.parameters = params.copy()
        elif isinstance(state, gvec.State):
            self.parameters = state.parameters.copy()
        else:
            raise ValueError(
                f"Neither valid 'params' nor a valid 'state' are provided! params: {params}, state: {state}"
            )
        self.parameters = cidict.recursive(self.parameters)

        # params: user provided parameters
        # self.parameters: global/fallback parameters
        # _state_parameters reflect the state of the Run, that is they change during each stage and/or restart
        #   after a stage they are reset to parameters. They are NOT state.parameters, as those reflect the state of a
        #   finished restart.
        if "ProjectName" not in self.parameters:
            self.parameters["ProjectName"] = "GVEC"
        if "picard_current" not in self.parameters:
            self.parameters["picard_current"] = "off"

        self.logger = logging.getLogger("gvec.run")
        self.filetype = parameter_format

        project_dir = Path(f"{self.parameters['ProjectName']}_gvec_stages")
        if project_dir.exists():
            self.logger.debug(f"Removing existing run directory {project_dir}")
            shutil.rmtree(project_dir)
        project_dir.mkdir()
        self.project_dir = project_dir

        self.nth_stage = 0
        self.nth_run = 0  # nth run in stage
        self.rundir: Path = None
        self.quiet = quiet
        self.diagnostics_run: xr.Dataset = None
        self.diagnostics_minimizer = None
        self.GVEC_iter_used = 0
        self.redirect_gvec_stdout = redirect_gvec_stdout

        self.totaliter = self.parameters.get("totaliter", DEFAULT_TOTALITER)
        if "maxIter" not in self.parameters:
            self.parameters["maxIter"] = self.totaliter
        if self.parameters["maxIter"] <= 0:
            raise ValueError(
                f"Value 'maxIter' must be greater than 0 (current value is {self.parameters['maxIter']})."
            )

        if "minimize_tol" not in self.parameters:
            self.parameters["minimize_tol"] = DEFAULT_MINIMIZE_TOL

        picard_current = self.parameters.get("picard_current", "off")

        if "I_tor" in self.parameters and picard_current != "off":
            self.curr_constraint = True
        elif "I_tor" in self.parameters and picard_current == "off":
            raise gvec.errors.MissingParameterError(
                "'I_tor' is provided but 'picard_current' is set to 'off' or not provided. Please provide a valid 'picard_current', e.g. 'auto'."
            )
        elif "I_tor" not in self.parameters and picard_current != "off":
            raise gvec.errors.MissingParameterError(
                "Expected 'I_tor' in the parameters since 'picard_current' is not 'off'."
                + " Please set 'picard_current' to 'off' if you want to use a fixed 'iota' profile or provide 'I_tor'."
            )
        else:
            self.curr_constraint = False

        # Automatically generate the stages for the current optimization
        if "stages" in self.parameters and picard_current == "auto":
            raise gvec.errors.InvalidParameterError(
                "Picard current is set to 'auto' but 'stages' is specified!"
            )
        if "stages" in self.parameters:
            self.stages = self.parameters["stages"]
            del self.parameters["stages"]
        else:
            self.parameters, self.stages = auto_generate_stages(
                self.parameters, self.curr_constraint
            )
            self.logger.info(f"... generated {len(self.stages)} stages.")

        # load restart state
        if isinstance(state, gvec.State):
            self.logger.info(f"reading restart state from {state.statefile}")
            self.state = state

            if self.parameters.get("which_hmap", 1) != self.state.parameters.get(
                "which_hmap", 1
            ):
                warnings.warn(
                    f"restarting with hmap={self.parameters.get('which_hmap', 1)} from hmap={self.state.parameters.get('which_hmap', 1)}."
                )

            # set initial iota profile
            if "iota" not in self.parameters:
                self.parameters["iota"] = self.state.parameters["iota"]

            # compute boundary perturbation relative to restart state
            base, perturbation = gvec.util.compute_boundary_perturbation(
                self.state.parameters, self.parameters
            )
            if perturbation:
                self.parameters |= base
                if "boundary_perturb_type" not in self.parameters:
                    self.parameters["boundary_perturb_type"] = "cosm"
                self.stages[0] |= perturbation
                self.stages[0]["boundary_perturb"] = True

        # load I_tor profile (and set initial iota if not provided)
        if self.curr_constraint:
            self._set_I_tor(self.parameters)
            self.iota_rms = None
            if "iota" not in self.parameters:
                self.parameters["iota"] = cidict({"type": "polynomial", "coefs": [0.0]})
        else:
            self.I_tor_target = None

        # account for the change in relative path of restart, hmap, etc. files
        for key, value in self.parameters.items():
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ]:
                self.parameters[key] = Path(value).resolve()
                try:
                    if self.parameters[key].exists():
                        Path(self.parameters[key].name).symlink_to(self.parameters[key])
                    else:
                        raise FileNotFoundError(
                            f"Could not find {key} at {self.parameters[key]}"
                        )
                except FileExistsError:
                    continue

        # count the number of runs in each stage, for dynamic progressbar during current optimization
        self.n_runs_in_stage = [0 for _ in self.stages]

        self._state_parameters = self.parameters.copy()

        # load restart state (2nd part)
        if not isinstance(state, gvec.State):  # Path, str or None
            if state is not None:
                warnings.warn("restarting from statefile without associated parameterfile")
            initial_param_file = self.project_dir / "parameters_initial.ini"
            gvec.util.write_parameters(self._state_parameters, initial_param_file)
            self.state = gvec.State(initial_param_file, state)

        parameters_stages = self.parameters.copy()
        parameters_stages["stages"] = self.stages
        parameters_stages_name = (
            f"parameter_{self.parameters['ProjectName']}.stages.{self.filetype}"
        )
        self.logger.info(f"... writing new parameters to '{parameters_stages_name}'")
        gvec.util.write_parameters(
            parameters_stages,
            project_dir / parameters_stages_name,
        )

    def _set_I_tor(self, params: Mapping):
        """Evaluate and set the target toroidal current profile at linearily spaced positions in rho.

        Raises
        ------
        gvec.errors.InvalidParameterError
            If an unknown profile type is provided.
        """
        match params["I_tor"].get("type", "polynomial"):
            case "polynomial":
                coefs = np.array(params["I_tor"]["coefs"])
                coefs *= params["I_tor"].get("scale", 1.0)
                self.I_tor = np.polynomial.Polynomial(coefs)
                if abs(coefs[0]) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis!  I_tor(rho=0): {coefs[0]}"
                    )
            case "bspline":
                from scipy.interpolate import BSpline

                coefs = np.array(params["I_tor"]["coefs"], dtype=float)
                coefs *= params["I_tor"].get("scale", 1.0)
                knots = np.array(params["I_tor"]["knots"], dtype=float)
                deg = np.sum(knots == knots[0]) - 1
                self.I_tor = BSpline(knots, coefs, deg)
                if abs(self.I_tor(0.0)) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis! I_tor(rho=0): {self.I_tor(0.0)}"
                    )
            case "interpolation":
                from scipy.interpolate import make_interp_spline

                y_vals = np.array(params["I_tor"]["vals"], dtype=float)
                rho2_vals = np.array(params["I_tor"]["rho2"], dtype=float)
                if min(np.sqrt(rho2_vals)) > 1e-4:
                    rho2_vals = np.append([0], rho2_vals)
                    y_vals = np.append([0], y_vals)
                self.I_tor = make_interp_spline(rho2_vals, y_vals)
                if abs(self.I_tor(0.0)) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis! I_tor(rho=0): {self.I_tor(0.0)}"
                    )
            case _:
                raise gvec.errors.InvalidParameterError(
                    f"Unknown Itor type: {params['I_tor']['type']}"
                )

    @property
    def iota_tol(self):
        return self._state_parameters["picard_current"]["iota_tol"]

    def run_single_minimization(self):
        """Run a single GVEC energy minimization using the current parameters. The run-state is updated after the run."""
        start_time = time.time()
        # find previous state
        if self.state.statefile:
            self.logger.debug(f"Restart from statefile {self.state.statefile}")

        # prepare the run directory
        self.rundir = self.project_dir / Path(f"{self.nth_stage:1d}-{self.nth_run:02d}")
        if self.rundir.exists():
            self.logger.debug(f"Removing existing run directory {self.rundir}")
            shutil.rmtree(self.rundir)
        self.rundir.mkdir()
        self.logger.debug(f"Created run directory {self.rundir}")

        # write parameterfile & run GVEC
        gvec.util.write_parameter_file_ini(
            gvec.util.flatten_parameters(self._state_parameters),
            self.rundir / "parameter.ini",
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage} run {self.nth_run})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )
        with gvec.util.chdir(self.rundir):
            fortran_run(
                "parameter.ini",
                self.state.statefile,
                stdout_path="stdout.txt" if self.redirect_gvec_stdout else None,
            )

        # postprocessing
        self.state = gvec.core.state.find_state(self.rundir)
        iterations = int(re.match(r".*State.*_(\d+)\.dat", self.state.statefile.name).group(1))
        iteration_offset = self.GVEC_iter_used
        self.GVEC_iter_used += iterations
        max_iterations = self._state_parameters["maxIter"]
        tolerance = self._state_parameters["minimize_tol"]
        self.logger.debug(f"Postprocessing statefile {self.state.statefile}")

        # update iota
        if self.curr_constraint and self._state_parameters["picard_current"] != "off":
            sgrid = self.state.get_radial_gridpoints()
            sgrid[0] = 0.0
            sgrid[-1] = 1.0
            Qs = ["iota", "iota_0", "iota_curr_0", "I_tor"]
            while True:
                knots = _linspace_sgrid(sgrid, self.knots_per_sgrid_element)
                knots = _add_knot_multiplicity(knots, IOTA_BSPLINE_DEGREE)
                rho = _linspace_sgrid(
                    sgrid, IOTA_SAMPLING_FACTOR * self.knots_per_sgrid_element
                )
                ev = _evaluate_extrapolate_in_rho2(
                    self.state, *Qs, rho=rho, theta="int", zeta="int"
                )[Qs]
                iota_target = ev.iota_0 + self.I_tor(rho**2) * ev.iota_curr_0

                # least-squares fit
                bspline = scipy.interpolate.make_lsq_spline(
                    rho**2, iota_target, knots**2, k=IOTA_BSPLINE_DEGREE
                )
                # enforce boundary conditions
                bspline.c[0] = iota_target[0]
                bspline.c[-1] = iota_target[-1]

                # check fit error
                delta_iota_fit = iota_target - bspline(rho**2)
                # ignore near-axis region as we CURRENTLY have accuracy problems there!!!
                delta_iota_fit = delta_iota_fit.sel(rho=slice(IOTA_MIN_RHO, None))
                self.rms_iota_fit = np.sqrt((delta_iota_fit**2).mean().item())
                if self.rms_iota_fit < self.iota_tol:
                    break
                if self.knots_per_sgrid_element >= IOTA_KNOTS_PER_ELEMENT_MAX:
                    self.logger.warning(
                        f"Could not fit iota profile up to target tolerance {self.iota_tol} "
                        f"with {self.knots_per_sgrid_element} knots per sgrid element. "
                        f"Current fit error: {self.rms_iota_fit}"
                    )
                    break
                if bspline.c.size * 2 > IOTA_COEFS_MAX:
                    self.logger.warning(
                        f"Could not fit iota profile up to target tolerance {self.iota_tol} "
                        f"with {bspline.c.size} coefficients (max is {IOTA_COEFS_MAX}). "
                        f"Current fit error: {self.rms_iota_fit}"
                    )
                    break
                self.logger.debug(
                    f"Increasing knots per sgrid element {self.knots_per_sgrid_element} "
                    f"-> {self.knots_per_sgrid_element * 2} to fit iota profile with tolerance {self.iota_tol}. "
                    f"Current fit error: {self.rms_iota_fit}"
                )
                self.knots_per_sgrid_element *= 2
            delta_iota = ev.iota - bspline(rho**2)
            # ignore near-axis region as we CURRENTLY have accuracy problems there!!!
            delta_iota = delta_iota.sel(rho=slice(IOTA_MIN_RHO, None))
            self.rms_iota = np.sqrt((delta_iota**2).mean().item())
            delta_Itor = ev.I_tor - self.I_tor(rho**2)
            # set new iota profile
            self._state_parameters["iota"] = {
                "type": "bspline",
                "coefs": bspline.c,
                "knots": bspline.t,
            }

        # diagnostics
        rho = np.linspace(0, 1, DIAGNOSTICS_NPOINTS)
        Qs = ["F_r_avg"]
        if self.curr_constraint:
            Qs += ["iota", "iota_curr_0", "iota_0", "I_tor"]
        ev = self.state.evaluate(*Qs, rho=rho, theta="int", zeta="int")[Qs]

        if self.curr_constraint and self._state_parameters["picard_current"] != "off":
            self.logger.info(
                f"max Δiota: {np.abs(delta_iota).max().item():.2e}, "
                f"rms Δiota: {self.rms_iota:.2e}, "
                f"iota_tol: {self.iota_tol:.2e}"
            )
            self.logger.info(f"max ΔItor: {np.abs(delta_Itor).max().item():.2e}")
            self.logger.debug(
                f"rms Δiota_fit: {self.rms_iota_fit:.2e}, "
                f"with {self.knots_per_sgrid_element} knots per sgrid element"
            )
        logfiles = sorted(self.rundir.glob("logMinimizer_*"))
        if logfiles:
            logfile = logfiles[-1]
            log_df = read_csv(logfile, sep=",", header=0)

            diag_run = xr.Dataset(
                dict(
                    F_r_avg=ev.F_r_avg,
                    gvec_iterations=iterations,
                    gvec_max_iterations=max_iterations,
                    gvec_tolerance=tolerance,
                )
            )
            diag_minimizer = xr.Dataset(
                dict(
                    total_iteration=iteration_offset
                    + np.array(log_df["#iterations"], dtype=int),
                    force_X1=("total_iteration", np.array(log_df["normF_X1"])),
                    force_X2=("total_iteration", np.array(log_df["normF_X2"])),
                    force_LA=("total_iteration", np.array(log_df["normF_LA"])),
                    W_MHD=("total_iteration", np.array(log_df["W_MHD3D"])),
                )
            )
            if self.curr_constraint:
                diag_run["iota"] = ev.iota
                diag_run["I_tor"] = ev.I_tor
                diag_run["iota_delta"] = ev.iota - bspline(rho**2)
                diag_run["I_tor_delta"] = ev.I_tor - self.I_tor(rho**2)
                diag_run["rms_iota"] = self.rms_iota
                diag_run["rms_iota_fit"] = self.rms_iota_fit
                diag_run["max_Itor"] = np.abs(delta_Itor).max().item()
                diag_run["iota_dofs"] = bspline.c.size
            diag_run = diag_run.drop_vars(["pol_weight", "tor_weight"])
            if self.diagnostics_run is None:
                diag_run = diag_run.expand_dims(dict(run=[self.nth_run]))
                self.diagnostics_run = diag_run
            else:
                diag_run = diag_run.expand_dims(dict(run=[self.diagnostics_run.run.size]))
                self.diagnostics_run = xr.concat(
                    [self.diagnostics_run, diag_run],
                    dim="run",
                    join="outer",
                    coords="different",
                    compat="equals",
                )
            if self.diagnostics_minimizer is None:
                diag_minimizer.force_X1.attrs = dict(
                    long_name="absolute MHD force on X1", symbol=r"|F_{X^1}|"
                )
                diag_minimizer.force_X2.attrs = dict(
                    long_name="absolute MHD force on X2", symbol=r"|F_{X^2}|"
                )
                diag_minimizer.force_LA.attrs = dict(
                    long_name="absolute MHD force on lambda", symbol=r"|F_{\lambda}|"
                )
                diag_minimizer.W_MHD.attrs = dict(
                    long_name="total MHD energy", symbol=r"W_{MHD}"
                )
                self.diagnostics_minimizer = diag_minimizer
            else:
                self.diagnostics_minimizer = xr.concat(
                    [self.diagnostics_minimizer, diag_minimizer], dim="total_iteration"
                )

            forces = self.diagnostics_minimizer.isel(total_iteration=-1)[
                ["force_X1", "force_X2", "force_LA"]
            ]
            self.max_force = max([forces[force].item() for force in forces])

            end_time = time.time()

            self.logger.info(
                f"|force| = {self.max_force:.2e} (minimize_tol = {self._state_parameters['minimize_tol']:.2e})"
            )
        else:
            end_time = time.time()
            self.max_force = np.inf
        self.logger.info(
            f"GVEC run took {end_time - start_time:5.1f} seconds for {iterations} iterations. (max {max_iterations}, tol {tolerance:.1e})"
        )
        self.logger.info(
            f"GVEC iterations used in total: {self.GVEC_iter_used} / {self.totaliter}"
        )
        self.logger.info("-" * 40)

    def _reset_params_to_original(self):
        """
        Reset the parameters to the original values. Except for `iota` and `maxIter`, which is limited by `totaliter`.
        """
        params = self.parameters.copy()
        if "stages" in params:
            del params["stages"]
        params["iota"] = self._state_parameters["iota"]
        params["maxIter"] = min(
            self.totaliter - self.GVEC_iter_used,
            self.parameters.get("maxIter", self.totaliter),
        )
        self._state_parameters = params

    def _set_params_for_stage(self, stage: Mapping):
        """
        Set the run parameters to the values specified in the stage.

        Parameters
        ----------
        stage : Mapping
            Dictionary specifying which parameters are to be changed from the original parameter set.
        """
        stage = cidict(stage)
        set_I_tor = False

        for key, value in stage.lower_items():
            if key == "maxiter":
                self._state_parameters[key] = min(self.totaliter - self.GVEC_iter_used, value)

            if key == "picard_current":
                if value == "off":
                    pass
                elif (
                    value != "off"
                    and "I_tor" not in self._state_parameters
                    and "I_tor" not in stage
                ):
                    raise gvec.errors.MissingParameterError(
                        "Expected 'I_tor' in the parameters since 'picard_current' is not 'off'."
                        + " Please set 'picard_current' to 'off' if you want to use a fixed 'iota' profile or provide 'I_tor'."
                    )
                elif isinstance(value, Mapping):
                    if key not in self._state_parameters:
                        self._state_parameters[key] = cidict()
                    for subkey, subvalue in value.items():
                        self._state_parameters[key][subkey] = subvalue
                else:
                    raise gvec.errors.InvalidParameterError(
                        f"unknown picard_current value! {value}"
                    )

            if key in ["iota", "pres", "sgrid", "i_tor"]:
                if key not in self._state_parameters:
                    self._state_parameters[key] = cidict()
                if key == "i_tor":
                    set_I_tor = True
                for subkey, subvalue in value.items():
                    self._state_parameters[key][subkey] = subvalue
            if key in self._state_parameters and isinstance(value, Mapping):
                for subkey, subvalue in value.items():
                    self._state_parameters[key][subkey] = subvalue
            else:
                self._state_parameters[key] = value

        if set_I_tor:
            self._set_I_tor(self._state_parameters)

    def run(
        self,
        keep_intermediates: Literal["all", "stages"] | None = None,
    ):
        """Sequentially run the stages of the Run object.

        Parameters
        ----------
        keep_intermediates : Literal["all", "stages"] | None, optional
            Whether to keep intermediate results of GVEC. With `"all"`, all intermediate results are kept. With `"stages"`,
            only the final restarts from each stage are kept. With `None`, all intermediate results are deleted. The default is None.

        Returns
        -------
        tuple:
            The final state (gvec.State) and the diagnostics dataset (xr.DataSet)

        Raises
        ------
        gvec.errors.InvalidParameterError
            If keep_intermediates is not None, "stages" or "all"
        gvec.errors.InvalidParameterError
            If stages are set when 'picard_current="auto"'
        gvec.errors.MissingParameterError
            If 'iota_tol' is not specified when 'I_tor' is provided.
        gvec.errors.InvalidParameterError
            If 'picard_current.target' is not properly specified.
        """
        if keep_intermediates and keep_intermediates not in ["all", "stages"]:
            raise gvec.errors.InvalidParameterError(
                f"""'keep_intermediates' has to be either None, "stages" or "all" but is {keep_intermediates}"""
            )

        start_time = time.time()
        for s, stage in enumerate(self.stages):
            if self.GVEC_iter_used >= self.totaliter:
                warnings.warn("Maximum number of GVEC iterations reached!")
                break

            self.nth_stage = s
            self.nth_run = 0

            self._reset_params_to_original()
            self._set_params_for_stage(stage)
            self.knots_per_sgrid_element = IOTA_KNOTS_PER_ELEMENT_MIN

            M1, N1 = self._state_parameters["X1_mn_max"]
            M2, N2 = self._state_parameters["X2_mn_max"]
            Ml, Nl = self._state_parameters["LA_mn_max"]
            K = self._state_parameters["sgrid"]["nelems"]
            if (M1 == M2 == Ml) and (N1 == N2 == Nl):
                hM, hN, hK = compute_hfactor(M1, N1, K)
                self.logger.info(f"h-factors: pol={hM:.2e}, tor={hN:.2e}, rad={hK:.2e}")
            else:
                hM1, hN1, hK1 = compute_hfactor(M1, N1, K)
                hM2, hN2, _ = compute_hfactor(M2, N2, K)
                hMl, hNl, _ = compute_hfactor(Ml, Nl, K)
                self.logger.info(
                    f"h-factors: pol(X1)={hM1:.2e}, pol(X2)={hM2:.2e}, pol(LA)={hMl:.2e}, tor(X1)={hN1:.2e}, tor(X2)={hN2:.2e}, tor(LA)={hNl:.2e}, rad={hK1:.2e}"
                )

            # run the stage
            if self.curr_constraint:
                if self._state_parameters["picard_current"] == "auto":
                    raise gvec.errors.InvalidParameterError(
                        'Detected `picard_current="auto"` during stage evaluation. Auto mode has to be set outside of the stages.'
                    )
                if self._state_parameters["picard_current"] == "off":
                    self._run_stage_target_force(keep_intermediates=keep_intermediates)
                else:
                    if "iota_tol" not in self._state_parameters["picard_current"]:
                        raise gvec.errors.MissingParameterError(
                            f"During stage {s} 'iota_tol' is not specified."
                        )
                    target = self._state_parameters["picard_current"].get(
                        "target", "iota_and_force"
                    )
                    match target:
                        case "iota":
                            self._run_stage_target_iota(keep_intermediates=keep_intermediates)
                        case "iota_and_force":
                            self._run_stage_target_iota_and_force(
                                keep_intermediates=keep_intermediates
                            )
                        case _:
                            raise gvec.errors.InvalidParameterError(
                                f"Unknown picard_current target:{target}"
                            )
            else:
                self._run_stage_target_force(keep_intermediates=keep_intermediates)

        self.logger.info("Done.")
        end_time = time.time()

        final_iota_str = ""
        if self.curr_constraint:
            final_iota_str += f"\n and rms Δiota = {self.rms_iota:.2e}"
            # in case current constrained was turned off during the final stage
            if self._state_parameters["picard_current"] != "off":
                final_iota_str += (
                    f"(iota_tol={self._state_parameters['picard_current']['iota_tol']:.2e})"
                )

        final_message = (
            f"GVEC finished after {end_time - start_time:5.1f} seconds",
            f" using {self.GVEC_iter_used} iterations (totalIter = {self.totaliter})",
            f" with |force| = {self.max_force:.2e} (minimize_tol = {self._state_parameters['minimize_tol']:.2e})",
            final_iota_str,
        )
        if self.quiet:
            self.logger.info("".join(final_message))
        else:
            print(*final_message)
        final_statefile = Path(self._state_parameters["ProjectName"] + "_State_final.dat")
        final_parameter_file = Path(
            "parameter_" + self._state_parameters["ProjectName"] + "_final.ini"
        )

        shutil.copy(self.state.statefile, final_statefile)
        parameters_final = gvec.util.read_parameter_file_ini(
            self.state.statefile.parents[0] / "parameter.ini"
        )
        parameters_final["maxIter"] = -1
        for key in parameters_final:
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ]:
                parameters_final[key] = Path(self.parameters[key]).name
        gvec.util.write_parameter_file_ini(
            parameters=parameters_final,
            path=final_parameter_file,
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage} run {self.nth_run})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )
        self.state = gvec.State(final_parameter_file.absolute(), final_statefile.absolute())
        if keep_intermediates is None:
            shutil.rmtree(self.project_dir)
            self.project_dir = None

        diagnostics = xr.merge([self.diagnostics_run, self.diagnostics_minimizer])
        return (self.state, diagnostics)

    def _run_stage_target_iota(
        self,
        keep_intermediates: Literal["all", "stages"] | None = None,
    ):
        """
        Target only iota in the current optimization, ignoring the forces.

        This method typically performs many picard iterations with few GVEC iterations to find a suitable iota profile for the current optimization.

        Parameters
        ----------
        keep_intermediates : Literal["all", "stages"] | None, optional
            Whether to keep intermediate results of GVEC. With `"all"`, all intermediate results are kept. With `"stages"`,
            only the final restarts from each stage are kept. With `None`, all intermediate results are deleted. The default is None.

        Returns
        -------
        n_runs_in_stage: ArrayLike
            Updated number of runs completed in the current stage.
        """

        self.rms_iota = np.inf
        self.rms_iota_fit = 0.0
        self.nth_run = -1
        max_restarts = self._state_parameters["picard_current"].get(
            "maxRestarts", DEFAULT_PICARD_CURRENT_MAXRESTARTS
        )
        iota_tol = self._state_parameters["picard_current"]["iota_tol"]
        while (self.rms_iota > iota_tol) and (self.GVEC_iter_used < self.totaliter):
            if self.nth_run + 1 > max_restarts:
                warnings.warn(
                    f"Maximum number of restarts reached for stage {self.nth_stage}! Moving on to next stage."
                )
                break
            self._state_parameters["maxIter"] = min(
                self.totaliter - self.GVEC_iter_used, self._state_parameters["maxIter"]
            )
            self.nth_run += 1
            self.n_runs_in_stage[self.nth_stage] += 1
            self._print_progress()
            self.run_single_minimization()
            rm_dir = self.project_dir / Path(f"{self.nth_stage:1d}-{(self.nth_run - 1):02d}")
            if rm_dir.exists() and (
                keep_intermediates == "stages" or keep_intermediates is None
            ):
                shutil.rmtree(rm_dir)

        if self.rms_iota > iota_tol:
            warnings.warn(
                f"Targeted iota has not been reached during stage {self.nth_stage}!\n"
                f"iota_tol.: {iota_tol:.2e}, rms Δiota.: {self.rms_iota:.2e}, rms Δiota_fit: {self.rms_iota_fit:.2e}\n"
                f"GVEC iterations used: {self.GVEC_iter_used}"
            )

    def _run_stage_target_iota_and_force(
        self, keep_intermediates: Literal["all", "stages"] | None = None
    ):
        """
        Run GVEC until the force tolerance is reached and perform picard iterations until the iota tolerance is reached.
        The maximum number of total GVEC iterations and the maximum number of picard iterations are limited by totaliter.

        Parameters
        ----------
        keep_intermediates : Literal["all", "stages"] | None, optional
            Whether to keep intermediate results of GVEC. With `"all"`, all intermediate results are kept. With `"stages"`,
            only the final restarts from each stage are kept. With `None`, all intermediate results are deleted. The default is None.

        Returns
        -------
        n_runs_in_stage: ArrayLike
            Updated number of runs completed in the current stage.
        """
        self.rms_iota = np.inf
        self.rms_iota_fit = 0.0
        self.nth_run = -1
        max_restarts = self._state_parameters["picard_current"].get(
            "maxRestarts", DEFAULT_PICARD_CURRENT_MAXRESTARTS
        )
        self.logger.debug(f"maxRestarts: {max_restarts}")
        iota_tol = self._state_parameters["picard_current"]["iota_tol"]
        while (self.GVEC_iter_used < self.totaliter) and (self.rms_iota > iota_tol):
            self.logger.debug(f"nth run: {self.nth_run}")
            if self.nth_run + 1 > max_restarts:
                warnings.warn(
                    f"Maximum number of restarts reached for stage {self.nth_stage}! Moving on to next stage."
                )
                break
            self._state_parameters["maxIter"] = min(
                self.totaliter - self.GVEC_iter_used, self._state_parameters["maxIter"]
            )
            self.nth_run += 1
            self.n_runs_in_stage[self.nth_stage] += 1
            self._print_progress()
            self.run_single_minimization()
            rm_dir = self.project_dir / Path(f"{self.nth_stage:1d}-{(self.nth_run - 1):02d}")
            if rm_dir.exists() and (
                keep_intermediates == "stages" or keep_intermediates is None
            ):
                shutil.rmtree(rm_dir)
        if self.rms_iota > iota_tol:
            warnings.warn(
                f"Targeted iota has not been reached during stage {self.nth_stage}!\n"
                f"iota_tol.: {iota_tol:.2e}, rms Δiota.: {self.rms_iota:.2e}, rms Δiota_fit: {self.rms_iota_fit:.2e}\n"
                f"GVEC iterations used: {self.GVEC_iter_used}"
            )
        if self.max_force > self._state_parameters["minimize_tol"]:
            warnings.warn(
                f"Force tolerance was not reached in stage {self.nth_stage}!\n"
                f"max|force|: {self.max_force:.2e}, minimize_tol: {self._state_parameters['minimize_tol']:.2e}"
            )

    def _run_stage_target_force(
        self, keep_intermediates: Literal["all", "stages"] | None = None
    ):
        self._state_parameters["maxIter"] = min(
            self.totaliter - self.GVEC_iter_used,
            self._state_parameters["maxIter"],
        )
        self._print_progress()
        self.run_single_minimization()
        self.n_runs_in_stage[self.nth_stage] += 1
        rm_dir = self.project_dir / Path(f"{self.nth_stage:1d}-{(self.nth_run - 1):02d}")
        if (
            self.max_force > self._state_parameters["minimize_tol"]
            and self._state_parameters["maxIter"] > 0
        ):
            warnings.warn(
                f"Force tolerance was not reached! |force|: {self.max_force:.2e}, minimize_tol: {self._state_parameters['minimize_tol']:.2e}"
            )
        if rm_dir.exists() and (keep_intermediates == "stages" or keep_intermediates is None):
            shutil.rmtree(rm_dir)

    def _print_progress(self):
        """
        Evaluate and print a progress string for the current stage and run.
        """
        progressstr = "|"
        for i, ir in enumerate(self.n_runs_in_stage):
            if i < self.nth_stage:
                progressstr += "=" * ir + "|"
            elif i == self.nth_stage:
                progressstr += "=" * (ir - 1) + ">" + "|"
            else:
                progressstr += ".|"
        start_str = "GVEC"
        state_str = ""
        if len(self.stages) > 1:
            state_str += f" - completed {self.nth_stage}/{len(self.stages)} stages"
        restart_str = ": "
        if self.curr_constraint:
            restart_str = f", restarts in current stage - {self.nth_run}" + restart_str
        progressstr = start_str + state_str + restart_str + progressstr
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(progressstr)
        elif not self.quiet:
            print(progressstr, end="\r")

    def plot_diagnostics_run(self):
        import matplotlib.pyplot as plt

        diagnostics = self.diagnostics_run
        if self.curr_constraint:
            fig, axs = plt.subplots(1, 3, figsize=(12, 3), layout="constrained", sharex=True)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(5, 3), layout="constrained", sharex=True)
            axs = [ax]
        axs[0].plot(
            diagnostics.run,
            [
                np.sqrt((diagnostics.sel(run=run).F_r_avg ** 2).mean(dim="rad"))
                for run in diagnostics.run
            ],
            ".-",
        )
        axs[0].set(
            xlabel="restart number",
            ylabel=f"rms ${diagnostics.F_r_avg.attrs['symbol']}$",
            title=diagnostics.F_r_avg.attrs["long_name"],
            yscale="log",
        )
        if self.curr_constraint:
            axs[1].plot(diagnostics.run, diagnostics.rms_iota, ".-")
            axs[1].plot(diagnostics.run, diagnostics.rms_iota_fit, ".-")
            axs[1].set(
                xlabel="restart number",
                ylabel=r"$\sqrt{\sum \left(\Delta\iota\right)^2}$",
                title=f"Difference to target {diagnostics.iota.attrs['long_name']}\nroot mean square\nprofile error (blue) and fit error (orange)",
                yscale="log",
            )
            axs[2].plot(diagnostics.run, diagnostics.max_Itor, ".-")
            axs[2].set(
                xlabel="restart number",
                ylabel=r"$\max \left(\Delta I_{tor}\right)$",
                title=f"Difference to target\n{diagnostics.I_tor.attrs['long_name']}",
                yscale="log",
            )
        return fig

    def plot_diagnostics_current_profiles(self):
        import matplotlib.pyplot as plt

        diagnostics = self.diagnostics_run

        fig, axs = plt.subplots(2, 2, figsize=(15, 5), tight_layout=True, sharex=True)
        for r in diagnostics.run.data:
            if r == diagnostics.run.data[-1]:
                kwargs = dict(marker=".", color="C0", alpha=1.0, label="final profile")
            elif r == diagnostics.run.data[0]:
                kwargs = dict(color="black", linestyle="--", alpha=0.5, label="initial profile")
            else:
                kwargs = dict(color="black", alpha=0.2 + 0.3 * (r / diagnostics.run.data[-1]))
                if r == diagnostics.run.data[1]:
                    kwargs["label"] = "intermediate profiles"
            d = diagnostics.sel(run=r)
            axs[0, 0].plot(d.rho**2, d.iota, **kwargs)
            axs[1, 0].plot(d.rho**2, np.abs(d.iota_delta), **kwargs)
            axs[0, 1].plot(d.rho**2, d.I_tor, **kwargs)
            axs[1, 1].plot(d.rho**2, np.abs(d.I_tor_delta), **kwargs)
        for i, var in enumerate(["iota", "I_tor"]):
            axs[0, i].set(
                title=diagnostics[var].attrs["long_name"],
                ylabel=f"${diagnostics[var].attrs['symbol']}$",
            )
            axs[1, i].set(
                title=f"Difference to target {diagnostics[var].attrs['long_name']}",
                xlabel=r"$\rho^2$",
                ylabel=rf"$|\Delta {diagnostics[var].attrs['symbol']}|$",
                yscale="log",
            )
        axs[0, 0].legend()
        return fig

    def plot_diagnostics_minimization(self):
        import matplotlib.pyplot as plt

        diagnostics = self.diagnostics_run
        if self.curr_constraint:
            fig, axs = plt.subplots(3, 1, figsize=(10, 5), layout="constrained", sharex=True)
            axs[2].plot(
                np.cumsum(diagnostics.gvec_iterations),
                diagnostics.rms_iota,
                ".-",
            )
            axs[2].set(ylabel=r"$\Delta\iota_{rms}$", yscale="log")
        else:
            fig, axs = plt.subplots(2, 1, figsize=(10, 5), tight_layout=True, sharex=True)
        axs[0].plot(
            self.diagnostics_minimizer.total_iteration, self.diagnostics_minimizer.W_MHD, "-"
        )
        axs[0].set(
            ylabel=f"${self.diagnostics_minimizer.W_MHD.attrs['symbol']}$",
        )
        axf = axs[1]

        axf.plot(
            self.diagnostics_minimizer.total_iteration,
            self.diagnostics_minimizer.force_X1,
            color="C0",
            label=r"$X^1$",
        )
        axf.plot(
            self.diagnostics_minimizer.total_iteration,
            self.diagnostics_minimizer.force_X2,
            color="C1",
            label=r"$X^2$",
        )
        axf.plot(
            self.diagnostics_minimizer.total_iteration,
            self.diagnostics_minimizer.force_LA,
            color="C2",
            label=r"$\lambda$",
        )

        if len(self.stages) > 1:
            # stages vlines
            for ax in axs:
                # zeroth "initial" stage does not perform a run:
                # -> offset in n_runs_in_stage but NOT in gvec_iterations
                n_runs_till_stage = np.cumsum(self.n_runs_in_stage[1:-1])
                ax.vlines(
                    [np.sum(diagnostics.gvec_iterations[:i]) for i in n_runs_till_stage],
                    *ax.get_ylim(),
                    colors="grey",
                    linestyle="solid",
                    alpha=0.6,
                    zorder=-1000,
                    label="stages",
                )

        # runs vlines
        if np.sum(np.array(self.n_runs_in_stage) > 0) >= 2:
            axf.vlines(
                np.cumsum(diagnostics.gvec_iterations[:-1]),
                *axf.get_ylim(),
                colors="k",
                linestyle="dashed",
                alpha=0.6,
                zorder=-900,
                label="restarts",
            )

        axf.set(ylabel="|Force|", yscale="log")
        axf.yaxis.grid(True, linestyle="--", alpha=0.5)
        axf.legend(bbox_to_anchor=(1.0, 0.5), loc="center left")
        axs[-1].set(xlabel="GVEC iteration")

        return fig


def auto_generate_stages(parameters: cidict, Itor: bool):
    target_minimize_tol = parameters.get("minimize_tol", DEFAULT_MINIMIZE_TOL)
    stages = []

    if Itor:
        parameters["picard_current"] = dict(
            iota_tol=AUTO_PICARD_CURRENT_IOTA_TOL,
        )

    # initial stage
    stages.append(cidict(maxIter=0, init_LA=parameters.get("init_LA", True)))
    parameters["init_LA"] = False
    # move boundary perturbation to initial stage
    if parameters.get("boundary_perturb", False):
        stages[0]["boundary_perturb"] = True
        parameters["boundary_perturb"] = False
    # move init from vmec to initial stage
    if parameters.get("whichInitEquilibrium", 0) == 1:
        stages[0]["whichInitEquilibrium"] = 1
        parameters["whichInitEquilibrium"] = 0
    if Itor:
        stages[0]["picard_current"] = "off"

    # iota initialization stage
    if Itor:
        stages.append(
            cidict(
                minimize_tol=AUTO_MINIMIZE_TOL_MIN,
                maxIter=AUTO_STAGE_IOTA_INIT_MAXITER,
                picard_current=cidict(
                    target="iota",
                    iota_tol=AUTO_STAGE_IOTA_INIT_IOTA_TOL,
                ),
            )
        )

    # ramp up minimize_tol
    minimize_tol_seq = 10 ** np.arange(
        np.log10(AUTO_MINIMIZE_TOL_MIN), np.log10(target_minimize_tol), -1
    )
    for minimize_tol in minimize_tol_seq:
        stages.append(
            cidict(
                minimize_tol=minimize_tol,
            )
        )
    # final minimize_tol
    stages.append(cidict())

    # set iota_tol in final stage
    if Itor:
        stages[-1]["picard_current"] = cidict(
            iota_tol=AUTO_PICARD_CURRENT_IOTA_TOL_TARGET,
        )
    return parameters, stages


def _linspace_sgrid(sgrid, n_points):
    """Return a refined radial grid built from a coarse `sgrid` partition.

    The interval between each pair of neighboring grid points is sampled with
    ``n_points`` uniformly spaced points. The right endpoint of each interval is
    excluded so that adjacent intervals do not duplicate shared boundaries. A
    final point at ``1`` is appended to close the grid.

    Parameters
    ----------
    sgrid : array_like
        Monotonically increasing radial grid points in the closed interval ``[0, 1]``.
    n_points : int
        Number of equally spaced samples to generate inside each grid interval.

    Returns
    -------
    numpy.ndarray
        One-dimensional array containing the refined grid points in ascending
        order.

    Notes
    -----
    This helper preserves the original cell boundaries from ``sgrid`` while
    increasing the resolution inside each cell. It is primarily used to sample
    and fit radial profiles on a finer mesh than the staging grid.
    """
    return np.concatenate(
        [
            np.linspace(left, right, n_points, endpoint=False)
            for left, right in zip(sgrid[:-1], sgrid[1:])
        ]
        + [np.ones(1)]
    )


def _add_knot_multiplicity(knots, degree):
    """Pad a knot vector with boundary multiplicity for B-spline fitting.

    The input knot vector is assumed to contain only the interior knots. This
    helper prepends and appends ``degree`` copies of ``0`` and ``1`` so that
    the resulting knot vector has the clamped boundary multiplicity required
    by B-Splines.

    Parameters
    ----------
    knots : array_like
        Knot locations without multiplicity, in the closed interval ``[0, 1]``.
    degree : int
        B-spline degree. The same value is used as the multiplicity of both
        boundary knots.

    Returns
    -------
    numpy.ndarray
        Knot vector with repeated boundary knots added at both ends.

    Notes
    -----
    The returned knot vector is suitable for clamped least-squares spline fits
    where the end conditions are controlled by boundary knot repetition.
    """
    return np.concatenate([np.zeros(degree)] + [knots] + [np.ones(degree)])


def _evaluate_extrapolate_in_rho2(state, *Qs, rho, theta, zeta):
    rho_aux = np.sqrt([1e-8, 2e-8, 3e-8])
    ev_aux = state.evaluate(*Qs, rho=rho_aux, theta=theta, zeta=zeta, extrapolate=False)
    r1 = ev_aux.isel(rad=0)
    r2 = ev_aux.isel(rad=1)
    r3 = ev_aux.isel(rad=2)
    on_axis = 3 * (r1 - r2) + r3
    on_axis.rho[...] = 0.0

    ev = state.evaluate(*Qs, rho=rho, theta=theta, zeta=zeta, extrapolate=False)
    for q in ev.data_vars:
        if "rad" in ev[q].dims:
            ev[q].loc[dict(rho=0.0)] = on_axis[q]
    return ev


def compute_hfactor(M, N, K):
    hM = 2 * np.pi / (2 * M + 1)
    hN = 2 * np.pi / (2 * N + 1)
    hK = 1.0 / K
    return hM, hN, hK


def fortran_run(
    parameterfile: str | Path,
    restartfile: str | Path | None = None,
    MPIcomm: int | None = None,
    stdout_path: str | Path | None = "stdout.txt",
):
    """
    Run gvec from python

    Parameters
    ----------
    parameterfile : str
        Path to / name of parameter file
    restartfile : str
        Path to / name of GVEC restart file, optional
    MPIcomm : int
        MPI communicator, optional (default in GVEC (if compiled with MPI) is MPI_COMM_WORLD)
    stdout_path : str
        Path to / name of file to redirect the standard output of GVEC. Optional, default is "stdout.txt".
        If set to None, stdout is not redirected
    """
    logger = logging.getLogger("gvec.run")
    logger.debug(f"Running GVEC with parameter file: {parameterfile}")
    if gvec.core.state._binding.bound_state is not None:
        gvec.core.state._binding.bound_state.unbind()

    _binding.redirect_abort()
    if stdout_path is not None:
        _binding.redirect_stdout(str(stdout_path))
        logger.debug(f"Redirecting GVEC stdout to {stdout_path}")

    if not Path(parameterfile).exists():
        raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
    MAXLEN = _globals.maxlen
    with open(parameterfile, "r") as f:
        content = f.readlines()
        if any([len(line) > MAXLEN for line in content]):
            raise ValueError(
                f"Parameter file {parameterfile} contains lines longer than {MAXLEN} characters."
            )
    if restartfile is not None:
        if not Path(restartfile).exists():
            raise FileNotFoundError(f"Restart file {restartfile} does not exist.")

    try:
        with catch_gvec_errors():
            _run.start_rungvec(str(parameterfile), restartfile_in=restartfile, comm_in=MPIcomm)
    except:
        logger.info("attempting cleanup")
        _run.cleanup()
        logger.debug("cleanup done")
        raise
