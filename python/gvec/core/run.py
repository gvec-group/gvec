# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""run gvec from python"""

import logging
import shutil
import time
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Literal

import numpy as np
import xarray as xr
from pandas import read_csv

import gvec
from gvec.core.state import State
from gvec.core.run_state import RunState
from gvec.util import CaseInsensitiveDict as cidict


DEFAULT_MINIMIZE_TOL = 1e-6  # different to default `minimize_tol` in fortran
DEFAULT_TOTALITER = 10**5  # different to default `maxIter` in fortran
AUTO_IOTA_TARGET = 1e-10


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

        # params: user provided parameters
        # self.parameters: global/fallback parameters
        # _state_parameters reflect the state of the Run, that is they change during each stage and/or restart
        #   after a stage they are reset to parameters. They are NOT state.parameters, as those reflect the state of a
        #   finished restart.
        if "ProjectName" not in self.parameters:
            self.parameters["ProjectName"] = "GVEC"

        self.logger = logging.getLogger("gvec.run")

        if "picard_current" not in self.parameters:
            self.parameters["picard_current"] = "off"
            self.curr_constraint = False

        match self.parameters["picard_current"]:
            case "off":
                self.curr_constraint = False
                self.parameters["picard_current"] = cidict()
            case "auto":
                self.curr_constraint = True
                self.parameters["picard_current"] = cidict()
            case cidict():
                self.curr_constraint = True

        self.filetype = parameter_format

        project_dir = Path(f"{self.parameters['ProjectName']}_gvec_stages")
        if project_dir.exists():
            self.logger.debug(f"Removing existing run directory {project_dir}")
            shutil.rmtree(project_dir)
        project_dir.mkdir()
        self.project_dir = project_dir

        self.nth_stage = 0
        self.rundir: Path = None
        self.quiet = quiet
        self.diagnostics_stage: xr.Dataset = None
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

        if "I_tor" in self.parameters:
            self.curr_constraint = True
        else:
            self.curr_constraint = False

        if "stages" in self.parameters:
            self.stages = self.parameters["stages"]
            del self.parameters["stages"]
        else:
            # initial stage
            self.stages = [cidict()]
            self.stages[0]["maxIter"] = 0
            self.stages[0]["init_LA"] = self.parameters.get("init_LA", True)
            self.parameters["init_LA"] = False
            # move boundary perturbation to initial stage
            if self.parameters.get("boundary_perturb", False):
                self.stages[0]["boundary_perturb"] = True
                self.parameters["boundary_perturb"] = False
            # move init from vmec to initial stage
            if self.parameters.get("whichInitEquilibrium", 0) == 1:
                self.stages[0]["whichInitEquilibrium"] = 1
                self.parameters["whichInitEquilibrium"] = 0

            # ensure at least one stage beyond initial stage
            self.stages += [cidict()]

        parameters_stages = self.parameters.copy()
        parameters_stages["stages"] = self.stages
        parameters_stages_name = (
            f"parameter_{self.parameters['ProjectName']}.stages.{self.filetype}"
        )
        self.logger.info(f"... generated {len(self.stages)} stages.")

        # load restart state
        if isinstance(state, gvec.State):
            self.logger.info(f"reading restart state from {state.statefile}")
            self.state = state

            if self.parameters.get("which_hmap", 1) != self.state.parameters.get(
                "which_hmap", 1
            ):
                self.logger.warning(
                    f"restarting with hmap={self.parameters.get('which_hmap', 1)} from hmap={self.state.parameters.get('which_hmap', 1)}."
                )

            # set initial iota profile
            if "iota" not in self.parameters:
                if "I_tor" in self.state.parameters:
                    ev = self.state.evaluate("iota", rho=np.linspace(0, 1, 51))
                    self.parameters["iota"] = {
                        "type": "interpolation",
                        "vals": ev.iota.data,
                        "rho2": (ev.rho**2).data,
                    }
                else:
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
            self._set_I_tor_target(self.parameters)
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

        self._state_parameters = self.parameters.copy()

        # load restart state (2nd part)
        if not isinstance(state, gvec.State):  # Path, str or None
            if state is not None:
                self.logger.warning(
                    "restarting from statefile without associated parameterfile"
                )
            initial_param_file = self.project_dir / "parameters_initial.ini"
            gvec.util.write_parameters(self._state_parameters, initial_param_file)
            self.state = gvec.State(initial_param_file, state)

        self.logger.info(f"... writing new parameters to '{parameters_stages_name}'")
        gvec.util.write_parameters(
            parameters_stages,
            project_dir / parameters_stages_name,
        )

    def _set_I_tor_target(self, params: Mapping):
        """Evaluate and set the target toroidal current profile at linearly spaced positions in rho.

        Raises
        ------
        gvec.errors.InvalidParameterError
            If an unknown profile type is provided.
        """

        match params["I_tor"].get("type", "polynomial"):
            case "polynomial":
                coefs = np.array(params["I_tor"]["coefs"][::-1])
                self.logger.debug(f"polynomial coefs: {coefs}")
                coefs *= params["I_tor"].get("scale", 1.0)
                self.I_tor_target = np.poly1d(coefs)
                if (
                    abs(coefs[-1]) > 1e-8
                ):  # poly1d is reverse to GVEC, e.g. coefs is ordered x²+x+1
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis!  I_tor(rho=0): {coefs[-1]}"
                    )

            case "bspline":
                from scipy.interpolate import BSpline

                coefs = np.array(params["I_tor"]["coefs"], dtype=float)
                coefs *= params["I_tor"].get("scale", 1.0)
                knots = np.array(params["I_tor"]["knots"], dtype=float)
                deg = np.sum(knots == knots[0]) - 1
                I_tor_bspl = BSpline(knots, coefs, deg)
                self.I_tor_target = I_tor_bspl
                if abs(I_tor_bspl(0.0)) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis! I_tor(rho=0): {I_tor_bspl(0.0)}"
                    )

            case "interpolation":
                from scipy.interpolate import make_splrep

                y_vals = np.array(params["I_tor"]["vals"], dtype=float)
                rho2_vals = np.array(params["I_tor"]["rho2"], dtype=float)
                if min(np.sqrt(rho2_vals)) > 1e-4:
                    rho2_vals = np.append([0], rho2_vals)
                    y_vals = np.append([0], y_vals)
                I_tor_bspl = make_splrep(rho2_vals, y_vals)
                self.I_tor_target = I_tor_bspl
                if abs(I_tor_bspl(0.0)) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis! I_tor(rho=0): {I_tor_bspl(0.0)}"
                    )

            case _:
                raise gvec.errors.InvalidParameterError(
                    f"Unknown Itor type: {params['I_tor']['type']}"
                )

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
                match value:
                    case "auto":
                        self._state_parameters[key] = cidict()
                        self.curr_constraint = True
                        continue
                    case "off":
                        self.curr_constraint = False
                        self._state_parameters[key] = cidict()
                        continue
                    case Mapping():
                        self.curr_constraint = True

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
            self._set_I_tor_target(self._state_parameters)

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
                self.logger.warning("Maximum number of GVEC iterations reached!")
                break
            self.logger.debug(f"Entering stage {s}.")
            self.nth_stage = s

            self._reset_params_to_original()
            self._set_params_for_stage(stage)

            # run the stage
            self._minimize_stage()

        self.logger.info("Done.")
        end_time = time.time()

        final_iota_str = ""
        if self.curr_constraint:
            final_iota_str += f"\n and rms Δiota = {self.state.rms_iota:.2e}"

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
        parameters_final = gvec.util.read_parameters(
            self.state.statefile.parents[0] / "parameter.ini"
        )
        # # make sure that the correct iota is in the final parameter file
        if self.curr_constraint:
            parameters_final["iota"] = self.state._parameter_dict["iota"]
        parameters_final["maxIter"] = -1
        for key in parameters_final:
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ]:
                parameters_final[key] = Path(self.parameters[key]).name
        parameters_final = gvec.util.flatten_parameters(parameters_final)
        gvec.util.write_parameter_file_ini(
            parameters=parameters_final,
            path=final_parameter_file,
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )
        self.state = gvec.State(final_parameter_file.absolute(), final_statefile.absolute())
        if keep_intermediates is None:
            shutil.rmtree(self.project_dir)
            self.project_dir = None

        diagnostics = (self.diagnostics_stage, self.diagnostics_minimizer)
        return (self.state, diagnostics)

    def _minimize_stage(self):
        self._state_parameters["maxIter"] = min(
            self.totaliter - self.GVEC_iter_used,
            self._state_parameters["maxIter"],
        )
        self._print_progress()

        self._run_single_minimization()
        if self.curr_constraint:
            self._state_parameters["iota"] = self.state._parameter_dict["iota"]
        if (
            self.max_force > self._state_parameters["minimize_tol"]
            and self._state_parameters["maxIter"] > 0
        ):
            self.logger.warning(
                f"Force tolerance was not reached! |force|: {self.max_force:.2e}, minimize_tol: {self._state_parameters['minimize_tol']:.2e}"
            )

    def _print_progress(self):
        """
        Evaluate and print a progress string for the current stage and run.
        """

        progressstr = "|"
        progressstr += "=|" * self.nth_stage
        progressstr += ">|" + ".|" * (len(self.stages) - self.nth_stage - 1)
        start_str = "GVEC"
        state_str = ""
        if len(self.stages) > 1:
            state_str += f" - completed {self.nth_stage}/{len(self.stages)} stages"
        progressstr = start_str + state_str + progressstr
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(progressstr)
        elif not self.quiet:
            print(progressstr, end="\r")

    def _run_single_minimization(self):
        """Run a single GVEC energy minimization using the current parameters, updating the run-state."""
        start_time = time.time()
        # find previous state
        if self.state.statefile:
            self.logger.debug(f"Restart from statefile {self.state.statefile}")

        # prepare the run directory
        self.rundir = self.project_dir / Path(f"stage-{self.nth_stage:1d}")
        if self.rundir.exists():
            self.logger.debug(f"Removing existing run directory {self.rundir}")
            shutil.rmtree(self.rundir)
        self.rundir.mkdir()
        self.logger.debug(f"Created run directory {self.rundir}")

        gvec.util.write_parameter_file_ini(
            gvec.util.flatten_parameters(self._state_parameters),
            self.rundir / "parameter.ini",
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )

        with gvec.util.chdir(self.rundir):
            # write parameterfile & run GVEC
            self.state = RunState(
                parameterfile="parameter.ini",
                statefile=self.state.statefile,
                I_tor=self.I_tor_target if self.curr_constraint else None,
                reset_LA=self._state_parameters["picard_current"].get("reset_LA", None),
                n_iota_rhos=self._state_parameters["picard_current"].get("nPoints", 101),
                stdout_path=Path("stdout.txt") if self.redirect_gvec_stdout else None,
            )
            iota_tol = self._state_parameters["picard_current"].get(
                "iota_tol", AUTO_IOTA_TARGET
            )
            minimize_target = self._state_parameters.get("minimize_tol", DEFAULT_MINIMIZE_TOL)
            self.state.minimize(
                iota_target=iota_tol,
                minimize_tol=minimize_target,
                maxIter=self._state_parameters["maxIter"],
                maxUpdates=self._state_parameters["picard_current"].get("maxRestarts", 30),
                target=self._state_parameters["picard_current"].get("target", "auto"),
            )

        # postprocessing
        statefiles_found = sorted(self.rundir.glob("*State*.dat"))
        if len(statefiles_found) > 0:
            self.state.statefile = Path(statefiles_found[-1]).absolute()

        iterations = self.state.get_iterations()
        iteration_offset = self.GVEC_iter_used
        self.GVEC_iter_used += iterations
        max_iterations = self._state_parameters["maxIter"]
        tolerance = self._state_parameters["minimize_tol"]
        self.logger.debug(f"Postprocessing statefile {self.state.statefile}")

        quantities = ["F_r_avg"]

        # diagnostics
        if self.curr_constraint:
            quantities += ["iota", "I_tor"]

            self.logger.info(
                f"rms Δiota: {self.state.rms_iota:.2e}"
                + (
                    f", iota_tol: {iota_tol:.2e}"
                    if self._state_parameters["picard_current"] != "off"
                    else ""
                )
            )
        if self.state.rho2_iota is not None:
            self.rho = np.sqrt(self.state.rho2_iota)
        else:
            self.rho = np.linspace(0, 1, 101)

        if self.curr_constraint:
            rho_eval = np.concatenate(
                [[np.sqrt(1e-8), np.sqrt(2e-8), np.sqrt(3e-8)], self.rho[1:]]
            )
        else:
            rho_eval = self.rho

        ev = self.state.evaluate(*quantities, rho=rho_eval)
        ev = ev[quantities]

        if self.curr_constraint:
            # extrapolate ev dataset, from evaluations at s=1e-8,2e-8,3e-8 to s=0, quadratically. Only keep s=0 in dataset.
            r1 = ev.isel(rad=0)
            r2 = ev.isel(rad=1)
            r3 = ev.isel(rad=2)
            ev = ev.isel(rad=slice(2, None))

            # workaround to modify xarray coordinate & index (with pandas >=3.0)
            rho = ev.rho.data.copy()
            rho[0] = 0.0
            ev = ev.assign_coords(rho=("rad", rho)).set_xindex("rho")

            for var in ev.data_vars:
                ev[var].data[0] = 3 * (r1[var].data - r2[var].data) + r3[var].data

        logfiles = sorted(self.rundir.glob("logMinimizer_*"))
        if logfiles:
            logfile = logfiles[-1]
            log_df = read_csv(logfile, sep=",", header=0)

            diag_stage = xr.Dataset(
                data_vars=dict(
                    F_r_avg=(("stage", "rho"), ev.F_r_avg.data.reshape(1, -1)),
                    gvec_tolerance=("stage", [tolerance]),
                ),
                coords=dict(
                    iteration=("iteration", [iterations]),
                    stage=("stage", [self.nth_stage]),
                    rho=("rad", ev.rho.data),
                ),
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
                diag_stage["iota"] = ev.iota
                diag_stage["I_tor"] = ev.I_tor
                diag_stage["n_updates"] = self.state.nth_update
                diag_stage["iota_delta"] = self.state.delta_iota
                diag_stage["I_tor_delta"] = ev.I_tor - self.state.I_tor
                diag_stage.reindex(iteration=self.state.iters_per_run, fill_value=np.nan)
                diag_stage["iota_delta_picard"] = xr.DataArray(
                    data=self.state.delta_iota_over_run,
                    dims=["picard_iteration"],
                    coords={"picard_iteration": self.state.iters_per_run + iteration_offset},
                )

            if self.diagnostics_stage is None:
                self.diagnostics_stage = xr.DataTree()
                self.diagnostics_stage[f"stage_{self.nth_stage}"] = diag_stage
            else:
                self.diagnostics_stage[f"stage_{self.nth_stage}"] = diag_stage

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

    def plot_diagnostics_run(self):
        import matplotlib.pyplot as plt

        diagnostics = self.diagnostics_stage
        if self.curr_constraint:
            fig, axs = plt.subplots(1, 2, figsize=(10, 3), tight_layout=True)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(5, 3), tight_layout=True)
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
            axs[1].plot(diagnostics.run, np.sqrt((diagnostics.iota_delta**2).mean("rad")), ".-")
            axs[1].set(
                xlabel="restart number",
                ylabel=r"$\sqrt{\sum \left(\Delta\iota\right)^2}$",
                title=f"Difference to target {diagnostics.iota.attrs['long_name']}\nroot mean square",
                yscale="log",
            )
        return fig

    def plot_diagnostics_current_profiles(self):
        import matplotlib.pyplot as plt

        diagnostics = self.diagnostics_stage

        fig, axs = plt.subplots(2, 2, figsize=(15, 5), tight_layout=True, sharex=True)
        for r in diagnostics.stage.data:
            if r == diagnostics.stage.data[-1]:
                kwargs = dict(marker=".", color="C0", alpha=1.0, label="final profile")
            elif r == diagnostics.stage.data[0]:
                kwargs = dict(color="black", linestyle="--", alpha=0.5, label="initial profile")
            else:
                kwargs = dict(color="black", alpha=0.2 + 0.3 * (r / diagnostics.stage.data[-1]))
                if r == diagnostics.stage.data[1]:
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

        diagnostics = self.diagnostics_stage
        if self.curr_constraint:
            fig, axs = plt.subplots(3, 1, figsize=(10, 5), layout="constrained", sharex=True)
            for stage in diagnostics:
                ds_aux = diagnostics[stage]
                axs[2].plot(
                    ds_aux.picard_iteration,
                    ds_aux.iota_delta_picard,
                    ".-",
                    label=stage.replace("_", " "),
                )
            axs[2].set(ylabel=r"$\Delta\iota_{rms}$", yscale="log")
            axs[2].legend(bbox_to_anchor=(1.0, 0.5), loc="center left")
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

        axf.set(ylabel="|Force|", yscale="log")
        axf.yaxis.grid(True, linestyle="--", alpha=0.5)
        axf.legend(bbox_to_anchor=(1.0, 0.5), loc="center left")
        axs[-1].set(xlabel="GVEC iteration")

        return fig
