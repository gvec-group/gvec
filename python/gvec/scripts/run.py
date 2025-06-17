# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC run script for running GVEC using stages and current constraints."""

import argparse
from datetime import datetime
import logging
import re
import shutil
import time
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike
import xarray as xr
from pandas import read_csv
import gvec

# === Argument Parser === #

parser = argparse.ArgumentParser(
    prog="pygvec-run",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Run GVEC with a given parameterfile, optionally restarting from an existing statefile.\n\n"
    "When given an INI parameterfile, GVEC is called directly.\n"
    "With YAML and TOML parameterfiles, GVEC can be run in several stages and a current constraint with picard iterations can be performed.",
)
parser.add_argument("parameterfile", type=Path, help="input GVEC parameterfile")
parser.add_argument(
    "restartfile",
    type=Path,
    help="GVEC statefile to restart from (optional)",
    nargs="?",
)
param_type = parser.add_mutually_exclusive_group()
param_type.add_argument(
    "--ini",
    action="store_const",
    const="ini",
    dest="param_type",
    help="interpret GVEC parameterfile classicly (INI)",
)
param_type.add_argument(
    "--yaml",
    action="store_const",
    const="yaml",
    dest="param_type",
    help="interpret GVEC parameterfile as YAML",
)
param_type.add_argument(
    "--toml",
    action="store_const",
    const="toml",
    dest="param_type",
    help="interpret GVEC parameterfile as TOML",
)
verbosity = parser.add_mutually_exclusive_group()
verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug, -vvv for GVEC output",
)
verbosity.add_argument("-q", "--quiet", action="store_true", help="suppress output")
parser.add_argument(
    "-d",
    "--diagnostics",
    type=Path,
    default=None,
    help="output netCDF file for diagnostics",
)
parser.add_argument("-p", "--plots", action="store_true", help="plot diagnostics")

# === Script === #


class RunWithStages:
    def __init__(
        self,
        params: Mapping,
        statefile: Path | None = None,
        logger: logging.Logger | None = None,
        diagnostics: xr.Dataset | None = None,
        redirect_gvec_stdout: bool = True,
        diagnosticfile: Path | None = None,
        progressbar: bool = True,
        plots: bool = False,
    ):
        """
        State of a GVEC run during a stage, e.g. a picard iteration during a current constraint run.

        Parameters
        ----------
        params : Mapping
            GVEC parameter dictionary.
        statefile : Path | None, optional
            Statefile to restart from. The default is None.
        logger : logging.Logger | None:
            Logger to write to. The default is None.
        diagnostics : xr.Dataset | None, optional
            NetCDF file for logging diagnostic quantities during the stages and runs. The default is None.
        redirect_gvec_stdout : bool, optional
            Whether to redirect GVEC's stdout. The default is True.
        diagnosticfile : Path | None, optional
            File to write diagnostic output to. The default is default None.
        progressbar : bool, optional
            Whether to show a progress bar. The default is True.
        plots : bool, optional
            Whether to plot diagnostics. The default is False.
        """
        self.statefile = statefile
        self.original_params = params.copy()
        self.params = params.copy()
        if "ProjectName" not in self.params:
            self.params["ProjectName"] = "GVEC"
            self.original_params["ProjectName"] = "GVEC"
        if "picard_current" not in self.original_params:
            self.original_params["picard_current"] = "off"

        self.logger = logging.getLogger("pyGVEC.script")

        project_dir = Path(f"{self.params['ProjectName']}_gvec_stages")
        if project_dir.exists():
            self.logger.debug(f"Removing existing run directory {project_dir}")
            shutil.rmtree(project_dir)
        project_dir.mkdir()
        self.project_dir = project_dir

        self.nth_stage = 0
        self.nth_run = 0
        self.rundir = None
        self.progressbar = progressbar
        self.diagnostics = diagnostics
        self.diagnosticfile = diagnosticfile
        self.GVEC_iter_used = 0
        self.redirect_gvec_stdout = redirect_gvec_stdout

        self.log_dia = None
        self.plots = plots

        self.stages = params.get("stages", [{}])

        self.totaliter = params.get("totaliter", int(1e5))
        if "MaxIter" not in params:
            self.params["MaxIter"] = self.totaliter
        sgrid = gvec.util.CaseInsensitiveDict(self.params["sgrid"])
        self.has_Itor = "I_tor" in params
        self.has_iota = "iota" in params

        rho = np.sqrt(np.linspace(0, 1, max(101, 2 * sgrid["nElems"])))
        rho[0] = 1e-4
        self.rho = rho

        picard_current = params.get("picard_current", "off")
        if self.has_Itor and (picard_current == "off"):
            raise KeyError(
                "'I_tor' is provided but 'picard_current' is set to 'off' or not provided. Please provide a valid 'picard_current', e.g. 'auto'."
            )
        if not self.has_Itor and (picard_current != "off"):
            raise KeyError(
                "Expected 'I_tor' in the parameters since 'picard_current' is not 'off'."
                + " Please set 'picard_current' to 'off' if you want to use a fixed 'iota' profile or provide 'I_tor'."
            )

        # Automatically generate the stages for the current constraint
        if picard_current == "auto":
            self.logger.info("Using `picard_current` automatic mode. Generating stages ...")
            if self.stages != [{}]:
                raise ValueError(
                    "Picard current is set to 'auto' but 'stages' is specified!"
                    + " Please remove `stages` from the parameter file to use `picard_current` in automatic mode."
                )

            minimize_tol = params.get("minimize_tol", 1e-6)
            self.stages = _auto_generate_stages(minimize_tol, 1e-10)
            paramters_stages_toml = params.copy()
            paramters_stages_toml["stages"] = self.stages
            self.logger.info(f"... generated {len(self.stages)} stages.")
            self.logger.info(
                f"... writing new parameters to 'parameter_{self.params['ProjectName']}.stages.toml'"
            )
            gvec.util.write_parameters(
                paramters_stages_toml,
                project_dir / f"parameter_{self.params['ProjectName']}.stages.toml",
            )
            self.params["picard_current"] = gvec.util.CaseInsensitiveDict()

        if self.has_Itor:
            match params["I_tor"].get("type", "polynomial"):
                case "polynomial":
                    coefs = np.array(params["I_tor"]["coefs"][::-1])
                    coefs *= params["I_tor"].get("scale", 1.0)
                    I_tor_target = np.poly1d(coefs)(rho**2)
                    try:
                        assert abs(np.poly1d(coefs)(0.0)) < 1e-9
                    except AssertionError:
                        self.logger.warning(
                            f"WARNING: Toroidal current profile not zero at magnetic axis!  I_tor(rho=0):{np.poly1d(coefs)(0.0)}"
                        )
                case "bspline":
                    from scipy.interpolate import BSpline

                    coefs = np.array(params["I_tor"]["coefs"], dtype=float)
                    coefs *= params["I_tor"].get("scale", 1.0)
                    knots = np.array(params["I_tor"]["knots"], dtype=float)
                    deg = np.sum(knots == knots[0]) - 1
                    I_tor_bspl = BSpline(knots, coefs, deg)
                    I_tor_target = I_tor_bspl(rho**2)
                    try:
                        assert abs(I_tor_bspl(0.0)) < 1e-9
                    except AssertionError:
                        self.logger.warning(
                            f"WARNING: Toroidal current profile not zero at magnetic axis! I_tor(rho=0):{I_tor_bspl(0.0)}"
                        )
                case "interpolation":
                    I_tor_target = np.array(params["I_tor"]["vals"], dtype=float)
                    rho = np.sqrt(np.array(params["I_tor"]["rho2"], dtype=float))
                    self.rho = rho
                    if min(rho) < 1e-6:
                        try:
                            assert min(abs(I_tor_target)) < 1e-9
                        except AssertionError:
                            self.logger.warning(
                                "WARNING: Toroidal current profile not zero at magnetic axis!"
                            )

                case _:
                    raise ValueError(f"Unknown Itor type: {params['Itor']['type']}")

            if not self.has_iota:
                self.params["iota"] = {"type": "polynomial", "coefs": [0.0]}
        else:
            self.I_tor_target = None

        if "I_tor" in params:
            self.I_tor_target = I_tor_target
            self.iota_rms = None
            self.curr_constraint = True
        else:
            self.curr_constraint = False

        # account for the change in relative path of restart, hmap, etc. files
        for key, value in self.params.items():
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ] and not value.startswith("/"):
                self.params[key] = f"../../{value}"

        # count the number of runs in each stage, for dynamic progressbar during current constraint
        self.n_runs_in_stage = [0 for _ in self.stages]

    def single_energy_minimization(self):
        """Run a single GVEC energy minimization using the current parameters. The run-state is updated after the run."""
        start_time = time.time()
        # find previous state
        if self.statefile:
            self.logger.debug(f"Restart from statefile {self.statefile}")
            self.params["init_LA"] = False

        # prepare the run directory
        self.rundir = self.project_dir / Path(f"{self.nth_stage:1d}-{self.nth_run:02d}")
        if self.rundir.exists():
            self.logger.debug(f"Removing existing run directory {self.rundir}")
            shutil.rmtree(self.rundir)
        self.rundir.mkdir()
        self.logger.debug(f"Created run directory {self.rundir}")

        # write parameterfile & run GVEC
        gvec.util.write_parameter_file_ini(
            gvec.util.flatten_parameters(self.params),
            self.rundir / "parameter.ini",
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage} run {self.nth_run})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )
        with gvec.util.chdir(self.rundir):
            gvec.run(
                "parameter.ini",
                "../../" / self.statefile if self.statefile else None,
                stdout_path="stdout.txt" if self.redirect_gvec_stdout else None,
            )

        # postprocessing
        self.statefile = sorted(self.rundir.glob("*State*.dat"))[-1]
        iterations = int(re.match(r".*State.*_(\d+)\.dat", self.statefile.name).group(1))
        iteration_offset = self.GVEC_iter_used
        self.GVEC_iter_used += iterations
        max_iterations = self.params.get("MaxIter")
        tolerance = self.params.get("minimize_tol")
        self.logger.debug(f"Postprocessing statefile {self.statefile}")

        with gvec.State(
            self.rundir / "parameter.ini",
            self.statefile,
            redirect_stdout=self.redirect_gvec_stdout,
        ) as state:
            ev = gvec.Evaluations(rho=self.rho, theta="int", zeta="int", state=state)
            state.compute(ev, "W_MHD", "N_FP")
            if self.curr_constraint:
                state.compute(ev, "iota", "iota_curr_0", "iota_0", "I_tor")

        # update iota
        if self.curr_constraint:
            iota_values = ev.iota_0 + self.I_tor_target * ev.iota_curr_0
            self.params["iota"] = {
                "type": "interpolation",
                "vals": iota_values.data,
                "rho2": (ev.rho**2).data,
            }

        # diagnostics
        if self.curr_constraint:
            iota_delta = ev.iota - iota_values
            self.rms_iota = np.sqrt((iota_delta**2).mean("rad"))
            self.logger.info(f"max Δiota: {np.abs(iota_delta).max().item():.2e}")
            self.logger.info(
                f"rms Δiota: {self.rms_iota.item():.2e}, iota_tol: {self.params['picard_current']['iota_tol']:.2e}"
            )
            self.logger.info(
                f"max ΔItor: {np.abs(ev.I_tor - self.I_tor_target).max().item():.2e}"
            )
        logfile = sorted(self.rundir.glob("logMinimizer_*"))[-1]
        log_df = read_csv(logfile, sep=",", header=0)

        self.logger.info(f"W_MHD: {ev.W_MHD.item():.2e}")

        d = xr.Dataset(
            dict(
                W_MHD=ev.W_MHD,
                gvec_iterations=iterations,
                gvec_max_iterations=max_iterations,
                gvec_tolerance=tolerance,
            )
        )
        d2 = xr.Dataset(
            dict(
                total_iteration=iteration_offset + np.array(log_df["#iterations"], dtype=int),
                force_X1=(["total_iteration"], np.array(log_df["normF_X1"])),
                force_X2=(["total_iteration"], np.array(log_df["normF_X2"])),
                force_LA=(["total_iteration"], np.array(log_df["normF_LA"])),
            )
        )
        if self.curr_constraint:
            d["iota"] = ev.iota
            d["I_tor"] = ev.I_tor
            d["iota_delta"] = iota_delta
            d["I_tor_delta"] = ev.I_tor - self.I_tor_target
        d = d.drop_vars(["pol_weight", "tor_weight"])
        if self.diagnostics is None:
            d = d.expand_dims(dict(run=[self.nth_run]))
            self.diagnostics = d
        else:
            d = d.expand_dims(dict(run=[self.diagnostics.run.size]))
            self.diagnostics = xr.concat([self.diagnostics, d], dim="run")
        if self.log_dia is None:
            self.log_dia = d2
        else:
            self.log_dia = xr.concat([self.log_dia, d2], dim="total_iteration")
        if self.diagnosticfile:
            xr.merge([self.diagnostics, self.log_dia]).to_netcdf(self.diagnosticfile)

        end_time = time.time()
        self.logger.info(
            f"GVEC run took {end_time - start_time:5.1f} seconds for {iterations} iterations. (max {max_iterations}, tol {tolerance:.1e})"
        )
        self.logger.info(
            f"GVEC iterations used in total: {self.GVEC_iter_used} / {self.totaliter}"
        )
        self.logger.info("-" * 40)

    def _reset_params_to_original(self):
        """
        Reset the parameters to the original values. Except for `iota` and `MaxIter`, which is limited by `totaliter`.
        """
        params = self.original_params.copy()
        params["iota"] = self.params["iota"]
        params["MaxIter"] = min(
            self.totaliter - self.GVEC_iter_used, self.original_params["MaxIter"]
        )
        if params["picard_current"] == "auto":
            params["picard_current"] = gvec.util.CaseInsensitiveDict()

        # account for the change in relative path of restart, hmap, etc. files
        for key, value in params.items():
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ] and not value.startswith("/"):
                params[key] = f"../../{value}"

        if self.has_Itor and (self.nth_run > 0 or self.nth_stage > 0):
            params["init_with_profile_iota"] = True

        self.params = params

    def _set_params_for_stage(self, stage: Mapping):
        """
        Set the run parameters to the values specified in the stage.

        Parameters
        ----------
        stage : Mapping
            Dictionary specifying which parameters are to be changed from the original parameter set.
        """
        for key, value in stage.items():
            if key == "MaxIter":
                self.params[key] = min(self.totaliter - self.GVEC_iter_used, value)

            if key == "picard_current" and value == "off":
                self.curr_constraint = False
            elif key == "picard_current" and value != "off" and not self.has_Itor:
                raise KeyError(
                    "Expected 'I_tor' in the parameters since 'picard_current' is not 'off'."
                    + " Please set 'picard_current' to 'off' if you want to use a fixed 'iota' profile or provide 'I_tor'."
                )
            elif key == "picard_current" and not isinstance(value, str):
                if key not in self.params:
                    self.params[key] = gvec.util.CaseInsensitiveDict()
                for subkey, subvalue in value.items():
                    self.params[key][subkey] = subvalue

            if key in ["iota", "pres", "sgrid"]:
                if key not in self.params:
                    self.params[key] = gvec.util.CaseInsensitiveDict()
                for subkey, subvalue in value.items():
                    self.params[key][subkey] = subvalue
            if key in self.params and isinstance(value, Mapping):
                for subkey, subvalue in value.items():
                    self.params[key][subkey] = subvalue
            else:
                self.params[key] = value

    def run_stages(self, return_output: bool = False):
        """Sequentially run the stages of the RunWithStages object.

        Parameters
        ----------
        return_output : bool, optional
            Flag to return the rundirectory, the final statefile and the diagnostics dataset, by default False.

        Returns
        -------
        tuple:
            Rundirectory (Path), the final statefile (Path) and the diagnostics dataset (xr.DataSet)

        Raises
        ------
        ValueError
            If stages are set when 'picard_current="auto"'
        KeyError
            If 'iota_tol' is not specified when 'I_tor' is provided.
        ValueError
            If 'picard_current.target' is not properly specified.
        """
        for s, stage in enumerate(self.stages):
            if self.GVEC_iter_used >= self.totaliter:
                self.logger.warning(
                    "WARNING: Maximum number of GVEC iterations reached. Aborting stages!"
                )
                break

            self.nth_stage = s
            self.nth_run = 0

            self._reset_params_to_original()
            self._set_params_for_stage(stage)

            # run the stage
            if self.curr_constraint:
                if self.params["picard_current"] == "auto":
                    raise ValueError(
                        'Detected `picard_current="auto"` during stage evaluation. Auto mode has to be set outside of the stages.'
                    )
                if "iota_tol" not in self.params["picard_current"]:
                    raise KeyError(f"During stage {s} 'iota_tol' is not specified.")

                iota_tol = self.params["picard_current"].get("iota_tol")
                target = self.params["picard_current"].get("target", "iota_and_force")
                match target:
                    case "iota":
                        self.n_runs_in_stage = self._current_constraint_target_iota(
                            self.totaliter, iota_tol, self.n_runs_in_stage
                        )
                    case "iota_and_force":
                        self.nth_run = -1
                        self.n_runs_in_stage = self._current_constraint_target_iota_and_force(
                            self.totaliter, iota_tol, self.n_runs_in_stage
                        )
                    case _:
                        raise ValueError(f"Unknown picard_current target:{target}")
            else:
                self.params["MaxIter"] = min(
                    self.totaliter - self.GVEC_iter_used,
                    self.params["MaxIter"],
                )
                self._eval_progressstr(self.n_runs_in_stage)
                self.single_energy_minimization()
                self.n_runs_in_stage[self.nth_stage] += 1

        self.logger.info("Done.")
        final_state = Path(self.params["ProjectName"] + "_State_final.dat")
        parameter_final = Path("parameter_" + self.params["ProjectName"] + "_final.ini")

        shutil.copy(self.statefile, final_state)
        parameters_final = gvec.util.read_parameter_file_ini(
            self.statefile.parents[0] / "parameter.ini"
        )
        parameters_final["MaxIter"] = -1
        for key in parameters_final:
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ]:
                parameters_final[key] = self.original_params[key]
        gvec.util.write_parameter_file_ini(
            parameters=parameters_final,
            path=parameter_final,
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage} run {self.nth_run})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )

        diagnostics = xr.merge([self.diagnostics, self.log_dia])

        if self.diagnosticfile:
            diagnostics.to_netcdf(self.diagnosticfile)
        if self.plots:
            self.plot_diagnostics(save_figs=True)

        if return_output:
            return (self.rundir, final_state, diagnostics)

    def _current_constraint_target_iota(
        self,
        totaliter: int,
        iota_tol: float,
        n_runs_in_stage: ArrayLike,
    ):
        """
        Target only iota in the current constraint, ignoring the forces.

        This method typically performs many picard iterations with few GVEC iterations to find a suitable iota profile for the current constraint.

        Parameters
        ----------
        totaliter : int
            The maximum total number of GVEC iterations allowed.
        iota_tol : float
            The tolerance for the iota convergence: abort criterion for the picard iterations.
        n_runs_in_stage : ArrayLike
            An ArrayLike object tracking the number of runs completed in each stage. Used for progressbar.

        Returns
        -------
        n_runs_in_stage: ArrayLike
            Updated number of runs completed in the current stage.
        """

        self.rms_iota = 1e6
        self.nth_run = -1
        if "maxRestarts" in self.params["picard_current"]:
            max_restarts = self.params["picard_current"]["maxRestarts"]
        else:
            max_restarts = 30
        self.logger.debug(f"maxRestarts:{max_restarts}")
        while (self.rms_iota > iota_tol) and (self.GVEC_iter_used < totaliter):
            self.logger.debug(f"nth run:{self.nth_run}")
            if self.nth_run + 1 > max_restarts:
                self.logger.warning(
                    "WARNING: Maximum number of restarts reached for this stage! Moving on to next stage."
                )
                break
            self.params["MaxIter"] = min(
                totaliter - self.GVEC_iter_used, self.params["MaxIter"]
            )
            self.nth_run += 1
            n_runs_in_stage[self.nth_stage] += 1
            self._eval_progressstr(n_runs_in_stage)
            self.single_energy_minimization()

        if self.rms_iota > iota_tol:
            self.logger.warning(
                f"WARNING: targeted iota has not been reached during stage {self.nth_stage}!\n"
                + f"target tol.: {iota_tol:.2e}, achieved tol.: {self.rms_iota.data:.2e}\n"
                + f"GVEC iterations used: {self.GVEC_iter_used}"
            )
        return n_runs_in_stage

    def _current_constraint_target_iota_and_force(
        self,
        totaliter: int,
        iota_tol: int,
        n_runs_in_stage: ArrayLike,
    ):
        """
        Run GVEC until the force tolerance is reached and perform picrad iterations until the iota tolerance is reached.
        The maximum number of total GVEC iterations and the maximum number of picard iterations are limited by totaliter.

        Parameters
        ----------
        totaliter : int
            The maximum total number of GVEC iterations allowed.
        iota_tol : int
            The tolerance for the iota convergence: abort criterion for the picard iterations.
        n_runs_in_stage : ArrayLike
            An ArrayLike object tracking the number of runs completed in each stage. Used for progressbar.

        Returns
        -------
        n_runs_in_stage: ArrayLike
            Updated number of runs completed in the current stage.
        """
        self.rms_iota = 1e6
        if "maxRestarts" in self.params["picard_current"]:
            max_restarts = self.params["picard_current"]["maxRestarts"]
        else:
            max_restarts = 30
        self.logger.debug(f"maxRestarts:{max_restarts}")
        while (self.GVEC_iter_used < totaliter) and (self.rms_iota > iota_tol):
            self.logger.debug(f"nth run:{self.nth_run}")
            if self.nth_run + 1 > max_restarts:
                self.logger.warning(
                    "WARNING: Maximum number of restarts reached for this stage! Moving on to next stage."
                )
                break
            self.params["MaxIter"] = min(
                totaliter - self.GVEC_iter_used, self.params["MaxIter"]
            )
            self.nth_run += 1
            n_runs_in_stage[self.nth_stage] += 1
            self._eval_progressstr(n_runs_in_stage)
            self.single_energy_minimization()
        if self.rms_iota > iota_tol:
            self.logger.warning(
                f"WARNING: targeted iota has not been reached during stage {self.nth_stage}!\n"
                + f"target tol.: {iota_tol:.2e}, achieved tol.: {self.rms_iota.data:.2e}\n"
                + f"GVEC iterations used: {self.GVEC_iter_used}"
            )
        return n_runs_in_stage

    def _eval_progressstr(self, n_runs_in_stage):
        """
        Evaluate and print a progress string for the current stage and run.

        Parameters
        ----------
        n_runs_in_stage : ArrayLike
            An ArrayLike object tracking the number of runs completed in each stage.

        Returns
        -------
        None
        """
        progressstr = "|"
        for i, ir in enumerate(n_runs_in_stage):
            if i < self.nth_stage:
                progressstr += "=" * ir + "|"
            elif i == self.nth_stage:
                progressstr += "=" * (ir - 1) + ">" + "|"
            else:
                progressstr += ".|"
        if self.progressbar:
            print(
                f"GVEC stage {self.nth_stage} run {self.nth_run}: {progressstr}",
                end="\r",
            )
        self.logger.info(f"GVEC stage {self.nth_stage} run {self.nth_run}: {progressstr}")

    def plot_diagnostics(self, save_figs: bool = False):
        """
        Plot diagnostic quantities from the GVEC runs.

        This method creates two figures: 'iterations.png' and 'profiles.png'.
        The first figure shows the evolution of the MHD energy and the root mean square difference
        to the target iota profile as a function of the run number.
        The second figure shows the profiles of the iota and I_tor values till the last run
        and the difference to the target profiles.

        Parameters
        ----------
        save_figs: bool
            Flag to save the figures

        Returns
        -------
        None
        """
        try:
            import matplotlib.pyplot as plt
        except ModuleNotFoundError:
            self.logger.warning(
                "WARNING: matplotlib not found! Diagnostic plots can not be generated."
            )
            return

        self.logger.debug("Plotting diagnostics...")
        diagnostics = self.diagnostics
        if self.curr_constraint:
            fig, axs = plt.subplots(1, 2, figsize=(10, 3), tight_layout=True)
        else:
            fig, ax = plt.subplots(1, 1, figsize=(5, 3), tight_layout=True)
            axs = [ax]
        axs[0].plot(diagnostics.run, diagnostics.W_MHD, ".-")
        axs[0].set(
            xlabel="run number",
            ylabel=f"${diagnostics.W_MHD.attrs['symbol']}$",
            title=diagnostics.W_MHD.attrs["long_name"],
        )
        if self.curr_constraint:
            axs[1].plot(diagnostics.run, np.sqrt((diagnostics.iota_delta**2).mean("rad")), ".-")
            axs[1].set(
                xlabel="run number",
                ylabel=r"$\sqrt{\sum \left(\Delta\iota\right)^2}$",
                title=f"Difference to target {diagnostics.iota.attrs['long_name']}\nroot mean square",
                yscale="log",
            )
        if save_figs:
            fig.savefig(f"{self.params['projectName']}_runs.png")
        else:
            fig.show()

        if self.curr_constraint:
            fig, axs = plt.subplots(2, 2, figsize=(15, 5), tight_layout=True, sharex=True)
            for r in diagnostics.run.data:
                if r == diagnostics.run.data[-1]:
                    kwargs = dict(marker=".", color="C0", alpha=1.0)
                else:
                    kwargs = dict(
                        color="black", alpha=0.2 + 0.3 * (r / diagnostics.run.data[-1])
                    )
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
            if save_figs:
                fig.savefig(f"{self.params['projectName']}_profiles.png")
            else:
                fig.show()

        if self.curr_constraint:
            fig_f, axs = plt.subplots(3, 1, figsize=(10, 5), tight_layout=True, sharex=True)
            axs[2].plot(
                np.cumsum(diagnostics.gvec_iterations),
                np.sqrt((diagnostics.iota_delta**2).mean("rad")),
                ".-",
            )
            axs[2].set(ylabel=r"$\Delta\iota_{rms}$", yscale="log")
        else:
            fig_f, axs = plt.subplots(2, 1, figsize=(10, 5), tight_layout=True, sharex=True)
        axs[0].plot(np.cumsum(diagnostics.gvec_iterations), diagnostics.W_MHD, ".-")
        axs[0].set(
            ylabel=f"${diagnostics.W_MHD.attrs['symbol']}$",
        )
        axf = axs[1]

        axf.plot(
            self.log_dia.total_iteration,
            self.log_dia.force_X1,
            color="C0",
            label=r"$X_1$",
        )
        axf.plot(
            self.log_dia.total_iteration,
            self.log_dia.force_X2,
            color="C1",
            label=r"$X_2$",
        )
        axf.plot(
            self.log_dia.total_iteration,
            self.log_dia.force_LA,
            color="C2",
            label=r"$\lambda$",
        )

        # stages vlines
        for ax in axs:
            n_runs_till_stage = np.cumsum(self.n_runs_in_stage)
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
        axf.vlines(
            np.cumsum(diagnostics.gvec_iterations),
            *axf.get_ylim(),
            colors="k",
            linestyle="dashed",
            alpha=0.6,
            zorder=-900,
            label="runs",
        )

        axf.set(ylabel="|Force|", yscale="log")
        axf.legend(bbox_to_anchor=(1.15, 1))
        # axf.legend(["stages","runs", r"$X_1$", r"$X_2$", r"$\lambda$"], bbox_to_anchor=(1.15, 1))
        axs[-1].set(xlabel="GVEC iteration")
        if save_figs:
            fig_f.savefig(f"{self.params['projectName']}_iterations.png")
        else:
            fig_f.show()


def _auto_generate_stages(minimize_target: float, iota_target: float):
    """Generate stages for 'picard_current' by ramping 'minimize_tol' and 'iota_target'.
    The first stage always targets 'iota', the other stages target 'iota_and_force'

    Parameters
    ----------
    minimize_target : float
        Final 'minimize_tol', i.e. the MHD force tolerance.
    iota_target : float
        Final 'iota_tol', i.e. the rms. tolerance on the targeted 'I_tor' profile.

    Returns
    -------
    stages: list
        List of dicts containing the changed parameters for each stage.
    """
    log_minimize_target = np.log10(minimize_target)
    log_iota_target = np.log10(iota_target)
    n_stages = max(int(max(-2 - log_minimize_target, log_minimize_target)), 1)
    minimize_tols = np.logspace(max(-3, log_minimize_target), log_minimize_target, n_stages)
    iota_tols = np.logspace(max(-3, log_iota_target), log_iota_target, n_stages)
    stages = [
        {
            "minimize_tol": minimize_tols[0],
            "MaxIter": 10,
            "picard_current": {"iota_tol": iota_tols[0], "target": "iota"},
        }
    ]
    for s, minimize_tol in enumerate(minimize_tols):
        iota_tol = iota_tols[s]
        stage = {
            "minimize_tol": minimize_tol,
            "picard_current": {"iota_tol": iota_tol, "target": "iota_and_force"},
        }
        stages.append(stage)
    return stages


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)
    if args.param_type is None:
        args.param_type = args.parameterfile.suffix[1:]

    if args.param_type == "ini":
        gvec.run(
            args.parameterfile,
            args.restartfile,
            stdout_path="stdout.txt" if args.quiet else None,
        )
    elif args.param_type in ["yaml", "toml"]:
        parameters = gvec.util.read_parameters(args.parameterfile)
        picard_mode = parameters.get("picard_current", "off")
        if "stages" not in parameters and picard_mode == "off":
            parameters = gvec.util.flatten_parameters(parameters)
            parameterfile = f"{args.parameterfile.name}.ini"
            gvec.util.write_parameter_file_ini(
                parameters,
                parameterfile,
                header=f"!Auto-generated from {args.parameterfile.name} with `pygvec run`\n!Created at {datetime.now().isoformat()}\n!pyGVEC v{gvec.__version__}\n",
            )
            gvec.run(
                parameterfile,
                args.restartfile,
                stdout_path="stdout.txt" if args.quiet else None,
            )
        else:
            logging.basicConfig(level=logging.WARNING)  # show warnings and above as normal
            logger = logging.getLogger(
                "pyGVEC.script"
            )  # show info/debug messages for this script
            logger.propagate = False
            loghandler = logging.StreamHandler()
            logformatter = logging.Formatter("{message}", style="{")
            loghandler.setFormatter(logformatter)
            logger.addHandler(loghandler)
            if args.verbose == 1:
                logger.setLevel(logging.INFO)
            elif args.verbose >= 2:
                logger.setLevel(logging.DEBUG)
            run_with_stages = RunWithStages(
                parameters,
                args.restartfile,
                logger=logger,
                progressbar=not args.quiet and not args.verbose,
                redirect_gvec_stdout=args.verbose < 3,
                diagnosticfile=args.diagnostics,
                plots=args.plots,
            )
            run_with_stages.run_stages()
    else:
        raise ValueError("Cannot determine parameterfile type")


if __name__ == "__main__":
    main()
