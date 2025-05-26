# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC run script for running GVEC using stages and current constraints."""

import argparse
import copy
from datetime import datetime
import logging
import re
import shutil
import time
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import xarray as xr

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


def run_stages(
    parameters: Mapping,
    statefile: Path | None = None,
    progressbar: bool = False,
    redirect_gvec_stdout: bool = True,
    diagnosticfile: Path | None = None,
    plots: bool = False,
    init_LA: bool = False,
) -> tuple[Path, Path, xr.Dataset]:
    """Run GVEC with several stages (assuming hierarchical parameters)"""
    logger = logging.getLogger("pyGVEC.script")
    diagnostics: xr.Dataset | None = None
    rho = np.sqrt(np.linspace(0, 1, 101))
    rho[0] = 1e-4

    if "Itor" in parameters:
        match parameters["Itor"].get("type", "polynomial"):
            case "polynomial":
                coefs = np.array(parameters["Itor"]["coefs"][::-1])
                coefs *= parameters["Itor"].get("scale", 1.0)
                I_tor_target = np.poly1d(coefs)(rho**2)
            case "bspline":
                from scipy.interpolate import BSpline

                coefs = np.array(parameters["Itor"]["coefs"], dtype=float)
                coefs *= parameters["Itor"].get("scale", 1.0)
                knots = np.array(parameters["Itor"]["knots"], dtype=float)
                I_tor_target = BSpline(knots, coefs)(rho**2)
            case "interpolation":
                I_tor_target = np.array(parameters["Itor"]["vals"], dtype=float)
                rho = np.sqrt(np.array(parameters["Itor"]["rho2"], dtype=float))
            case _:
                raise ValueError(f"Unknown Itor type: {parameters['Itor']['type']}")

    stages = parameters.get("stages", [{}])

    # prepare the run directory
    project_dir = Path(f"{parameters['ProjectName']}_gvec_stages")
    if project_dir.exists():
        logger.debug(f"Removing existing run directory {project_dir}")
        shutil.rmtree(project_dir)
    project_dir.mkdir()

    run_params = gvec.util.CaseInsensitiveDict(copy.deepcopy(parameters))
    # add additional directory for path parameters
    for key, value in run_params.items():
        if key.lower() in [
            "vmecwoutfile",
            "boundary_filename",
            "hmap_ncfile",
        ] and not value.startswith("/"):
            run_params[key] = f"../{value}"
    for s, stage in enumerate(stages):
        # adapt parameters for this stage
        for key in ["stages", "Itor"]:
            if key in run_params:
                del run_params[key]
        for key, value in stage.items():
            if key in ["runs"]:
                continue
            if key in ["iota", "pres", "sgrid"]:
                if key not in run_params:
                    run_params[key] = {}
                for subkey, subvalue in value.items():
                    run_params[key][subkey] = subvalue
            if key in run_params and isinstance(value, Mapping):
                for subkey, subvalue in value.items():
                    run_params[key][subkey] = subvalue
            else:
                run_params[key] = value

        # run the stage
        runs = range(stage.get("runs", 1))
        for r in runs:
            progressstr = (
                "".join("|" + "=" * st.get("runs", 1) for st in stages[:s])
                + "|"
                + "=" * r
                + ">"
                + "." * (stage.get("runs", 1) - r - 1)
                + "|"
                + "".join("." * st.get("runs", 1) + "|" for st in stages[s + 1 :])
            )
            if progressbar:
                print(f"GVEC stage {s} run {r}: {progressstr}", end="\r")
            logger.info(f"GVEC stage {s} run {r}: {progressstr}")
            start_time = time.time()
            # find previous state
            if statefile:
                logger.debug(f"Restart from statefile {statefile}")
                run_params["init_LA"] = init_LA

            # prepare the run directory
            rundir = project_dir / Path(f"{s:1d}-{r:02d}")
            if rundir.exists():
                logger.debug(f"Removing existing run directory {rundir}")
                shutil.rmtree(rundir)
            rundir.mkdir()
            logger.debug(f"Created run directory {rundir}")

            # write parameterfile & run GVEC
            gvec.util.write_parameter_file(
                gvec.util.flatten_parameters(run_params),
                rundir / "parameter.ini",
                header=f"!Auto-generated with `pygvec run` (stage {s} run {r})\n"
                f"!Created at {datetime.now().isoformat()}\n"
                f"!pyGVEC v{gvec.__version__}\n",
            )
            with gvec.util.chdir(rundir):
                gvec.run(
                    "parameter.ini",
                    "../../" / statefile if statefile else None,
                    stdout_path="stdout.txt" if redirect_gvec_stdout else None,
                )

            # postprocessing
            statefile = sorted(rundir.glob("*State*.dat"))[-1]
            iterations = int(re.match(r".*State.*_(\d+)\.dat", statefile.name).group(1))
            max_iterations = run_params.get("maxiter")
            tolerance = run_params.get("minimize_tol")
            logger.debug(f"Postprocessing statefile {statefile}")

            with gvec.State(
                rundir / "parameter.ini",
                statefile,
                redirect_stdout=redirect_gvec_stdout,
            ) as state:
                ev = gvec.Evaluations(rho=rho, theta="int", zeta="int", state=state)
                state.compute(ev, "W_MHD", "N_FP")
                if "Itor" in parameters:
                    state.compute(ev, "iota", "iota_curr_0", "iota_0", "I_tor")

            if "Itor" in parameters:
                iota_values = ev.iota_0 + I_tor_target * ev.iota_curr_0
                run_params["iota"] = {
                    "type": "interpolation",
                    "vals": iota_values.data,
                    "rho2": (ev.rho**2).data,
                }

            # diagnostics
            # ToDo: possible early stop condition

            if "Itor" in parameters:
                iota_delta = ev.iota - iota_values
                rms_iota = np.sqrt((iota_delta**2).mean("rad"))
                logger.info(f"max Δiota: {np.abs(iota_delta).max().item():.2e}")
                logger.info(f"rms Δiota: {rms_iota.item():.2e}")
                logger.info(
                    f"max ΔItor: {np.abs(ev.I_tor - I_tor_target).max().item():.2e}"
                )

            logger.info(f"W_MHD: {ev.W_MHD.item():.2e}")

            d = xr.Dataset(
                dict(
                    W_MHD=ev.W_MHD,
                    gvec_iterations=iterations,
                    gvec_max_iterations=max_iterations,
                    gvec_tolerance=tolerance,
                )
            )
            if "Itor" in parameters:
                d["iota"] = ev.iota
                d["I_tor"] = ev.I_tor
                d["iota_delta"] = iota_delta
                d["I_tor_delta"] = ev.I_tor - I_tor_target
            d = d.drop_vars(["pol_weight", "tor_weight"])
            if diagnostics is None:
                d = d.expand_dims(dict(run=[r]))
                diagnostics = d
            else:
                d = d.expand_dims(dict(run=[diagnostics.run.size]))
                diagnostics = xr.concat([diagnostics, d], dim="run")
            if diagnosticfile:
                diagnostics.to_netcdf(diagnosticfile)

            end_time = time.time()
            logger.info(
                f"GVEC run took {end_time - start_time:5.1f} seconds for {iterations} iterations. (max {max_iterations}, tol {tolerance:.1e})"
            )
            logger.info("-" * 40)

            if "boundary_perturb" in run_params:
                run_params["boundary_perturb"] = False

    if plots:
        import matplotlib.pyplot as plt

        logger.debug("Plotting diagnostics...")

        if "Itor" in parameters:
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
        if "Itor" in parameters:
            axs[1].plot(
                diagnostics.run, np.sqrt((diagnostics.iota_delta**2).mean("rad")), ".-"
            )
            axs[1].set(
                xlabel="run number",
                ylabel=r"$\sqrt{\sum \left(\Delta\iota\right)^2}$",
                title=f"Difference to target {diagnostics.iota.attrs['long_name']}\nroot mean square",
                yscale="log",
            )
        fig.savefig("iterations.png")

        if "Itor" in parameters:
            fig, axs = plt.subplots(
                2, 2, figsize=(15, 5), tight_layout=True, sharex=True
            )
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
            fig.savefig("profiles.png")

    logger.info("Done.")
    final_state = Path(parameters["ProjectName"] + "_State_final.dat")
    parameter_final = Path("parameter_" + parameters["ProjectName"] + "_final.ini")

    shutil.copy(statefile, final_state)
    shutil.copy(statefile.parents[0] / "parameter.ini", parameter_final)
    return rundir, final_state, diagnostics


class StagesState:
    def __init__(
        self,
        params: Mapping,
        project_dir: Path,
        statefile: Path | None = None,
        logger: logging.Logger | None = None,
        stage: int = 0,
        nth_run: int = 0,
        I_tor_target: np.ndarray | None = None,
        rho: np.ndarray | None = None,
        diagnostics: xr.Dataset | None = None,
        GVEC_iter_used: int = 0,
        redirect_gvec_stdout: bool = True,
        diagnosticfile: Path | None = None,
        progressbar=True,
    ):
        self.statefile = statefile
        self.params = params
        self.project_dir = project_dir
        self.logger = logger
        self.stage = stage
        self.nth_run = nth_run
        self.rundir = None
        self.progressbar = progressbar
        if "Itor" in params:
            self.I_tor_target = I_tor_target
            self.iota_rms = None
            self.curr_constraint = True
        else:
            self.curr_constraint = False
        self.rho = rho
        self.diagnostics = diagnostics
        self.diagnosticfile = diagnosticfile
        self.GVEC_iter_used = GVEC_iter_used
        self.redirect_gvec_stdout = redirect_gvec_stdout

    def run(self):
        start_time = time.time()
        # find previous state
        if self.statefile:
            self.logger.debug(f"Restart from statefile {self.statefile}")
            self.params["init_LA"] = False

        # prepare the run directory
        self.rundir = self.project_dir / Path(f"{self.stage:1d}-{self.nth_run:02d}")
        if self.rundir.exists():
            self.logger.debug(f"Removing existing run directory {self.rundir}")
            shutil.rmtree(self.rundir)
        self.rundir.mkdir()
        self.logger.debug(f"Created run directory {self.rundir}")

        # write parameterfile & run GVEC
        gvec.util.write_parameter_file(
            gvec.util.flatten_parameters(self.params),
            self.rundir / "parameter.ini",
            header=f"!Auto-generated with `pygvec run` (stage {self.stage} run {self.nth_run})\n"
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
        iterations = int(
            re.match(r".*State.*_(\d+)\.dat", self.statefile.name).group(1)
        )
        self.GVEC_iter_used += iterations
        max_iterations = self.params.get("maxiter")
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
                f"rms Δiota: {self.rms_iota.item():.2e}, iota_tol: {self.params['iota_tol']:.2e}"
            )
            self.logger.info(
                f"max ΔItor: {np.abs(ev.I_tor - self.I_tor_target).max().item():.2e}"
            )

        self.logger.info(f"W_MHD: {ev.W_MHD.item():.2e}")

        d = xr.Dataset(
            dict(
                W_MHD=ev.W_MHD,
                gvec_iterations=iterations,
                gvec_max_iterations=max_iterations,
                gvec_tolerance=tolerance,
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
        if self.diagnosticfile:
            self.diagnostics.to_netcdf(self.diagnosticfile)

        end_time = time.time()
        self.logger.info(
            f"GVEC run took {end_time - start_time:5.1f} seconds for {iterations} iterations. (max {max_iterations}, tol {tolerance:.1e})"
        )
        self.logger.info("-" * 40)

        if "boundary_perturb" in self.params:
            self.params["boundary_perturb"] = False

    def _initialize_current_constraint(
        self, stage, maxIter_total, iota_tol, runs, n_runs_in_stage
    ):
        self.rms_iota = 1e6
        self.nth_run = -1
        while (
            (self.nth_run + 1 < runs)
            and (self.rms_iota > iota_tol)
            and (self.GVEC_iter_used < maxIter_total)
        ):
            self.nth_run += 1
            n_runs_in_stage[self.stage] += 1
            if self.progressbar:
                self._eval_progressstr(n_runs_in_stage)
            self.run()
            return n_runs_in_stage

    def _run_current_constraint(
        self, stage, maxIter_total, iota_tol, maxiter_per_run, runs, n_runs_in_stage
    ):
        self.rms_iota = 1e6
        while (
            (self.GVEC_iter_used < maxIter_total)
            and (self.rms_iota > iota_tol)
            and (self.nth_run + 1 < runs)
        ):
            self.params["maxIter"] = min(
                maxIter_total - self.GVEC_iter_used, maxiter_per_run
            )
            self.nth_run += 1
            n_runs_in_stage[self.stage] += 1
            if self.progressbar:
                self._eval_progressstr(n_runs_in_stage)
            self.run()
        return n_runs_in_stage

    def _eval_progressstr(self, n_runs_in_stage):
        progressstr = "|"
        for i, ir in enumerate(n_runs_in_stage):
            if i < self.stage:
                progressstr += "=" * ir + "|"
            elif i == self.stage:
                progressstr += "=" * (ir - 1) + ">" + "|"
            else:
                progressstr += ".|"
        print(f"GVEC stage {self.stage} run {self.nth_run}: {progressstr}", end="\r")
        self.logger.info(f"GVEC stage {self.stage} run {self.nth_run}: {progressstr}")

    def plot_diagnostics(self):
        import matplotlib.pyplot as plt

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
            axs[1].plot(
                diagnostics.run, np.sqrt((diagnostics.iota_delta**2).mean("rad")), ".-"
            )
            axs[1].set(
                xlabel="run number",
                ylabel=r"$\sqrt{\sum \left(\Delta\iota\right)^2}$",
                title=f"Difference to target {diagnostics.iota.attrs['long_name']}\nroot mean square",
                yscale="log",
            )
        fig.savefig("iterations.png")

        if self.curr_constraint:
            fig, axs = plt.subplots(
                2, 2, figsize=(15, 5), tight_layout=True, sharex=True
            )
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
            fig.savefig("profiles.png")


def run_with_class(
    parameters: Mapping,
    statefile: Path | None = None,
    redirect_gvec_stdout: bool = True,
    diagnosticfile: Path | None = None,
    progressbar: bool = True,
    plots: bool = False,
) -> tuple[Path, Path, xr.Dataset]:
    logger = logging.getLogger("pyGVEC.script")
    rho = np.sqrt(np.linspace(0, 1, 101))
    rho[0] = 1e-4
    I_tor_target = None
    if "Itor" in parameters:
        match parameters["Itor"].get("type", "polynomial"):
            case "polynomial":
                coefs = np.array(parameters["Itor"]["coefs"][::-1])
                coefs *= parameters["Itor"].get("scale", 1.0)
                I_tor_target = np.poly1d(coefs)(rho**2)
            case "bspline":
                from scipy.interpolate import BSpline

                coefs = np.array(parameters["Itor"]["coefs"], dtype=float)
                coefs *= parameters["Itor"].get("scale", 1.0)
                knots = np.array(parameters["Itor"]["knots"], dtype=float)
                I_tor_target = BSpline(knots, coefs)(rho**2)
            case "interpolation":
                I_tor_target = np.array(parameters["Itor"]["vals"], dtype=float)
                rho = np.sqrt(np.array(parameters["Itor"]["rho2"], dtype=float))
            case _:
                raise ValueError(f"Unknown Itor type: {parameters['Itor']['type']}")

    stages = parameters.get("stages", [{}])
    # prepare the run directory
    project_dir = Path(f"{parameters['ProjectName']}_gvec_stages")
    if project_dir.exists():
        logger.debug(f"Removing existing run directory {project_dir}")
        shutil.rmtree(project_dir)
    project_dir.mkdir()
    run_params = gvec.util.CaseInsensitiveDict(copy.deepcopy(parameters))
    maxIter_total = parameters.get("maxIter", int(1e5))
    state_of_stage = StagesState(
        params=run_params,
        project_dir=project_dir,
        statefile=statefile,
        I_tor_target=I_tor_target,
        rho=rho,
        logger=logger,
        diagnosticfile=diagnosticfile,
        redirect_gvec_stdout=redirect_gvec_stdout,
        progressbar=progressbar,
    )
    for key, value in state_of_stage.params.items():
        if key.lower() in [
            "vmecwoutfile",
            "boundary_filename",
            "hmap_ncfile",
        ] and not value.startswith("/"):
            state_of_stage.params[key] = f"../{value}"
    n_runs_in_stage = [0 for _ in stages]
    progressstr = "".join(["|" + "." * max(1, ir) for ir in n_runs_in_stage]) + "|"
    for s, stage in enumerate(stages):
        # adapt parameters for this stage
        for key in ["stages", "Itor"]:
            if key in state_of_stage.params:
                del state_of_stage.params[key]
        for key, value in stage.items():
            if key in ["runs"]:
                continue
            if key in ["iota", "pres", "sgrid"]:
                if key not in state_of_stage.params:
                    state_of_stage.params[key] = {}
                for subkey, subvalue in value.items():
                    state_of_stage.params[key][subkey] = subvalue
            if key in state_of_stage.params and isinstance(value, Mapping):
                for subkey, subvalue in value.items():
                    state_of_stage.params[key][subkey] = subvalue
            else:
                state_of_stage.params[key] = value
        maxIter_total = stage.get("maxIter", maxIter_total)
        state_of_stage.stage = s
        # run the stage
        state_of_stage.nth_run = 0
        if "minimize_tol" not in state_of_stage.params:
            print("... no minimize_tol found in parameters, setting to 1e-7! \r")
            state_of_stage.params["minimize_tol"] = 1e-7

        if state_of_stage.curr_constraint:
            iota_tol = stage.get("iota_tol")
            if iota_tol is None:
                raise ValueError(
                    "iota_tol must be set in stage for toroidal current constraint!"
                )
            init_iota = stage.get("init_iota", False)
            runs = stage.get("runs", 100)
            if init_iota:
                maxIter_per_run = stage.get("maxiter_per_run", 10)
                state_of_stage.params["maxIter"] = maxIter_per_run
                n_runs_in_stage = state_of_stage._initialize_current_constraint(
                    stage, maxIter_total, iota_tol, runs, n_runs_in_stage
                )
            else:
                maxIter_per_run = stage.get("maxiter_per_run", int(1e6))
                state_of_stage.nth_run = -1
            n_runs_in_stage = state_of_stage._run_current_constraint(
                stage, maxIter_total, iota_tol, maxIter_per_run, runs, n_runs_in_stage
            )
        else:
            runs = stage.get("runs", 1)
            for r in range(runs):
                state_of_stage.nth_run = r
                progressstr = (
                    "".join("|" + "=" * st.get("runs", 1) for st in stages[:s])
                    + "|"
                    + "=" * r
                    + ">"
                    + "." * (stage.get("runs", 1) - r - 1)
                    + "|"
                    + "".join("." * st.get("runs", 1) + "|" for st in stages[s + 1 :])
                )
                if progressbar:
                    print(f"GVEC stage {s} run {r}: {progressstr}", end="\r")
                state_of_stage.logger.info(f"GVEC stage {s} run {r}: {progressstr}")
                maxIter_per_run = stage.get("maxiter_per_run", int(1e6))
                state_of_stage.params["maxIter"] = min(
                    maxIter_total - state_of_stage.GVEC_iter_used, maxIter_per_run
                )
                state_of_stage.run()

    state_of_stage.logger.info("Done.")
    final_state = Path(parameters["ProjectName"] + "_State_final.dat")
    parameter_final = Path("parameter_" + parameters["ProjectName"] + "_final.ini")

    shutil.copy(state_of_stage.statefile, final_state)
    shutil.copy(state_of_stage.statefile.parents[0] / "parameter.ini", parameter_final)
    if diagnosticfile:
        state_of_stage.diagnostics.to_netcdf(diagnosticfile)
    if plots:
        state_of_stage.plot_diagnostics()
    return state_of_stage.rundir, final_state, state_of_stage.diagnostics


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
        if "stages" not in parameters:
            parameters = gvec.util.flatten_parameters(parameters)
            parameterfile = f"{args.parameterfile.name}.ini"
            gvec.util.write_parameter_file(
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
            logging.basicConfig(
                level=logging.WARNING
            )  # show warnings and above as normal
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
            # run_stages(
            #     parameters,
            #     args.restartfile,
            #     progressbar=not args.quiet and not args.verbose,
            #     redirect_gvec_stdout=args.verbose < 3,
            #     diagnosticfile=args.diagnostics,
            #     plots=args.plots,
            # )
            run_with_class(
                parameters,
                args.restartfile,
                progressbar=not args.quiet and not args.verbose,
                redirect_gvec_stdout=args.verbose < 3,
                diagnosticfile=args.diagnostics,
                plots=args.plots,
            )
    else:
        raise ValueError("Cannot determine parameterfile type")


if __name__ == "__main__":
    main()
