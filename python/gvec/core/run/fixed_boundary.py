# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
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
from .abstract import (
    AbstractRun,
    DEFAULT_MINIMIZE_TOL,
    AUTO_PICARD_CURRENT_IOTA_TOL_TARGET,
    DEFAULT_PICARD_CURRENT_MAXRESTARTS,
    DIAGNOSTICS_NPOINTS,
)


class FixedBoundaryRun(AbstractRun):
    def __init__(
        self,
        params: Mapping,
        state: Path | State | None = None,
        redirect_gvec_stdout: bool = True,
        quiet: bool = False,
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
        super().__init__(params, state, redirect_gvec_stdout, quiet, parameter_format)

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
            self.rundir / f"parameters_stage_{self.nth_stage:03d}_initial.ini",
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )

        with gvec.util.chdir(self.rundir):
            # write parameterfile & run GVEC
            self.state = State(
                parameterfile=f"parameters_stage_{self.nth_stage:03d}_initial.ini",
                statefile=self.state.statefile,
                init_for_minimization=True,
                stdout_path=Path("stdout.txt") if self.redirect_gvec_stdout else None,
            )
            iota_tol = self._state_parameters["picard_current"].get(
                "iota_tol", AUTO_PICARD_CURRENT_IOTA_TOL_TARGET
            )
            minimize_target = self._state_parameters.get("minimize_tol", DEFAULT_MINIMIZE_TOL)
            if self.save_intermediate:
                self.state.save(
                    f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                    f"!Created at {datetime.now().isoformat()}\n"
                    f"!pyGVEC v{gvec.__version__}\n",
                )
            self.state.minimize(
                I_tor=self.I_tor if self.curr_constraint else None,
                reset_LA=self._state_parameters["picard_current"].get("reset_LA", False),
                iota_tol=iota_tol,
                minimize_tol=minimize_target,
                maxIter=self._state_parameters["maxIter"],
                maxUpdates=self._state_parameters["picard_current"].get(
                    "maxRestarts", DEFAULT_PICARD_CURRENT_MAXRESTARTS
                ),
                save_intermediate=self.save_intermediate,
                picard_target=self._state_parameters["picard_current"].get("target", "auto"),
                save_header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n",
            )
            self.state.save(
                f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                f"!Created at {datetime.now().isoformat()}\n"
                f"!pyGVEC v{gvec.__version__}\n",
            )
        # postprocessing
        statefiles_found = sorted(self.rundir.glob("*State*.dat"))
        if len(statefiles_found) > 0:
            FileNotFoundError(f"Could not find any State files in rundir {self.rundir}")
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
        rho_eval = np.linspace(0, 1, DIAGNOSTICS_NPOINTS)

        ev = self.state.evaluate(*quantities, rho=rho_eval)
        ev = ev[quantities]

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
                diag_stage["iota_delta"] = ev.iota - self.state.evaluate_profile(
                    "iota", rho_eval
                )
                diag_stage["I_tor_delta"] = ev.I_tor - self.state.I_tor(rho_eval**2)
                diag_stage["rms_iota"] = self.state.rms_iota
                diag_stage["rms_iota_fit"] = self.state.rms_iota_fit_over_run[-1]
                diag_stage["max_Itor"] = self.state.delta_I_tor_max_over_run[-1]
                diag_stage["iota_dofs"] = (
                    self.state.iota.c.size if self.state.iota is not None else np.nan
                )
                diag_stage.reindex(iteration=self.state.iters_per_run, fill_value=np.nan)
                diag_stage["iota_delta_picard"] = xr.DataArray(
                    data=self.state.rms_iota_over_run,
                    dims=["picard_iteration"],
                    coords={
                        "picard_iteration": np.array(self.state.iters_per_run)
                        + iteration_offset
                    },
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
