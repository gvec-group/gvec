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
from scipy.interpolate import BSpline

import gvec
from gvec.core.state import State
from gvec.core.minimize import MinimizationState
from .abstract import (
    AbstractRun,
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

    def _run_stage(self):
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
        self.state = MinimizationState.new(
            self._state_parameters,
            rundir=self.rundir,
            statefile=self.state.statefile,
            stdout_path=self.rundir / "stdout.txt",
            redirect_stdout=self.redirect_gvec_stdout,
            header=[
                "GVEC parameter file",
                f"Auto-generated with `pygvec run` (stage {self.nth_stage})",
                f"Created at {datetime.now().isoformat()}",
                f"pyGVEC v{gvec.__version__}",
            ],
        )

        with gvec.util.chdir(self.rundir):
            # write parameterfile & run GVEC
            if self.save_intermediate:
                self.state.save(
                    f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                    f"!Created at {datetime.now().isoformat()}\n"
                    f"!pyGVEC v{gvec.__version__}\n",
                )
            self.state.solve(
                save_intermediate=self.save_intermediate,
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
        if self.state.I_tor:
            quantities += ["iota", "I_tor"]

            self.logger.info(
                f"rms Δiota: {self.state.rms_iota:.2e}"
                f", iota_tol: {self.state.picard_current.iota_tol:.2e}"
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

            if self.state.I_tor:
                diag_stage["iota"] = ev.iota
                diag_stage["I_tor"] = ev.I_tor
                diag_stage["n_updates"] = self.state.nth_update
                diag_stage["iota_delta"] = ev.iota - self.state.evaluate_profile(
                    "iota", rho_eval
                )
                diag_stage["I_tor_delta"] = ev.I_tor - self.state.I_tor(rho_eval**2)

                diag_stage["rms_iota"] = self.state.rms_iota
                diag_stage["rms_iota_evolution"] = (
                    ("iota_update"),
                    self.state.rms_iota_over_run,
                )

                diag_stage["rms_iota_fit"] = self.state.rms_iota_fit_over_run[-1]
                diag_stage["rms_iota_fit_evolution"] = (
                    ("iota_update"),
                    self.state.rms_iota_fit_over_run,
                )

                diag_stage["max_Itor"] = self.state.delta_I_tor_max_over_run[-1]
                diag_stage["max_Itor_evolution"] = (
                    ("iota_update"),
                    self.state.delta_I_tor_max_over_run,
                )

                diag_stage["iota_dofs"] = (
                    self.state.iota.c.size if isinstance(self.state.iota, BSpline) else np.nan
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

    def plot_diagnostics_current_profiles(self):
        import matplotlib.pyplot as plt
        from gvec.plotting.plots1d import _add_rationals_to_iota_plot

        diagnostics = self.diagnostics_stage

        fig, axs = plt.subplots(2, 2, sharex="col", layout="constrained", figsize=(12, 8))
        updates_offset = 0
        for i, stage in enumerate(diagnostics):
            ds = diagnostics[stage]

            # delta iota
            axs[0, 0].semilogy(
                ds.iota_update + updates_offset,
                ds.rms_iota_fit_evolution,
                "x--",
                color=f"C{i}",
                label=r"$\Delta \iota_\text{fit}$ " + stage.replace("_", " "),
            )
            axs[0, 0].semilogy(
                ds.iota_update + updates_offset,
                ds.rms_iota_evolution,
                "d--",
                color=f"C{i}",
                label=r"$\Delta \iota_\text{rms}$ " + stage.replace("_", " "),
            )

            # delta Itor
            axs[1, 0].semilogy(
                ds.iota_update + updates_offset,
                ds.max_Itor_evolution,
                "x--",
                label=stage.replace("_", " "),
            )
            updates_offset += ds.n_updates

            # iota profile
            axs[0, 1].plot(ds.iota.rho, ds.iota, color=f"C{i}", label=stage.replace("_", " "))
            _add_rationals_to_iota_plot(axs[0, 1])
            axs[0, 1].set(ylabel=r"$\iota$")

            # Itor profile
            axs[1, 1].plot(ds.rho, ds.I_tor, color=f"C{i}", label=stage.replace("_", " "))
            axs[1, 1].set(ylabel=r"$I_\text{tor}$ [A]", xlabel=r"$\rho$")

        for ax in axs.flat:
            ax.grid("on")
            ax.legend(
                bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
                loc="lower left",
                ncols=2,
                mode="expand",
                borderaxespad=0.0,
            )

        axs[1, 0].set(
            ylabel=r"$\max\left(\Delta I_\text{tor}\right)$ [A]", xlabel=r"$\iota$ update"
        )
        return fig

    def plot_diagnostics_minimization(self):
        import matplotlib.pyplot as plt

        diagnostics = self.diagnostics_stage
        if self.state.I_tor:
            fig, axs = plt.subplots(3, 1, figsize=(10, 5), layout="constrained", sharex=True)
        else:
            fig, axs = plt.subplots(2, 1, figsize=(10, 5), tight_layout=True, sharex=True)

        # Energy
        axs[0].plot(
            self.diagnostics_minimizer.total_iteration, self.diagnostics_minimizer.W_MHD, "-"
        )
        axs[0].set(
            ylabel=f"${self.diagnostics_minimizer.W_MHD.attrs['symbol']}$",
        )

        # Forces
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

        if self.state.I_tor:
            # rms iota
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

            # vlines
            # after all the plots to get proper ylims
            for ax in axs:
                ylims = ax.get_ylim()
                for stage in diagnostics:
                    ds_aux = diagnostics[stage]
                    ax.vlines(
                        ds_aux.picard_iteration,
                        ylims[0],
                        ylims[1],
                        color="grey",
                        alpha=0.4,
                        zorder=-100,
                    )

        return fig
