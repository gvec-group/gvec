# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
import numpy as np
import xarray as xr
from collections.abc import Mapping
from scipy.constants import mu_0
from pandas import read_csv
from tqdm import tqdm
import time
import logging
import shutil
from datetime import datetime
from typing import Literal
from pathlib import Path

# GVEC
import gvec
from gvec.coils import CoilSet, Coil
from gvec.core.state import State
from gvec.core.minimize import MinimizationState
from .abstract import (
    AbstractRun,
    DEFAULT_MINIMIZE_TOL,
    DEFAULT_TOTALITER,
    DIAGNOSTICS_NPOINTS,
)
from gvec.util import write_parameter_file_ini, flatten_parameters, chdir
from gvec.core.compute import ev2ft
from gvec import fourier
from gvec.lib import _state


DEFAULT_FORCE_RATIO = 0.01
DEFAULT_FORCE_TOL_START = 1e-3
DEFAULT_ALPHA = 10.0
DEFAULT_BLEND_DEPTH = 0.6
DEFAULT_GTOL = 1e-8
DEFAULT_XTOL = 1e-5
DEFAULT_DIGITS_BIEST = 8
DEFAULT_FB_TOLERANCE = 1e-8
FORCE_RATIO_MINIMUM_FORCE = 1e-12
DEFAULT_FTOL = 1e-5
MOVING_AVERAGE_WINDOW = 10

logger = logging.getLogger("gvec.run.free_boundary")


class FreeBoundaryRun(AbstractRun):
    def __init__(
        self,
        params: Mapping,
        coil_set: CoilSet | Coil,
        state: Path | State | None = None,
        redirect_gvec_stdout: bool = True,
        quiet: bool = False,
        parameter_format: Literal["toml", "yaml"] = "toml",
        save_fb_updates: bool = False,
    ):
        params["FreeBoundary"] = True
        super().__init__(params, state, redirect_gvec_stdout, quiet, parameter_format)
        self.coil_set = coil_set
        self.save_fb_updates = save_fb_updates
        self.iters_used = 0
        self._resolved_force_ratio = 1.0

        self.m_max, self.n_max = self.state.get_mn_max()

        self.target_mavg_last = np.nan
        self.target_mavg = np.nan

        self.xnorm = np.nan
        self.gnorm = np.nan
        self.fdelta = np.nan

        self.xnorm_mavg = np.nan
        self.gnorm_mavg = np.nan

        if self._has_initial_stage:
            self.stages[0]["free_boundary"] = {"updates": 0}

    def _allocate_diagnostics(self, track_boundary_fields: bool):
        """Setup the diagnostics dataset."""
        diagnostics = self.state.evaluate("xyz", "theta", "zeta", rho=1.0).isel(rad=0)
        fb_updates = max(self._state_parameters["free_boundary"]["updates"], 1)
        diagnostics["updates"] = ("updates", np.arange(fb_updates))
        diagnostics["force_balance"] = (
            ("tor", "pol", "updates"),
            np.full(shape=(self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan),
        )
        diagnostics.force_balance.attrs = dict(
            long_name="free-boundary force balance on the LCFS.", symbol=r"$F_{b}}$"
        )
        if track_boundary_fields:
            for key, long_name, symbol in [
                ("B_jump", "jump of the magnetic field on the LCFS", r"$B_{\text{jump}}$"),
                ("B_int", "internal magnetic field on the LCFS", r"$B_{\text{int}}$"),
                ("B_ext", "external magnetic field on the LCFS", r"$B_{\text{ext}}$"),
                ("K", "surface current on the LCFS", r"$K$"),
            ]:
                diagnostics[key] = (
                    ("xyz", "tor", "pol", "updates"),
                    np.full(
                        shape=(3, self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan
                    ),
                )
                diagnostics[key].attrs = dict(long_name=long_name, symbol=symbol)

        for key in [
            "free_boundary_target",
            "free_boundary_target_mavg",
            "fdelta",
            "xnorm",
            "xnorm_mavg",
            "gnorm",
            "gnorm_mavg",
            "Fnorm_X1_b_IP",
            "Fnorm_X2_b_IP",
            "Fnorm_X1_b_modes",
            "Fnorm_X2_b_modes",
            "Fnorm_X1_b",
            "Fnorm_X2_b",
            "spectral_width",
            "dt",
            "resolved_force_ratio",
        ]:
            diagnostics[key] = (
                "updates",
                np.full(shape=fb_updates, fill_value=np.nan),
            )

        m, n, F = self.state._get_boundary_forces()
        diagnostics["m_modes"] = ("mode", m)
        diagnostics["n_modes"] = ("mode", n)
        diagnostics["boundary_forces"] = (
            ("modes", "updates"),
            np.full(shape=(len(m), fb_updates), fill_value=np.nan),
        )

        if self.state.I_tor:
            diagnostics["iota"] = (
                ("rho_iota", "updates"),
                np.full(shape=(DIAGNOSTICS_NPOINTS, fb_updates), fill_value=np.nan),
            )
            diagnostics["rho_iota"] = np.linspace(0, 1, DIAGNOSTICS_NPOINTS)

        if self.diagnostics_stage is None:
            self.diagnostics_stage = xr.DataTree()
        self.track_boundary_fields = track_boundary_fields
        self.diagnostics_stage[f"stage_{self.nth_stage}"] = diagnostics

    def _update_force_ratio(self, diag_aux, update):
        force_balance_ds = diag_aux.dataset.isel(updates=update)
        force_balance_ds["N_FP"] = self.state.nfp
        fourier_force_balance = ev2ft(force_balance_ds[["N_FP", "force_balance"]])
        cut_fb = fourier_force_balance.where(
            (fourier_force_balance.m <= self.m_max), drop=True
        ).where(np.abs(fourier_force_balance.n) <= self.n_max, drop=True)
        full_force = max(
            (fourier_force_balance.force_balance_mns**2).sum().item()
            + (fourier_force_balance.force_balance_mnc**2).sum().item(),
            FORCE_RATIO_MINIMUM_FORCE,
        )
        cut_force = (cut_fb.force_balance_mns**2).sum().item() + (
            cut_fb.force_balance_mnc**2
        ).sum().item()
        self._resolved_force_ratio = cut_force / full_force

    def _update_diagnostics(
        self,
        B_ext_2: np.ndarray,
        B_coils: np.ndarray,
        B_plas_ext: np.ndarray,
        update: int,
        dt: float,
    ):
        """Helper function to fill the diagnostics dataset

        Parameters
        ----------
        B_ext_2 : np.ndarray
            Squared external magnetic field to evaluate the free-boundary force balance
        B_coils : np.ndarray
            Coil field at the boundary of shape (3,tor,pol)
        B_plas_ext : np.ndarray
            Plasma contributionn to the external field at the boundary (BIEST) of shape (3,tor,pol)
        update : int
            Determines which entry in the allocated dataset is filled.
        """

        free_boundary_target, force_balance, ev_gvec = self.state._free_boundary_target(B_ext_2)
        B_ext = B_coils + B_plas_ext
        B_jump = ev_gvec.B - B_ext

        diag_aux = self.diagnostics_stage[f"stage_{self.nth_stage}"]
        diag_aux.force_balance[:, :, update] = force_balance
        diag_aux.free_boundary_target[update] = free_boundary_target
        Fnorm = self.state.get_Fnorm_boundary()
        if self.track_boundary_fields:
            # mu_0*K = n x [[B]]
            diag_aux.K[:, :, :, update] = np.cross(ev_gvec.normal, B_jump, axis=0) / mu_0
            diag_aux.B_ext[:, :, :, update] = B_ext
            diag_aux.B_jump[:, :, :, update] = B_jump.data
            diag_aux.B_int[:, :, :, update] = ev_gvec.B.data

        diag_aux.Fnorm_X1_b_IP[update] = Fnorm[0]
        diag_aux.Fnorm_X1_b_modes[update] = Fnorm[1]
        diag_aux.Fnorm_X2_b_IP[update] = Fnorm[2]
        diag_aux.Fnorm_X2_b_modes[update] = Fnorm[3]
        diag_aux.Fnorm_X1_b[update] = Fnorm[4]
        diag_aux.Fnorm_X2_b[update] = Fnorm[5]

        diag_aux.spectral_width[update] = _state.get_spectral_width()

        diag_aux.dt[update] = dt
        self.free_boundary_target = free_boundary_target.item()

        self._update_force_ratio(diag_aux, update)
        diag_aux.resolved_force_ratio[update] = self._resolved_force_ratio

        diag_aux.boundary_forces[:, update] = self.state._get_boundary_forces()[2]

        diag_aux.xnorm[update] = self.xnorm
        diag_aux.gnorm[update] = self.gnorm

        # moving averages
        diag_aux.free_boundary_target_mavg[update] = np.mean(
            diag_aux.free_boundary_target[
                update - min(MOVING_AVERAGE_WINDOW, update) : update + 1
            ]
        )
        diag_aux.fdelta[update] = (
            np.std(
                diag_aux.free_boundary_target_mavg[
                    update - min(MOVING_AVERAGE_WINDOW, update) : update + 1
                ]
            )
            if update != 0
            else np.nan
        )

        diag_aux.xnorm_mavg[update] = np.mean(
            diag_aux.xnorm[update - min(MOVING_AVERAGE_WINDOW, update) : update + 1].dropna(
                dim="updates"
            )
        )
        diag_aux.gnorm_mavg[update] = np.mean(
            diag_aux.gnorm[update - min(MOVING_AVERAGE_WINDOW, update) : update + 1].dropna(
                dim="updates"
            )
        )

        self.fdelta = diag_aux.fdelta[update].item()
        self.xnorm_mavg = diag_aux.xnorm_mavg[update].item()
        self.gnorm_mavg = diag_aux.gnorm_mavg[update].item()

        if np.isnan(self.target_mavg_last):
            self.target_mavg_last = diag_aux.free_boundary_target_mavg[update]
        else:
            self.target_mavg_last = diag_aux.free_boundary_target_mavg[update - 1]

        self.target_mavg = diag_aux.free_boundary_target_mavg[update]

        self.diagnostics_stage[f"stage_{self.nth_stage}"] = diag_aux

    def _run_stage(self):
        """Perform the free-boundary run.

        algorithm
        ---------
        1. start with a fixed-boundary minimization from the initial guess
        2. calculate the external field using BIEST and the coils
        3. update the boundary forces and push the boundary via boundary perturbations
        4. minimize the fixed-boundary equilibrium with the new boundary
        5. repeat from 2. until convergence  or until the budget is used up
        """

        # start with a fixed boundary equilibrium

        # how many free boundary updates,
        # NOTE: required parameter
        fb_updates = self._state_parameters["free_boundary"]["updates"]

        # target termination criterion
        free_boundary_tolerance = self._state_parameters["free_boundary"].get(
            "tolerance", DEFAULT_FB_TOLERANCE
        )

        # ftol termination criterion
        ftol = self._state_parameters["free_boundary"].get("ftol", DEFAULT_FTOL)

        minimizer_type = self._state_parameters["free_boundary"].get("minimizer_type", "GD")

        step_size_heuristic = False
        match minimizer_type:
            case "GD" | "mGD" | "lBFGS":
                step_size_heuristic = True
            case "BFGS" | "dBFGS" | "SR1" | "aGD" | "ADAM":
                step_size_heuristic = False
            case _:
                raise ValueError(f"Unknown free-boundary minimizer_type: {minimizer_type}")

        # If True, saves diagnostics after each free-b update
        save_diagnostics = self._state_parameters["free_boundary"].get(
            "save_diagnostics", False
        )

        # resolved force termination criterion
        force_ratio_limit = self._state_parameters["free_boundary"].get(
            "force_ratio_limit", DEFAULT_FORCE_RATIO
        )

        # If True, adaptively decrease the force tolerance according to the free-b target
        ramp_force_tol = self._state_parameters["free_boundary"].get("ramp_force_tol", True)

        track_boundary_fields = self._state_parameters["free_boundary"].get(
            "track_boundary_fields", False
        )

        # starting tolerance for the ramping
        ramp_force_tol_start = self._state_parameters["free_boundary"].get(
            "ramp_force_tol_start", DEFAULT_FORCE_TOL_START
        )

        # ramping factor for force tolerance: force_tolerance = free_boundary_target / alpha
        alpha = self._state_parameters["free_boundary"].get("alpha", DEFAULT_ALPHA)

        # how far the boundary is blend into an equilibrium after a free-boundary push
        blend_depth = self._state_parameters["free_boundary"].get(
            "blend_depth", DEFAULT_BLEND_DEPTH
        )

        # gradient termination criterion
        gtol = self._state_parameters["free_boundary"].get("gtol", DEFAULT_GTOL)

        # DOF change termination criterion
        # NOTE: fourier coefficients are weighted according to expected spectral decay
        xtol = self._state_parameters["free_boundary"].get("xtol", DEFAULT_XTOL)

        # how accurate BIEST solves for the exterior field
        digits_BIEST = self._state_parameters["free_boundary"].get(
            "digits_BIEST", DEFAULT_DIGITS_BIEST
        )

        # option to limit the number of updates per free-boundary update
        # TODO: Can this  be replaced by maxIter?
        max_iter_per_update = self._state_parameters["free_boundary"].get(
            "max_iter_per_update", DEFAULT_TOTALITER
        )

        # targeted minimize_tol
        minimize_tol_min = self._state_parameters.get("minimize_tol", DEFAULT_MINIMIZE_TOL)

        # minimum stepsize
        min_dt = self._state_parameters["free_boundary"].get("min_dt", 1e-4)

        # Flag to restart the free boundary minimizer a step-size reduction.
        reset_minimizer = False

        if ramp_force_tol:
            minimize_tol = ramp_force_tol_start
            current_lowest_min_tol = ramp_force_tol_start
        else:
            minimize_tol = minimize_tol_min

            # find previous state
        if self.state.statefile:
            logger.debug(f"Restart from statefile {self.state.statefile}")

        # prepare the run directory
        self.rundir = self.project_dir / Path(f"stage-{self.nth_stage:1d}")
        if self.rundir.exists():
            logger.debug(f"Removing existing run directory {self.rundir}")
            shutil.rmtree(self.rundir)
        self.rundir.mkdir()
        logger.debug(f"Created run directory {self.rundir}")

        self.state = MinimizationState.new(
            self._state_parameters,
            rundir=self.rundir,
            statefile=self.state.statefile,
            redirect_stdout=self.redirect_gvec_stdout,
            stdout_path=self.rundir / "stdout.txt",
            header=[
                "GVEC parameter file",
                f"Auto-generated with `pygvec run` (stage {self.nth_stage})",
                f"Created at {datetime.now().isoformat()}",
                f"pyGVEC v{gvec.__version__}",
            ],
        )
        with chdir(self.rundir):
            # write parameterfile & run GVEC
            self.m_max, self.n_max = self.state.get_mn_max()

            # initial stepsize for the free-b minimization
            dt = self.state.dt_free
            max_dt = self.state.max_dt_free

            if self.save_intermediate:
                self.state.save(
                    f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                    f"!Created at {datetime.now().isoformat()}\n"
                    f"!pyGVEC v{gvec.__version__}\n",
                )

            self.state.solve(
                minimize_tol=minimize_tol,
                maxIter=self._state_parameters["maxIter"],
                save_intermediate=self.save_intermediate,
                save_header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n",
            )
            self._allocate_diagnostics(track_boundary_fields)
            self.iters_used = self.state.get_iterations()

            # NOTE: decoupled in current constrain mode for future optimization
            if self.state.I_tor:
                self.I_tor_boundary = self.state.I_tor(1.0)
            else:
                self.I_tor_boundary = (
                    self.state.evaluate("I_tor", rho=1.0).squeeze()["I_tor"].item()
                )

            # Call to BIEST
            (B2, B_plas_ext, B_coils, B_ext_2) = _get_B2_fp(
                state=self.state,
                coil_set=self.coil_set,
                ntheta=self.state.t_n,
                nzeta=self.state.z_n,
                I_tor=self.I_tor_boundary,
                diagnostic_output=True,
                digits=digits_BIEST,
            )

            self._update_diagnostics(B_ext_2, B_coils, B_plas_ext, 0, dt)

            fb_update = 0
            try:  # allow clean abort with KeybordInterrupt
                for fb_update in range(1, fb_updates):
                    logger.info(
                        f"============= free-boundary update {fb_update} / {fb_updates} ============="
                    )
                    # abort if budget is used up
                    if self.iters_used >= self.totaliter:
                        logger.warning(
                            "Maximum number of iterations reached! aborting free-boundary iterations."
                        )
                        break

                    # ensure that the maximum number of iterations is kept
                    max_iter_per_update_aux = min(
                        self.totaliter - self.iters_used, max_iter_per_update
                    )

                    if ramp_force_tol:
                        # minimize_tol is ahead of free boundary target
                        minimize_tol = np.max(
                            [
                                np.min(
                                    [self.free_boundary_target / alpha, ramp_force_tol_start]
                                ),
                                minimize_tol_min,
                            ]
                        )
                        # this is to ensure that minimize tol never increases
                        if minimize_tol < current_lowest_min_tol:
                            current_lowest_min_tol = minimize_tol
                        minimize_tol = np.min([minimize_tol, current_lowest_min_tol])
                        logger.info(f"Setting gvec minimize_tol: {minimize_tol:.2e}")

                    fb_target_reached = self.free_boundary_target <= free_boundary_tolerance
                    ftol_reached = self.fdelta <= ftol
                    xtol_reached = self.xnorm_mavg <= xtol
                    gtol_reached = self.gnorm_mavg <= gtol
                    force_ratio_limit_reached = self._resolved_force_ratio <= force_ratio_limit

                    if fb_target_reached:
                        logger.info(
                            "Reached free boundary target tolerance! Terminating free boundary updates."
                        )
                        break
                    elif ftol_reached:
                        logger.info("Reached ftol! Terminating free boundary updates.")
                        break
                    elif xtol_reached:
                        logger.info("Reached xtol! Terminating free boundary updates.")
                        break
                    elif gtol_reached:
                        logger.info("Reached gtol! Terminating free boundary updates.")
                        break
                    elif force_ratio_limit_reached:
                        logger.info(
                            "Force ratio below specified limit! Terminating free boundary updates."
                        )
                        break
                    else:
                        logger.info(f"Δf: {self.fdelta:.2e}, ftol: {ftol:.2e}.")
                        logger.info(
                            f"xnorm: {self.xnorm:.2e}, m.avg. xnorm: {self.xnorm_mavg:.2e}, xtol: {xtol:.2e}."
                        )
                        logger.info(
                            f"gnorm: {self.gnorm:.2e}, m.avg. gnorm: {self.gnorm_mavg:.2e}, gtol: {gtol:.2e}."
                        )
                        logger.info(
                            f"force ratio {self._resolved_force_ratio:.2e}, limit {force_ratio_limit:.2e}"
                        )
                        t_fb_0 = time.time()

                        # Step-size dt is reduced if the moving average of the target increases until the minimum step-size is reached,
                        # otherwise it is increased. However, some algorithms such as BFGS also perform step-size estimations, which overwrite
                        # the provided step-size, if the estimated step-size is SMALLER.
                        if (
                            self.target_mavg_last < self.target_mavg or dt < min_dt
                        ) and step_size_heuristic:
                            dt = max(0.5 * dt, min_dt)
                        elif step_size_heuristic:
                            dt = min(1.1 * dt, max_dt)

                        # this is relevant to whipe the memory of methods like BFGS and aGD
                        if dt <= min_dt and reset_minimizer and step_size_heuristic:
                            reset_minimizer = True
                            logger.warning(
                                "Reached minimum step-size! Resetting free boundary minimizer."
                            )
                        elif step_size_heuristic:
                            reset_minimizer = False

                        logger.info(f"Free boundary dt: {dt:.2e}")

                        # pushes boundary and does a fixed-b run
                        self.state._free_boundary_update(
                            maxIter=min(
                                max_iter_per_update_aux, self._state_parameters["maxIter"]
                            ),
                            B2=B2,
                            dt=dt,
                            reset_minimizer=reset_minimizer,
                            blend_depth=blend_depth,
                            minimize_tol=minimize_tol,
                            save_header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                            + f"!free-boundary update {fb_update}\n",
                            save_intermediate=self.save_intermediate,
                        )
                        dt = self.state.dt_free

                        if self.save_fb_updates:
                            self.state.save(
                                f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                                f"!free-boundary update {fb_update}\n"
                                f"!Created at {datetime.now().isoformat()}\n"
                                f"!pyGVEC v{gvec.__version__}\n"
                            )

                        # abort criteria
                        self.xnorm = self.state.get_xnorm() / dt
                        self.gnorm = self.state.get_gnorm()
                        # NOTE: decoupled in current constrain mode for future optimization
                        if self.state.I_tor:
                            self.I_tor_boundary = self.state.I_tor(1.0)
                        else:
                            self.I_tor_boundary = (
                                self.state.evaluate("I_tor", rho=1.0).squeeze()["I_tor"].item()
                            )

                        # Call to BIEST
                        (B2, B_plas_ext, B_coils, B_ext_2) = _get_B2_fp(
                            state=self.state,
                            coil_set=self.coil_set,
                            ntheta=self.state.t_n,
                            nzeta=self.state.z_n,
                            I_tor=self.I_tor_boundary,
                            diagnostic_output=True,
                            digits=digits_BIEST,
                        )

                        self._update_diagnostics(B_ext_2, B_coils, B_plas_ext, fb_update, dt)

                        t_fb_1 = time.time()
                        logger.info(f"Time for free boundary update: {t_fb_1 - t_fb_0:.4}s")
                        reset_minimizer = False
                        self.iters_used = self.state.get_iterations()
                        if self.state.I_tor is None:
                            self.I_tor_boundary = self.state.evaluate(
                                "I_tor", rho=1.0
                            ).I_tor.item()
                        else:
                            self.diagnostics_stage[f"stage_{self.nth_stage}"].iota[
                                :, fb_update
                            ] = self.state.evaluate_profile(
                                "iota",
                                self.diagnostics_stage[f"stage_{self.nth_stage}"]["rho_iota"],
                            )

                    logger.info(
                        f"free-boundary - target: {self.free_boundary_target:.3e}, tolerance: {free_boundary_tolerance}, iters: {self.iters_used} / {self.totaliter} "
                    )
                    if save_diagnostics:
                        try:
                            self.diagnostics_stage[f"stage_{self.nth_stage}"].to_netcdf(
                                "diagnostics_aux.nc"
                            )
                        except PermissionError:
                            logger.warning(
                                "Permission error when attemting to write diagnostics to netcdf! Maybe the data-set is opened somewhere else?"
                            )

                # ensure that the final result is in fixed-boundary equilibrium with specified tolerance
                logger.info("=================================================================")
                logger.info("FINAL free-boundary optimization criteria")
                logger.info(
                    f"free-boundary - target: {self.free_boundary_target:.3e}, tolerance: {free_boundary_tolerance}, iters: {self.iters_used} / {self.totaliter} "
                )
                logger.info(f"Δf: {self.fdelta:.2e}, ftol: {ftol}.")
                logger.info(f"xnorm: {self.xnorm:.2e}, xtol: {xtol:.2e}.")
                logger.info(f"gnorm: {self.gnorm:.2e}, gtol: {gtol:.2e}.")
                logger.info(
                    f"resolved force-balance ratio: {self._resolved_force_ratio:.2e}, specified limit: {force_ratio_limit:.2e}."
                )
                logger.info("=================================================================")
                if self._state_parameters["maxIter"] > 0:
                    logger.info(
                        "Minimizing fixed-boundary equilibrium for the final boundary ..."
                    )
                    self.state.solve(
                        minimize_tol=minimize_tol_min,
                        maxIter=self._state_parameters["maxIter"],
                        save_header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n",
                    )

            except KeyboardInterrupt:
                logger.warning("Aborting free-boundary due to keyboard interrupt!")
            if (
                self._resolved_force_ratio < force_ratio_limit
                and self.free_boundary_target > free_boundary_tolerance
            ):
                fb_aux = (
                    self.diagnostics_stage[f"stage_{self.nth_stage}"]
                    .force_balance.isel(updates=fb_update - 1)
                    .squeeze()
                    .transpose("pol", "tor")
                )
                M_proposed, N_proposed = minimal_modes_force_balance(
                    fb_aux,
                    zeta0=fb_aux.zeta[0].item(),
                    theta0=fb_aux.theta[0].item(),
                    tolerance=free_boundary_tolerance,
                )
                logger.info(
                    "Not all modes in the free-boundary force  balance can be resolved with the provided angular resolution.\n"
                    + "Try reducing minimze_tol or adapt the angular resolution.\n"
                    + "Proposed angular resolution increase for the solution variables X1 / X2 / LA:\n"
                    + f"    m_max - current {self.m_max} -> proposed {M_proposed}\n"
                    + f"    n_max - current {self.n_max} -> proposed {N_proposed}"
                )

            if self.free_boundary_target > free_boundary_tolerance and self.nth_stage > 0:
                logger.warning(
                    f"Free-boundary target not reached! target: {self.free_boundary_target:.3e}, specified tol.: {free_boundary_tolerance:.3e}"
                )

            self.state.save(
                f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
                f"!free-boundary update {fb_update}\n"
                f"!Created at {datetime.now().isoformat()}\n"
                f"!pyGVEC v{gvec.__version__}\n"
            )

        iterations = self.state.get_iterations()
        iteration_offset = self.GVEC_iter_used
        self.GVEC_iter_used += iterations

        # postprocessing
        statefiles_found = sorted(self.rundir.glob("*State*.dat"))
        if len(statefiles_found) > 0:
            FileNotFoundError(f"Could not find any State files in rundir {self.rundir}")
        self.state.statefile = Path(statefiles_found[-1]).absolute()

        logfiles = sorted(self.rundir.glob("logMinimizer_*"))
        if logfiles:
            logfile = logfiles[-1]
            log_df = read_csv(logfile, sep=",", header=0)

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
        else:
            self.max_force = np.inf

    def plot_diagnostics_minimization(self):
        import matplotlib.pyplot as plt
        from gvec.plotting.plots1d import _add_rationals_to_iota_plot

        fig_fb, axs_fb = plt.subplots(5, 1, figsize=(10, 10), layout="constrained", sharex=True)
        quantities = [
            "free_boundary_target",
            "free_boundary_target_mavg",
            "Fnorm_X1_b",
            "Fnorm_X2_b",
            "resolved_force_ratio",
            "xnorm",
            "xnorm_mavg",
            "gnorm",
            "gnorm_mavg",
        ]
        if self.state.I_tor:
            fig_iota, ax_iota = plt.subplots()
            quantities.append("iota")
        diagnostics = self.diagnostics_stage

        updates_offset = 0
        for ith_stage, stage in enumerate(diagnostics):
            ds = xr.Dataset(diagnostics[stage])[quantities].dropna(dim="updates")
            # skip the setup stage with zero free boundary updates.
            if len(ds.updates) == 0:
                continue

            axs_fb[0].semilogy(
                ds.updates + updates_offset,
                ds.free_boundary_target,
                color=f"C{ith_stage}",
                label="value " + stage.replace("_", " "),
            )
            axs_fb[0].semilogy(
                ds.updates + updates_offset,
                ds.free_boundary_target_mavg,
                color=f"C{ith_stage}",
                linestyle="-.",
                label="mov. avg. " + stage.replace("_", " "),
            )

            axs_fb[1].semilogy(
                ds.updates + updates_offset,
                ds.Fnorm_X1_b,
                color=f"C{ith_stage}",
                label=r"$X^1$ " + stage.replace("_", " "),
            )
            axs_fb[1].semilogy(
                ds.updates + updates_offset,
                ds.Fnorm_X2_b,
                "--",
                color=f"C{ith_stage}",
                label=r"$X^2$ " + stage.replace("_", " "),
            )
            axs_fb[2].semilogy(
                ds.updates + updates_offset,
                ds.resolved_force_ratio,
                color=f"C{ith_stage}",
                label=stage.replace("_", " "),
            )

            axs_fb[3].semilogy(
                ds.updates + updates_offset,
                ds.xnorm,
                color=f"C{ith_stage}",
                label="value " + stage.replace("_", " "),
            )
            axs_fb[3].semilogy(
                ds.updates + updates_offset,
                ds.xnorm_mavg,
                color=f"C{ith_stage}",
                linestyle="-.",
                label="mov. avg. " + stage.replace("_", " "),
            )

            axs_fb[4].semilogy(
                ds.updates + updates_offset,
                ds.gnorm,
                color=f"C{ith_stage}",
                label="value " + stage.replace("_", " "),
            )
            axs_fb[4].semilogy(
                ds.updates + updates_offset,
                ds.gnorm_mavg,
                color=f"C{ith_stage}",
                linestyle="-.",
                label="mov. avg. " + stage.replace("_", " "),
            )

            updates_offset += ds.updates.max()
            if self.state.I_tor:
                max_iter = ds.updates.max()
                for update in ds.updates:
                    ds_aux = ds.sel(updates=update)
                    ax_iota.plot(
                        ds_aux.iota.rho_iota,
                        ds_aux.iota,
                        "-.",
                        color=f"C{ith_stage}",
                        alpha=min(0.01 + (update / (5 * max_iter)).item(), 1),
                    )
                ax_iota.plot(
                    ds_aux.iota.rho_iota,
                    ds_aux.iota,
                    color=f"C{ith_stage}",
                    label=r"final $\iota$ " + stage.replace("_", " "),
                )
        if self.state.I_tor:
            ax_iota.set(xlabel=r"$\rho$", ylabel=r"$\iota$")
            _add_rationals_to_iota_plot(ax_iota, nfp=self.state.nfp)
            ax_iota.legend()
        axs_fb[0].set(ylabel="free boundary target")
        axs_fb[1].set(ylabel=r"$\left\|\hat{\mathbf{F}}_\text{b}\right\|^\text{max}$")
        axs_fb[2].set(ylabel="resolved force ratio")
        axs_fb[3].set(ylabel="xnorm")
        axs_fb[4].set(ylabel="gnorm", xlabel="free boundary update")

        for ax in axs_fb:
            ax.grid("on")
            ax.legend(ncol=2)
        out = ((fig_fb, axs_fb),)
        if self.state.I_tor:
            out += ((fig_iota, ax_iota),)
        return out

    def plot_diagnostics_fixed_boundary(self):
        import matplotlib.pyplot as plt

        dm = self.diagnostics_minimizer
        fig, axs = plt.subplots(2, 1, tight_layout=True, sharex=True)
        axs[0].plot(dm.total_iteration, dm.W_MHD)
        axs[1].semilogy(dm.total_iteration, dm.force_X1, label=r"$X^1$")
        axs[1].semilogy(dm.total_iteration, dm.force_X2, label=r"$X^2$")
        axs[1].semilogy(dm.total_iteration, dm.force_LA, label=r"$\lambda$")
        axs[0].set(ylabel=r"$W_\text{MHD}$")
        axs[1].set(
            xlabel="fixed boundary iteration",
            ylabel=r"$\left\|\hat{\mathbf{F}}_\text{i}\right\|^\text{max}$",
        )
        axs[1].legend()
        return (fig, axs)


# ============= BIEST CALL ===========================
def _get_B2_fp(
    state: State,
    coil_set: CoilSet | Coil,
    ntheta: int,
    nzeta: int,
    I_tor: float,
    diagnostic_output: bool,
    digits: int,
):
    """Setup BIEST with points corresponding to one field period and return the flattened squared exterior field for free boundary gvec.

    Parameters
    ----------
    state : State
        GVEC state object. This is used to extract the boundary points
    coil_set : CoilSet | Coil
        Free-boundary coil set. The coil magnetic field is evaluated at the boundary-integration points.
    ntheta : int
        Number of theta integration points.
    nzeta : int
        Number of zeta integration points.
    I_tor : float
        Toroidal current at the boundary.
    diagnostic_output : bool, optional
        Debug-flag to enable the return of the intermediate fields, by default False.
        If True B2_flat, B_plas_ext, B_coils, B2 for one field period are returned

    Returns
    -------
    np.ndarray(B2_ext)
        Squared exterior field at the GVEC integration points flattened in GVEC order.
    """

    try:
        from fbtools.exterior_vacuum_field import get_B_ext_plasma
    except ImportError as e:
        raise ImportError(
            "Could not import fbtools. Was free-boundary gvec (`gvec[free_boundary]`) installed?"
        ) from e

    theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    theta += 0.5 * theta[1]
    zeta = np.linspace(0, 2 * np.pi / state.nfp, nzeta, endpoint=False)
    zeta += 0.5 * zeta[1]

    # get the boundary for the full device (=all field periods)
    boundary = state.evaluate("pos", rho=1.0, zeta=zeta, theta=theta)
    boundary = boundary.sel(rho=1).transpose("xyz", "tor", "pol")

    # coil-field at the boundary
    ev_coils = coil_set.eval_B(boundary.pos)
    B_coils = ev_coils.B.transpose("xyz", "tor", "pol").values
    if state.clockwise_zeta:
        boundary_biest = np.flip(boundary.pos.values, axis=1)
        B_coils_biest = np.flip(B_coils.copy(), axis=1)
    else:
        boundary_biest = boundary.pos.values
        B_coils_biest = B_coils.copy()

    # main BIEST calculation
    time_0 = time.time()
    B_plas_ext = get_B_ext_plasma(
        boundary_biest,
        state.nfp,
        B_coils_biest,
        I_tor * mu_0,
        digits=digits,
    )
    time_1 = time.time()
    logger.info(f"time for BIEST: {time_1 - time_0:.3}s for {ntheta * nzeta} points")
    if state.clockwise_zeta:
        B_plas_ext = np.flip(B_plas_ext, axis=1)

    # Bext^2
    B_external = B_plas_ext + B_coils

    B2 = np.sum(B_external**2, axis=0)

    B2_flat = np.asfortranarray(B2.flatten(order="C"))
    if diagnostic_output:
        return (
            B2_flat,
            B_plas_ext,
            B_coils,
            B2,
        )
    else:
        return B2_flat


# ========= RESOLVED BOUNDARY FORCES ========================
def minimal_modes_force_balance(X, tolerance=1e-8, **kwargs_fft2d):
    """
    Find the minimal maximum mode numbers (M, N) such that the error is below the tolerance.
    First dimension of X to be theta,
    second dimension is assumed to be zeta.
    """
    Xcos, Xsin = fourier.fft2d(X, **kwargs_fft2d)

    M, N = Xcos.shape[0] - 1, Xcos.shape[1] // 2

    m, n = fourier.fft2d_modes(M, N, grid=True)
    Mrange, Nrange = np.arange(M + 1), np.arange(0, N + 1)
    error = np.full((M + 1, N + 1), np.inf)
    norm = np.sqrt(np.sum(Xcos**2 + Xsin**2))
    for Mnew in Mrange:
        for Nnew in Nrange:
            # sum magnitudes of all modes above the cutoff
            mask = (m > Mnew) | (n > Nnew) | (n < -Nnew)
            err = Xcos[mask] ** 2 + Xsin[mask] ** 2
            error[Mnew, Nnew] = np.sqrt(np.sum(err)) / norm

    # select candidates with error below the tolerance
    mcan, ncan = np.meshgrid(Mrange, Nrange, indexing="ij")
    mask = error < tolerance
    # restrict candidates to those with minimum DoFs
    dofs = ncan + 1 + mcan * (2 * ncan + 1)
    mask &= dofs == dofs[mask].min()
    # select candidate with minimum error
    mask &= error == error[mask].min()

    return mcan[mask].item(), ncan[mask].item()
