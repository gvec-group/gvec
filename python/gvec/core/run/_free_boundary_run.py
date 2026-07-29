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
from gvec.core.state import State, DEFAULT_N_IOTA_RHOS
from ._abstract_run import (
    AbstractRun,
    DEFAULT_MINIMIZE_TOL,
    DEFAULT_TOTALITER,
    AUTO_IOTA_TARGET,
    DEFAULT_PICARD_CURRENT_MAXRESTARTS,
)
from gvec.util import write_parameter_file_ini, flatten_parameters, chdir
from gvec.core.compute import ev2ft
from gvec import fourier
from gvec.lib import _state

try:
    import fbtools

    has_fbtools = True
except ImportError:
    has_fbtools = False


DEFAULT_FORCE_RATIO = 0.01
DEFAULT_FORCE_TOL_START = 1e-3
DEFAULT_ALPHA = 10.0
DEFAULT_BLEND_DEPTH = 0.6
DEFAULT_GTOL = 1e-8
DEFAULT_XTOL = 1e-8
DEFAULT_DIGITS_BIEST = 8
DEFAULT_FB_TOLERANCE = 1e-8

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

        self.stages[0]["free_boundary"] = {"updates": 0}

    def _allocate_diagnostics(self, track_boundary_fields: bool):
        """Setup the diagnostics dataset."""
        diagnostics = self.state.evaluate("xyz", "theta", "zeta", rho=1.0).isel(rad=0)
        fb_updates = max(self._state_parameters["free_boundary"]["updates"], 1)
        diagnostics["force_balance"] = (
            ("tor", "pol", "iteration"),
            np.full(shape=(self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan),
        )
        diagnostics.force_balance.attrs = dict(
            long_name="free-boundary force balance on the LCFS.", symbol=r"$F_{b}}$"
        )
        if track_boundary_fields:
            diagnostics["B_jump"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(
                    shape=(3, self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan
                ),
            )
            diagnostics.B_jump.attrs = dict(
                long_name="jump of the magnetic field on the LCFS", symbol=r"$B_{\text{jump}}$"
            )

            diagnostics["B_int"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(
                    shape=(3, self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan
                ),
            )
            diagnostics.B_int.attrs = dict(
                long_name="internal magnetic field on the LCFS", symbol=r"$B_{\text{int}}$"
            )

            diagnostics["B_ext"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(
                    shape=(3, self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan
                ),
            )
            diagnostics.B_ext.attrs = dict(
                long_name="external magnetic field on the LCFS", symbol=r"$B_{\text{ext}}$"
            )

            diagnostics["K"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(
                    shape=(3, self.state.z_n, self.state.t_n, fb_updates), fill_value=np.nan
                ),
            )
            diagnostics.K.attrs = dict(long_name="surface current on the LCFS", symbol=r"$K$")

        diagnostics["free_boundary_target"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )

        diagnostics["Fnorm_X1_b_IP"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X2_b_IP"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X1_b_modes"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X2_b_modes"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )

        diagnostics["spectral_width"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )

        diagnostics["Fnorm_X2_b"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X1_b"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )
        diagnostics["dt"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )
        diagnostics["resolved_force_ratio"] = (
            "iteration",
            np.full(shape=fb_updates, fill_value=np.nan),
        )

        m, n, F = self.state._get_boundary_forces()
        diagnostics["m_modes"] = ("mode", m)
        diagnostics["n_modes"] = ("mode", n)
        diagnostics["boundary_forces"] = (
            ("modes", "iteration"),
            np.full(shape=(len(m), fb_updates), fill_value=np.nan),
        )

        if self.curr_constraint:
            diagnostics["iota"] = (
                ("rho_iota", "iteration"),
                np.full(shape=(len(self.state.rho2_iota), fb_updates), fill_value=np.nan),
            )
            diagnostics["rho_iota"] = np.sqrt(self.state.rho2_iota)

        if self.diagnostics_stage is None:
            self.diagnostics_stage = xr.DataTree()
        self.track_boundary_fields = track_boundary_fields
        self.diagnostics_stage[f"stage_{self.nth_stage}"] = diagnostics

    def _update_force_ratio(self, iteration):
        force_balance_ds = self.diagnostics_stage[f"stage_{self.nth_stage}"].dataset.isel(
            iteration=iteration
        )
        force_balance_ds["N_FP"] = self.state.nfp
        fourier_force_balance = ev2ft(force_balance_ds[["N_FP", "force_balance"]])
        cut_fb = fourier_force_balance.where(
            (fourier_force_balance.m <= self.m_max), drop=True
        ).where(np.abs(fourier_force_balance.n) <= self.n_max, drop=True)
        full_force = max(
            (fourier_force_balance.force_balance_mns**2).sum().item()
            + (fourier_force_balance.force_balance_mnc**2).sum().item(),
            1e-12,
        )
        cut_force = (cut_fb.force_balance_mns**2).sum().item() + (
            cut_fb.force_balance_mnc**2
        ).sum().item()
        self._resolved_force_ratio = cut_force / full_force
        self.diagnostics_stage[f"stage_{self.nth_stage}"] = force_balance_ds

    def _update_diagnostics(
        self,
        B_ext_2: np.ndarray,
        B_coils: np.ndarray,
        B_plas_ext: np.ndarray,
        iteration: int,
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
        iteration : int
            Determines which entry in the allocated dataset is filled.
        """

        free_boundary_target, force_balance, ev_gvec = self.state._free_boundary_target(B_ext_2)
        B_ext = B_coils + B_plas_ext
        B_jump = ev_gvec.B - B_ext

        diag_aux = self.diagnostics_stage[f"stage_{self.nth_stage}"]
        diag_aux.force_balance[:, :, iteration] = force_balance
        diag_aux.free_boundary_target[iteration] = free_boundary_target
        Fnorm = self.state.get_Fnorm_boundary()
        if self.track_boundary_fields:
            # mu_0*K = n x [[B]]
            diag_aux.K[:, :, :, iteration] = np.cross(ev_gvec.normal, B_jump, axis=0) / mu_0
            diag_aux.B_ext[:, :, :, iteration] = B_ext
            diag_aux.B_jump[:, :, :, iteration] = B_jump.data
            diag_aux.B_int[:, :, :, iteration] = ev_gvec.B.data

        diag_aux.Fnorm_X1_b_IP[iteration] = Fnorm[0]
        diag_aux.Fnorm_X1_b_modes[iteration] = Fnorm[1]
        diag_aux.Fnorm_X2_b_IP[iteration] = Fnorm[2]
        diag_aux.Fnorm_X2_b_modes[iteration] = Fnorm[3]
        diag_aux.Fnorm_X1_b[iteration] = Fnorm[4]
        diag_aux.Fnorm_X2_b[iteration] = Fnorm[4]

        diag_aux.spectral_width[iteration] = _state.get_spectral_width()

        diag_aux.dt[iteration] = dt
        self.free_boundary_target = free_boundary_target.item()
        self._update_force_ratio(iteration)
        diag_aux.resolved_force_ratio[iteration] = self._resolved_force_ratio

        diag_aux.boundary_forces[:, iteration] = self.state._get_boundary_forces()[2]
        self.diagnostics_stage[f"stage_{self.nth_stage}"] = diag_aux

    def _run_single_minimization(self):
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

        # ftol termination criterion
        free_boundary_tolerance = self._state_parameters["free_boundary"].get(
            "tolerance", DEFAULT_FB_TOLERANCE
        )

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

        # Flag to restart the free boundary minimizer after a stage.
        reset_minimizer = self._state_parameters["free_boundary"].get("reset_minimizer", False)

        # initial stepsize for the free-b minimization
        dt = self._state_parameters["free_boundary"].get("dt", 1e-3)

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

        write_parameter_file_ini(
            flatten_parameters(self._state_parameters),
            self.rundir / "parameter.ini",
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )
        with chdir(self.rundir):
            # write parameterfile & run GVEC
            self.state = State(
                parameterfile="parameter.ini",
                statefile=self.state.statefile,
                stdout_path=Path("stdout.txt") if self.redirect_gvec_stdout else None,
                init_for_minimization=True,
            )
            iota_tol = self._state_parameters["picard_current"].get(
                "iota_tol", AUTO_IOTA_TARGET
            )

            self.state.minimize(
                iota_target=iota_tol,
                minimize_tol=minimize_tol,
                maxIter=self._state_parameters["maxIter"],
                maxUpdates=self._state_parameters["picard_current"].get(
                    "maxRestarts", DEFAULT_PICARD_CURRENT_MAXRESTARTS
                ),
                target=self._state_parameters["picard_current"].get("target", "auto"),
                I_tor=self.I_tor_target if self.curr_constraint else None,
                reset_LA=self._state_parameters["picard_current"].get("reset_LA", None),
                n_iota_rhos=self._state_parameters["picard_current"].get(
                    "nPoints", DEFAULT_N_IOTA_RHOS
                ),
            )
            self._allocate_diagnostics(track_boundary_fields)
            self.iters_used = self.state.get_iterations()

            # NOTE: always use self-consistent I_tor values
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

            try:  # allow clean abort with KeybordInterrupt
                for i in range(1, fb_updates):
                    logger.info(
                        f"============= free-boundary step {i} / {fb_updates} ============="
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

                    # abort criteria
                    xnorm = self.state.get_xnorm()
                    gnorm = self.state.get_gnorm()

                    fb_target_reached = self.free_boundary_target <= free_boundary_tolerance
                    xtol_reached = xnorm <= xtol
                    gtol_reached = gnorm <= gtol
                    force_ratio_limit_reached = self._resolved_force_ratio <= force_ratio_limit

                    if (
                        fb_target_reached
                        or xtol_reached
                        or gtol_reached
                        or force_ratio_limit_reached
                    ):
                        logger.info(
                            "Reached free boundary stopping criterion!\n"
                            + f"    free-boundary target: {self.free_boundary_target:.3e}, specified tol.: {free_boundary_tolerance:.3e}.\n"
                            + f"    xnorm: {xnorm:.2e}, xtol: {xtol:.2e}.\n"
                            + f"    gnorm: {gnorm:.2e}, gtol: {gtol:.2e}.\n"
                            + f"    resolved force-balance ratio: {self._resolved_force_ratio:.2e} %, specified limit: {force_ratio_limit:.2e}.\n"
                        )
                        break
                    else:
                        logger.info(
                            f"\n    xnorm: {xnorm:.2e}, xtol: {xtol:.2e}.\n"
                            + f"    gnorm: {gnorm:.2e}, gtol: {gtol:.2e}.\n"
                            + f"    force ratio {self._resolved_force_ratio:.2}, limit {force_ratio_limit:.2}"
                        )
                        t_fb_0 = time.time()
                        # pushes boundary and does a fixed-b run
                        dt = self.state._free_boundary_update(
                            maxIter=min(
                                max_iter_per_update_aux, self._state_parameters["maxIter"]
                            ),
                            B2=B2,
                            dt=dt,
                            reset_minimizer=reset_minimizer,
                            blend_depth=blend_depth,
                            iota_target=iota_tol,
                            minimize_tol=minimize_tol,
                            maxUpdates=self._state_parameters["picard_current"].get(
                                "maxRestarts", DEFAULT_PICARD_CURRENT_MAXRESTARTS
                            ),
                            target=self._state_parameters["picard_current"].get(
                                "target", "auto"
                            ),
                            I_tor=self.I_tor_target if self.curr_constraint else None,
                            reset_LA=self._state_parameters["picard_current"].get(
                                "reset_LA", None
                            ),
                            n_iota_rhos=self._state_parameters["picard_current"].get(
                                "nPoints", DEFAULT_N_IOTA_RHOS
                            ),
                        )

                        if self.save_fb_updates:
                            self.state.save()
                        # NOTE: always use self consistent I_tor
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

                        self._update_diagnostics(B_ext_2, B_coils, B_plas_ext, i, dt)

                        t_fb_1 = time.time()
                        logger.info(f"Time for free boundary update: {t_fb_1 - t_fb_0:.4}s")
                        reset_minimizer = False
                        self.iters_used = self.state.get_iterations()
                        if not self.curr_constraint:
                            self.I_tor_boundary = self.state.evaluate(
                                "I_tor", rho=1.0
                            ).I_tor.item()
                        else:
                            self.diagnostics_stage[f"stage_{self.nth_stage}"].iota[:, i] = (
                                self.state.iota.data
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
                if self._state_parameters["maxIter"] > 0:
                    self.state.minimize(
                        iota_target=iota_tol,
                        minimize_tol=minimize_tol_min,
                        maxIter=self._state_parameters["maxIter"],
                        maxUpdates=self._state_parameters["picard_current"].get(
                            "maxRestarts", 30
                        ),
                        target=self._state_parameters["picard_current"].get("target", "auto"),
                    )

            except KeyboardInterrupt:
                logger.warning("Aborting free-boundary due to keyboard interrupt!")
            if (
                self._resolved_force_ratio < force_ratio_limit
                and self.free_boundary_target > free_boundary_tolerance
            ):
                fb_aux = (
                    self.diagnostics_stage.force_balance.isel(iteration=i)
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

        iterations = self.state.get_iterations()
        iteration_offset = self.GVEC_iter_used
        self.GVEC_iter_used += iterations

        # postprocessing
        statefiles_found = sorted(self.rundir.glob("*State*.dat"))
        if len(statefiles_found) > 0:
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


# ============= BIEST CALLS ===========================
def _get_B2_all_fp(
    state: State,
    coil_set: CoilSet | Coil,
    ntheta: int,
    nzeta: int,
    I_tor: float,
    diagnostic_output=False,
    fac: int = 1,
):
    """Setup BIEST with points corresponding to all field periods and return the flattened squared exterior field for free boundary gvec.
    All field periods are used to avoid automatic copying of from inside BIEST which may expect different orientation of the points.

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
    fac: int
        Factor which increases the number of source points with respect to the GVEC integration points.
        For example, if fac = 2, twice as many points in theta (=2*ntheta) and zeta (=2*NFP*nzeta) will be used to evaluate the field,
        but still only the GVEC integration points will be returned. That is B2_ext will be of size(ntheta*nzeta).
        Defaults to 1.

    Returns
    -------
    np.ndarray(B2_ext)
        Squared exterior field at the GVEC integration points flattened in GVEC order.
    """

    if has_fbtools:
        from fbtools.exterior_vacuum_field import get_B_ext_plasma
    else:
        logger.warning(
            "Could not import 'fbtools'! Please make sure the package is installed. Aborting Bext calulation."
        )
        return

    # upsample the gvec internal integration points
    ntheta_aux = fac * ntheta
    nzeta_aux = fac * nzeta

    theta = np.linspace(0, 2 * np.pi, ntheta_aux, endpoint=False)
    theta += 0.5 * theta[1]
    N_FP = state.evaluate("N_FP").N_FP
    zeta = np.linspace(0, 2 * np.pi, nzeta_aux * N_FP.item(), endpoint=False)
    zeta += 0.5 * zeta[1]

    # get the boundary for the full device (=all field periods)
    boundary = state.evaluate("pos", "N_FP", rho=1.0, theta=theta, zeta=zeta)
    boundary = boundary.sel(rho=1).transpose("xyz", "tor", "pol")

    # coil-field at the boundary
    ev_coils = coil_set.eval_B(boundary.pos)
    B_coils = ev_coils.B.transpose("xyz", "tor", "pol").values

    # main BIEST calculation
    time_0 = time.time()
    B_plas_ext = get_B_ext_plasma(
        boundary.pos.values,
        1,  # one FP!
        B_coils,
        I_tor * mu_0,
        digits=10,
    )
    time_1 = time.time()
    logger.info(
        f"time for BIEST: {time_1 - time_0:.3}s for {ntheta_aux * nzeta_aux * N_FP.item()} points"
    )

    # Bext^2
    B_external = B_plas_ext + B_coils
    B2 = np.sum(B_external**2, axis=0)

    B2 = B2[:nzeta_aux:fac, ::fac]
    B2_flat = np.ones(ntheta * nzeta) * 1e6
    counter = 0

    # flatten according to internal GVEC order
    for nIP in range(nzeta):
        for mIP in range(ntheta):
            B2_flat[counter] = B2[nIP, mIP]
            counter += 1
    B2_flat = np.asfortranarray(B2_flat)
    if diagnostic_output:
        return (
            B2_flat,
            B_plas_ext[:, :nzeta_aux:fac, ::fac],
            B_coils[:, :nzeta_aux:fac, ::fac],
            B2,
        )
    else:
        return B2_flat


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

    if has_fbtools:
        from fbtools.exterior_vacuum_field import get_B_ext_plasma
    else:
        logger.warning(
            "Could not import 'fbtools'! Please make sure the package is installed. Aborting Bext calulation."
        )
        return

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

    B2_flat = np.ones(ntheta * nzeta) * 1e6
    counter = 0

    # flatten according to internal GVEC order
    for nIP in range(nzeta):
        for mIP in range(ntheta):
            B2_flat[counter] = B2[nIP, mIP]
            counter += 1
    B2_flat = np.asfortranarray(B2_flat)
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
