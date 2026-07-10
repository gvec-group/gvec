import numpy as np
from gvec import State, fluxsurface_integral
from gvec.coils import CoilSet, Coil
from gvec.util import read_parameters, write_parameters, CaseInsensitiveDict, chdir
from gvec.core.compute import ev2ft
from gvec import fourier
from gvec.core.state import with_binding
from scipy.constants import mu_0
import shutil
from tqdm import tqdm
from typing import Callable

from gvec.lib import modgvec_py_state as _state

import time
import logging
from pathlib import Path

try:
    import fbtools

    has_fbtools = True
except ImportError:
    has_fbtools = False

logger = logging.getLogger("gvec.free-boundary")


class FreeBoundaryState(State):
    def __init__(
        self,
        parameterfile: str | Path,
        coil_set: CoilSet | Coil,
        fb_updates: int,
        totalIter: int,
        max_iter_per_update: int = int(1e12),
        free_boundary_tolerance: float = 1e-8,
        I_tor: float | Callable | None = None,
        iter_B_update: int = 1,
        dt: float = 0.001,
        statefile: str | Path | None = None,
        reset_LA: bool = False,
        track_boundary_fields: bool = False,
    ):
        """Free-boundary state object.

        Parameters
        ----------
        parameterfile : str | Path
            Parameter.ini file used for initialization.
        coil_set : CoilSet | Coil
            Coils that are used for the free-boundary equilibrium.
        fb_updates : int
            Maximum number of free-boundary updates.
        totalIter : int
            Maximum number of GVEC iterations.
        max_iter_per_update : int, optional
            Maximum number of GVEC iterations per free-boundary update, by default int(1e12). The default tries to converge all fixed-boundary steps.
        free_boundary_tolerance : float, optional
            Targeted tolerance on the free-boundary target, by default 1e-8
        I_tor : float | Callable | None, optional
            Toroidal current profile, by default None. If None, fixed-iota runs are preformed.
        iter_B_update : int, optional
            Number of iterations after which BIEST is called again, by default 1. If >1, the external magnetic field is linearly extrapolated between updates.
        dt : float, optional
            Stepsize for the free-boundary push, by default 0.001
        statefile : str | Path | None, optional
            Statefile to restart from, by default None
        reset_LA: bool, optional
            Whether to keep the previous lambda after a boundary push.
            If True, lambda is either set to zero or reinitialized depending on `init_LA`.
            The default is False.
        track_boundary_fields : bool, optional
            Whether to track relevant fields on the boundary in the diagnostics data-set.
            The default is False. Note that tracking the fields can result in comparatively large data-sets.
        """
        super().__init__(parameterfile, statefile, True)
        self.coil_set = coil_set

        # these parameters are modified during current constraint!
        self._parameter_dict = read_parameters(parameterfile)

        self.fb_updates = fb_updates
        self.totalIter = totalIter

        self.reset_LA = reset_LA

        if max_iter_per_update > totalIter:
            logger.warning(
                "max_iter_per_update > totalIter! Setting max_iter_per_update = totalIter"
            )
            max_iter_per_update = totalIter
        self.max_iter_per_update = max_iter_per_update

        self.free_boundary_tolerance = free_boundary_tolerance
        self.iter_B_update = iter_B_update
        self.dt = dt

        if I_tor is not None:
            self.current_constraint = True
            self.rho_iota = np.linspace(0, 1, 101)
            if callable(I_tor):
                self.I_tor = I_tor(self.rho_iota)
            else:
                self.I_tor = I_tor * np.ones_like(self.rho_iota)
            self.I_tor_boundary = self.I_tor[-1]
        else:
            self.current_constraint = False
            self.rho_iota = None
            self.I_tor = None
            self.I_tor_boundary = self.evaluate("I_tor", rho=1.0).I_tor.item()

        # initialize diagnostics dataset for quantities at the boundary for each free-boundary update

        r_GP, r_w, t_n, t_w, z_n, z_w = self.get_integration_points()

        self.t_n = t_n
        self.z_n = z_n

        self.m_max, self.n_max = self.get_mn_max()
        self.track_boundary_fields = track_boundary_fields
        self._allocate_diagnostics()

        self.iter_B2_history = np.array([0, 1])
        self.B2_history = np.empty([t_n * z_n, 2])
        self.iters_used = 0
        self.free_boundary_target = 1e30
        self._resolved_force_ratio = 1.0

        self.clockwise_zeta = self._determine_clockwise_zeta()

        if self.current_constraint:
            self.diagnostics["iota"] = (
                ("rho_iota", "iteration"),
                np.full(shape=(len(self.rho_iota), self.fb_updates), fill_value=np.nan),
            )
            self.diagnostics["rho_iota"] = self.rho_iota

    @with_binding
    def _minimize(self, dt=None, logfinaliter=True, maxIter=None, minimize_tol=None):
        """Minimize the MHD functional for the current state."""
        logger.debug("Minimizing MHD functional.")
        with chdir(self.rundir):
            _state.minimize(dt, logfinaliter, maxIter, minimize_tol)

    @with_binding
    def set_Bext_squared(self, B2):
        logger.debug("Updating exterior field")
        _state.set_bext_squared(B2)

    @with_binding
    def save_statefile(self):
        logger.debug("Writing state to file")
        _state.save_statefile()

    @with_binding
    def update_iota_for_current_constraint(self, rho2, iota, reset_LA):
        logger.debug("Updating iota")
        _state.update_iota_for_current_constraint(rho2, iota, reset_LA)

    @with_binding
    def _get_boundary_forces(self):
        n_X1 = self._get_X1_n_modes()
        n_X2 = self._get_X2_n_modes()
        forces_out = np.asfortranarray(np.zeros((3, n_X1 + n_X2)))
        _state.get_boundary_forces(n_X1, n_X2, forces_out)
        m = forces_out[0, :]
        n = forces_out[1, :]
        F = forces_out[2, :]
        return (m, n, F)

    @with_binding
    def get_Fnorm_boundary(self):
        logger.debug("Extracting Fnorm at the boundary")
        return _state.get_fnorm_boundary()

    @with_binding
    def get_outputlevel(self):
        logger.debug("Extracting outputlevel")
        outid = _state.get_outputlevel()
        return outid

    @with_binding
    def get_iterations(self):
        logger.debug("Extracting iterations")
        iters = _state.get_iterations()
        return iters

    @with_binding
    def get_minimize_tol(self):
        logger.debug("Extracting minimize tol")
        tol = _state.get_minimize_tol()
        return tol

    @with_binding
    def get_gnorm(self):
        logger.debug("Extracting gnorm")
        gnorm = _state.get_gnorm()
        return gnorm

    @with_binding
    def get_xnorm(self):
        logger.debug("Extracting xnorm")
        xnorm = _state.get_xnorm()
        return xnorm

    @with_binding
    def reset_descent(self):
        logger.debug("Resetting accelerated descent")
        _state.reset_descent()

    @with_binding
    def update_boundary(self, dt: float, reset_minimizer: bool, blend_depth: float = 0.6):
        logger.debug("Stepping in free boundary direction")
        dt = _state.update_boundary(dt, reset_minimizer, blend_depth=blend_depth)
        return dt

    @with_binding
    def save(self):
        self.save_statefile()
        iters = self.get_iterations()
        outid = self.get_outputlevel()
        write_parameters(
            self._parameter_dict,
            self._parameter_dict["ProjectName"] + f"_parameters_{outid:04}_{iters:08}.ini",
        )

    @with_binding
    def _determine_clockwise_zeta(self):
        """Rotates a point of the axis with the Rodrigues equation to see if zeta runs clockwise or counter-clockwise.

        Returns
        -------
        bool
            Whether zeta runs clockwise ort not.

        Raises
        ------
        ValueError
            If there is an implementation error that causes this approach to fail.
        """

        zeta = 2 * np.pi / self.nfp
        ev = self.evaluate("pos", rho=0.0, theta=0.0, zeta=[0, 2 * np.pi / self.nfp])
        v0 = ev.pos.isel(pol=0, tor=0, rad=0).data
        v1 = ev.pos.isel(pol=0, tor=1, rad=0).data
        k = np.array([0, 0, 1])
        vrot = (
            v0 * np.cos(zeta)
            + np.sin(zeta) * np.cross(k, v0)
            + (1 - np.cos(zeta)) * np.dot(k, v0) * k
        )
        zeta_inv = 2 * np.pi - zeta
        vrot_inv = (
            v0 * np.cos(zeta_inv)
            + np.sin(zeta_inv) * np.cross(k, v0)
            + (1 - np.cos(zeta_inv)) * np.dot(k, v0) * k
        )

        if all(abs(vrot - v1) < 1e-12):
            clockwise_zeta = False
        elif all(abs(vrot_inv - v1) < 1e-12):
            clockwise_zeta = True
        else:
            raise ValueError("Issue when determining zeta-rotation.")
        return clockwise_zeta

    def _allocate_diagnostics(self):
        """Setup the diagnostics dataset."""
        diagnostics = self.evaluate("xyz", "theta", "zeta", rho=1.0).isel(rad=0)
        diagnostics["force_balance"] = (
            ("tor", "pol", "iteration"),
            np.full(shape=(self.z_n, self.t_n, self.fb_updates), fill_value=np.nan),
        )
        diagnostics.force_balance.attrs = dict(
            long_name="free-boundary force balance on the LCFS.", symbol=r"$F_{b}}$"
        )
        if self.track_boundary_fields:
            diagnostics["B_jump"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(shape=(3, self.z_n, self.t_n, self.fb_updates), fill_value=np.nan),
            )
            diagnostics.B_jump.attrs = dict(
                long_name="jump of the magnetic field on the LCFS", symbol=r"$B_{\text{jump}}$"
            )

            diagnostics["B_int"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(shape=(3, self.z_n, self.t_n, self.fb_updates), fill_value=np.nan),
            )
            diagnostics.B_int.attrs = dict(
                long_name="internal magnetic field on the LCFS", symbol=r"$B_{\text{int}}$"
            )

            diagnostics["B_ext"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(shape=(3, self.z_n, self.t_n, self.fb_updates), fill_value=np.nan),
            )
            diagnostics.B_ext.attrs = dict(
                long_name="external magnetic field on the LCFS", symbol=r"$B_{\text{ext}}$"
            )

            diagnostics["K"] = (
                ("xyz", "tor", "pol", "iteration"),
                np.full(shape=(3, self.z_n, self.t_n, self.fb_updates), fill_value=np.nan),
            )
            diagnostics.K.attrs = dict(long_name="surface current on the LCFS", symbol=r"$K$")

        diagnostics["free_boundary_target"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )

        diagnostics["Fnorm_X1_b_IP"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X2_b_IP"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X1_b_modes"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X2_b_modes"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )

        diagnostics["spectral_width"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )

        diagnostics["Fnorm_X2_b"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )
        diagnostics["Fnorm_X1_b"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )
        diagnostics["dt"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )
        diagnostics["resolved_force_ratio"] = (
            "iteration",
            np.full(shape=self.fb_updates, fill_value=np.nan),
        )

        m, n, F = self._get_boundary_forces()
        diagnostics["m_modes"] = ("mode", m)
        diagnostics["n_modes"] = ("mode", n)
        diagnostics["boundary_forces"] = (
            ("modes", "iteration"),
            np.full(shape=(len(m), self.fb_updates), fill_value=np.nan),
        )

        self.diagnostics = diagnostics

    @with_binding
    def _get_new_iota(self):
        """Calculate a new iota profile. Extracted from gvec-run, to be used for free-boundary current constraint

        Returns
        -------
        (xr.DataArray, xr.DataArray, xr.DataArray)
            new iota values, rho^2 interpolation points, root mean square iota error
        """
        rho_eval = np.concatenate(
            [[np.sqrt(1e-8), np.sqrt(2e-8), np.sqrt(3e-8)], self.rho_iota[1:]]
        )
        quantites = ["iota", "iota_curr_0", "iota_0", "I_tor"]
        ev = self.evaluate(*quantites, rho=rho_eval)
        ev = ev[quantites]
        r1 = ev.isel(rad=0)
        r2 = ev.isel(rad=1)
        r3 = ev.isel(rad=2)
        ev = ev.isel(rad=slice(2, None))

        # workaround to modify xarray coordinate & index (with pandas >=3.0)
        rho = ev.rho.data.copy()
        rho[0] = 0.0
        ev = ev.assign_coords(rho=("rad", rho)).set_xindex("rho")

        # extrapolation to axis
        for var in ev.data_vars:
            ev[var].data[0] = 3 * (r1[var].data - r2[var].data) + r3[var].data

        iota_values = ev.iota_0 + self.I_tor * ev.iota_curr_0
        delta_iota = ev.iota - iota_values
        rms_iota = np.sqrt((delta_iota**2).mean("rad"))
        return (iota_values, ev.rho**2, rms_iota)

    @with_binding
    def _minimize_single_current_constraint(
        self,
        maxIter: int,
        maxIter_per_update: int = int(1e12),
        iota_tol: float = 1e-10,
        minimize_tol: float | None = None,
        maxUpdates: int = 30,
        reset_LA: bool = False,
    ):
        """Single current optimization stage.

        Parameters
        ----------
        maxIter : int
            Maximum number of fixed-boundary iterations for this stage.
        maxIter_per_update : int, optional
            Maximum number of fixed-boundary iterations per iota update, by default int(1e12).
        iota_tol : float, optional
            Tolerance on the RMS of the iota update, by default 1e-10.
        minimize_tol : float | None, optional
            Tolerance for the norm of the forces, by default None
        maxUpdates : int, optional
            Maximum number of iota updates, by default 30

        Returns
        -------
        float
            RMS of the iota update.
        """
        rms_iota = 1e6
        nth_update = -1
        iters_used_start = self.get_iterations()
        self.reset_descent()
        iters_used = 0

        while (rms_iota > iota_tol) and (iters_used < maxIter):
            if nth_update + 1 > maxUpdates:
                logger.warning("Maximum number of iota updates reached reached!")
                break
            self._minimize(
                logfinaliter=False,
                maxIter=min(maxIter - iters_used, maxIter_per_update),
                minimize_tol=minimize_tol,
            )
            new_iota, new_rho2, rms_iota = self._get_new_iota()

            # only reset LA at the start
            if iters_used > 0:
                reset_LA = False

            self.update_iota_for_current_constraint(new_rho2, new_iota, reset_LA)
            self._parameter_dict["iota"] = CaseInsensitiveDict(
                type="interpolation",
                rho2=new_rho2.data,
                vals=new_iota.data,
            )
            self.iota = new_iota

            iters_used = self.get_iterations() - iters_used_start

            nth_update += 1

        if rms_iota > iota_tol:
            logger.warning(
                "Targeted iota has not been reached!\n"
                + f"iota_tol.: {iota_tol:.2e}, achieved rms Δiota.: {rms_iota:.2e}"
            )
        logger.debug(
            f"Terminated current minimization after {nth_update} updates. Δiota = {rms_iota:.3e}, target Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
        )
        return rms_iota

    @with_binding
    def _minimize_current_constraint(
        self,
        iota_target: float,
        minimize_target: float,
        maxIter: int,
        maxUpdates: int = 30,
        reset_LA: bool = False,
    ):
        """Fixed-boundary run with prescribed current profile.

        Parameters
        ----------
        iota_target : float
            Tolerance on the RMS of the iota update.
        minimize_target : float
            Tolerance for the norm of the forces.
        maxIter : int
            Maximum number of iterations for the fixed boundary run.
        maxUpdates : int, optional
            Maximum number of iota-updates per stage, by default 30
        """
        t0 = time.time()
        log_minimize_target = np.log10(minimize_target)
        log_iota_target = np.log10(iota_target)
        n_stages = max(int(max(-2 - log_minimize_target, log_minimize_target)), 1)
        minimize_tols = np.logspace(max(-4, log_minimize_target), log_minimize_target, n_stages)
        if n_stages > 1:
            iota_tols = np.logspace(max(-3, log_iota_target), log_iota_target, n_stages)
        else:
            iota_tols = [1e-10]

        iota_tol = iota_tols[0]

        iters_used_start = self.get_iterations()
        # target iota only
        rms_iota = self._minimize_single_current_constraint(
            maxIter=maxIter,
            maxIter_per_update=min(10, maxIter),
            iota_tol=iota_tol,
            minimize_tol=minimize_tols[0],
            maxUpdates=maxUpdates,
            reset_LA=reset_LA,
        )

        iters_used = self.get_iterations() - iters_used_start
        if iters_used < maxIter:
            for minimize_tol, iota_tol in zip(minimize_tols, iota_tols):
                logger.debug(f"minimize_tol:{minimize_tol}")
                logger.debug(f"iota_tol:{minimize_tol}")
                logger.debug(f"iters_used_start:{iters_used_start}")
                logger.debug(f"iters_used:{iters_used}")
                logger.debug(f"maxIter:{maxIter - iters_used}")
                rms_iota = self._minimize_single_current_constraint(
                    maxIter=maxIter - iters_used,
                    iota_tol=iota_tol,
                    minimize_tol=minimize_tol,
                    maxUpdates=maxUpdates,
                )
                iters_used = self.get_iterations() - iters_used_start
        t1 = time.time()
        logger.info(
            f"Terminated current minimization after {t1 - t0:.2}s with Δiota = {rms_iota:.3e}, Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
        )
        self.save()

    @with_binding
    def _free_boundary_target(self, B_ext_2: np.ndarray):
        """Evaluate the free-boundary target function based on the force-balance.

        Parameters
        ----------
        B_ext_2 : np.ndarray
            Squared external magnetic field

        Returns
        -------
        float
            free-boundary target
        """
        ev_gvec = (
            self.evaluate("B", "p", "normal", "Phi_edge", rho=1.0)
            .sel(rho=1.0)
            .transpose("xyz", "tor", "pol")
        )
        B_gvec_2 = (ev_gvec.B**2).sum(dim="xyz")
        force_balance = (B_gvec_2 - B_ext_2) / 2 + mu_0 * ev_gvec.p

        B_avg = (B_gvec_2 + B_ext_2) / 2
        # flux surface averaged normalized force balance
        # <|[[B²/2+mu_0*P]] / B² |>
        free_boundary_target = fluxsurface_integral(
            abs(force_balance / B_avg) * ev_gvec.Jac
        ) / fluxsurface_integral(ev_gvec.Jac)
        return free_boundary_target, force_balance, ev_gvec

    @with_binding
    def _update_force_ratio(self, iteration):
        force_balance_ds = self.diagnostics.isel(iteration=iteration)
        force_balance_ds["N_FP"] = self.nfp
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

    @with_binding
    def _update_diagnostics(
        self, B_ext_2: np.ndarray, B_coils: np.ndarray, B_plas_ext: np.ndarray, iteration: int
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

        free_boundary_target, force_balance, ev_gvec = self._free_boundary_target(B_ext_2)
        B_ext = B_coils + B_plas_ext
        B_jump = ev_gvec.B - B_ext

        self.diagnostics.force_balance[:, :, iteration] = force_balance
        self.diagnostics.free_boundary_target[iteration] = free_boundary_target
        Fnorm = self.get_Fnorm_boundary()
        if self.track_boundary_fields:
            # mu_0*K = n x [[B]]
            self.diagnostics.K[:, :, :, iteration] = (
                np.cross(ev_gvec.normal, B_jump, axis=0) / mu_0
            )
            self.diagnostics.B_ext[:, :, :, iteration] = B_ext
            self.diagnostics.B_jump[:, :, :, iteration] = B_jump.data
            self.diagnostics.B_int[:, :, :, iteration] = ev_gvec.B.data

        self.diagnostics.Fnorm_X1_b_IP[iteration] = Fnorm[0]
        self.diagnostics.Fnorm_X1_b_modes[iteration] = Fnorm[1]
        self.diagnostics.Fnorm_X2_b_IP[iteration] = Fnorm[2]
        self.diagnostics.Fnorm_X2_b_modes[iteration] = Fnorm[3]
        self.diagnostics.Fnorm_X1_b[iteration] = Fnorm[4]
        self.diagnostics.Fnorm_X2_b[iteration] = Fnorm[4]

        self.diagnostics.spectral_width[iteration] = _state.get_spectral_width()

        self.diagnostics.dt[iteration] = self.dt
        self.free_boundary_target = free_boundary_target.item()
        self._update_force_ratio(iteration)
        self.diagnostics.resolved_force_ratio[iteration] = self._resolved_force_ratio

        self.diagnostics.boundary_forces[:, iteration] = self._get_boundary_forces()[2]

    @with_binding
    def _free_boundary_update(
        self,
        maxIter: int,
        B2: np.ndarray,
        tol: float,
        copy_stdout: bool = True,
        reset_minimizer: bool = False,
        blend_depth: float = 0.6,
    ):
        """Single free boundary push.

        Updates the FORTRAN B_ext^2, resets the fixed-boundary minimizer,
        pushes the boundary and minimizes the fixed-boundary problem

        Parameters
        ----------
        maxIter : int
            Maximum number of iterations for the fixed-boundary minimization.
        B2 : np.ndarray
            External magnetic field, flattened in GVEC order (zeta->theta)
        tol : float
            Tolerance on fixed boundary forces.
        copy_stdout : bool, optional
            Whether to save stdout to stoudt.txt, by default True
        reset_minimizer: bool, optional
            Whether to reset the free boundary minimizers to their initial states.
        """
        # set Bext
        self.set_Bext_squared(B2)
        self.reset_descent()
        self.dt = self.update_boundary(self.dt, reset_minimizer, blend_depth)

        # get minimizer tolerance
        if self.current_constraint:
            self._minimize_current_constraint(
                iota_target=1e-10,
                minimize_target=tol,
                maxIter=maxIter,
                reset_LA=self.reset_LA,
            )
        else:
            # fixed iota
            self._minimize(logfinaliter=False, maxIter=maxIter)

        if copy_stdout:
            shutil.copyfile(self._stdout.name, "stdout.txt")

    @with_binding
    def run_free_boundary(
        self,
        save_diagnostics: bool = False,
        force_ratio_limit: float = 1e-1,
        ramp_force_tol=True,
        ramp_force_tol_start: float = 1e-3,
        alpha: float = 10.0,
        blend_depth: float = 0.6,
        gtol: float = 1e-8,
        xtol: float = 1e-8,
        digits_BIEST: int = 8,
    ):
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
        minimize_tol_min = self.get_minimize_tol()
        if ramp_force_tol:
            minimize_tol = ramp_force_tol_start
            current_lowest_min_tol = ramp_force_tol_start
        else:
            minimize_tol = minimize_tol_min
        reset_minimizer = False  # Flag to restart minimizers after a stage.

        if self.current_constraint:
            self._minimize_current_constraint(
                iota_target=1e-10,
                minimize_target=minimize_tol,
                maxIter=self.max_iter_per_update,
            )
        else:
            # fixed iota
            self._minimize(logfinaliter=False, maxIter=self.max_iter_per_update)
        self.iters_used = self.get_iterations()
        try:  # allow clean abort with KeybordInterrupt
            for i in tqdm(range(self.fb_updates)):
                # abort if budget is used up
                if self.iters_used >= self.totalIter:
                    logger.warning(
                        "Maximum number of iterations reached! aborting free-boundary iterations."
                    )
                    break

                # ensure that the maximum number of iterations is kept
                max_iter_per_update_aux = min(
                    self.totalIter - self.iters_used, self.max_iter_per_update
                )

                # TODO: use the consistent current or the prescribed current?
                # if self.current_constraint:
                #     self.I_tor_boundary = self.evaluate("I_tor", rho=1.0).squeeze().I_tor.item()

                # whether to interpolate the external field or not
                if i < 2:
                    (B2, B_plas_ext, B_coils, B_ext_2) = _get_B2_fp(
                        state=self,
                        coil_set=self.coil_set,
                        ntheta=self.t_n,
                        nzeta=self.z_n,
                        I_tor=self.I_tor_boundary,
                        diagnostic_output=True,
                        digits=digits_BIEST,
                    )
                    self.B2_history[:, i] = B2

                    self._update_diagnostics(B_ext_2, B_coils, B_plas_ext, i)

                else:
                    if i % self.iter_B_update == 0:
                        (B2, B_plas_ext, B_coils, B_ext_2) = _get_B2_fp(
                            state=self,
                            coil_set=self.coil_set,
                            ntheta=self.t_n,
                            nzeta=self.z_n,
                            I_tor=self.I_tor_boundary,
                            diagnostic_output=True,
                            digits=digits_BIEST,
                        )
                        self.B2_history[:, 0] = self.B2_history[:, 1]
                        self.B2_history[:, 1] = B2
                        self.iter_B2_history[0] = self.iter_B2_history[1]
                        self.iter_B2_history[1] = i

                        self._update_diagnostics(B_ext_2, B_coils, B_plas_ext, i)
                    else:
                        # linear interpolation of the last two B2 solutions along the iterations
                        B2 = self.B2_history[:, 0] + (i - self.iter_B2_history[0]) * (
                            self.B2_history[:, 1] - self.B2_history[:, 0]
                        ) / (self.iter_B2_history[1] - self.iter_B2_history[0])

                if ramp_force_tol:
                    # minimize_tol is ahead of free boundary target
                    minimize_tol = np.max(
                        [
                            np.min([self.free_boundary_target / alpha, ramp_force_tol_start]),
                            minimize_tol_min,
                        ]
                    )
                    # this is to ensure that minimize tol never increases
                    if minimize_tol < current_lowest_min_tol:
                        current_lowest_min_tol = minimize_tol
                    minimize_tol = np.min([minimize_tol, current_lowest_min_tol])
                    logger.info(f"Setting gvec minimize_tol: {minimize_tol:.2e}")

                # abort criteria
                xnorm = self.get_xnorm()
                gnorm = self.get_gnorm()
                fb_target_reached = self.free_boundary_target <= self.free_boundary_tolerance
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
                        + f"    free-boundary target: {self.free_boundary_target:.3e}, specified tol.: {self.free_boundary_tolerance:.3e}.\n"
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
                    self._free_boundary_update(
                        maxIter=max_iter_per_update_aux,
                        B2=B2,
                        tol=minimize_tol,
                        reset_minimizer=reset_minimizer,
                        blend_depth=blend_depth,
                    )
                    t_fb_1 = time.time()
                    logger.info(f"Time for free boundary update: {t_fb_1 - t_fb_0:.4}s")
                    reset_minimizer = False
                    self.iters_used = self.get_iterations()
                    if not self.current_constraint:
                        self.I_tor_boundary = self.evaluate("I_tor", rho=1.0).I_tor.item()
                    else:
                        self.diagnostics.iota[:, i] = self.iota.data

                logger.info(
                    f"free-boundary - target: {self.free_boundary_target:.3e}, tolerance: {self.free_boundary_tolerance}, iters: {self.iters_used} / {self.totalIter} "
                )
                if save_diagnostics:
                    try:
                        self.diagnostics.to_netcdf("diagnostics_aux.nc")
                    except PermissionError:
                        logger.warning(
                            "Permission error when attemting to write diagnostics to netcdf! Maybe the data-set is opened somewhere else?"
                        )
            # ensure that the final result is in fixed-boundary equilibrium
            self._minimize(logfinaliter=False, maxIter=(self.totalIter - self.iters_used))
        except KeyboardInterrupt:
            logger.warning("Aborting free-boundary due to keyboard interrupt!")
        if (
            self._resolved_force_ratio < force_ratio_limit
            and self.free_boundary_target > self.free_boundary_tolerance
        ):
            fb_aux = (
                self.diagnostics.force_balance.isel(iteration=i)
                .squeeze()
                .transpose("pol", "tor")
            )
            M_proposed, N_proposed = minimal_modes_force_balance(
                fb_aux,
                zeta0=fb_aux.zeta[0].item(),
                theta0=fb_aux.theta[0].item(),
                tolerance=self.free_boundary_tolerance,
            )
            logger.info(
                "Not all modes in the free-boundary force  balance can be resolved with the provided angular resolution.\n"
                + "Try reducing minimze_tol or adapt the angular resolution.\n"
                + "Proposed angular resolution increase for the solution variables X1 / X2 / LA:\n"
                + f"    m_max - current {self.m_max} -> proposed {M_proposed}\n"
                + f"    n_max - current {self.n_max} -> proposed {N_proposed}"
            )

        if self.free_boundary_target > self.free_boundary_tolerance:
            logger.warning(
                f"Free-boundary target not reached! target: {self.free_boundary_target:.3e}, specified tol.: {self.free_boundary_tolerance:.3e}"
            )


def _get_B2_all_fp(
    state: State | FreeBoundaryState,
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
    state: State | FreeBoundaryState,
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


def minimal_modes_force_balance(X, tolerance=1e-8, **kwargs):
    """
    Find the minimal maximum mode numbers (M, N) such that the error is below the tolerance.
    First dimension of X to be theta,
    second dimension is assumed to be zeta.
    """
    Xcos, Xsin = fourier.fft2d(X, **kwargs)

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
