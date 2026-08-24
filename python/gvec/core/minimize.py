# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.minimize - contains the MinimizeState class"""

# === Imports === #

from .state import State
from .state._binding import with_binding, bind
from gvec.lib import _state

import numpy as np

import scipy
from scipy.interpolate import BSpline
from scipy.constants import mu_0

from gvec.util import (
    write_parameter_file_ini,
    CaseInsensitiveDict,
    chdir,
    flatten_parameters,
    read_parameters,
)
from typing import Literal, Callable
from gvec.core.compute import fluxsurface_integral
from gvec import __version__
from datetime import datetime

import time
import logging

logger = logging.getLogger("gvec.state")

AUTO_STAGE_IOTA_INIT_MAXITER = 10
AUTO_STAGE_IOTA_INIT_IOTA_TOL = 1e-3
AUTO_PICARD_CURRENT_IOTA_TOL = 1e-5
AUTO_MINIMIZE_TOL_MIN = 1e-3
IOTA_BSPLINE_DEGREE = 5
IOTA_SAMPLING_FACTOR = 4
IOTA_KNOTS_PER_ELEMENT_MIN = 2
IOTA_KNOTS_PER_ELEMENT_MAX = 256
IOTA_COEFS_MAX = 768  # at 20 characters per coefficient: 15360 characters + overhead < 2^14
IOTA_MIN_RHO = 1e-2  # ignore near-axis region for Δiota & Δfit, TODO: fix accuracy problems
IOTA_MAXITER_PER_UPDATE = 1e12  # maximum number of fixed-boundary iterations per iota update


class MinimizationState(State):
    def __init__(self, parameterfile, statefile=None, stdout_path=None):
        super().__init__(parameterfile, statefile, stdout_path)

        # these parameters are modified during current constraint!
        self._parameter_dict = read_parameters(parameterfile)

        r_GP, r_w, t_n, t_w, z_n, z_w = self.get_integration_points()

        self.t_n = t_n
        self.z_n = z_n

        self.iters_used = 0
        self.rms_iota = None
        self.delta_iota = None
        self.iota = None

    def bind(self, force: bool = False):
        bind(self, force=force, init_for_minimize=True)

    logger = logging.getLogger("gvec.state")

    @with_binding
    def _minimize(
        self,
        dt: float | None = None,
        logfinaliter: bool = True,
        maxIter: int | None = None,
        minimize_tol: float | None = None,
    ):
        """Interface to FORTRAN fixed boundary minimize. Propagates the  fixed-bondary DOFs.
        If None is given for any of the  arguments, the FORTRAN allocated variables will be used.

        Parameters
        ----------
        dt : float | None, optional
            Initial step-size of  the minimizer, by default None
        logfinaliter : bool, optional
            Whether or not to log and write the final state after the minimization, by default True
        maxIter : int | None, optional
            Maximum number of iterations used for this minimization call, by default None
        minimize_tol : float | None, optional
            Termination criterion on the norm of the forces X1,X2 and lambda, by default None
        """
        logger.debug("Minimizing MHD functional.")

        if dt is not None and not isinstance(dt, float):
            raise TypeError(f"Specified step-size  dt is not float or None but {type(dt)}")

        if not isinstance(logfinaliter, bool):
            raise TypeError(f"Specified logfinaliter is not bool but {type(logfinaliter)}")

        if maxIter is not None and not isinstance(maxIter, int):
            raise TypeError(f"Specified maxIter is not int or None but {type(maxIter)}")

        if minimize_tol is not None and not isinstance(minimize_tol, float):
            raise TypeError(
                f"Specified minimize_tol is not float or None but {type(minimize_tol)}"
            )

        with chdir(self.rundir):
            _state.minimize(dt, logfinaliter, maxIter, minimize_tol)

    @with_binding
    def save_statefile(self):
        logger.debug("Writing state to file.")
        _state.save_statefile()

    @with_binding
    def update_iota_for_current_constraint(self, iota_bspl: BSpline, reset_LA: bool):
        """update the iota and chi-prime B-spline in FORTRAN.

        Parameters
        ----------
        iota_bspl: BSpline
            Scipy basis spline object used for updating the iota profile
        reset_LA : bool
            Whether to reset the DOFs for lambda.
        """
        logger.debug("Updating iota")
        knots = iota_bspl.t
        coefs = iota_bspl.c

        if not isinstance(knots, np.ndarray) and not np.issubdtype(knots.dtype, np.floating):
            raise TypeError("Specified knots are not a numpy array of type float.")

        if knots.ndim != 1:
            raise ValueError(f"Specified knots.ndim is not 1 but {knots.ndim}.")

        if not isinstance(coefs, np.ndarray) and not np.issubdtype(coefs.dtype, np.floating):
            raise TypeError("Specified coefs are not a numpy array of type float.")

        if coefs.ndim != 1:
            raise ValueError(f"Specified coefs.ndim is not 1 but {coefs.ndim}.")

        if not isinstance(reset_LA, bool):
            raise TypeError(f"Specified reset_LA is not bool but {type(reset_LA)}")

        _state.update_iota_for_current_constraint(knots, coefs, reset_LA)

        # forget momentum and adapted time step after a iota update
        self.reset_descent()
        self._parameter_dict["iota"] = CaseInsensitiveDict(
            type="bspline",
            coefs=iota_bspl.c,
            knots=iota_bspl.t,
        )
        self.iota = iota_bspl

    @with_binding
    def reset_descent(self):
        logger.debug("Resetting accelerated descent")
        _state.reset_descent()

    def save(self, header: str = ""):
        """Write a parameter and statefile for the current state.
        The files are numbered according to the outputlevel (=number of restarts) and number of iterations used during this restart.

        Parameters
        ----------
        header : str, optional
            Header string for the parameter file, by default ""
        """
        logger.debug("Saving state")
        self.save_statefile()
        iters = self.get_iterations()
        outid = self.get_outputlevel()
        write_parameter_file_ini(
            flatten_parameters(self._parameter_dict),
            self._parameter_dict["ProjectName"] + f"_parameters_{outid:04}_{iters:08}.ini",
            header=header,
        )

    def _get_new_iota(self, iota_tol: float):
        """Fit an iota B-spline according to iota_geom and iota_curr for the current constraint mode

        Parameters
        ----------
        iota_tol : float
            Fit tolerance

        Returns
        -------
        tuple(scipy.bspl, float, xr.DataArray, float, float)
            (iota B-spline, root mean squared iota error, maximum difference to target Itor, root mean squared iota fit error)
        """
        sgrid = self.get_radial_gridpoints()
        sgrid[0] = 0.0
        sgrid[-1] = 1.0
        Qs = ["iota", "iota_geom", "iota_curr_bar", "I_tor"]
        while True:
            knots = _linspace_sgrid(sgrid, self.knots_per_sgrid_element)
            knots = _add_knot_multiplicity(knots, IOTA_BSPLINE_DEGREE)
            rho = _linspace_sgrid(sgrid, IOTA_SAMPLING_FACTOR * self.knots_per_sgrid_element)
            ev = self._evaluate_extrapolate_in_rho2(*Qs, rho=rho, theta="int", zeta="int")[Qs]
            iota_target = ev.iota_geom + self.I_tor(rho**2) * ev.iota_curr_bar

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
            rms_iota_fit = np.sqrt((delta_iota_fit**2).mean().item())
            if rms_iota_fit < iota_tol:
                break
            if self.knots_per_sgrid_element >= IOTA_KNOTS_PER_ELEMENT_MAX:
                logger.warning(
                    f"Could not fit iota profile up to target tolerance {iota_tol} "
                    f"with {self.knots_per_sgrid_element} knots per sgrid element. "
                    f"Current fit error: {rms_iota_fit}"
                )
                break
            if bspline.c.size * 2 > IOTA_COEFS_MAX:
                logger.warning(
                    f"Could not fit iota profile up to target tolerance {iota_tol} "
                    f"with {bspline.c.size} coefficients (max is {IOTA_COEFS_MAX}). "
                    f"Current fit error: {rms_iota_fit}"
                )
                break
            logger.debug(
                f"Increasing knots per sgrid element {self.knots_per_sgrid_element} "
                f"-> {self.knots_per_sgrid_element * 2} to fit iota profile with tolerance {iota_tol}. "
                f"Current fit error: {rms_iota_fit}"
            )
            self.knots_per_sgrid_element *= 2
        delta_iota = ev.iota - bspline(rho**2)

        # ignore near-axis region as we CURRENTLY have accuracy problems there!!!
        delta_iota = delta_iota.sel(rho=slice(IOTA_MIN_RHO, None))
        rms_iota = np.sqrt((delta_iota**2).mean().item())
        delta_Itor = ev.I_tor - self.I_tor(rho**2)
        delta_I_tor_max = np.abs(delta_Itor).max().item()
        return (bspline, rms_iota, delta_iota, delta_I_tor_max, rms_iota_fit)

    def _minimize_single_current_constraint(
        self,
        maxIter: int,
        maxIter_per_update: int,
        iota_tol: float,
        minimize_tol: float | None,
        maxUpdates: int,
        reset_LA: bool,
        save_intermediate: bool,
        save_header: str,
        **kwargs_minimize,
    ):
        """Runs a single current-minimization stage (e.g, one minimize tolerance and one iota tolerance)

        Parameters
        ----------
        maxIter : int
            Maximum number of fixed-boundary iterations
        maxIter_per_update : int
            Maximum number of fixed-boundary iterations per iota update
        iota_tol : float
            Tolerance on the iota-discrepancy.
        minimize_tol : float | None
            Tolerance on the norm of the forces for X1, X2 and lambda.
        maxUpdates : int
            Maximum number of iota updates (Picard iterations).
        reset_LA : bool
            If True, lambda is reset after an iota update.
        save_intermediate : bool
            If True, the state after each iota update is saved.
        save_header : str
            Header string for the saved parameter files.

        Returns
        -------
        tuple(float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int)
            (
            root mean squared iota error,
            maximum iota error,
            evolution of the rms_iota error over iota updates,
            evolution of the iota fit error over iota updates,
            number of used iota updates,
            )
        """
        self.rms_iota = np.inf
        self.delta_iota_max = np.nan
        delta_iota = None

        nth_update = 0
        iters_used_start = self.get_iterations()
        self.reset_descent()
        iters_used = 0
        if "logfinaliter" in kwargs_minimize:
            logfinaliter = kwargs_minimize["logfinaliter"]
            del kwargs_minimize["logfinaliter"]
        else:
            logfinaliter = False

        iota_bspl = self.iota

        while (self.rms_iota > iota_tol) and (iters_used < maxIter):
            if nth_update >= maxUpdates:
                logger.warning("Maximum number of iota updates reached!")
                break

            # only reset LA at the start
            if iters_used > 0:
                reset_LA = False

            if iota_bspl is not None:
                self.update_iota_for_current_constraint(iota_bspl, reset_LA)

            self._minimize(
                logfinaliter=logfinaliter,
                maxIter=min(maxIter - iters_used, maxIter_per_update),
                minimize_tol=minimize_tol,
                **kwargs_minimize,
            )

            (iota_bspl, self.rms_iota, delta_iota, delta_I_tor_max, rms_iota_fit) = (
                self._get_new_iota(iota_tol=iota_tol)
            )

            iters = self.get_iterations()
            iters_used = iters - iters_used_start

            nth_update += 1
            self.rms_iota_over_run.append(self.rms_iota)
            self.iters_per_run.append(iters)
            self.delta_I_tor_max_over_run.append(delta_I_tor_max)
            self.rms_iota_fit_over_run.append(rms_iota_fit)

            if save_intermediate:
                self.save(
                    header=save_header + f"!Iota update {nth_update + self.nth_update}\n"
                    f"!Created at {datetime.now().isoformat()}\n"
                    f"!pyGVEC v{__version__}\n",
                )

        if self.rms_iota > iota_tol and maxIter > 0:
            logger.warning(
                "Targeted iota has not been reached!\n"
                + f"iota_tol.: {iota_tol:.2e}, achieved rms Δiota.: {self.rms_iota:.2e}"
            )
        logger.debug(
            f"Terminated current minimization after {nth_update} iota updates. Δiota = {self.rms_iota:.3e}, target Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
        )
        if delta_iota is not None:
            self.delta_iota_max = np.abs(delta_iota).max().item()
        self.nth_update += nth_update

    def _minimize_current_constraint(
        self,
        iota_tol: float,
        minimize_tol: float,
        maxIter: int,
        maxUpdates: int,
        reset_LA: bool,
        target: Literal["iota", "iota_and_force", "auto"],
        save_intermediate: bool,
        save_header: str,
        **kwargs_minimize,
    ):
        """Fixed-boundary run with prescribed current profile.

        Parameters
        ----------
        iota_tol : float
            Tolerance on the iota-discrepancy.
        minimize_tol : float | None
            Tolerance on the norm of the forces for X1, X2 and lambda.
        minimize_tol : float | None
            Tolerance on the norm of the forces for X1, X2 and lambda.
        maxUpdates : int
            Maximum number of Picard iterations.
        reset_LA : bool
            If True, lambda is reset after an iota update.
        target : Literal["iota", "iota_and_force", "auto"]
            Optimization strategy for the Picard algorithm:
            "iota" = aggressively optimize for iota regardless of the force tolerance;
            "iota_and_force" = optimize for both iota and force tolerances;
            "auto" = start with target "iota" and then switch to target "iota_and_force";
        save_intermediate : bool
            If True, the state after a Picard iteration is saved.
        save_header : str
            Header string for the saved parameter files.
        kwargs_minimize : dict
            Additional keyword arguments passed to the energy minimization routine.
        """

        t0 = time.time()
        log_minimize_target = np.log10(minimize_tol)
        n_stages = max(int(max(-2 - log_minimize_target, log_minimize_target)), 1)
        minimize_tols = np.logspace(
            max(np.log10(AUTO_MINIMIZE_TOL_MIN), log_minimize_target),
            log_minimize_target,
            n_stages,
        )
        self.knots_per_sgrid_element = IOTA_KNOTS_PER_ELEMENT_MIN

        if n_stages > 1:
            iota_tols = np.ones(n_stages) * AUTO_PICARD_CURRENT_IOTA_TOL
            iota_tols[0] = AUTO_STAGE_IOTA_INIT_IOTA_TOL
            iota_tols[-1] = iota_tol
        else:
            iota_tols = [iota_tol]
        iota_tol = iota_tols[0]

        iters_used_start = self.get_iterations()

        self.nth_update = 0
        self.rms_iota_over_run = []
        self.iters_per_run = []
        self.delta_I_tor_max_over_run = []
        self.rms_iota_fit_over_run = []

        # target iota only
        if target == "iota" or target == "auto":
            self._minimize_single_current_constraint(
                maxIter=maxIter,
                maxIter_per_update=min(AUTO_STAGE_IOTA_INIT_MAXITER, maxIter),
                iota_tol=iota_tol,
                minimize_tol=minimize_tols[-1],
                maxUpdates=maxUpdates,
                reset_LA=reset_LA,
                save_intermediate=save_intermediate,
                save_header=save_header,
                **kwargs_minimize,
            )

        iters_used = self.get_iterations() - iters_used_start

        # target iota and forces
        if target == "iota_and_force" or target == "auto":
            if iters_used < maxIter:
                for minimize_tol, iota_tol in zip(minimize_tols, iota_tols):
                    logger.debug(f"minimize_tol:{minimize_tol}")
                    logger.debug(f"iota_tol:{minimize_tol}")
                    logger.debug(f"iters_used_start:{iters_used_start}")
                    logger.debug(f"iters_used:{iters_used}")
                    logger.debug(f"maxIter:{maxIter - iters_used}")
                    self.knots_per_sgrid_element = IOTA_KNOTS_PER_ELEMENT_MIN
                    self._minimize_single_current_constraint(
                        maxIter=maxIter - iters_used,
                        maxIter_per_update=maxIter - iters_used,
                        iota_tol=iota_tol,
                        minimize_tol=minimize_tol,
                        maxUpdates=maxUpdates,
                        reset_LA=reset_LA,
                        save_intermediate=save_intermediate,
                        save_header=save_header,
                        **kwargs_minimize,
                    )
                    iters_used = self.get_iterations() - iters_used_start
        t1 = time.time()

        logger.info(f"Terminated current minimization after {t1 - t0:.2}s.")
        logger.info(f"  used {iters_used} / {maxIter} iterations,")
        logger.info(f"  max Δiota: {self.delta_iota_max:.2e},")
        logger.info(f"  rms Δiota: {self.rms_iota:.2e},")
        logger.info(f"  iota_tol: {iota_tol:.2e}")

        if maxIter > 0:
            logger.info(f"  max ΔItor: {self.delta_I_tor_max_over_run[-1]:.2e}")
            logger.debug(
                f"  rms Δiota_fit: {self.rms_iota_fit_over_run[-1]:.2e}, "
                f"  with {self.knots_per_sgrid_element} knots per sgrid element"
            )

    def minimize(
        self,
        iota_tol: float,
        maxUpdates: int,
        picard_target: Literal["iota", "iota_and_force", "auto"],
        I_tor: Callable | None = None,
        reset_LA: bool = False,
        save_intermediate: bool = False,
        save_header: str = "",
        **kwargs,
    ):
        """Propagate the current states fixed-boundary DOFs.

        Parameters
        ----------
        iota_tol : float
            Tolerance on the iota-discrepancy. Only relevant during current optimization.
        maxUpdates : int
            Maximum number of Picard iterations. Only relevant during current optimization
        target : Literal["iota", "iota_and_force", "auto"]
            Optimization strategy for the Picard algorithm:
            "iota" = aggressively optimize for iota regardless of the force tolerance;
            "iota_and_force" = optimize for both iota and force tolerances;
            "auto" = start with target "iota" and then switch to target "iota_and_force";
        I_tor : Callable | None, optional
            Function to evaluate Itor on rho**2. If None, a fixed iota run is performed.
        reset_LA : bool, optional
            If True lambda is reset, by default False
        save_intermediate : bool
            If True, the state after a Picard iteration is saved.
        save_header : str
            Header string for the saved parameter files.

        Raises
        ------
        ValueError
            If iota_target is unknown.
        """
        if I_tor is not None:
            self.current_constraint = True
            self.I_tor = I_tor
        else:
            self.current_constraint = False
            self.I_tor = None

        if picard_target not in ["iota", "iota_and_force", "auto"]:
            raise ValueError(
                f"""Unknown picard_current["target"]: {picard_target}! target has to be "iota","iota_and_force" or "auto"."""
            )

        if self.current_constraint:
            logger.debug("minimize with current optimization")
            self._minimize_current_constraint(
                iota_tol=iota_tol,
                maxUpdates=maxUpdates,
                reset_LA=reset_LA,
                target=picard_target,
                save_intermediate=save_intermediate,
                save_header=save_header,
                **kwargs,
            )
        else:
            logger.debug("minimize with fixed iota")
            self._minimize(**kwargs)

    # ==== Free Boundary ==== #

    @with_binding
    def set_Bext_squared(self, B2: np.ndarray):
        """Update the FORTRAN exterior field array

        Parameters
        ----------
        B2 : np.nd.array
            Squared exterior field  used for free-boundary
        """
        logger.debug("Updating exterior field")

        if not isinstance(B2, np.ndarray) and not np.issubdtype(B2.dtype, np.floating):
            raise TypeError("Specified B2 is not a numpy array of type float.")

        if B2.ndim != 1:
            raise ValueError(f"Specified B2.ndim is not 1 but {B2.ndim}.")

        _state.set_bext_squared(B2)

    @with_binding
    def _get_boundary_forces(self):
        """Extract the forces on the boundary DOFs.

        Returns
        -------
        tuple(np.ndarray, np.ndarray, np.ndarray)
            (poloidal mode numbers, toroidal mode numbers, forces)
        """
        n_X1 = self._get_X1_n_modes()
        n_X2 = self._get_X2_n_modes()

        forces_out = np.asfortranarray(np.zeros((n_X1 + n_X2)))
        m = np.asfortranarray(np.zeros((n_X1 + n_X2), dtype=np.int32))
        n = np.asfortranarray(np.zeros_like(m, dtype=np.int32))

        _state.get_boundary_forces(m, n, forces_out)

        return (m, n, forces_out)

    @with_binding
    def get_Fnorm_boundary(self):
        logger.debug("Extracting Fnorm at the boundary")
        return _state.get_fnorm_boundary()

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
    def update_boundary(self, dt: float, reset_minimizer: bool, blend_depth: float = 0.6):
        logger.debug("Stepping in free boundary direction")
        dt = _state.update_boundary(dt, reset_minimizer, blend_depth=blend_depth)
        return dt

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
            .squeeze()
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
    def _free_boundary_update(
        self,
        maxIter: int,
        B2: np.ndarray,
        dt: float,
        reset_minimizer: bool = False,
        blend_depth: float = 0.6,
        **minimize_kwargs,
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
        reset_minimizer: bool, optional
            Whether to reset the free boundary minimizers to their initial states.
        blend_depth: float
            How far in rho boundary perturbations are blend into the interior. The default  is 0.6.

        Returns
        -------
        float
            New stepsize
        """
        # set Bext
        self.set_Bext_squared(B2)
        self.reset_descent()
        dt = self.update_boundary(dt, reset_minimizer, blend_depth)

        self.minimize(logfinaliter=False, maxIter=maxIter, **minimize_kwargs)
        return dt

    def _evaluate_extrapolate_in_rho2(self, *Qs, rho, theta, zeta):
        rho_aux = np.sqrt([1e-8, 2e-8, 3e-8])
        ev_aux = self.evaluate(*Qs, rho=rho_aux, theta=theta, zeta=zeta, extrapolate=False)
        r1 = ev_aux.isel(rad=0)
        r2 = ev_aux.isel(rad=1)
        r3 = ev_aux.isel(rad=2)
        on_axis = 3 * (r1 - r2) + r3
        on_axis.rho[...] = 0.0

        ev = self.evaluate(*Qs, rho=rho, theta=theta, zeta=zeta, extrapolate=False)
        for q in ev.data_vars:
            if "rad" in ev[q].dims:
                ev[q].loc[dict(rho=0.0)] = on_axis[q]
        return ev

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
