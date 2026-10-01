# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.minimize - contains the MinimizeState class"""

# === Imports === #

import time
import logging
from datetime import datetime
from typing import Literal, Callable
from pathlib import Path
import warnings
from collections.abc import Mapping

import numpy as np
import scipy
from scipy.interpolate import BSpline
from scipy.constants import mu_0
from numpy.polynomial import Polynomial

from gvec import __version__
from gvec.util import (
    write_parameter_file_ini,
    chdir,
    flatten_parameters,
    profile_to_parameters,
    profile_from_parameters,
)
from gvec.core.compute import fluxsurface_integral
from gvec.core.state import State
from gvec.core.state._binding import with_binding
from gvec.config import ConfigABC
from gvec.errors import InvalidParameterError
from gvec.lib import _state, c_rProfile, t_rProfile_poly, t_rProfile_bspl


# === Globals === #


logger = logging.getLogger("gvec.state")

AUTO_STAGE_IOTA_INIT_MAXITER = 10
AUTO_STAGE_IOTA_INIT_IOTA_TOL = 1e-3  # iota tolerance during the initial stage
AUTO_PICARD_CURRENT_IOTA_TOL = 1e-5  # iota tolerance during the intermediate stage
AUTO_MINIMIZE_TOL_MIN = 1e-3
IOTA_BSPLINE_DEGREE = 5
IOTA_SAMPLING_FACTOR = 4
IOTA_KNOTS_PER_ELEMENT_MIN = 2
IOTA_KNOTS_PER_ELEMENT_MAX = 256
IOTA_COEFS_MAX = 768  # at 20 characters per coefficient: 15360 characters + overhead < 2^14
IOTA_MIN_RHO = 1e-2  # ignore near-axis region for Δiota & Δfit, TODO: fix accuracy problems
IOTA_MAXITER_PER_UPDATE = 1e12  # maximum number of fixed-boundary iterations per iota update


# === Config class for picard_current === #


class PicardCurrentConfig(ConfigABC):
    aliases = {"maxRestarts": "max_updates"}

    def __init__(
        self,
        target: Literal["iota", "iota_and_force", "auto"] = "auto",
        iota_tol: float = 1e-8,  # least-squares fit: limited to single precision!
        max_updates: int = 50,
        reset_LA: bool = False,
    ):
        self.target = target
        self.iota_tol = iota_tol
        self.max_updates = max_updates
        self.reset_LA = reset_LA

    @classmethod
    def from_dict(cls, config: Mapping | Literal["auto", "off"] | None):
        if config == "off":
            warnings.warn(
                "picard_current='off' is deprecated. Specify I_tor=None (default) to use a fixed iota profile.",
                DeprecationWarning,
            )
            return cls()
        elif config == "auto":
            warnings.warn(
                "picard_current='auto' is deprecated. Leave picard_current=None for default behavior.",
                DeprecationWarning,
            )
            return cls()
        return super().from_dict(config)


# === MinimizationState class === #


class MinimizationState(State):
    def __init__(
        self,
        parameterfile,
        statefile=None,
        stdout_path=None,
        redirect_stdout=True,
        I_tor: Mapping | None = None,
        picard_current: Mapping | Literal["auto", "off"] | None = None,
    ):
        if picard_current == "off" and I_tor is not None:
            raise InvalidParameterError(
                "picard_current='off' is deprecated and incompatible with I_tor. Specify I_tor=None (default) to use a fixed iota profile."
            )
        elif picard_current == "auto" and I_tor is None:
            raise InvalidParameterError(
                "picard_current='auto' is deprecated and incompatible with I_tor=None. Specify a current profile I_tor to use picard_current='auto'."
            )

        super().__init__(parameterfile, statefile, stdout_path, redirect_stdout)
        r_GP, r_w, t_n, t_w, z_n, z_w = self.get_integration_points()

        self.t_n = t_n
        self.z_n = z_n

        if I_tor is not None:
            self.I_tor = profile_from_parameters(I_tor)
            if abs(value := self.I_tor(0.0)) > 1e-8:
                raise InvalidParameterError(
                    f"Toroidal current profile is not zero at magnetic axis! I_tor(rho=0): {value}"
                )
        else:
            self.I_tor = None
        self.picard_current = PicardCurrentConfig.from_dict(picard_current)
        self.iters_used = 0
        self.rms_iota = None
        self.delta_iota = None
        self.dt_fixed_initial = self.dt_fixed
        self._iota: BSpline | Polynomial = profile_from_parameters(self.parameters["iota"])

    def bind(self, force: bool = False):
        super().bind(force=force, init_for_minimize=True)

    @classmethod
    def new(cls, parameters: Mapping, rundir: str | Path | None = None, **kwargs):
        picard_current = parameters.get("picard_current", None)
        I_tor = parameters.get("I_tor", None)
        return super().new(
            parameters, rundir=rundir, picard_current=picard_current, I_tor=I_tor, **kwargs
        )

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

    @property
    def iota(self) -> BSpline | Polynomial:
        """The rotational transform profile as a function of rho^2."""
        return self._iota

    @iota.setter
    @with_binding
    def iota(self, value: BSpline | Polynomial):
        self.logger.debug("Setting new iota profile")
        match value:
            case BSpline():
                fobj = t_rProfile_bspl(value.t, value.c)
            case Polynomial():
                fobj = t_rProfile_poly(value.coef)
        self.parameters["iota"] = profile_to_parameters(value)
        self._iota = value
        _state.set_iota(fobj)
        del fobj  # deallocate
        self.reset_descent()

    @with_binding
    def _get_iota_fobj(self) -> c_rProfile:
        """Get the fortran object for the iota profile."""
        raise NotImplementedError(
            "_get_iota_fobj causes issues with some compilers and might be revisited when it is actually needed."
        )
        return _state.get_iota()

    @with_binding
    def reset_descent(self):
        logger.debug("Resetting accelerated descent")
        _state.reset_descent()
        self.dt_fixed = self.dt_fixed_initial

    @property
    @with_binding
    def dt_fixed(self):
        "The current time-step for the fixed-boundary minimizer"
        logger.debug("Extracting fixed boundary time-step...")
        dt = _state.get_dt()
        logger.debug(f"... dt = {dt:.2e}")
        return dt

    @dt_fixed.setter
    @with_binding
    def dt_fixed(self, dt):
        "Set the time-step for the fixed-boundary minimizer"
        if not isinstance(dt, float):
            raise TypeError(f"Expected type(dt) float but got {type(dt)} in _set_dt()")
        logger.debug(f"Setting fixed boundary time-step dt={dt:.2e}.")
        return _state.set_dt(dt)

    @property
    @with_binding
    def dt_free(self):
        "The current time-step for the free-boundary minimizer"
        logger.debug("Extracting free boundary time-step...")
        dt = _state.get_free_boundary_dt()
        logger.debug(f"... dt = {dt:.2e}")
        return dt

    @dt_free.setter
    @with_binding
    def dt_free(self, dt):
        "Set the time-step for the free-boundary minimizer"
        if not isinstance(dt, float):
            raise TypeError(
                f"Expected type(dt) float but got {type(dt)} in _set_free_boundary_dt()"
            )
        logger.debug(f"Setting free boundary time-step dt={dt:.2e}.")
        return _state.set_free_boundary_dt(dt)

    @property
    @with_binding
    def max_dt_free(self):
        "The maximum time-step for the free-boundary minimizer"
        logger.debug("Extracting free boundary max time-step...")
        dt = _state.get_free_boundary_max_dt()
        logger.debug(f"... max_dt = {dt:.2e}")
        return dt

    @max_dt_free.setter
    @with_binding
    def max_dt_free(self, dt):
        "Set the maximum time-step for the free-boundary minimizer"
        if not isinstance(dt, float):
            raise TypeError(
                f"Expected type(dt) float but got {type(dt)} in _set_free_boundary_dt()"
            )
        logger.debug(f"Setting free boundary max time-step max_dt={dt:.2e}.")
        return _state.set_free_boundary_max_dt(dt)

    @with_binding
    def reset_LA(self):
        logger.debug("Resetting lambda DoFs")
        _state.reset_LA()

    def save(self, header: str = ""):
        """Write a parameter and statefile for the current state.
        The files are numbered according to the outputlevel (=number of restarts) and number of iterations used during this restart.

        Parameters
        ----------
        header : str, optional
            Header string for the parameter file, by default ""
        """
        logger.debug("Saving state")
        iters = self.get_iterations()
        outid = self.get_outputlevel()

        self.save_statefile()
        # name hardcoded in pygvec/state.F90:save_statefile & state/restart.f90:WriteStateToASCII
        self.statefile = (
            Path(self.rundir)
            / f"{self.parameters['ProjectName']}_State_{outid:04}_{iters:08}.dat"
        )

        path = (
            Path(self.rundir)
            / f"{self.parameters['ProjectName']}_parameters_{outid:04}_{iters:08}.ini"
        )
        write_parameter_file_ini(
            flatten_parameters(self.parameters),
            path,
            header=header,
        )
        self.parameterfile = path

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

    def _optimize_current(
        self,
        maxIter: int,
        maxIter_per_update: int,
        iota_tol: float,
        minimize_tol: float | None,
        maxUpdates: int,
        reset_LA: bool,
        save_intermediate: bool,
        save_header: str,
        logfinaliter: bool = False,
        dt: float | None = None,
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
        logfinaliter : bool
            Whether or not to log and write the final state after the minimization, by default False. Passed to _minimize()
        dt: float | None
            Initial step-size of the minimizer, optional. Passed to _minimize()

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

        iota_bspl = self.iota

        while (self.rms_iota > iota_tol) and (iters_used < maxIter):
            if nth_update >= maxUpdates:
                logger.warning("Maximum number of iota updates reached!")
                break

            self.iota = iota_bspl

            # only reset LA at the start
            if reset_LA and iters_used == 0:
                self._state.reset_LA()

            self._minimize(
                dt=dt,
                logfinaliter=logfinaliter,
                maxIter=min(maxIter - iters_used, maxIter_per_update),
                minimize_tol=minimize_tol,
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

    def _solve_with_current(
        self,
        save_intermediate: bool,
        save_header: str,
        minimize_tol: float | None = None,
        maxIter: int | None = None,
        **kwargs_minimize,
    ):
        """Fixed-boundary run with prescribed current profile.

        Parameters
        ----------
        save_intermediate : bool
            If True, the state after a Picard iteration is saved.
        save_header : str
            Header string for the saved parameter files.
        maxIter : int | None, optional
            Maximum number of iterations used for this minimization call, None (default) reads the value from the parameter dictionary.
        minimize_tol : float | None, optional
            Tolerance for the energy minimization, None (default) reads the value from the parameter dictionary.
        kwargs_minimize : dict
            Additional keyword arguments passed to the energy minimization routine.
        """
        if maxIter is None:
            maxIter = self.parameters["maxIter"]
        if minimize_tol is None:
            minimize_tol = self.parameters["minimize_tol"]

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
            iota_tols[-1] = self.picard_current.iota_tol
        else:
            iota_tols = [self.picard_current.iota_tol]
        iota_tol = AUTO_STAGE_IOTA_INIT_IOTA_TOL

        iters_used_start = self.get_iterations()

        self.nth_update = 0
        self.rms_iota_over_run = []
        self.iters_per_run = []
        self.delta_I_tor_max_over_run = []
        self.rms_iota_fit_over_run = []

        # target iota only
        if self.picard_current.target == "iota" or self.picard_current.target == "auto":
            self._optimize_current(
                maxIter=maxIter,
                maxIter_per_update=min(AUTO_STAGE_IOTA_INIT_MAXITER, maxIter),
                iota_tol=iota_tol,
                minimize_tol=minimize_tols[0],
                maxUpdates=self.picard_current.max_updates,
                reset_LA=self.picard_current.reset_LA,
                save_intermediate=save_intermediate,
                save_header=save_header,
                **kwargs_minimize,
            )

        iters_used = self.get_iterations() - iters_used_start

        # target iota and forces
        logger.debug(f"minmize tols {minimize_tols}")
        logger.debug(f"iota tols {iota_tols}")
        if (
            self.picard_current.target == "iota_and_force"
            or self.picard_current.target == "auto"
        ):
            if iters_used < maxIter:
                for minimize_tol, iota_tol in zip(minimize_tols, iota_tols):
                    logger.debug(f"minimize_tol:{minimize_tol}")
                    logger.debug(f"iota_tol:{minimize_tol}")
                    logger.debug(f"iters_used_start:{iters_used_start}")
                    logger.debug(f"iters_used:{iters_used}")
                    logger.debug(f"maxIter:{maxIter - iters_used}")
                    self.knots_per_sgrid_element = IOTA_KNOTS_PER_ELEMENT_MIN
                    self._optimize_current(
                        maxIter=maxIter - iters_used,
                        maxIter_per_update=maxIter - iters_used,
                        iota_tol=iota_tol,
                        minimize_tol=minimize_tol,
                        maxUpdates=self.picard_current.max_updates,
                        reset_LA=self.picard_current.reset_LA,
                        save_intermediate=save_intermediate,
                        save_header=save_header,
                        **kwargs_minimize,
                    )
                    iters_used = self.get_iterations() - iters_used_start
        t1 = time.time()

        if maxIter > 0:
            logger.info(f"Terminated current minimization after {t1 - t0:.2}s.")
            logger.info(f"  used {iters_used} / {maxIter} iterations,")
            logger.info(f"  max Δiota: {self.delta_iota_max:.2e},")
            logger.info(f"  rms Δiota: {self.rms_iota:.2e},")
            logger.info(f"  iota_tol: {iota_tol:.2e}")

            logger.info(f"  max ΔItor: {self.delta_I_tor_max_over_run[-1]:.2e}")
            logger.debug(
                f"  rms Δiota_fit: {self.rms_iota_fit_over_run[-1]:.2e}, "
                f"  with {self.knots_per_sgrid_element} knots per sgrid element"
            )

    def solve(
        self,
        save_intermediate: bool = False,
        save_header: str = "",
        **kwargs,
    ):
        """Propagate the current states fixed-boundary DOFs.

        Parameters
        ----------
        save_intermediate : bool
            If True, the state after a Picard iteration is saved.
        save_header : str
            Header string for the saved parameter files.

        Raises
        ------
        ValueError
            If iota_target is unknown.
        """
        if self.I_tor is not None:
            logger.debug("minimize with current optimization")
            self._solve_with_current(
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
        _state.update_boundary(dt, reset_minimizer, blend_depth=blend_depth)

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

        self.solve(logfinaliter=False, maxIter=maxIter, **minimize_kwargs)
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
