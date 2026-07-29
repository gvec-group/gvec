# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._minimize - low-level python API for moving a gvec.State via some minimization"""

# === Imports === #
import numpy as np
from gvec.util import write_parameters, CaseInsensitiveDict, chdir
from ._binding import with_binding
from typing import Literal, Callable
from gvec.core.compute import fluxsurface_integral

from gvec.lib import _state

import time
import logging
from scipy.constants import mu_0

logger = logging.getLogger("gvec.run")

DEFAULT_N_IOTA_RHOS = 101


@with_binding
def _minimize(self, dt=None, logfinaliter=True, maxIter=None, minimize_tol=None):
    """Minimize the MHD functional for the current state."""
    logger.debug("Minimizing MHD functional.")
    with chdir(self.rundir):
        _state.minimize(dt, logfinaliter, maxIter, minimize_tol)


@with_binding
def save_statefile(self):
    logger.debug("Writing state to file")
    _state.save_statefile()


@with_binding
def update_iota_for_current_constraint(self, rho2, iota, reset_LA):
    logger.debug("Updating iota")
    _state.update_iota_for_current_constraint(rho2, iota, reset_LA)


@with_binding
def reset_descent(self):
    logger.debug("Resetting accelerated descent")
    _state.reset_descent()


def save(self):
    self.save_statefile()
    iters = self.get_iterations()
    outid = self.get_outputlevel()
    write_parameters(
        self._parameter_dict,
        self._parameter_dict["ProjectName"] + f"_parameters_{outid:04}_{iters:08}.ini",
    )


def _get_new_iota(self):
    """Calculate a new iota profile. Extracted from gvec-run, to be used for current constraint

    Returns
    -------
    (xr.DataArray, xr.DataArray, xr.DataArray)
        new iota values, rho^2 interpolation points, root mean square iota error
    """
    rho_eval = np.concatenate(
        [[np.sqrt(1e-8), np.sqrt(2e-8), np.sqrt(3e-8)], np.sqrt(self.rho2_iota[1:])]
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
    rms_iota = np.sqrt((delta_iota**2).mean(dim="rad")).item()
    return (iota_values, ev.rho**2, rms_iota, delta_iota)


def _minimize_single_current_constraint(
    self,
    maxIter: int,
    maxIter_per_update: int,
    iota_tol: float,
    minimize_tol: float | None,
    maxUpdates: int,
    reset_LA: bool,
    **kwargs_minimize,
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
        Maximum number of iota updates

    Returns
    -------
    float
        RMS of the iota update.
    """
    rms_iota = np.inf
    nth_update = -1
    iters_used_start = self.get_iterations()
    self.reset_descent()
    iters_used = 0
    delta_iota = np.inf
    delta_iota_over_run = []
    iters_per_run = []
    if "logfinaliter" in kwargs_minimize:
        logfinaliter = kwargs_minimize["logfinaliter"]
        del kwargs_minimize["logfinaliter"]
    else:
        logfinaliter = False

    while (rms_iota > iota_tol) and (iters_used < maxIter):
        if nth_update + 1 > maxUpdates:
            logger.warning("Maximum number of iota updates reached reached!")
            break
        self._minimize(
            logfinaliter=logfinaliter,
            maxIter=min(maxIter - iters_used, maxIter_per_update),
            minimize_tol=minimize_tol,
            **kwargs_minimize,
        )
        new_iota, new_rho2, rms_iota, delta_iota = _get_new_iota(self)

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

        iters = self.get_iterations()
        iters_used = iters - iters_used_start

        nth_update += 1
        delta_iota_over_run.append(rms_iota)
        iters_per_run.append(iters)

    if rms_iota > iota_tol and maxIter > 0:
        logger.warning(
            "Targeted iota has not been reached!\n"
            + f"iota_tol.: {iota_tol:.2e}, achieved rms Δiota.: {rms_iota:.2e}"
        )
    logger.debug(
        f"Terminated current minimization after {nth_update} updates. Δiota = {rms_iota:.3e}, target Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
    )
    return (
        rms_iota,
        delta_iota,
        np.array(delta_iota_over_run),
        np.array(iters_per_run),
        nth_update,
    )


def _minimize_current_constraint(
    self,
    iota_target: float,
    minimize_tol: float,
    maxIter: int,
    maxUpdates: int,
    reset_LA: bool,
    target: Literal["iota", "iota_and_force", "auto"],
    **kwargs,
):
    """Fixed-boundary run with prescribed current profile.

    Parameters
    ----------
    iota_target : float
        Tolerance on the RMS of the iota update.
    minimize_tol : float
        Tolerance for the norm of the forces.
    maxIter : int
        Maximum number of iterations for the fixed boundary run.
    maxUpdates : int, optional
        Maximum number of iota-updates per stage
    """

    t0 = time.time()
    self.iota = self._get_new_iota()[0]
    log_minimize_target = np.log10(minimize_tol)
    log_iota_target = np.log10(iota_target)
    n_stages = max(int(max(-2 - log_minimize_target, log_minimize_target)), 1)
    minimize_tols = np.logspace(max(-4, log_minimize_target), log_minimize_target, n_stages)
    if n_stages > 1:
        iota_tols = np.logspace(max(-3, log_iota_target), log_iota_target, n_stages)
    else:
        iota_tols = [1e-10]

    iota_tol = iota_tols[0]

    iters_used_start = self.get_iterations()

    self.nth_update = 0
    # target iota only
    if target == "iota" or target == "auto":
        (
            self.rms_iota,
            self.delta_iota,
            self.delta_iota_over_run,
            self.iters_per_run,
            nth_update,
        ) = self._minimize_single_current_constraint(
            maxIter=maxIter,
            maxIter_per_update=min(10, maxIter),
            iota_tol=iota_tol,
            minimize_tol=minimize_tols[0],
            maxUpdates=maxUpdates,
            reset_LA=reset_LA,
        )
        self.nth_update += nth_update
    else:
        self.delta_iota_over_run = np.array([])
        self.iters_per_run = np.array([])
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
                (
                    self.rms_iota,
                    self.delta_iota,
                    delta_iota_over_run,
                    iters_per_run,
                    nth_update,
                ) = self._minimize_single_current_constraint(
                    maxIter=maxIter - iters_used,
                    maxIter_per_update=1e12,
                    iota_tol=iota_tol,
                    minimize_tol=minimize_tol,
                    maxUpdates=maxUpdates,
                    reset_LA=reset_LA,
                    **kwargs,
                )
                self.nth_update += nth_update
                self.delta_iota_over_run = np.concat(
                    [self.delta_iota_over_run, delta_iota_over_run]
                )
                self.iters_per_run = np.concat([self.iters_per_run, iters_per_run])

                iters_used = self.get_iterations() - iters_used_start
    t1 = time.time()

    logger.info(
        f"Terminated current minimization after {t1 - t0:.2}s with Δiota = {self.rms_iota:.3e}, Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
    )
    self.save()


def minimize(
    self,
    iota_target: float,
    maxUpdates: int,
    target: Literal["iota", "iota_and_force", "auto"],
    I_tor: float | Callable | None = None,
    n_iota_rhos=DEFAULT_N_IOTA_RHOS,
    reset_LA: bool = False,
    **kwargs,
):
    if I_tor is not None:
        self.current_constraint = True
        self.rho2_iota = np.linspace(0, 1, n_iota_rhos) ** 2
        if callable(I_tor):
            self.I_tor = I_tor(self.rho2_iota)
        else:
            self.I_tor = I_tor * np.ones_like(self.rho2_iota)
    else:
        self.current_constraint = False
        self.rho2_iota = None
        self.I_tor = None

    if target not in ["iota", "iota_and_force", "auto"]:
        raise ValueError(
            f"""Unknown picard_current["target"]: {target}! target has to be "iota","iota_and_force" or "auto"."""
        )

    if self.current_constraint:
        logger.debug("minimize with current optimization")
        self._minimize_current_constraint(
            iota_target=iota_target,
            maxUpdates=maxUpdates,
            reset_LA=reset_LA,
            target=target,
            **kwargs,
        )
    else:
        logger.debug("minimize with fixed iota")
        self._minimize(**kwargs)


# ==== Free Boundary ==== #


@with_binding
def set_Bext_squared(self, B2):
    logger.debug("Updating exterior field")
    _state.set_bext_squared(B2)


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
    dt = self.update_boundary(dt, reset_minimizer, blend_depth)

    self.minimize(logfinaliter=False, maxIter=maxIter, **minimize_kwargs)
    return dt
