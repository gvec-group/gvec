import numpy as np
from gvec import State, fluxsurface_integral
from gvec.coils import CoilSet, Coil
from scipy.constants import mu_0
from tqdm import tqdm
from typing import Callable
import shutil
import time
import logging

try:
    import fbtools

    has_fbtools = True
except ImportError:
    has_fbtools = False

logger = logging.getLogger("gvec.free-boundary")


def _get_B2_one_fp(
    state: State,
    coil_set: CoilSet | Coil,
    ntheta: int,
    nzeta: int,
    I_tor: float,
    diagnostic_output=False,
    fac: int = 1,
):
    """Setup BIEST with points corresponding to one field period and return the flattened squared exterior field for free boundary gvec.
    One field period is used to avoid automatic copying of from inside beast which may switch directions.

    Parameters
    ----------
    state : State
        GVEC state object. This is used to extract the boundary integration points
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
        Squared exterior field at the GVEC integration points flattened in in GVEC order.
    """

    if has_fbtools:
        from fbtools.exterior_vacuum_field import get_B_ext_plasma
    else:
        logger.warning(
            "Could not import 'fbtools'! Please make sure the package is installed. Aborting Bext calulation."
        )
        return

    ntheta_aux = fac * ntheta
    nzeta_aux = fac * nzeta

    theta = np.linspace(0, 2 * np.pi, ntheta_aux, endpoint=False)
    theta += 0.5 * theta[1]
    N_FP = state.evaluate("N_FP").N_FP
    zeta = np.linspace(0, 2 * np.pi, nzeta_aux * N_FP.item(), endpoint=False)
    zeta += 0.5 * zeta[1]

    boundary = state.evaluate("pos", "N_FP", rho=1.0, theta=theta, zeta=zeta)
    boundary = boundary.sel(rho=1).transpose("xyz", "tor", "pol")

    ev_coils = coil_set.eval_B(boundary.pos)
    B_coils = ev_coils.B.transpose("xyz", "tor", "pol").values
    time_0 = time.time()
    B_plas_ext = get_B_ext_plasma(
        boundary.pos.values,
        1,  # boundary.N_FP,
        B_coils,
        I_tor * mu_0,
        digits=10,
    )
    time_1 = time.time()
    logger.debug(
        f"time for beast:{time_1 - time_0} for {nzeta_aux * nzeta_aux * N_FP.item()} points"
    )
    B_external = B_plas_ext + B_coils
    B2 = np.sum(B_external**2, axis=0)

    B2 = B2[:nzeta_aux:fac, ::fac]
    B2_flat = np.ones(ntheta * nzeta) * 1e6
    counter = 0

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


def _get_new_iota(state: State, I_tor_target: np.ndarray, rho: np.ndarray):
    """Calculate a new iota profile. Extracted from gvec-run, to be used for free-boundary current constraint

    Parameters
    ----------
    state : State
        Current gvec state object
    I_tor_target : np.ndarray
        Prescribed current profile at equally spaced radial rho points.

    Returns
    -------
    (xr.DataArray, xr.DataArray, xr.DataArray)
        new iota values, rho^2 interpolation points, root mean square iota error
    """
    rho_eval = np.concatenate([[np.sqrt(1e-8), np.sqrt(2e-8), np.sqrt(3e-8)], rho[1:]])
    quantites = ["iota", "iota_curr_0", "iota_0", "I_tor"]
    ev = state.evaluate(*quantites, rho=rho_eval)
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

    iota_values = ev.iota_0 + I_tor_target * ev.iota_curr_0
    delta_iota = ev.iota - iota_values
    rms_iota = np.sqrt((delta_iota**2).mean("rad"))
    return (iota_values, ev.rho**2, rms_iota)


def _free_boundary_update(
    state: State,
    maxIter: int,
    B2: np.ndarray,
    dt: float = 0.001,
    current_constraint: bool = False,
    I_tor: np.ndarray | None = None,
    rho: np.ndarray | None = None,
    copy_stdout: bool = True,
):
    # set Bext
    state.set_Bext_squared(B2)
    state.reset_descent()
    dt = state.update_boundary(dt)

    # get minimizer tolerance
    tol = state.get_minimize_tol()
    logger.debug("minimize tol:", tol)
    logger.debug("current constraint", current_constraint)
    if current_constraint:
        state = _minimize_current_constraint(
            state=state,
            iota_target=1e-10,
            minimize_target=tol,
            maxIter=maxIter,
            dt=dt,
            I_tor=I_tor,
            rho=rho,
        )
    else:
        # fixed iota
        state._minimize(logfinaliter=False, maxIter=maxIter)
    if copy_stdout:
        shutil.copyfile(state._stdout.name, "stdout.txt")
    return state, dt


def _minimize_current_constraint(
    state: State,
    iota_target: float,
    minimize_target: float,
    maxIter: int,
    dt: float,
    I_tor: np.ndarray,
    rho: np.ndarray,
):
    log_minimize_target = np.log10(minimize_target)
    log_iota_target = np.log10(iota_target)
    n_stages = max(int(max(-2 - log_minimize_target, log_minimize_target)), 1)
    minimize_tols = np.logspace(max(-3, log_minimize_target), log_minimize_target, n_stages)
    iota_tols = np.logspace(max(-3, log_iota_target), log_iota_target, n_stages)

    iters_used_start = state.get_iterations()
    # target iota only
    logger.debug(f"maxIter:{maxIter}")
    state, rms_iota = _minimize_single_current_constraint(
        state=state,
        iota_tol=iota_tols[0],
        maxIter=maxIter,
        I_tor=I_tor,
        rho=rho,
        minimize_tol=minimize_tols[0],
        maxIter_per_update=10,
        dt=dt,
    )

    iters_used = state.get_iterations() - iters_used_start

    for minimize_tol, iota_tol in zip(minimize_tols, iota_tols):
        logger.debug(f"iters_used_start:{iters_used_start}")
        logger.debug(f"iters_used:{iters_used}")
        logger.debug(f"maxIter:{maxIter - iters_used}")
        state, rms_iota = _minimize_single_current_constraint(
            state=state,
            iota_tol=iota_tol,
            maxIter=maxIter - iters_used,
            I_tor=I_tor,
            rho=rho,
            minimize_tol=minimize_tol,
            dt=dt,
        )
        iters_used = state.get_iterations() - iters_used_start
    logger.info(
        f"Terminated current minimization with Δiota = {rms_iota:.3e}, Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
    )
    return state


def _minimize_single_current_constraint(
    state: State,
    maxIter: int,
    I_tor: np.ndarray,
    rho: np.ndarray,
    iota_tol: float = 1e-10,
    minimize_tol: float | None = None,
    maxIter_per_update: int = int(
        1e12
    ),  # set this to a reasonable value to target iota rather then iota and forces
    maxUpdates: int = 30,
    dt: float | None = None,
):
    rms_iota = 1e6
    nth_update = -1
    iters_used_start = state.get_iterations()
    state.reset_descent()
    iters_used = 0
    while (rms_iota > iota_tol) and (iters_used < maxIter):
        if nth_update + 1 > maxUpdates:
            logger.warning("Maximum number of iota updates reached reached!")
            break
        state._minimize(
            logfinaliter=False,
            maxIter=min(maxIter - iters_used, maxIter_per_update),
            minimize_tol=minimize_tol,
        )
        new_iota, new_rho2, rms_iota = _get_new_iota(state, I_tor, rho)
        state.update_iota_for_current_constraint(new_rho2, new_iota)
        iters_used = state.get_iterations() - iters_used_start
        nth_update += 1

    if rms_iota > iota_tol:
        logger.warning(
            "Targeted iota has not been reached!\n"
            + f"iota_tol.: {iota_tol:.2e}, achieved rms Δiota.: {rms_iota:.2e}"
        )
    logger.debug(
        f"Terminated current minimization after {nth_update} updates. Δiota = {rms_iota:.3e}, target Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
    )
    return state, rms_iota


def _free_boundary_target(state, B_ext_2):
    ev_gvec = state.evaluate("B", "p", rho=1.0).sel(rho=1.0).transpose("xyz", "tor", "pol")
    B_gvec_2 = (ev_gvec.B**2).sum(dim="xyz")
    B_jump = (B_gvec_2 - B_ext_2) / 2 + mu_0 * ev_gvec.p
    free_boundary_target = abs(fluxsurface_integral(B_jump) / fluxsurface_integral(B_gvec_2))
    return free_boundary_target


def run_free_boundary(
    state: State,
    coil_set: CoilSet | Coil,
    fb_updates: int,
    totalIter: int,
    max_iter_per_update: int = 1,
    free_boundary_tolerance=1e-8,
    I_tor: float | Callable | None = None,
    fac: int = 1,
    iter_B_update: int = 1,
    dt: float = 0.001,
    **kwargs,
):
    if I_tor is not None:
        current_constraint = True
        rho = np.linspace(0, 1, 101)
        if callable(I_tor):
            I_tor_at_rho = I_tor(rho)
        else:
            I_tor_at_rho = I_tor * np.ones_like(rho)
        I_tor_boundary = I_tor_at_rho[-1]
    else:
        current_constraint = False
        rho = None
        I_tor_at_rho = None
        I_tor_boundary = state.evaluate("I_tor", rho=1.0).I_tor.item()

    r_GP, r_w, t_n, t_w, z_n, z_w = state.get_integration_points()

    iter_B2_history = np.array([0, 1])
    B2_history = np.empty([t_n * z_n, 2])
    iters_used = 0
    free_boundary_target = 1e30
    for i in tqdm(range(fb_updates)):
        # abort if budget is used up
        if iters_used >= totalIter:
            return state, free_boundary_target

        # ensure that the maximum number of iterations is kept
        max_iter_per_update_aux = min(totalIter - iters_used, max_iter_per_update)

        # whether to interpolate the external field or not
        if i < 2:
            (B2, B_plas_ext, B_coils, B_ext_2) = _get_B2_one_fp(
                state=state,
                coil_set=coil_set,
                ntheta=t_n,
                nzeta=z_n,
                fac=fac,
                I_tor=I_tor_boundary,
                diagnostic_output=True,
            )
            B2_history[:, i] = B2

            free_boundary_target = _free_boundary_target(state, B_ext_2)
        else:
            if i % iter_B_update == 0:
                (B2, B_plas_ext, B_coils, B_ext_2) = _get_B2_one_fp(
                    state=state,
                    coil_set=coil_set,
                    ntheta=t_n,
                    nzeta=z_n,
                    fac=fac,
                    I_tor=I_tor_boundary,
                    diagnostic_output=True,
                )
                B2_history[:, 0] = B2_history[:, 1]
                B2_history[:, 1] = B2
                iter_B2_history[0] = iter_B2_history[1]
                iter_B2_history[1] = i
                free_boundary_target = _free_boundary_target(state, B_ext_2)
            else:
                # linear interpolation of the last two B2 solutions along the iterations
                B2 = B2_history[:, 0] + (i - iter_B2_history[0]) * (
                    B2_history[:, 1] - B2_history[:, 0]
                ) / (iter_B2_history[1] - iter_B2_history[0])
        if free_boundary_target > free_boundary_tolerance:
            state, dt = _free_boundary_update(
                state=state,
                maxIter=max_iter_per_update_aux,
                B2=B2,
                current_constraint=current_constraint,
                I_tor=I_tor_at_rho,
                rho=rho,
                dt=dt,
                **kwargs,
            )
            iters_used = state.get_iterations()
            if not current_constraint:
                I_tor_boundary = state.evaluate("I_tor", rho=1.0).I_tor.item()
        else:
            logger.info(
                f"Reached free boundary tolerance. free-boundary target: {free_boundary_target.item():3e}, specified tol.: {free_boundary_tolerance:3e}"
            )
            break
        logger.info(
            f"free-boundary target: {free_boundary_target.item()}, free-boundary target tolerance: {free_boundary_tolerance}"
        )

    state._minimize(logfinaliter=False, maxIter=(totalIter - iters_used))
    if free_boundary_target > free_boundary_tolerance:
        logger.warning(
            f"Free- boundary target not reached! free-boundary target: {free_boundary_target.item():3e}, specified tol.: {free_boundary_tolerance:3e}"
        )
    return state, free_boundary_target
