import numpy as np
from gvec import State
from gvec.coils import CoilSet, Coil
from warnings import warn
from scipy.constants import mu_0
from tqdm import tqdm
from typing import Callable

try:
    import fbtools

    has_fbtools = True
except ImportError:
    has_fbtools = False


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
        warn(
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
    B_plas_ext = get_B_ext_plasma(
        boundary.pos.values,
        1,  # boundary.N_FP,
        B_coils,
        I_tor * mu_0,
        digits=10,
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
    (xr.DataArray, xr.DataArray)
        new iota values, rho^2 interpolation points
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
    return (iota_values, ev.rho**2)


def _free_boundary_update(
    state: State,
    maxIter: int,
    B2: np.ndarray,
    current_constraint: bool = False,
    I_tor: np.ndarray | None = None,
    rho: np.ndarray | None = None,
):
    # set Bext
    state.set_Bext_squared(B2)

    # minimize
    state._minimize(logfinaliter=False, maxIter=maxIter)
    if current_constraint:
        new_iota, new_rho = _get_new_iota(state, I_tor, rho)
        state.update_iota_for_current_constraint(new_rho, new_iota)
    return state


def run_free_boundary(
    state: State,
    coil_set: CoilSet | Coil,
    fb_updates: int,
    max_iter_per_update: int = 1,
    I_tor: float | Callable | None = None,
    fac: int = 1,
    iter_B_update: int = 1,
    **kwargs,
):
    if I_tor:
        current_constraint = True
        rho = np.linspace(0, 1, 101)
        if callable(I_tor):
            I_tor_at_rho = I_tor(rho)
        else:
            I_tor_at_rho = I_tor
        I_tor_boundary = I_tor_at_rho[-1]
    else:
        current_constraint = False
        rho = None
        I_tor_at_rho = None
        I_tor_boundary = state.evaluate("I_tor", rho=1.0).I_tor.item()

    r_GP, r_w, t_n, t_w, z_n, z_w = state.get_integration_points()

    iter_B2_history = np.array([0, 1])
    B2_history = np.empty([t_n * z_n, 2])

    for i in tqdm(range(fb_updates)):
        if i < 2:
            B2 = _get_B2_one_fp(
                state=state,
                coil_set=coil_set,
                ntheta=t_n,
                nzeta=z_n,
                fac=fac,
                I_tor=I_tor_boundary,
            )
            B2_history[:, i] = B2
        else:
            if i % iter_B_update == 0:
                B2 = _get_B2_one_fp(
                    state=state,
                    coil_set=coil_set,
                    ntheta=t_n,
                    nzeta=z_n,
                    fac=fac,
                    I_tor=I_tor_boundary,
                )
                B2_history[:, 0] = B2_history[:, 1]
                B2_history[:, 1] = B2
                iter_B2_history[0] = iter_B2_history[1]
                iter_B2_history[1] = i
            else:
                # linear interpolation of the last two B2 solutions along the iterations
                B2 = B2_history[:, 0] + (i - iter_B2_history[0]) * (
                    B2_history[:, 1] - B2_history[:, 0]
                ) / (iter_B2_history[1] - iter_B2_history[0])

        state = _free_boundary_update(
            state=state,
            maxIter=max_iter_per_update,
            B2=B2,
            current_constraint=current_constraint,
            I_tor=I_tor_at_rho,
            rho=rho,
            **kwargs,
        )
        if not current_constraint:
            I_tor_boundary = state.evaluate("I_tor", rho=1.0).I_tor.item()
    return state
