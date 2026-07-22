# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._evaluate - low-level python API for evaluating a gvec.State"""

# === Imports === #

import inspect
import re
from collections.abc import Iterable

import numpy as np

from gvec.lib import _state
from ._binding import with_binding

# === Decorators & Utility Functions === #


def _evaluate_1D_factory(
    func: callable, argnames: Iterable[str], n_out: int, vector_out: bool = False
):
    params = [inspect.Parameter("self", inspect.Parameter.POSITIONAL_OR_KEYWORD)] + [
        inspect.Parameter(name, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=np.ndarray)
        for name in argnames
    ]
    returns = tuple[tuple(np.ndarray for _ in range(n_out))]
    sig = inspect.Signature(params, return_annotation=returns)

    @with_binding
    def wrapper(self, *args, **kwargs):
        bound_args = sig.bind(self, *args, **kwargs)
        inputs = [
            np.asfortranarray(value, dtype=np.float64)
            for key, value in bound_args.arguments.items()
            if key != "self"
        ]
        n = inputs[0].size
        for value in inputs:
            if value.shape != (n,):
                raise ValueError("All arguments must be 1D arrays of the same size.")

        if vector_out:
            outputs = [np.zeros((3, n), dtype=np.float64, order="F") for _ in range(n_out)]
        else:
            outputs = [np.zeros(n, dtype=np.float64) for _ in range(n_out)]
        func(n, *inputs, *outputs)
        return outputs

    wrapper.__signature__ = sig
    wrapper.__name__ = func.__name__
    return wrapper


# === Evaluation Methods === #


@with_binding
def evaluate_base_tens(
    self,
    quantity: str,
    derivs: str | None,
    rho: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")
    if derivs is not None:
        if not isinstance(derivs, str):
            raise ValueError("Derivatives must be a string.")
        if m := re.match(r"(r{0,2})(t{0,2}|z{0,2}|tz)$", derivs):
            sel_derivs = m.groups()
        else:
            raise ValueError(f"Unknown derivative: {derivs}")
    else:
        sel_derivs = ("", "")

    rho = np.asfortranarray(rho, dtype=np.float64)
    theta = np.asfortranarray(theta, dtype=np.float64)
    zeta = np.asfortranarray(zeta, dtype=np.float64)
    if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
        raise ValueError("rho, theta, and zeta must be 1D arrays.")
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")

    result = np.zeros((rho.size, theta.size, zeta.size), dtype=np.float64, order="F")
    _state.evaluate_base_tens(rho, theta, zeta, quantity, *sel_derivs, result)
    return result


@with_binding
def evaluate_base_list_tz(
    self, quantity: str, derivs: str | None, rho: np.ndarray, thetazeta: np.ndarray
):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")
    if derivs is not None:
        if not isinstance(derivs, str):
            raise ValueError("Derivatives must be a string.")
        if m := re.match(r"(r{0,2})(t{0,2}|z{0,2}|tz)$", derivs):
            sel_derivs = m.groups()
        else:
            raise ValueError(f"Unknown derivative: {derivs}")
    else:
        sel_derivs = ("", "")

    rho = np.asfortranarray(rho, dtype=np.float64)
    thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
    if rho.ndim != 1:
        raise ValueError("rho must be a 1D array.")
    if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
        raise ValueError("thetazeta must be a 2D array with shape (2, n).")
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")

    result = np.zeros((rho.size, thetazeta.shape[1]), dtype=np.float64, order="F")
    _state.evaluate_base_list_tz(
        rho.size, thetazeta.shape[1], rho, thetazeta, quantity, *sel_derivs, result
    )
    return result


@with_binding
def evaluate_base_list_tz_all(self, quantity: str, rho: np.ndarray, thetazeta: np.ndarray):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")

    rho = np.asfortranarray(rho, dtype=np.float64)
    thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
    if rho.ndim != 1:
        raise ValueError("rho must be a 1D array.")
    if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
        raise ValueError("thetazeta must be a 2D array with shape (2, n).")
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")

    # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
    outputs = [
        np.zeros((rho.size, thetazeta.shape[1]), dtype=np.float64, order="F") for _ in range(10)
    ]

    _state.evaluate_base_list_tz_all(
        rho.size, thetazeta.shape[1], rho, thetazeta, quantity, *outputs
    )
    return outputs


@with_binding
def evaluate_base_list_rtz_all(self, quantity: str, rhothetazeta: np.ndarray):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")

    rhothetazeta = np.asarray(rhothetazeta)
    if rhothetazeta.ndim != 2 or rhothetazeta.shape[0] != 3:
        raise ValueError(
            f"rhothetazeta must be a 2D array with shape (3, n), but has shape {rhothetazeta.shape}"
        )
    rho = np.asfortranarray(rhothetazeta[0, :], dtype=np.float64)
    thetazeta = np.asfortranarray(rhothetazeta[1:, :], dtype=np.float64)
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")

    # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
    outputs = [
        np.zeros((rhothetazeta.shape[1]), dtype=np.float64, order="F") for _ in range(10)
    ]

    _state.evaluate_base_list_stz_all(rhothetazeta.shape[1], rho, thetazeta, quantity, *outputs)
    return outputs


@with_binding
def evaluate_base_tens_all(
    self, quantity: str, rho: np.ndarray, theta: np.ndarray, zeta: np.ndarray
):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")

    rho = np.asfortranarray(rho, dtype=np.float64)
    theta = np.asfortranarray(theta, dtype=np.float64)
    zeta = np.asfortranarray(zeta, dtype=np.float64)
    if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
        raise ValueError("rho, theta, and zeta must be 1D arrays.")
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")

    # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
    outputs = [
        np.zeros((rho.size, theta.size, zeta.size), dtype=np.float64, order="F")
        for _ in range(10)
    ]

    _state.evaluate_base_tens_all(
        rho.size, theta.size, zeta.size, rho, theta, zeta, quantity, *outputs
    )
    return outputs


evaluate_hmap = _evaluate_1D_factory(
    _state.evaluate_hmap,
    [
        "X1",
        "X2",
        "zeta",
    ]
    + [f"d{Q}_d{i}" for i in "rtz" for Q in ["X1", "X2"]],
    4,
    True,
)  # -> pos, e_rho, e_theta, e_zeta

evaluate_hmap_only = _evaluate_1D_factory(
    _state.evaluate_hmap_only, ["X1", "X2", "zeta"], 4, True
)  # -> pos, e_q1, e_q2, e_q3

evaluate_hmap_derivs = _evaluate_1D_factory(
    _state.evaluate_hmap_derivs, ["X1", "X2", "zeta"], 6, True
)  # -> k_q1q1, k_q1q2, k_q1q3, k_q2q2, k_q2q3, k_q3q3

evaluate_metric_derivs = _evaluate_1D_factory(
    _state.evaluate_metric_derivs,
    [
        "X1",
        "X2",
        "zeta",
    ]
    + [f"d{Q}_d{i}" for i in "r t z rr rt rz tt tz zz".split() for Q in ["X1", "X2"]],
    18,
)  # -> dg_rr_dr, dg_rt_dr ... dg_zz_dz

evaluate_jac_h_derivs = _evaluate_1D_factory(
    _state.evaluate_jac_h_derivs,
    [
        "X1",
        "X2",
        "zeta",
        "dX1_dr",
        "dX2_dr",
        "dX1_dt",
        "dX2_dt",
        "dX1_dz",
        "dX2_dz",
    ],
    3,
)  # ->  dJac_h_dr, dJac_h_dt, dJac_h_dz


@with_binding
def evaluate_profile(self, quantity: str, rho: np.ndarray, deriv: int = 0):
    """Evaluate 1D profiles at the provided positions of the radial coordinate ``rho``.

    Parameters
    ----------
    quantity : str
        name of the profile. Has to be either ``"iota"`` (rotational transform), ``"p"`` (pressure), ``"chi"`` (poloidal magn. flux), ``"Phi"`` (toroidal magn. flux)
    rho : ndarray
        Positions at the radial flux coordinate ``rho``.
    deriv : int, optional
        Order of the derivative in ``rho``. Note that for some quantities not all derivatives can be calculated, e.g. for rotational transform and pressure the maximum is ``deriv=4``. Defaults to ``0``.

    Raises
    ------
    ValueError
        If ``quantity`` is not a string, or if an invalid quantity is provided or if `rho`` is not a 1D array or if ``rho`` is not in ``[0, 1]``.
    NotImplementedError
        If ``deriv > 1`` for ``quantity="chi"``.

    Returns
    -------
    result : ndarray
        profile values at ``rho``.
    """
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["iota", "p", "chi", "Phi"]:
        raise ValueError(f"Unknown quantity: {quantity}")

    rho = np.asfortranarray(rho, dtype=np.float64)
    if rho.ndim != 1:
        raise ValueError("rho must be a 1D array.")
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")

    result = np.zeros(rho.size, dtype=np.float64, order="F")

    _state.evaluate_profile(rho.size, rho, deriv, quantity, result)
    return result


@with_binding
def evaluate_rho2_profile(self, quantity: str, rho2: np.ndarray, deriv: int = 0):
    r"""Evaluate 1D profiles at the provided positions of the radial coordinate ``rho2``:math:`=s=\rho^2`.
    Note: Use this routine to obtain derivarives with respect to :math:`s` coordinate, else use ``evaluate_profile``.

    Parameters
    ----------
    quantity : str
        name of the profile. Has to be either ``"iota"`` (rotational transform), ``"p"`` (pressure), ``"chi"`` (poloidal magn. flux), ``"Phi"`` (toroidal magn. flux)
    rho2 : ndarray
        Positions at the radial flux coordinate :math:`s=\rho^2`.
    deriv : int, optional
        Order of the derivative, in coordinate :math:`s=rho^2` (!). Defaults to ``0``.

    Raises
    ------
    ValueError
        If ``quantity`` is not a string, or if an invalid quantity is provided or if `rho`` is not a 1D array or if ``rho`` is not in ``[0, 1]``.

    Returns
    -------
    result : ndarray
        Profile values at ``rho2``.
    """
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["iota", "p", "chi", "Phi"]:
        raise ValueError(f"Unknown quantity: {quantity}")

    rho2 = np.asfortranarray(rho2, dtype=np.float64)
    if rho2.ndim != 1:
        raise ValueError("rho2 must be a 1D array.")
    if rho2.max() > 1.0 or rho2.min() < 0.0:
        raise ValueError("rho2 must be in the range [0, 1].")

    result = np.zeros(rho2.size, dtype=np.float64, order="F")

    _state.evaluate_rho2_profile(rho2.size, rho2, deriv, quantity, result)
    return result
