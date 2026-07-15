# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._sfl - low-level python API for finding and evaluating straight-fieldline coordinates"""

# === Imports === #

import logging
from typing import Literal

import numpy as np

import gvec.lib
from gvec.lib import modgvec_py_state as _state
from ._binding import with_binding

# === Globals === #

logger = logging.getLogger("gvec.state")

# === Boozer Potential === #


@with_binding
def get_boozer(
    self,
    rho: np.ndarray,
    MNfactor: int = 5,
    *,
    M: int | None = None,
    N: int | None = None,
    M_nyq: int | None = None,
    N_nyq: int | None = None,
    sincos: Literal["sin", "cos", "sincos"] = "sin",
    recompute_lambda: bool = True,
):
    r"""
    Initialize a new Boozer potential with M poloidal and N toroidal nodes for all fluxsurfaces given by rho.

    Parameters
    ----------
    rho
        Array of (radius-like) flux surface labels.
    MNfactor
        Multiplication factor between the maximum M, N of the equilbrium and the maximum M, N of the Boozer potential.
        Only used if M, N are not explicitly given.
    M
        Number of poloidal nodes of the Boozer potential :math:`\nu_B`. Defaults to the maximum number of nodes of the basis.
    N
        Number of toroidal nodes of the Boozer potential :math:`\nu_B`. Defaults to the maximum number of nodes of the basis.

    Returns
    -------
    sfl_boozer
        Straight-fieldline Boozer object (wrapped Fortran object).
    """
    # --- Defaults --- #
    M_max, N_max = self.get_mn_max()
    _, M_nyq_LA, N_nyq_LA = _state.get_integration_points_num("LA")

    if M is None:
        M = MNfactor * M_max
    if N is None:
        N = MNfactor * N_max
    if M_nyq is None:
        M_nyq = max(4 * M + 1, M_nyq_LA)
    if N_nyq is None:
        N_nyq = max(4 * N + 1, N_nyq_LA)

    # --- Argument Handling --- #
    if not isinstance(M, int) or not isinstance(N, int) or M < 0 or N < 0:
        raise ValueError("M and N must be non-negative integers (or None).")
    if M < M_max or N < N_max:
        raise ValueError(
            f"The number of poloidal and toroidal nodes for the Boozer potential must be equal or larger to those of the original lambda: ({M}, {N}) < ({M_max}, {N_max})"
        )
    if (
        not isinstance(M_nyq, int)
        or not isinstance(N_nyq, int)
        or M_nyq < min(2 * M + 1, M_nyq_LA)
        or N_nyq < min(2 * N + 1, N_nyq_LA)
    ):
        raise ValueError(
            f"M_nyq and N_nyq must be integers larger than min({2 * M + 1=}, {M_nyq_LA=}) and min({2 * N + 1=}, {N_nyq_LA=}) (or None)."
        )

    rho = np.asfortranarray(rho, dtype=np.float64)
    if rho.ndim != 1 or rho.max() > 1.0 or rho.min() < 1e-4:
        raise ValueError("rho must be a 1D array in the range [1e-4, 1].")

    if sincos not in ["sin", "cos", "sincos"]:
        raise ValueError("sincos must be 'sin', 'cos', or 'sincos'.")
    sincos = {"sin": " _sin_", "cos": " _cos_", "sincos": "_sin_cos_"}[sincos]

    recompute_lambda = bool(recompute_lambda)

    # --- Create & compute Boozer potential --- #
    logger.debug("Initializing new Boozer potential.")
    sfl_boozer = _state.init_boozer(
        (M, N), (M_nyq, N_nyq), sincos, rho.size, rho, recompute_lambda
    )
    self._children.append(sfl_boozer)
    logger.debug(f"Computing Boozer potential {sfl_boozer!r}")
    _state.get_boozer(sfl_boozer)

    # ToDo: wrap sfl_boozer again to make it safer?
    return sfl_boozer


@with_binding
def get_boozer_angles(
    self,
    sfl_boozer: gvec.lib.Modgvec_Sfl_Boozer.t_sfl_boozer,
    tz_list: np.ndarray,
    rad: int | None = None,
):
    """
    Find the logical angles (theta, zeta) for the corresponding (theta_B, zeta_B) coordinates on the Boozer surface.

    Parameters
    ----------
    sfl_boozer : lib.Modgvec_Sfl_Boozer.t_sfl_boozer
        The Boozer potential object to use.
    tz_list : 2D array or Sequence of shape (2, n)
        The list of (theta_B, zeta_B) coordinates for which to find the logical angles.
        The first row contains theta_B and the second row contains zeta_B.
    rad : int, optional
        The (optional) radial index of the surface (with respect to the sfl_boozer object) on which to find the angles.

    Returns
    -------
    2D or 3D np.ndarray
        The logical angles (theta, zeta) corresponding to the input (theta_B, zeta_B) coordinates.
        If rad is None, returns a 3D array of shape (2, n, nrho), where nrho is the number of radial surfaces in sfl_boozer.
        If rad is specified, returns a 2D array of shape (2, n).
    """
    if not isinstance(sfl_boozer, gvec.lib.Modgvec_Sfl_Boozer.t_sfl_boozer):
        raise ValueError(f"Boozer object {sfl_boozer!r} must be of type `t_sfl_boozer`.")
    if sfl_boozer not in self._children:
        raise ValueError(f"Boozer object {sfl_boozer!r} is not known to the state {self!r}.")
    if not sfl_boozer.initialized:
        raise ValueError(f"Boozer object {sfl_boozer!r} is not initialized.")

    tz_list = np.asfortranarray(tz_list, dtype=np.float64)
    if tz_list.ndim != 2 or tz_list.shape[0] != 2:
        raise ValueError("thetazeta must be a 2D array with shape (2, n).")

    if rad is None:  # find angles on all surfaces
        tz_out = np.ndarray((2, tz_list.shape[1], sfl_boozer.nrho), dtype=np.float64, order="F")
        sfl_boozer.find_angles(tz_list.shape[1], tz_list, tz_out)
    else:  # find angles on a specific surface
        if rad not in range(sfl_boozer.nrho):
            raise ValueError(f"rad must be in the range [0, {sfl_boozer.nrho - 1}], got {rad}.")
        irho = rad + 1  # 1-indexed
        tz_out = np.ndarray((2, tz_list.shape[1]), dtype=np.float64, order="F")
        sfl_boozer.find_angles_irho(irho, tz_list.shape[1], tz_list, tz_out)
    return tz_out


@with_binding
def evaluate_boozer_list_tz_all(
    self,
    sfl_boozer: gvec.lib.Modgvec_Sfl_Boozer.t_sfl_boozer,
    quantity: str,
    rad: np.ndarray,
    thetazeta: np.ndarray,
):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["LA", "NU"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'LA', 'NU'.")
    if not isinstance(sfl_boozer, gvec.lib.Modgvec_Sfl_Boozer.t_sfl_boozer):
        raise ValueError(f"Boozer object {sfl_boozer!r} must be of type `t_sfl_boozer`.")
    if sfl_boozer not in self._children:
        raise ValueError(f"Boozer object {sfl_boozer!r} is not known to the state {self!r}.")
    if not sfl_boozer.initialized:
        raise ValueError(f"Boozer object {sfl_boozer!r} is not initialized.")

    rad = np.asfortranarray(rad, dtype=np.int64)
    thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
    if rad.ndim != 1:
        raise ValueError("rad must be a 1D array.")
    if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
        raise ValueError("thetazeta must be a 2D array with shape (2, n).")
    if rad.min() < 0:
        raise ValueError("rad must be a positive integer.")

    # Q, dQ_dtheta, dQ_dzeta, dQ_dtt, dQ_dtz, dQ_dzz
    outputs = [
        np.zeros((rad.size, thetazeta.shape[1]), dtype=np.float64, order="F") for _ in range(6)
    ]

    _state.evaluate_boozer_list_tz_all(
        sfl_boozer, rad.size, thetazeta.shape[1], rad, thetazeta, quantity, *outputs
    )
    return outputs


# === Straight-Fieldline PEST angles === #


@with_binding
def get_pest_angles(
    self,
    rho: np.ndarray,
    tz_list: np.ndarray,
):
    """
    Find the logical theta angle for the corresponding (theta_P, zeta) coordinates on the flux surface.

    Parameters
    ----------
    rho: 1D array
    tz_list : 2D array of shape (2, n) or 3D array of shape (2, n, rho.size)
        The list of (theta_P, zeta) coordinates for which to find the logical angles.
        The first row contains theta_P and the second row contains zeta.
        If a 2D array is given, the same postions are searched for on each surface.

    Returns
    -------
    2D np.ndarray of shape (n, rho.size)
        The logical theta angle corresponding to the input (theta_P, zeta) coordinates.
    """
    rho = np.asfortranarray(rho, dtype=np.float64)
    if rho.ndim != 1:
        raise ValueError(f"Expected a 1D array for rho, got shape {rho.shape}.")
    if rho.max() > 1.0 or rho.min() < 0.0:
        raise ValueError("rho must be in the range [0, 1].")
    tz_list = np.asfortranarray(tz_list, dtype=np.float64)
    match tz_list.shape:
        case (2, n):
            tz_out = np.zeros((2, n, rho.size), order="F")
            _state.find_pest_angles_2d(rho.size, rho, n, tz_list, tz_out)
            return tz_out[0, :, :]
        case (2, n, k) if k == rho.size:
            t_out = np.zeros((n, rho.size))
            tz_out = np.zeros((2, n, 1), order="F")
            for i, r in enumerate(rho):
                _state.find_pest_angles_2d(1, [r], n, tz_list[:, :, i], tz_out)
                t_out[:, i] = tz_out[0, :, 0]
            return t_out
        case _:
            raise ValueError(
                f"Expected 'tz_list' of shape (2, n) or (2, n, rho.size), got {tz_list.shape}."
            )
