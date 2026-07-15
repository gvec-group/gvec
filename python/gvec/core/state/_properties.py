# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._properties - low-level python API for accessing attributes of gvec.State"""

# === Imports === #

import logging

import numpy as np

from gvec.lib import _state
from ._binding import with_binding

# === Globals === #

logger = logging.getLogger("gvec.state")

# === Properties === #


@property
@with_binding
def nfp(self):
    return _state.nfp


@with_binding
def get_integration_points(self, quantity: str = "LA"):
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")

    r_n, t_n, z_n = _state.get_integration_points_num(quantity)
    r_GP, r_w = (np.zeros(r_n, dtype=np.float64) for _ in range(2))
    t_w, z_w = _state.get_integration_points(quantity, r_GP, r_w)
    return r_GP, r_w, t_n, t_w, z_n, z_w


@with_binding
def get_radial_gridpoints(self):
    nelems = _state.get_sgrid_nelems()
    s_pos = np.zeros(nelems + 1, dtype=np.float64, order="F")
    _state.get_sgrid_positions(s_pos)
    return s_pos


@with_binding
def get_mn_max(self, quantity: str = "all") -> tuple[int, int]:
    if not isinstance(quantity, str):
        raise ValueError("Quantity must be a string.")
    elif quantity not in ["X1", "X2", "LA", "all"]:
        raise ValueError(f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'.")

    if quantity == "all":
        m, n = zip(*[self.get_mn_max(q) for q in ["X1", "X2", "LA"]])
        return max(m), max(n)

    return _state.get_mn_max(quantity)
