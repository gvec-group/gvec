# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._properties - low-level python API for accessing attributes of gvec.State"""

# === Imports === #

import logging

import numpy as np

from gvec.lib import _state
from ._binding import with_binding
from gvec.core.compute import CoordinateSpec

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


@with_binding
def _get_angle_integration_points(self):
    """Extracts the theta and zeta integration points from fortran GVEC

    Returns
    -------
    (np.ndarray, np.ndarray)
        theta and zeta integration points.
    """
    r_GP, r_w, t_n, t_w, z_n, z_w = self.get_integration_points(quantity="X1")
    theta = np.zeros(t_n)
    zeta = np.zeros(z_n)
    _state.get_angle_integration_points(theta, t_n, zeta, z_n)
    return theta, zeta


@with_binding
def _get_X1_n_modes(self):
    return _state.get_x1_n_modes()


@with_binding
def _get_X2_n_modes(self):
    return _state.get_x2_n_modes()


@property
def winding_nr(self):
    """
    Determine the winding number of the axis around the z-axis.
    """
    z_n = self.get_integration_points()[4]
    nzetas = z_n * self.nfp
    zeta = np.linspace(0, 2 * np.pi, nzetas, endpoint=False)
    ev = self.evaluate("pos", "e_q3", rho=0.0, theta=0.0, zeta=zeta).squeeze()
    pos = ev.pos
    e_q3 = ev.e_q3
    denom = pos.sel(xyz="x") ** 2 + pos.sel(xyz="y") ** 2
    numer = pos.sel(xyz="x") * e_q3.sel(xyz="y") - pos.sel(xyz="y") * e_q3.sel(xyz="x")
    delta_zeta = zeta[1] - zeta[0]
    W = 1 / (2 * np.pi) * delta_zeta * (numer / denom).sum(dim="tor")
    return W.item()


@property
def clockwise_zeta(self):
    """
    Determines zeta direction.

    Returns
    -------
    bool
        clockwise_zeta

    Raises
    ------
    ValueError
        If the zeta direction cannot be determined.
    """
    if self.winding_nr > 0:
        return False
    elif self.winding_nr < 0:
        return True
    else:
        raise ValueError(
            f"Could not determine zeta direction! Winding number:{self.winding_nr}"
        )


def get_intersection_planes(self, zetas: CoordinateSpec = 3, min_box_size: float = 1e-3):
    """Create IntersectionPlane objects from a State object at specified zeta positions.

    Parameters
    ----------
    state : State
        GVEC state object to construct the planes from.
    zetas : CoordinateSpec, optional
        Zeta positions of the desired planes.
        If an integer is given linearly spaced zetas will be used. By default 3
    min_box_size : float, optional
        Minimum length of the bounding box, by default 1e-3.

    Returns
    -------
    List[IntersectionPlane]
        List of IntersectionPlane objects.
    """
    from gvec.coils import IntersectionPlane
    import xarray as xr

    boundary = self.evaluate("pos", zeta=zetas, theta=128, rho=1.0).pos

    planes = []
    for zeta in boundary.zeta:
        ev = self.evaluate(
            "pos",
            "e_q1",
            "e_q2",
            "e_q3",
            "X1",
            "X2",
            "grad_zeta",
            "mod_grad_zeta",
            rho=0.0,
            theta=0.0,
            zeta=zeta,
        ).squeeze()
        Jh = xr.dot(xr.cross(ev.e_q1, ev.e_q2, dim="xyz"), ev.e_q3, dim="xyz")
        surf_normal = (ev.grad_zeta / ev.mod_grad_zeta).data

        surf_point = ev.pos - ev.X1 * ev.e_q1 - ev.X2 * ev.e_q2
        grad_q1 = (xr.cross(ev.e_q2, ev.e_q3, dim="xyz") / Jh).data
        grad_q2 = (xr.cross(ev.e_q3, ev.e_q1, dim="xyz") / Jh).data

        boundary_x = boundary.sel(rho=1.0, xyz="x", zeta=zeta)
        boundary_y = boundary.sel(rho=1.0, xyz="y", zeta=zeta)
        boundary_z = boundary.sel(rho=1.0, xyz="z", zeta=zeta)

        xlim = np.array([boundary_x.min().data, boundary_x.max().data])
        ylim = np.array([boundary_y.min().data, boundary_y.max().data])
        zlim = np.array([boundary_z.min().data, boundary_z.max().data])
        dx = xlim[1] - xlim[0]
        dy = ylim[1] - ylim[0]
        dz = zlim[1] - zlim[0]

        if abs(dx) < min_box_size:
            xlim[0] -= min_box_size / 2
            xlim[1] += min_box_size / 2

        if abs(dy) < min_box_size:
            ylim[0] -= min_box_size / 2
            ylim[1] += min_box_size / 2

        if abs(dz) < min_box_size:
            zlim[0] -= min_box_size / 2
            zlim[1] += min_box_size / 2

        planes.append(
            IntersectionPlane(
                surf_normal=surf_normal,
                surf_point=surf_point,
                xlim=xlim,
                ylim=ylim,
                zlim=zlim,
                grad_q1=grad_q1,
                grad_q2=grad_q2,
            )
        )

    return planes
