# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""pyGVEC postprocessing - Surface representation

This module provides a Surface class for representing a flux surface in 3D.
"""

# === Imports === #

from typing import Iterable
import logging
import functools

import numpy as np
import xarray as xr

from . import fourier
from .comp import register, compute
from .quantities import latex_partial_smart, derivative_name_smart

# === Globals === #

QUANTITIES_SURFACE = {}
compute = functools.partial(compute, registry=QUANTITIES_SURFACE)
register = functools.partial(register, registry=QUANTITIES_SURFACE)

# === Surface === #


def init_surface(x: np.ndarray, y: np.ndarray, z: np.ndarray, nfp: int = 1):
    if x.shape != y.shape != z.shape or x.ndim not in [2, 3]:
        raise ValueError(
            "x, y, and z must have the same shape '(rad,pol,tor)' or '(pol,tor)'."
        )

    if x.ndim == 2:
        x = x[np.newaxis, :, :]
        y = y[np.newaxis, :, :]
        z = z[np.newaxis, :, :]

    theta1d = np.linspace(0, 2 * np.pi, x.shape[1], endpoint=False)
    zeta1d = np.linspace(0, 2 * np.pi / nfp, x.shape[2], endpoint=False)
    theta, zeta = np.meshgrid(theta1d, zeta1d, indexing="ij")
    surfs = []

    for radidx in range(x.shape[0]):
        xhat = np.cos(zeta) * x[radidx, ...] + np.sin(zeta) * y[radidx, ...]
        yhat = -np.sin(zeta) * x[radidx, ...] + np.cos(zeta) * y[radidx, ...]
        zhat = z[radidx, ...]

        # Ignore stellarator symmetry: will not store fourier coefficients
        # ToDo: performance impact?
        xhatc, xhats = fourier.fft2d(xhat)
        yhatc, yhats = fourier.fft2d(yhat)
        zhatc, zhats = fourier.fft2d(zhat)

        surf = xr.Dataset(
            coords=dict(
                theta=("pol", theta1d),
                zeta=("tor", zeta1d),
            )
        )

        for var, c, s, symbol, name in [
            ("xhat", xhatc, xhats, r"\hat{x}", "modified x"),
            ("yhat", yhatc, yhats, r"\hat{y}", "modified y"),
            ("zhat", zhatc, zhats, r"\hat{z}", "modified z"),
        ]:
            surf[var] = (("pol", "tor"), fourier.eval2d(c, s, theta, zeta, nfp=nfp))
            surf[var].attrs["long_name"] = f"{name}-coordinate"
            surf[var].attrs["symbol"] = symbol
            for deriv in ["t", "z", "tt", "tz", "zz"]:
                dvar = f"d{var}_d{deriv}"
                surf[dvar] = (
                    ("pol", "tor"),
                    fourier.eval2d(c, s, theta, zeta, deriv, nfp=nfp),
                )
                surf[dvar].attrs["long_name"] = derivative_name_smart(
                    f"{name}-coordinate", deriv
                )
                surf[dvar].attrs["symbol"] = latex_partial_smart(symbol, deriv)
        surfs.append(surf)

    if len(surfs) == 1:
        return surfs[0]
    else:
        return xr.concat(surfs, dim="rad")


# === Computable Quantities === #


@register(
    attrs=dict(long_name="cartesian vector components", symbol=r"\mathbf{x}"),
)
def xyz(ds: xr.Dataset):
    ds.coords["xyz"] = ("xyz", ["x", "y", "z"])


@register(
    requirements=["xhat", "yhat", "zhat", "zeta", "xyz"],
    attrs=dict(
        long_name="cartesian coordinates",
        symbol=r"\mathbf{x}",
    ),
)
def pos(ds: xr.Dataset):
    ds["pos"] = xr.concat(
        [
            ds.xhat * np.cos(ds.zeta) - ds.yhat * np.sin(ds.zeta),
            ds.xhat * np.sin(ds.zeta) + ds.yhat * np.cos(ds.zeta),
            ds.zhat,
        ],
        dim="xyz",
    )


@register(
    requirements=["dxhat_dt", "dyhat_dt", "dzhat_dt", "zeta", "xyz"],
    attrs=dict(
        long_name="poloidal tangent basis vector", symbol=r"\mathbf{e}_{\theta}"
    ),
)
def e_theta(ds: xr.Dataset):
    ds["e_theta"] = xr.concat(
        [
            ds.dxhat_dt * np.cos(ds.zeta) - ds.dyhat_dt * np.sin(ds.zeta),
            ds.dxhat_dt * np.sin(ds.zeta) + ds.dyhat_dt * np.cos(ds.zeta),
            ds.dzhat_dt,
        ],
        dim="xyz",
    )


@register(
    requirements=["dxhat_dz", "dyhat_dz", "dzhat_dz", "xhat", "yhat", "zeta", "xyz"],
    attrs=dict(long_name="toroidal tangent basis vector", symbol=r"\mathbf{e}_{\zeta}"),
)
def e_zeta(ds: xr.Dataset):
    ds["e_zeta"] = xr.concat(
        [
            (ds.dxhat_dz - ds.yhat) * np.cos(ds.zeta)
            - (ds.dyhat_dz + ds.xhat) * np.sin(ds.zeta),
            (ds.dxhat_dz - ds.yhat) * np.sin(ds.zeta)
            + (ds.dyhat_dz + ds.xhat) * np.cos(ds.zeta),
            ds.dzhat_dz,
        ],
        dim="xyz",
    )


@register(
    requirements=["e_theta", "e_zeta"],
    attrs=dict(long_name="surface normal vector", symbol=r"\mathbf{n}"),
)
def normal(ds: xr.Dataset):
    n = xr.cross(ds.e_theta, ds.e_zeta, dim="xyz")
    ds["normal"] = n / np.sqrt(xr.dot(n, n, dim="xyz"))


@register(
    requirements=["e_theta"],
    attrs=dict(
        long_name="poloidal component of the metric tensor / first fundamental form",
        symbol=r"g_{\theta\theta}",
    ),
)
def g_tt(ds: xr.Dataset):
    ds["g_tt"] = xr.dot(ds.e_theta, ds.e_theta, dim="xyz")


@register(
    requirements=["e_theta", "e_zeta"],
    attrs=dict(
        long_name="poloidal-toroidal component of the metric tensor / first fundamental form",
        symbol=r"g_{\theta\zeta}",
    ),
)
def g_tz(ds: xr.Dataset):
    ds["g_tz"] = xr.dot(ds.e_theta, ds.e_zeta, dim="xyz")


@register(
    requirements=["e_zeta"],
    attrs=dict(
        long_name="toroidal component of the metric tensor / first fundamental form",
        symbol=r"g_{\zeta\zeta}",
    ),
)
def g_zz(ds: xr.Dataset):
    ds["g_zz"] = xr.dot(ds.e_zeta, ds.e_zeta, dim="xyz")


@register(
    requirements=["dxhat_dtt", "dyhat_dtt", "zeta", "xyz"],
    attrs=dict(
        long_name="poloidal curvature vector", symbol=r"\mathbf{k}_{\theta\theta}"
    ),
)
def k_tt(ds: xr.Dataset):
    ds["k_tt"] = xr.concat(
        [
            ds.dxhat_dtt * np.cos(ds.zeta) - ds.dyhat_dtt * np.sin(ds.zeta),
            ds.dxhat_dtt * np.sin(ds.zeta) + ds.dyhat_dtt * np.cos(ds.zeta),
            ds.dzhat_dtt,
        ],
        dim="xyz",
    )


@register(
    requirements=[
        "dxhat_dtz",
        "dyhat_dtz",
        "dzhat_dtz",
        "dyhat_dt",
        "dxhat_dt",
        "zeta",
        "xyz",
    ],
    attrs=dict(
        long_name="poloidal-toroidal curvature vector",
        symbol=r"\mathbf{k}_{\theta\zeta}",
    ),
)
def k_tz(ds: xr.Dataset):
    ds["k_tz"] = xr.concat(
        [
            (ds.dxhat_dtz - ds.dyhat_dt) * np.cos(ds.zeta)
            - (ds.dyhat_dtz + ds.dxhat_dt) * np.sin(ds.zeta),
            (ds.dxhat_dtz - ds.dyhat_dt) * np.sin(ds.zeta)
            + (ds.dyhat_dtz + ds.dxhat_dt) * np.cos(ds.zeta),
            ds.dzhat_dtz,
        ],
        dim="xyz",
    )


@register(
    requirements=[
        "dxhat_dzz",
        "dyhat_dzz",
        "dzhat_dzz",
        "dxhat_dz",
        "dyhat_dz",
        "xhat",
        "yhat",
        "zeta",
        "xyz",
    ],
    attrs=dict(
        long_name="toroidal curvature vector", symbol=r"\mathbf{k}_{\zeta\zeta}"
    ),
)
def k_zz(ds: xr.Dataset):
    ds["k_zz"] = xr.concat(
        [
            (ds.dxhat_dzz - 2 * ds.dyhat_dz - ds.xhat) * np.cos(ds.zeta)
            - (ds.dyhat_dzz + 2 * ds.dxhat_dz - ds.yhat) * np.sin(ds.zeta),
            (ds.dxhat_dzz - 2 * ds.dyhat_dz - ds.xhat) * np.sin(ds.zeta)
            + (ds.dyhat_dzz + 2 * ds.dxhat_dz - ds.yhat) * np.cos(ds.zeta),
            ds.dzhat_dzz,
        ],
        dim="xyz",
    )


@register(
    requirements=["normal", "k_tt"],
    attrs=dict(
        long_name="poloidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta\theta}",
    ),
)
def II_tt(ds: xr.Dataset):
    ds["II_tt"] = xr.dot(ds.normal, ds.k_tt, dim="xyz")


@register(
    requirements=["normal", "k_tz"],
    attrs=dict(
        long_name="poloidal-toroidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta\zeta}",
    ),
)
def II_tz(ds: xr.Dataset):
    ds["II_tz"] = xr.dot(ds.normal, ds.k_tz, dim="xyz")


@register(
    requirements=["normal", "k_zz"],
    attrs=dict(
        long_name="toroidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\zeta\zeta}",
    ),
)
def II_zz(ds: xr.Dataset):
    ds["II_zz"] = xr.dot(ds.normal, ds.k_zz, dim="xyz")
