# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""pyGVEC postprocessing - Surface representation

This module provides a Surface class for representing a flux surface in 3D.
Currently this surface is always defined in terms of the Boozer angles (theta_B, zeta_B).
"""

# === Imports === #

from typing import Literal
import functools
import warnings

import numpy as np
import xarray as xr

from gvec import fourier
from gvec.core.compute import register, compute, latex_partial_smart, derivative_name_smart

# === Globals === #

QUANTITIES_SURFACE = {}
compute = functools.partial(compute, registry=QUANTITIES_SURFACE)
register = functools.partial(register, registry=QUANTITIES_SURFACE)

# === Surface === #


def init_surface(
    pos: xr.DataArray,
    nfp: int | xr.DataArray = 1,
    ift: Literal["fft", "eval"] | None = None,
    winding: int = 1,
) -> xr.Dataset:
    if set(pos.dims) > {"xyz", "rad", "pol", "tor"} or set(pos.dims) < {
        "xyz",
        "pol",
        "tor",
    }:
        raise ValueError(
            f"expected pos to be a DataArray with dimensions ('xyz', 'rad', 'pol', 'tor') or ('xyz', 'pol', 'tor'), not {pos.dims}"
        )

    if "rad" not in pos.dims:
        pos = pos.expand_dims("rad")

    if isinstance(nfp, xr.DataArray):
        nfp = nfp.item()

    if ift == "fft" or ift is None:
        use_fft = False
        if (
            pos.theta_B.ndim == 1
            and pos.theta_B.size % 2 == 1
            and pos.zeta_B.ndim == 1
            and pos.zeta_B.size % 2 == 1
        ):
            theta1d = np.linspace(0, 2 * np.pi, pos.theta_B.size, endpoint=False)
            zeta1d = np.linspace(0, 2 * np.pi / nfp, pos.zeta_B.size, endpoint=False)
            if np.allclose(theta1d, pos.theta_B) and np.allclose(zeta1d, pos.zeta_B):
                use_fft = True
        if not use_fft:
            if ift == "fft":
                raise ValueError(
                    "Unaligned boozer angles: cannot use `ift=fft`; use `ift=eval` instead"
                )
            ift = "eval"
    if ift == "eval":
        theta1d = pos.theta_B
        zeta1d = pos.zeta_B
        theta2d, zeta2d = xr.broadcast(theta1d, zeta1d)
    if ift not in ["fft", "eval"]:
        raise ValueError(f"expected ift to be 'fft', 'eval' or None, not {ift}")

    surf = xr.Dataset(coords=pos.coords)
    surf["winding"] = ((), winding)
    surf["winding"].attrs = dict(
        long_name="signed winding number",
        symbol=r"\mathrm{winding}",
    )

    for r in pos.rad:
        x = pos.sel(xyz="x", rad=r)
        y = pos.sel(xyz="y", rad=r)
        xhat = np.cos(winding * zeta1d) * x + np.sin(winding * zeta1d) * y
        yhat = -np.sin(winding * zeta1d) * x + np.cos(winding * zeta1d) * y
        zhat = pos.sel(xyz="z", rad=r)

        # Ignore stellarator symmetry: will not store fourier coefficients - performance impact is negligible
        xhatc, xhats = fourier.fft2d(xhat.transpose("pol", "tor").data)
        yhatc, yhats = fourier.fft2d(yhat.transpose("pol", "tor").data)
        zhatc, zhats = fourier.fft2d(zhat.transpose("pol", "tor").data)

        for var, c, s, symbol, name in [
            ("xhat", xhatc, xhats, r"\hat{x}", "field periodic x"),
            ("yhat", yhatc, yhats, r"\hat{y}", "field periodic y"),
            ("zhat", zhatc, zhats, r"\hat{z}", "field periodic z"),
        ]:
            if var not in surf:
                surf[var] = (
                    ("rad", "pol", "tor"),
                    np.zeros((pos.rad.size, pos.pol.size, pos.tor.size)),
                )
            if ift == "fft":
                surf[var][r, :, :] = fourier.ifft2d(c, s)
            else:
                surf[var][r, :, :] = fourier.eval2d(c, s, theta2d.data, zeta2d.data, nfp=nfp)
            surf[var].attrs["long_name"] = f"{name}-coordinate"
            surf[var].attrs["symbol"] = symbol
            for deriv in ["t", "z", "tt", "tz", "zz"]:
                dvar = f"d{var}_d{deriv}"
                if dvar not in surf:
                    surf[dvar] = (
                        ("rad", "pol", "tor"),
                        np.zeros((pos.rad.size, pos.pol.size, pos.tor.size)),
                    )
                if ift == "fft":
                    surf[dvar][r, :, :] = fourier.ifft2d(c, s, deriv, nfp)
                else:
                    surf[dvar][r, :, :] = fourier.eval2d(
                        c, s, theta2d.data, zeta2d.data, deriv, nfp=nfp
                    )
                surf[dvar].attrs["long_name"] = derivative_name_smart(
                    f"{name}-coordinate", deriv
                )
                surf[dvar].attrs["symbol"] = latex_partial_smart(symbol, deriv)
    return surf


def get_xyz_hat(ds: xr.Dataset, winding: int | None = None):
    """Compute the xhat, yhat, zhat coordiantes for a dataset with pointwise x, y, z, zeta coordinates."""
    if winding is None:
        if "winding" in ds:
            winding = ds.winding.item()
        else:
            raise ValueError("winding must be specified or present in the dataset!")
    else:
        if "winding" in ds and winding != ds.winding.item():
            raise ValueError(
                "winding is specified but also present in the dataset with a different value!"
            )
        elif "winding" not in ds:
            ds["winding"] = xr.DataArray(
                winding,
                dims=(),
                attrs=dict(long_name="signed winding number", symbol=r"\mathrm{winding}"),
            )

    x = ds.pos.sel(xyz="x")
    y = ds.pos.sel(xyz="y")
    ds["xhat"] = np.cos(winding * ds.zeta) * x + np.sin(winding * ds.zeta) * y
    ds["yhat"] = -np.sin(winding * ds.zeta) * x + np.cos(winding * ds.zeta) * y
    ds["zhat"] = ds.pos.sel(xyz="z")

    for i in "xyz":
        ds[f"{i}hat"].attrs = dict(
            long_name=f"field periodic {i}-coordinate",
            symbol=rf"\hat{{{i}}}",
        )


# === Computable Quantities === #


@register(
    attrs=dict(long_name="cartesian vector components", symbol=r"\mathbf{x}"),
)
def xyz(ds: xr.Dataset):
    ds.coords["xyz"] = ("xyz", ["x", "y", "z"])


@register(
    requirements=["xhat", "yhat", "zhat", "zeta_B", "xyz", "winding"],
    attrs=dict(
        long_name="cartesian coordinates",
        symbol=r"\mathbf{x}",
    ),
)
def pos(ds: xr.Dataset):
    ds["pos"] = xr.concat(
        [
            ds.xhat * np.cos(ds.winding * ds.zeta_B) - ds.yhat * np.sin(ds.winding * ds.zeta_B),
            ds.xhat * np.sin(ds.winding * ds.zeta_B) + ds.yhat * np.cos(ds.winding * ds.zeta_B),
            ds.zhat,
        ],
        dim="xyz",
    )


@register(
    requirements=["dxhat_dt", "dyhat_dt", "dzhat_dt", "zeta_B", "xyz", "winding"],
    attrs=dict(long_name="poloidal tangent basis vector", symbol=r"\mathbf{e}_{\theta_B}"),
)
def e_theta_B(ds: xr.Dataset):
    ds["e_theta_B"] = xr.concat(
        [
            (
                ds.dxhat_dt * np.cos(ds.winding * ds.zeta_B)
                - ds.dyhat_dt * np.sin(ds.winding * ds.zeta_B)
            ),
            (
                ds.dxhat_dt * np.sin(ds.winding * ds.zeta_B)
                + ds.dyhat_dt * np.cos(ds.winding * ds.zeta_B)
            ),
            ds.dzhat_dt,
        ],
        dim="xyz",
    )


@register(
    requirements=[
        "dxhat_dz",
        "dyhat_dz",
        "dzhat_dz",
        "xhat",
        "yhat",
        "zeta_B",
        "xyz",
        "winding",
    ],
    attrs=dict(long_name="toroidal tangent basis vector", symbol=r"\mathbf{e}_{\zeta_B}"),
)
def e_zeta_B(ds: xr.Dataset):
    ds["e_zeta_B"] = xr.concat(
        [
            (ds.dxhat_dz - ds.yhat * ds.winding) * np.cos(ds.winding * ds.zeta_B)
            - (ds.dyhat_dz + ds.xhat * ds.winding) * np.sin(ds.winding * ds.zeta_B),
            (ds.dxhat_dz - ds.yhat * ds.winding) * np.sin(ds.winding * ds.zeta_B)
            + (ds.dyhat_dz + ds.xhat * ds.winding) * np.cos(ds.winding * ds.zeta_B),
            ds.dzhat_dz,
        ],
        dim="xyz",
    )


@register(
    requirements=["e_theta_B", "e_zeta_B"],
    attrs=dict(long_name="surface normal vector", symbol=r"\mathbf{n}"),
)
def normal(ds: xr.Dataset):
    n = xr.cross(ds.e_theta_B, ds.e_zeta_B, dim="xyz")
    ds["normal"] = n / np.sqrt(xr.dot(n, n, dim="xyz"))


@register(
    requirements=["e_theta_B"],
    attrs=dict(
        long_name="poloidal component of the metric tensor / first fundamental form",
        symbol=r"g_{\theta_B\theta_B}",
    ),
)
def g_tt_B(ds: xr.Dataset):
    ds["g_tt_B"] = xr.dot(ds.e_theta_B, ds.e_theta_B, dim="xyz")


@register(
    requirements=["e_theta_B", "e_zeta_B"],
    attrs=dict(
        long_name="poloidal-toroidal component of the metric tensor / first fundamental form",
        symbol=r"g_{\theta_B\zeta_B}",
    ),
)
def g_tz_B(ds: xr.Dataset):
    ds["g_tz_B"] = xr.dot(ds.e_theta_B, ds.e_zeta_B, dim="xyz")


@register(
    requirements=["e_zeta_B"],
    attrs=dict(
        long_name="toroidal component of the metric tensor / first fundamental form",
        symbol=r"g_{\zeta_B\zeta_B}",
    ),
)
def g_zz_B(ds: xr.Dataset):
    ds["g_zz_B"] = xr.dot(ds.e_zeta_B, ds.e_zeta_B, dim="xyz")


@register(
    requirements=["dxhat_dtt", "dyhat_dtt", "zeta_B", "xyz", "winding"],
    attrs=dict(long_name="poloidal curvature vector", symbol=r"\mathbf{k}_{\theta_B\theta_B}"),
)
def k_tt_B(ds: xr.Dataset):
    ds["k_tt_B"] = xr.concat(
        [
            (
                ds.dxhat_dtt * np.cos(ds.winding * ds.zeta_B)
                - ds.dyhat_dtt * np.sin(ds.winding * ds.zeta_B)
            ),
            (
                ds.dxhat_dtt * np.sin(ds.winding * ds.zeta_B)
                + ds.dyhat_dtt * np.cos(ds.winding * ds.zeta_B)
            ),
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
        "zeta_B",
        "xyz",
        "winding",
    ],
    attrs=dict(
        long_name="poloidal-toroidal curvature vector",
        symbol=r"\mathbf{k}_{\theta_B\zeta_B}",
    ),
)
def k_tz_B(ds: xr.Dataset):
    ds["k_tz_B"] = xr.concat(
        [
            (ds.dxhat_dtz - ds.dyhat_dt * ds.winding) * np.cos(ds.winding * ds.zeta_B)
            - (ds.dyhat_dtz + ds.dxhat_dt * ds.winding) * np.sin(ds.winding * ds.zeta_B),
            (ds.dxhat_dtz - ds.dyhat_dt * ds.winding) * np.sin(ds.winding * ds.zeta_B)
            + (ds.dyhat_dtz + ds.dxhat_dt * ds.winding) * np.cos(ds.winding * ds.zeta_B),
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
        "zeta_B",
        "xyz",
        "winding",
    ],
    attrs=dict(long_name="toroidal curvature vector", symbol=r"\mathbf{k}_{\zeta_B\zeta_B}"),
)
def k_zz_B(ds: xr.Dataset):
    ds["k_zz_B"] = xr.concat(
        [
            (ds.dxhat_dzz - 2 * ds.dyhat_dz * ds.winding - ds.xhat * ds.winding**2)
            * np.cos(ds.winding * ds.zeta_B)
            - (ds.dyhat_dzz + 2 * ds.dxhat_dz * ds.winding - ds.yhat * ds.winding**2)
            * np.sin(ds.winding * ds.zeta_B),
            (ds.dxhat_dzz - 2 * ds.dyhat_dz * ds.winding - ds.xhat * ds.winding**2)
            * np.sin(ds.winding * ds.zeta_B)
            + (ds.dyhat_dzz + 2 * ds.dxhat_dz * ds.winding - ds.yhat * ds.winding**2)
            * np.cos(ds.winding * ds.zeta_B),
            ds.dzhat_dzz,
        ],
        dim="xyz",
    )


@register(
    requirements=["normal", "k_tt_B"],
    attrs=dict(
        long_name="poloidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta_B\theta_B}",
    ),
)
def II_tt_B(ds: xr.Dataset):
    ds["II_tt_B"] = xr.dot(ds.normal, ds.k_tt_B, dim="xyz")


@register(
    requirements=["normal", "k_tz_B"],
    attrs=dict(
        long_name="poloidal-toroidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta_B\zeta_B}",
    ),
)
def II_tz_B(ds: xr.Dataset):
    ds["II_tz_B"] = xr.dot(ds.normal, ds.k_tz_B, dim="xyz")


@register(
    requirements=["normal", "k_zz_B"],
    attrs=dict(
        long_name="toroidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\zeta_B\zeta_B}",
    ),
)
def II_zz_B(ds: xr.Dataset):
    ds["II_zz_B"] = xr.dot(ds.normal, ds.k_zz_B, dim="xyz")
