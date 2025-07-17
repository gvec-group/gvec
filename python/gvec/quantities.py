# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""
GVEC Postprocessing - computable quantities

This module defines various quantities and their computation functions for the GVEC package.

The module contains functions that compute different physical quantities such as rotational transform and pressure profiles,
coordinate mappings and their derivatives, magnetic field, current density and more.

These functions are registered with `compute.QUANTITIES` using the `@register` function decorator.
"""

import logging

import xarray as xr
import numpy as np

from gvec.core.state import State
from gvec.core.compute import (
    QUANTITIES,
    register,
    radial_integral,
    fluxsurface_integral,
    volume_integral,
    rtz_directions,
    rtz_symbols,
    derivative_name_smart,
    latex_partial,
    latex_partial_smart,
)


# === special ========================================================================== #


@register(
    attrs=dict(long_name="magnetic constant", symbol=r"\mu_0"),
)
def mu0(ds: xr.Dataset):
    from scipy.constants import mu_0

    ds["mu0"] = mu_0


@register(
    attrs=dict(long_name="adiabatic index", symbol=r"\gamma"),
)
def gamma(ds: xr.Dataset):
    # only gamma=0 is supported by gvec currently, in order to prescibe pressure profile directly.
    ds["gamma"] = 0.0


@register(
    attrs=dict(long_name="cartesian vector components", symbol=r"(x,y,z)"),
)
def xyz(ds: xr.Dataset):
    ds.coords["xyz"] = ("xyz", ["x", "y", "z"])


# === profiles ========================================================================= #


def _profile(var, evalvar, deriv, long_name, symbol):
    """Factory function for profile quantities."""

    @register(
        quantities=var,
        attrs=dict(long_name=long_name, symbol=symbol),
    )
    def profile(ds: xr.Dataset, state: State):
        if "rho" not in ds:
            raise KeyError("Evaluation of profiles requires the radial coordinate 'rho'.")
        if ds.rho.dims == ("rad",):
            ds[var] = ("rad", state.evaluate_profile(evalvar, ds.rho, deriv=deriv))
        else:
            rho = ds.rho.data.flatten()
            output = state.evaluate_profile(evalvar, rho, deriv=deriv)
            ds[var] = (ds.rho.dims, output.reshape(ds.rho.shape))

    return profile


for var, name, symbol in [
    ("iota", "rotational transform", r"\iota"),
    ("p", "pressure", r"p"),
    ("chi", "poloidal magnetic flux", r"\chi"),
    ("Phi", "toroidal magnetic flux", r"\Phi"),
]:
    globals()[var] = _profile(var, var, 0, name, symbol)
    globals()[f"d{var}_dr"] = _profile(
        f"d{var}_dr", var, 1, f"{name} gradient", f"\\frac{{d{symbol}}}{{d\\rho}}"
    )
    globals()[f"d{var}_drr"] = _profile(
        f"d{var}_drr", var, 2, f"{name} curvature", f"\\frac{{d^2{symbol}}}{{d\\rho^2}}"
    )


@register(
    attrs=dict(long_name="toroidal magnetic flux at the edge", symbol=r"\Phi_0"),
)
def Phi_edge(ds: xr.Dataset, state: State):
    ds["Phi_edge"] = ((), state.evaluate_profile("Phi", [1.0])[0])


# === base ============================================================================= #


def _base(var, long_name, symbol):
    """Factory function for base quantities."""

    @register(
        quantities=[var] + [f"d{var}_d{i}" for i in "r t z rr rt rz tt tz zz".split()],
        attrs={var: dict(long_name=long_name, symbol=symbol)}
        | {
            f"d{var}_d{i}": dict(
                long_name=derivative_name_smart(long_name, i),
                symbol=latex_partial_smart(symbol, i),
            )
            for i in ("r", "t", "z", "rr", "rt", "rz", "tt", "tz", "zz")
        },
    )
    def base(ds: xr.Dataset, state: State):
        if "rho" not in ds or "theta" not in ds or "zeta" not in ds:
            raise KeyError(
                "Evaluation of base variables requires 'rho', 'theta', 'zeta' to be defined."
            )

        # mesh in logical coordinates (rho, theta, zeta) -> rho(rad), theta(pol), zeta(tor)
        if ds.rho.dims == ("rad",) and ds.theta.dims == ("pol",) and ds.zeta.dims == ("tor",):
            outputs = state.evaluate_base_tens_all(var, ds.rho, ds.theta, ds.zeta)
            for key, value in zip(base.quantities, outputs):
                ds[key] = (("rad", "pol", "tor"), value)

        # mesh in other flux aligned coordinates e.g. (rho, theta_B, zeta_B) -> rho(rad), theta(rad, ...), zeta(rad, ...)
        elif ds.rho.dims == ("rad",):
            if "rad" in ds.theta.dims or "rad" in ds.zeta.dims:
                theta, zeta = xr.broadcast(ds.theta, ds.zeta)
                theta = theta.transpose("rad", ...)
                zeta = zeta.transpose("rad", ...)
                assert theta.dims == zeta.dims
                output_dims = theta.dims[1:]
                theta = theta.values.reshape(ds.rad.size, -1)
                zeta = zeta.values.reshape(ds.rad.size, -1)

                # Compute base on each radial position
                outputs = []
                for r, rho in enumerate(ds.rho.data):
                    thetazeta = np.stack([theta[r, :], zeta[r, :]], axis=0)
                    outputs.append(state.evaluate_base_list_tz_all(var, [rho], thetazeta))
                outputs = [np.stack(value) for value in zip(*outputs)]

            else:
                theta, zeta = xr.broadcast(ds.theta, ds.zeta)
                assert theta.dims == zeta.dims
                output_dims = theta.dims
                theta = theta.values.flatten()
                zeta = zeta.values.flatten()

                # Compute base on each radial position
                thetazeta = np.stack([theta, zeta], axis=0)
                outputs = state.evaluate_base_list_tz_all(var, ds.rho, thetazeta)

            # Write to dataset
            output_shape = [ds[dim].size for dim in output_dims]
            for key, value in zip(base.quantities, outputs):
                value = value.reshape(ds.rad.size, *output_shape)
                ds[key] = (("rad", *output_dims), value)

        # mesh in other coordinates
        else:
            rho, theta, zeta = xr.broadcast(ds.rho, ds.theta, ds.zeta)
            output_dims = rho.dims
            assert theta.dims == zeta.dims == output_dims
            rho = rho.values.flatten()
            theta = theta.values.flatten()
            zeta = zeta.values.flatten()
            rhothetazeta = np.stack([rho, theta, zeta], axis=0)

            # Compute base on each point individually
            outputs = state.evaluate_base_list_rtz_all(var, rhothetazeta)

            # Write to dataset
            output_shape = [ds[dim].size for dim in output_dims]
            for key, value in zip(base.quantities, outputs):
                value = value.reshape(output_shape)
                ds[key] = (output_dims, value)

    return base


for var, long_name, symbol in [
    ("X1", "first reference coordinate", r"X^1"),
    ("X2", "second reference coordinate", r"X^2"),
    ("LA", "straight field line potential", r"\lambda"),
]:
    globals()[var] = _base(var, long_name, symbol)


@register(
    attrs=dict(long_name="number of field periods", symbol=r"N_{FP}"),
)
def N_FP(ds: xr.Dataset, state: State):
    ds["N_FP"] = state.nfp


# === mapping ========================================================================== #


@register(
    quantities=("pos", "e_q1", "e_q2", "e_q3"),
    requirements=("xyz", "X1", "X2", "zeta"),
    attrs=dict(
        pos=dict(long_name="position vector", symbol=r"\mathbf{x}"),
        e_q1=dict(long_name="first reference tangent basis vector", symbol=r"\mathbf{e}_{q^1}"),
        e_q2=dict(
            long_name="second reference tangent basis vector",
            symbol=r"\mathbf{e}_{q^2}",
        ),
        e_q3=dict(
            long_name="toroidal reference tangent basis vector",
            symbol=r"\mathbf{e}_{q^3}",
        ),
    ),
)
def hmap(ds: xr.Dataset, state: State):
    X1, X2, zeta = xr.broadcast(ds.X1, ds.X2, ds.zeta)
    outputs = state.evaluate_hmap_only(*[v.values.flatten() for v in (X1, X2, zeta)])
    for key, value in zip(hmap.quantities, outputs):
        ds[key] = (
            ("xyz", *X1.dims),
            value.reshape(3, *X1.shape),
        )


# === metric =========================================================================== #
@register(
    quantities=["g_rr"],
    requirements=["e_rho"],
    attrs={
        "g_rr": dict(
            long_name="rr component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols['r'] + rtz_symbols['r']}}}",
        ),
    },
)
def g_rr(ds: xr.Dataset):
    ds["g_rr"] = xr.dot(ds.e_rho, ds.e_rho, dim="xyz")


@register(
    quantities=["g_rt"],
    requirements=["e_rho", "e_theta"],
    attrs={
        "g_rt": dict(
            long_name="rt component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols['r'] + rtz_symbols['t']}}}",
        ),
    },
)
def g_rt(ds: xr.Dataset):
    ds["g_rt"] = xr.dot(ds.e_rho, ds.e_theta, dim="xyz")


@register(
    quantities=["g_rz"],
    requirements=["e_rho", "e_zeta"],
    attrs={
        "g_rz": dict(
            long_name="rz component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols['r'] + rtz_symbols['z']}}}",
        ),
    },
)
def g_rz(ds: xr.Dataset):
    ds["g_rz"] = xr.dot(ds.e_rho, ds.e_zeta, dim="xyz")


@register(
    quantities=["g_tt"],
    requirements=["e_theta"],
    attrs={
        "g_tt": dict(
            long_name="tt component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols['t'] + rtz_symbols['t']}}}",
        ),
    },
)
def g_tt(ds: xr.Dataset):
    ds["g_tt"] = xr.dot(ds.e_theta, ds.e_theta, dim="xyz")


@register(
    quantities=["g_tz"],
    requirements=["e_theta", "e_zeta"],
    attrs={
        "g_tz": dict(
            long_name="tz component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols['t'] + rtz_symbols['z']}}}",
        ),
    },
)
def g_tz(ds: xr.Dataset):
    ds["g_tz"] = xr.dot(ds.e_theta, ds.e_zeta, dim="xyz")


@register(
    quantities=["g_zz"],
    requirements=["e_zeta"],
    attrs={
        "g_zz": dict(
            long_name="zz component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols['z'] + rtz_symbols['z']}}}",
        ),
    },
)
def g_zz(ds: xr.Dataset):
    ds["g_zz"] = xr.dot(ds.e_zeta, ds.e_zeta, dim="xyz")


@register(
    quantities=[
        pattern.format(ij=ij)
        for pattern in ("dg_{ij}_dr", "dg_{ij}_dt", "dg_{ij}_dz")
        for ij in ("rr", "rt", "rz", "tt", "tz", "zz")
    ],
    requirements=["X1", "X2", "zeta"]
    + [
        f"d{Xi}_d{j}"
        for j in ("r", "t", "z", "rr", "rt", "rz", "tt", "tz", "zz")
        for Xi in ("X1", "X2")
    ],
    attrs={
        f"dg_{ij}_d{k}": dict(
            long_name=derivative_name_smart(f"{ij} component of the metric tensor", k),
            symbol=latex_partial(f"g_{{{rtz_symbols[ij[0]] + rtz_symbols[ij[1]]}}}", k),
        )
        for ij in ("rr", "rt", "rz", "tt", "tz", "zz")
        for k in "rtz"
    },
)
def metric(ds: xr.Dataset, state: State):
    outputs = state.evaluate_metric_derivs(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in metric.requirements]
    )
    for key, value in zip(metric.quantities, outputs):
        ds[key] = (
            ds.X1.dims,
            value.reshape(ds.X1.shape),
        )


# === jacobian determinant ============================================================= #


@register(
    quantities=("Jac_h"),
    requirements=[
        "e_q1",
        "e_q2",
        "e_q3",
    ],
    attrs={
        "Jac_h": dict(long_name="reference Jacobian determinant", symbol=r"\mathcal{J}_h"),
    },
)
def Jac_h(ds: xr.Dataset):
    ds["Jac_h"] = xr.dot(ds.e_q1, xr.cross(ds.e_q2, ds.e_q3, dim="xyz"), dim="xyz")


@register(
    quantities=(*(f"dJac_h_d{i}" for i in "rtz"),),
    requirements=[
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
    attrs={
        f"dJac_h_d{i}": dict(
            long_name=derivative_name_smart("reference Jacobian determinant", i),
            symbol=latex_partial_smart(r"\mathcal{J}_h", i),
        )
        for i in "rtz"
    },
)
def Jac_h_derivs(ds: xr.Dataset, state: State):
    outputs = state.evaluate_jac_h_derivs(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in Jac_h_derivs.requirements]
    )
    for key, value in zip(Jac_h_derivs.quantities, outputs):
        ds[key] = (
            ds.X1.dims,
            value.reshape(ds.X1.shape),
        )


@register(
    quantities=(
        "Jac",
        "Jac_l",
    ),
    requirements=(
        "Jac_h",
        *(f"d{Q}_d{i}" for Q in "X1 X2".split() for i in "r t".split()),
    ),
    attrs={
        "Jac": dict(long_name="Jacobian determinant", symbol=r"\mathcal{J}"),
        "Jac_l": dict(long_name="logical Jacobian determinant", symbol=r"\mathcal{J}_l"),
    },
)
def Jac(ds: xr.Dataset):
    ds["Jac_l"] = ds.dX1_dr * ds.dX2_dt - ds.dX1_dt * ds.dX2_dr
    ds["Jac"] = ds.Jac_h * ds.Jac_l


@register(
    quantities=(*(f"dJac{suf}_d{i}" for suf in ["", "_l"] for i in "rtz"),),
    requirements=(
        "Jac_h",
        "Jac_l",
        *(f"dJac_h_d{i}" for i in "rtz"),
        *(f"d{Q}_d{i}" for Q in "X1 X2".split() for i in "r t z rr rt rz tt tz zz".split()),
    ),
    attrs={
        f"dJac_d{i}": dict(
            long_name=derivative_name_smart("Jacobian determinant", i),
            symbol=latex_partial_smart(r"\mathcal{J}", i),
        )
        for i in "rtz"
    }
    | {
        f"dJac_l_d{i}": dict(
            long_name=derivative_name_smart("logical Jacobian determinant", i),
            symbol=latex_partial_smart(r"\mathcal{J}_l", i),
        )
        for i in "rtz"
    },
)
def Jac_derivs(ds: xr.Dataset):
    ds["dJac_l_dr"] = (
        ds.dX1_drr * ds.dX2_dt
        + ds.dX1_dr * ds.dX2_drt
        - ds.dX1_drt * ds.dX2_dr
        - ds.dX1_dt * ds.dX2_drr
    )
    ds["dJac_l_dt"] = (
        ds.dX1_drt * ds.dX2_dt
        + ds.dX1_dr * ds.dX2_dtt
        - ds.dX1_dtt * ds.dX2_dr
        - ds.dX1_dt * ds.dX2_drt
    )
    ds["dJac_l_dz"] = (
        ds.dX1_drz * ds.dX2_dt
        + ds.dX1_dr * ds.dX2_dtz
        - ds.dX1_dtz * ds.dX2_dr
        - ds.dX1_dt * ds.dX2_drz
    )
    ds["dJac_dr"] = ds.dJac_h_dr * ds.Jac_l + ds.Jac_h * ds.dJac_l_dr
    ds["dJac_dt"] = ds.dJac_h_dt * ds.Jac_l + ds.Jac_h * ds.dJac_l_dt
    ds["dJac_dz"] = ds.dJac_h_dz * ds.Jac_l + ds.Jac_h * ds.dJac_l_dz


# === straight field line coordinates - PEST =========================================== #


@register(
    requirements=("LA",),
    attrs=dict(long_name="poloidal angle in PEST coordinates", symbol=r"\theta_P"),
)
def theta_P(ds: xr.Dataset):
    ds["theta_P"] = ds.theta + ds.LA


@register(
    requirements=("xyz", "theta_sfl", "dLA_dr", "dLA_dt", "dLA_dz"),
    attrs=dict(
        long_name="poloidal reciprocal basis vector in PEST coordinates",
        symbol=r"\nabla \theta_P",
    ),
)
def grad_theta_P(ds: xr.Dataset):
    ds["grad_theta_P"] = (
        ds.grad_theta * (1 + ds.dLA_dt) + ds.grad_rho * ds.dLA_dr + ds.grad_zeta * ds.dLA_dz
    )


@register(
    requirements=("Jac", "dLA_dt"),
    attrs=dict(long_name="Jacobian determinant in PEST coordinates", symbol=r"\mathcal{J}_P"),
)
def Jac_P(ds: xr.Dataset):
    ds["Jac_P"] = ds.Jac / (1 + ds.dLA_dt)


# === derived ========================================================================== #


@register(
    requirements=("xyz", "e_q1", "e_q2", "dX1_dr", "dX2_dr"),
    attrs=dict(long_name="radial tangent basis vector", symbol=r"\mathbf{e}_\rho"),
)
def e_rho(ds: xr.Dataset):
    ds["e_rho"] = ds.e_q1 * ds.dX1_dr + ds.e_q2 * ds.dX2_dr


@register(
    requirements=("xyz", "e_q1", "e_q2", "dX1_dt", "dX2_dt"),
    attrs=dict(long_name="poloidal tangent basis vector", symbol=r"\mathbf{e}_\theta"),
)
def e_theta(ds: xr.Dataset):
    ds["e_theta"] = ds.e_q1 * ds.dX1_dt + ds.e_q2 * ds.dX2_dt


@register(
    requirements=("xyz", "e_q1", "e_q2", "e_q3", "dX1_dz", "dX2_dz"),
    attrs=dict(long_name="toroidal tangent basis vector", symbol=r"\mathbf{e}_\zeta"),
)
def e_zeta(ds: xr.Dataset):
    ds["e_zeta"] = ds.e_q1 * ds.dX1_dz + ds.e_q2 * ds.dX2_dz + ds.e_q3


@register(
    requirements=("xyz", "Jac", "e_theta", "e_zeta"),
    attrs=dict(long_name="radial reciprocal basis vector", symbol=r"\nabla\rho"),
)
def grad_rho(ds: xr.Dataset):
    ds["grad_rho"] = xr.cross(ds.e_theta, ds.e_zeta, dim="xyz") / ds.Jac


@register(
    requirements=("xyz", "Jac", "e_rho", "e_zeta"),
    attrs=dict(long_name="poloidal reciprocal basis vector", symbol=r"\nabla\theta"),
)
def grad_theta(ds: xr.Dataset):
    ds["grad_theta"] = xr.cross(ds.e_zeta, ds.e_rho, dim="xyz") / ds.Jac


@register(
    requirements=("xyz", "Jac", "e_rho", "e_theta"),
    attrs=dict(long_name="toroidal reciprocal basis vector", symbol=r"\nabla\zeta"),
)
def grad_zeta(ds: xr.Dataset):
    ds["grad_zeta"] = xr.cross(ds.e_rho, ds.e_theta, dim="xyz") / ds.Jac


@register(
    quantities=("B", "B_contra_t", "B_contra_z"),
    requirements=(
        "xyz",
        "iota",
        "dLA_dt",
        "dLA_dz",
        "dPhi_dr",
        "Jac",
        "e_theta",
        "e_zeta",
    ),
    attrs=dict(
        B=dict(long_name="magnetic field", symbol=r"\mathbf{B}"),
        B_contra_t=dict(long_name="poloidal magnetic field", symbol=r"B^\theta"),
        B_contra_z=dict(long_name="toroidal magnetic field", symbol=r"B^\zeta"),
    ),
)
def B(ds: xr.Dataset):
    ds["B_contra_t"] = (ds.iota - ds.dLA_dz) * ds.dPhi_dr / ds.Jac
    ds["B_contra_z"] = (1 + ds.dLA_dt) * ds.dPhi_dr / ds.Jac
    ds["B"] = ds.B_contra_t * ds.e_theta + ds.B_contra_z * ds.e_zeta


@register(
    quantities=[f"dB_contra_{i}_d{j}" for i in "tz" for j in "rtz"],
    requirements=[
        "Jac",
        "dPhi_dr",
        "dPhi_drr",
        "iota",
        "diota_dr",
    ]
    + [f"dJac_d{i}" for i in "r t z".split()]
    + [f"dLA_d{i}" for i in "t z rt rz tt tz zz".split()],
    attrs={
        f"dB_contra_{i}_d{j}": dict(
            long_name=derivative_name_smart(f"{rtz_directions[i]} magnetic field", j),
            symbol=latex_partial(f"B^{rtz_symbols[i]}", j),
        )
        for i in "tz"
        for j in "rtz"
    },
)
def dB(ds: xr.Dataset):
    ds["dB_contra_t_dr"] = -ds.dPhi_dr / ds.Jac * (
        ds.dJac_dr / ds.Jac * (ds.iota - ds.dLA_dz) + ds.dLA_drz - ds.diota_dr
    ) + ds.dPhi_drr / ds.Jac * (ds.iota - ds.dLA_dz)
    ds["dB_contra_t_dt"] = -(ds.dPhi_dr / ds.Jac) * (
        ds.dJac_dt / ds.Jac * (ds.iota - ds.dLA_dz) + ds.dLA_dtz
    )
    ds["dB_contra_t_dz"] = -(ds.dPhi_dr / ds.Jac) * (
        ds.dJac_dz / ds.Jac * (ds.iota - ds.dLA_dz) + ds.dLA_dzz
    )
    ds["dB_contra_z_dr"] = -ds.dPhi_dr / ds.Jac * (
        ds.dJac_dr / ds.Jac * (1 + ds.dLA_dt) - ds.dLA_drt
    ) + ds.dPhi_drr / ds.Jac * (1 + ds.dLA_dt)
    ds["dB_contra_z_dt"] = (
        -ds.dPhi_dr / ds.Jac * (ds.dJac_dt / ds.Jac * (1 + ds.dLA_dt) - ds.dLA_dtt)
    )
    ds["dB_contra_z_dz"] = (
        -ds.dPhi_dr / ds.Jac * (ds.dJac_dz / ds.Jac * (1 + ds.dLA_dt) - ds.dLA_dtz)
    )


@register(
    quantities=["J", "J_contra_r", "J_contra_t", "J_contra_z"],
    requirements=[
        "B_contra_t",
        "B_contra_z",
        "Jac",
        "mu0",
    ]
    + [f"g_{ij}" for ij in "rt rz tt tz zz".split()]
    + [f"dg_{ij}_d{k}" for ij in "rt rz tt tz zz".split() for k in "rtz"]
    + [f"dB_contra_{i}_d{j}" for i in "tz" for j in "rtz"]
    + [f"e_{i}" for i in "rho theta zeta".split()],
    attrs={
        "J": dict(long_name="current density", symbol=r"\mathbf{J}"),
    }
    | {
        f"J_contra_{i}": dict(
            long_name=f"contravariant {rtz_directions[i]} current density",
            symbol=rf"J^{{{rtz_symbols[i]}}}",
        )
        for i in "rtz"
    },
)
def J(ds: xr.Dataset):
    def ij(i, j):
        if i < j:
            return i + j
        return j + i

    dB_co = {}
    for i in "rtz":
        for j in "rtz":
            if i == j:
                continue
            dB_co[i, j] = sum(
                ds[f"dg_{ij(i, k)}_d{j}"] * ds[f"B_contra_{k}"]
                + ds[f"g_{ij(i, k)}"] * ds[f"dB_contra_{k}_d{j}"]
                for k in "tz"
            )
    ds["J_contra_r"] = (dB_co["z", "t"] - dB_co["t", "z"]) / (ds.mu0 * ds.Jac)
    ds["J_contra_t"] = (dB_co["r", "z"] - dB_co["z", "r"]) / (ds.mu0 * ds.Jac)
    ds["J_contra_z"] = (dB_co["t", "r"] - dB_co["r", "t"]) / (ds.mu0 * ds.Jac)
    ds["J"] = ds.J_contra_r * ds.e_rho + ds.J_contra_t * ds.e_theta + ds.J_contra_z * ds.e_zeta


@register(
    requirements=("J", "B", "dp_dr", "grad_rho"),
    attrs=dict(
        long_name="MHD force",
        symbol=r"F",
    ),
)
def F(ds: xr.Dataset):
    ds["F"] = xr.cross(ds.J, ds.B, dim="xyz") - ds.dp_dr * ds.grad_rho


@register(
    requirements=("F", "e_rho", "Jac"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="radial force balance",
        symbol=r"\overline{F_\rho}",
    ),
)
def F_r_avg(ds: xr.Dataset):
    ds["F_r_avg"] = fluxsurface_integral(
        xr.dot(ds.F, ds.e_rho, dim="xyz") * ds.Jac
    ) / fluxsurface_integral(ds.Jac)


def _modulus(v):
    """Factory function for modulus (absolute value) quantities."""

    @register(
        quantities=f"mod_{v}",
        requirements=(v,),
        attrs=dict(
            long_name=f"modulus of the {QUANTITIES[v].attrs[v]['long_name']}",
            symbol=rf"\left|{QUANTITIES[v].attrs[v]['symbol']}\right|",
        ),
    )
    def mod_v(ds: xr.Dataset):
        ds[f"mod_{v}"] = np.sqrt(xr.dot(ds[v], ds[v], dim="xyz"))

    return mod_v


for v in [
    "e_rho",
    "e_theta",
    "e_zeta",
    "grad_rho",
    "grad_theta",
    "grad_zeta",
    "B",
    "J",
    "F",
]:
    globals()[v] = _modulus(v)


# === Straight Field Line Coordinates - Boozer ========================================= #
# dNU_B_dt and dNU_B_dz are overwritten when a boozer transform is performed!


@register(
    requirements=("B", "e_theta", "dLA_dt", "iota", "B_theta_avg", "B_zeta_avg"),
    attrs=dict(
        long_name="poloidal derivative of the Boozer potential computed from the magnetic field",
        symbol=r"\left." + latex_partial(r"\nu_B", "t") + r"\right|",
    ),
)
def dNU_B_dt(ds: xr.Dataset):
    Bt = xr.dot(ds.B, ds.e_theta, dim="xyz")
    ds["dNU_B_dt"] = (Bt - ds.B_theta_avg * (1 + ds.dLA_dt)) / (
        ds.iota * ds.B_theta_avg + ds.B_zeta_avg
    )


@register(
    requirements=("B", "e_zeta", "dLA_dz", "iota", "B_theta_avg", "B_zeta_avg"),
    attrs=dict(
        long_name="toroidal derivative of the Boozer potential computed from the magnetic field",
        symbol=r"\left." + latex_partial(r"\nu_B", "z") + r"\right|",
    ),
)
def dNU_B_dz(ds: xr.Dataset):
    Bz = xr.dot(ds.B, ds.e_zeta, dim="xyz")
    ds["dNU_B_dz"] = (Bz - ds.B_theta_avg * ds.dLA_dz - ds.B_zeta_avg) / (
        ds.iota * ds.B_theta_avg + ds.B_zeta_avg
    )


@register(
    requirements=("Jac", "dLA_dt", "dLA_dz", "dNU_B_dt", "dNU_B_dz", "iota"),
    attrs=dict(
        long_name="Jacobian determinant in Boozer coordinates",
        symbol=r"\mathcal{J}_B",
    ),
)
def Jac_B(ds: xr.Dataset):
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    ds["Jac_B"] = ds.Jac / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)


@register(
    requirements=("dPhi_dr", "iota", "Jac_B"),
    attrs=dict(
        long_name="contravariant poloidal magnetic field in Boozer coordinates",
        symbol=r"B^{\theta_B}",
    ),
)
def B_contra_t_B(ds: xr.Dataset):
    ds["B_contra_t_B"] = ds.dPhi_dr * ds.iota / ds.Jac_B


@register(
    requirements=("dPhi_dr", "Jac_B"),
    attrs=dict(
        long_name="contravariant toroidal magnetic field in Boozer coordinates",
        symbol=r"B^{\zeta_B}",
    ),
)
def B_contra_z_B(ds: xr.Dataset):
    ds["B_contra_z_B"] = ds.dPhi_dr / ds.Jac_B


@register(
    requirements=(
        "iota",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dz",
        "dNU_B_dt",
        "e_theta",
        "e_zeta",
    ),
    attrs=dict(
        long_name="poloidal tangent basis vector in Boozer coordinates",
        symbol=r"\mathbf{e}_{\theta_B}",
    ),
)
def e_theta_B(ds: xr.Dataset):
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    ds["e_theta_B"] = (dzB_dz * ds.e_theta - dzB_dt * ds.e_zeta) / (
        dtB_dt * dzB_dz - dtB_dz * dzB_dt
    )


@register(
    requirements=(
        "iota",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dz",
        "dNU_B_dt",
        "e_theta",
        "e_zeta",
    ),
    attrs=dict(
        long_name="toroidal tangent basis vector in Boozer coordinates",
        symbol=r"\mathbf{e}_{\zeta_B}",
    ),
)
def e_zeta_B(ds: xr.Dataset):
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    ds["e_zeta_B"] = (dtB_dt * ds.e_zeta - dtB_dz * ds.e_theta) / (
        dtB_dt * dzB_dz - dtB_dz * dzB_dt
    )


# === integrals ======================================================================== #


@register(
    requirements=("Jac",),
    integration=("rho", "theta", "zeta"),
    attrs=dict(long_name="plasma volume", symbol=r"V"),
)
def V(ds: xr.Dataset):
    ds["V"] = volume_integral(ds.Jac)


@register(
    requirements=("Jac", "Phi_edge", "dPhi_dr"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="derivative of the plasma volume w.r.t. normalized toroidal magnetic flux",
        symbol=r"\frac{dV}{d\Phi_n}",
    ),
)
def dV_dPhi_n(ds: xr.Dataset):
    # d/dPhi_n = dr/dPhi_n * d/dr = Phi_0 / dPhi_dr * d/dr
    ds["dV_dPhi_n"] = fluxsurface_integral(ds.Jac) * ds.Phi_edge / ds.dPhi_dr


@register(
    requirements=("Jac", "dJac_dr", "Phi_edge", "dPhi_dr", "dPhi_drr"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="second derivative of the plasma volume w.r.t. normalized toroidal magnetic flux",
        symbol=r"\frac{d^2V}{d\Phi_n^2}",
    ),
)
def dV_dPhi_n2(ds: xr.Dataset):
    # d/dPhi_n = dr/dPhi_n * d/dr = Phi_0 / dPhi_dr * d/dr
    # d/dr 1/dPhi_dr = -1/dPhi_dr**2 * dPhi_drr
    ds["dV_dPhi_n2"] = (
        fluxsurface_integral(ds.dJac_dr) * (ds.Phi_edge / ds.dPhi_dr) ** 2
        - fluxsurface_integral(ds.Jac) * ds.Phi_edge**2 / ds.dPhi_dr**3 * ds.dPhi_drr
    )


@register(
    quantities=("minor_radius", "major_radius"),
    requirements=("V", "Jac_l"),
    integration=("rho", "theta", "zeta"),
    attrs=dict(
        minor_radius=dict(long_name="minor radius", symbol=r"r_{min}"),
        major_radius=dict(long_name="major radius", symbol=r"r_{maj}"),
    ),
)
def minor_major_radius(ds: xr.Dataset):
    surface_average = volume_integral(ds.Jac_l) / (2 * np.pi)
    ds["minor_radius"] = np.sqrt(surface_average / np.pi)
    ds["major_radius"] = np.sqrt(ds.V / (2 * np.pi * surface_average))


@register(
    requirements=("iota",),
    integration=("rho",),
    attrs=dict(long_name="average rotational transform", symbol=r"\bar{\iota}"),
)
def iota_avg(ds: xr.Dataset):
    ds["iota_avg"] = radial_integral(ds.iota)


@register(
    requirements=("mu0", "dPhi_dr", "Jac", "g_tt"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="factor to the toroidal current contribution to the rotational transform",
        description="iota = iota_0 + iota_curr, iota_curr = I_tor * iota_curr_0",
        symbol=r"\iota_{curr,0}",
    ),
)
def iota_curr_0(ds: xr.Dataset):
    Gamma_t = fluxsurface_integral(ds.g_tt / ds.Jac)
    ds["iota_curr_0"] = 2 * np.pi * ds.mu0 / ds.dPhi_dr / Gamma_t


@register(
    requirements=("I_tor", "iota_curr_0"),
    attrs=dict(
        long_name="toroidal current contribution to the rotational transform",
        description="iota = iota_0 + iota_curr, iota_curr = I_tor * iota_curr_0",
        symbol=r"\iota_{curr}",
    ),
)
def iota_curr(ds: xr.Dataset):
    ds["iota_curr"] = ds.I_tor * ds.iota_curr_0
    # = 2 * np.pi * ds.mu0 * ds.dI_tor_dr / (ds.dPhi_drr * Gamma_t + ds.dPhi_dr * dGamma_t_dr)
    # = 2 * np.pi * ds.mu0 * ds.dI_tor_drr / (ds.dPhi_drr * dGamma_t_dr + ds.dPhi_dr * dGamma_t_drr)
    # = 2 * np.pi * ds.mu0 * ds.dI_tor_drr / (ds.dPhi_drr * dGamma_t_dr)

    # Gamma_t_dr = fluxsurface_integral(ds.dg_tt_dr / ds.Jac - ds.g_tt / ds.Jac**2 * ds.dJac_dr)
    # = ds.dg_tt_drr / ds.dJac_dr - ds.dg_tt_dr / (2 * ds.Jac * ds.dJac_dr) * ds.dJac_dr


@register(
    requirements=("g_tt", "g_tz", "Jac", "dLA_dt", "dLA_dz"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="geometric contribution to the rotational transform",
        description="iota = iota_0 + iota_curr",
        symbol=r"\iota_0",
    ),
)
def iota_0(ds: xr.Dataset):
    ds["iota_0"] = (
        fluxsurface_integral(ds.g_tt / ds.Jac * ds.dLA_dz)
        - fluxsurface_integral(ds.g_tz / ds.Jac)
        - fluxsurface_integral(ds.g_tz / ds.Jac * ds.dLA_dt)
    ) / fluxsurface_integral(ds.g_tt / ds.Jac)


@register(
    requirements=("iota", "diota_dr"),
    attrs=dict(long_name="global magnetic shear", symbol=r"s_g"),
)
def shear(ds: xr.Dataset):
    ds["shear"] = -ds.rho / ds.iota * ds.diota_dr


@register(
    requirements=("B", "e_theta"),
    integration=("theta", "zeta"),
    attrs=dict(long_name="average poloidal magnetic field", symbol=r"\overline{B_\theta}"),
)
def B_theta_avg(ds: xr.Dataset):
    ds["B_theta_avg"] = fluxsurface_integral(xr.dot(ds.B, ds.e_theta, dim="xyz")) / (
        4 * np.pi**2
    )


@register(
    requirements=("B_theta_avg", "mu0"),
    attrs=dict(long_name="toroidal current", symbol=r"I_{tor}"),
)
def I_tor(ds: xr.Dataset):
    ds["I_tor"] = ds.B_theta_avg * 2 * np.pi / ds.mu0


@register(
    requirements=("B", "e_zeta"),
    integration=("theta", "zeta"),
    attrs=dict(long_name="average toroidal magnetic field", symbol=r"\overline{B_\zeta}"),
)
def B_zeta_avg(ds: xr.Dataset):
    ds["B_zeta_avg"] = fluxsurface_integral(xr.dot(ds.B, ds.e_zeta, dim="xyz")) / (4 * np.pi**2)


@register(
    requirements=("B_zeta_avg", "mu0"),
    integration=("theta", "zeta"),
    attrs=dict(long_name="poloidal current", symbol=r"I_{pol}"),
)
def I_pol(ds: xr.Dataset):
    ds["I_pol"] = ds.B_zeta_avg * 2 * np.pi / ds.mu0
    ds["I_pol"] = ds.I_pol.sel(rho=0, method="nearest") - ds.I_pol
    if not np.isclose(ds.rho.sel(rho=0, method="nearest"), 0):
        logging.warning(
            f"Computation of `I_pol` uses `rho={ds.rho[0].item():e}` instead of the magnetic axis."
        )


@register(
    requirements=("mu0", "gamma", "mod_B", "p", "Jac"),
    integration=("rho", "theta", "zeta"),
    attrs=dict(
        long_name="total MHD energy",
        symbol=r"W_{MHD}",
    ),
)
def W_MHD(ds: xr.Dataset):
    ds["W_MHD"] = volume_integral((0.5 * ds.mod_B**2 + (ds.gamma - 1) * ds.mu0 * ds.p) * ds.Jac)
