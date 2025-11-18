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
    poloidal_integral,
    toroidal_integral,
    fluxsurface_integral,
    volume_integral,
    rtz_directions,
    rtz_symbols,
    rtz_variables,
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


def _profile_factory(var, evalvar, deriv, long_name, symbol):
    """Factory function for profile quantities."""

    @register(
        quantities=var,
        attrs=dict(long_name=long_name, symbol=symbol),
    )
    def _profile(ds: xr.Dataset, state: State):
        if "rho" not in ds:
            raise KeyError("Evaluation of profiles requires the radial coordinate 'rho'.")
        if ds.rho.dims == ("rad",):
            ds[var] = ("rad", state.evaluate_profile(evalvar, ds.rho, deriv=deriv))
        else:
            rho = ds.rho.data.flatten()
            output = state.evaluate_profile(evalvar, rho, deriv=deriv)
            ds[var] = (ds.rho.dims, output.reshape(ds.rho.shape))

    return _profile


# generate functions from factory function
for var, name, symbol in [
    ("iota", "rotational transform", r"\iota"),
    ("p", "pressure", r"p"),
    ("chi", "poloidal magnetic flux", r"\chi"),
    ("Phi", "toroidal magnetic flux", r"\Phi"),
]:
    globals()[var] = _profile_factory(var, var, 0, name, symbol)
    globals()[f"d{var}_dr"] = _profile_factory(
        f"d{var}_dr", var, 1, f"{name} gradient", f"\\frac{{d{symbol}}}{{d\\rho}}"
    )
    globals()[f"d{var}_drr"] = _profile_factory(
        f"d{var}_drr", var, 2, f"{name} curvature", f"\\frac{{d^2{symbol}}}{{d\\rho^2}}"
    )


@register(
    attrs=dict(long_name="toroidal magnetic flux at the edge", symbol=r"\Phi_0"),
)
def Phi_edge(ds: xr.Dataset, state: State):
    ds["Phi_edge"] = ((), state.evaluate_profile("Phi", [1.0])[0])


# === base ============================================================================= #


def _base_factory(var, long_name, symbol):
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
    def _base(ds: xr.Dataset, state: State):
        if "rho" not in ds or "theta" not in ds or "zeta" not in ds:
            raise KeyError(
                "Evaluation of base variables requires 'rho', 'theta', 'zeta' to be defined."
            )

        # mesh in logical coordinates (rho, theta, zeta) -> rho(rad), theta(pol), zeta(tor)
        if ds.rho.dims == ("rad",) and ds.theta.dims == ("pol",) and ds.zeta.dims == ("tor",):
            outputs = state.evaluate_base_tens_all(var, ds.rho, ds.theta, ds.zeta)
            for key, value in zip(_base.quantities, outputs):
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
            for key, value in zip(_base.quantities, outputs):
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
            for key, value in zip(_base.quantities, outputs):
                value = value.reshape(output_shape)
                ds[key] = (output_dims, value)

    return _base


# generate functions from factory function
for var, long_name, symbol in [
    ("X1", "first reference coordinate", r"X^1"),
    ("X2", "second reference coordinate", r"X^2"),
    ("LA", "straight field line potential", r"\lambda"),
]:
    globals()[var] = _base_factory(var, long_name, symbol)


@register(
    attrs=dict(long_name="number of field periods", symbol=r"N_\text{FP}"),
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
def _hmap(ds: xr.Dataset, state: State):
    X1, X2, zeta = xr.broadcast(ds.X1, ds.X2, ds.zeta)
    outputs = state.evaluate_hmap_only(*[v.values.flatten() for v in (X1, X2, zeta)])
    for key, value in zip(_hmap.quantities, outputs):
        ds[key] = (
            ("xyz", *X1.dims),
            value.reshape(3, *X1.shape),
        )


@register(
    quantities=("k_q1q1", "k_q1q2", "k_q1q3", "k_q2q2", "k_q2q3", "k_q3q3"),
    requirements=("xyz", "X1", "X2", "zeta"),
    attrs={
        f"k_q{i}q{j}": dict(
            long_name=f"q{i}-q{j} reference curvature vector", symbol=f"k_{{q^{i}q^{j}}}"
        )
        for i, j in ("11", "12", "13", "22", "23", "33")
    },
)
def _hmap_derivs(ds: xr.Dataset, state: State):
    X1, X2, zeta = xr.broadcast(ds.X1, ds.X2, ds.zeta)
    outputs = state.evaluate_hmap_derivs(
        X1=X1.values.flatten(), X2=X2.values.flatten(), zeta=zeta.values.flatten()
    )
    for key, value in zip(_hmap_derivs.quantities, outputs):
        ds[key] = (
            ("xyz", *X1.dims),
            value.reshape(3, *X1.shape),
        )


def _k_ab_factory(a, b):
    r"""Factory function for logical curvature vectors.

    The curvature vector is computed in cartesian space, with

    k_{\alpha\beta} = \frac{\partial}{\partial \alpha} \left(\frac{\partial x}{\partial \beta}\right)
                    = \sum_i \frac{\partial^2 q^i}{\partial \alpha \partial \beta} \mathbf{e}_q^i
                    + \sum_{ij} \frac{\partial q^i}{\partial \alpha} \frac{\partial q^j}{\partial \beta} \mathbf{k}_{q^iq^j}

    Sums in $i,j$ are in $1,2,3$,  and $\alpha$ and $\beta$ can be chosen as $\rho$, $\vartheta$ or $\zeta$.
    Note that since $q^3=\zeta$, some terms become active only if $\alpha$ and/or $\beta$ are equal $\zeta$.
    """

    @register(
        quantities=f"k_{a}{b}",
        requirements=sum(
            ([f"dX{i}_d{a}{b}", f"dX{i}_d{a}", f"dX{i}_d{b}"] for i in "12"), start=[]
        )
        + ["e_q1", "e_q2", "k_q1q1", "k_q1q2", "k_q2q2"]
        + (["k_q1q3", "k_q2q3"] if a == "z" or b == "z" else [])
        + (["k_q3q3"] if a == "z" and b == "z" else []),
        attrs=dict(
            long_name=f"{a}{b} logical curvature vector",
            symbol=rf"\mathbf{{k}}_{{{rtz_symbols[a] + rtz_symbols[b]}}}",
        ),
    )
    def _k_ab(ds: xr.Dataset):
        da = ds[f"dX1_d{a}{b}"] * ds.e_q1 + ds[f"dX2_d{a}{b}"] * ds.e_q2
        da += ds[f"dX1_d{a}"] * ds[f"dX1_d{b}"] * ds.k_q1q1
        da += (
            ds[f"dX1_d{a}"] * ds[f"dX2_d{b}"] + ds[f"dX2_d{a}"] * ds[f"dX1_d{b}"]
        ) * ds.k_q1q2
        da += ds[f"dX2_d{a}"] * ds[f"dX2_d{b}"] * ds.k_q2q2
        if b == "z":
            da += ds[f"dX1_d{a}"] * ds.k_q1q3 + ds[f"dX2_d{a}"] * ds.k_q2q3
        if a == "z":
            da += ds[f"dX1_d{b}"] * ds.k_q1q3 + ds[f"dX2_d{b}"] * ds.k_q2q3
        if a == "z" and b == "z":
            da += ds.k_q3q3
        ds[f"k_{a}{b}"] = da

    return _k_ab


# generate functions from factory function
for a, b in ["rr", "rt", "rz", "tt", "tz", "zz"]:
    globals()[f"k_{a}{b}"] = _k_ab_factory(a, b)


@register(
    requirements=("normal", "k_tt"),
    attrs=dict(
        long_name="poloidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta\theta}",
    ),
)
def II_tt(ds: xr.Dataset):
    ds["II_tt"] = xr.dot(ds.normal, ds.k_tt, dim="xyz")


@register(
    requirements=("normal", "k_tz"),
    attrs=dict(
        long_name="poloidal-toroidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta\zeta}",
    ),
)
def II_tz(ds: xr.Dataset):
    ds["II_tz"] = xr.dot(ds.normal, ds.k_tz, dim="xyz")


@register(
    requirements=("normal", "k_zz"),
    attrs=dict(
        long_name="toroidal component of the second fundamental form",
        symbol=r"\mathrm{II}_{\zeta\zeta}",
    ),
)
def II_zz(ds: xr.Dataset):
    ds["II_zz"] = xr.dot(ds.normal, ds.k_zz, dim="xyz")


# === metric =========================================================================== #


def _g_ij_factory(i, j):
    """Factory function for metric tensor components."""

    e_i = f"e_{rtz_variables[i]}"
    e_j = f"e_{rtz_variables[j]}"

    @register(
        quantities=f"g_{i}{j}",
        requirements={e_i, e_j},
        attrs=dict(
            long_name=f"{i}{j} component of the metric tensor",
            symbol=rf"g_{{{rtz_symbols[i] + rtz_symbols[j]}}}",
        ),
    )
    def _g_ij(ds: xr.Dataset):
        ds[f"g_{i}{j}"] = xr.dot(ds[e_i], ds[e_j], dim="xyz")

    return _g_ij


# generate functions from factory function
for i, j in ["rr", "rt", "rz", "tt", "tz", "zz"]:
    globals()[f"g_{i}{j}"] = _g_ij_factory(i, j)


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
def _metric(ds: xr.Dataset, state: State):
    outputs = state.evaluate_metric_derivs(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in _metric.requirements]
    )
    for key, value in zip(_metric.quantities, outputs):
        ds[key] = (
            ds.X1.dims,
            value.reshape(ds.X1.shape),
        )


# === jacobian determinant ============================================================= #


@register(
    requirements=(
        "e_q1",
        "e_q2",
        "e_q3",
    ),
    attrs={
        "Jac_h": dict(long_name="reference Jacobian determinant", symbol=r"\mathcal{J}_h"),
    },
)
def Jac_h(ds: xr.Dataset):
    ds["Jac_h"] = xr.dot(ds.e_q1, xr.cross(ds.e_q2, ds.e_q3, dim="xyz"), dim="xyz")


@register(
    quantities=[f"dJac_h_d{i}" for i in "rtz"],
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
def _Jac_h_derivs(ds: xr.Dataset, state: State):
    outputs = state.evaluate_jac_h_derivs(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in _Jac_h_derivs.requirements]
    )
    for key, value in zip(_Jac_h_derivs.quantities, outputs):
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
    quantities=[f"dJac{suf}_d{i}" for suf in ["", "_l"] for i in "rtz"],
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
def _Jac_derivs(ds: xr.Dataset):
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
    requirements=("grad_rho",),
    attrs=dict(long_name="surface normal", symbol=r"\mathbf{n}"),
)
def normal(ds: xr.Dataset):
    ds["normal"] = ds.grad_rho / np.sqrt(xr.dot(ds.grad_rho, ds.grad_rho, dim="xyz"))


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
        B_contra_t=dict(
            long_name="poloidal component of the magnetic field", symbol=r"B^\theta"
        ),
        B_contra_z=dict(
            long_name="toroidal component of the magnetic field", symbol=r"B^\zeta"
        ),
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
def _dB(ds: xr.Dataset):
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
    requirements=(
        "B_contra_t",
        "B_contra_z",
        "dB_contra_t_dr",
        "dB_contra_z_dr",
        "e_theta",
        "e_zeta",
        "k_rt",
        "k_rz",
    ),
    attrs=dict(
        long_name="radial derivative of the magnetic field",
        symbol=r"\frac{\partial \mathbf{B}}{\partial \rho}",
    ),
)
def dB_dr(ds: xr.Dataset):
    ds["dB_dr"] = (
        ds.dB_contra_t_dr * ds.e_theta
        + ds.B_contra_t * ds.k_rt
        + ds.dB_contra_z_dr * ds.e_zeta
        + ds.B_contra_z * ds.k_rz
    )


@register(
    requirements=(
        "B_contra_t",
        "B_contra_z",
        "dB_contra_t_dt",
        "dB_contra_z_dt",
        "e_theta",
        "e_zeta",
        "k_tt",
        "k_tz",
    ),
    attrs=dict(
        long_name="poloidal derivative of the magnetic field",
        symbol=r"\frac{\partial \mathbf{B}}{\partial \theta}",
    ),
)
def dB_dt(ds: xr.Dataset):
    ds["dB_dt"] = (
        ds.dB_contra_t_dt * ds.e_theta
        + ds.B_contra_t * ds.k_tt
        + ds.dB_contra_z_dt * ds.e_zeta
        + ds.B_contra_z * ds.k_tz
    )


@register(
    requirements=(
        "B_contra_t",
        "B_contra_z",
        "dB_contra_t_dz",
        "dB_contra_z_dz",
        "e_theta",
        "e_zeta",
        "k_tz",
        "k_zz",
    ),
    attrs=dict(
        long_name="toroidal derivative of the magnetic field",
        symbol=r"\frac{\partial \mathbf{B}}{\partial \zeta}",
    ),
)
def dB_dz(ds: xr.Dataset):
    ds["dB_dz"] = (
        ds.dB_contra_t_dz * ds.e_theta
        + ds.B_contra_t * ds.k_tz
        + ds.dB_contra_z_dz * ds.e_zeta
        + ds.B_contra_z * ds.k_zz
    )


def _dmod_B_factory(a):
    @register(
        quantities=[f"dmod_B_d{a}"],
        requirements=[
            "mod_B",
            "B_contra_t",
            "B_contra_z",
            f"dB_contra_t_d{a}",
            f"dB_contra_z_d{a}",
            "g_tt",
            "g_tz",
            "g_zz",
        ]
        + [f"dg_{ij}_d{a}" for ij in ["tt", "tz", "zz"]],
        attrs=dict(
            long_name=derivative_name_smart("modulus of the magnetic field", a),
            symbol=latex_partial(r"\left|\mathbf{B}\right|", a),
        ),
    )
    def _dmod_B(ds: xr.Dataset):
        dmod_B2_da = (
            2 * ds[f"dB_contra_t_d{a}"] * ds.B_contra_t * ds.g_tt
            + ds.B_contra_t**2 * ds[f"dg_tt_d{a}"]
            + 2 * ds[f"dB_contra_z_d{a}"] * ds.B_contra_t * ds.g_tz
            + 2 * ds[f"dB_contra_t_d{a}"] * ds.B_contra_z * ds.g_tz
            + 2 * ds.B_contra_t * ds.B_contra_z * ds[f"dg_tz_d{a}"]
            + 2 * ds[f"dB_contra_z_d{a}"] * ds.B_contra_z * ds.g_zz
            + ds.B_contra_z**2 * ds[f"dg_zz_d{a}"]
        )
        ds[f"dmod_B_d{a}"] = dmod_B2_da / (2 * ds.mod_B)


# generate functions from factory function
for a in "rtz":
    globals()[f"dmod_B_d{a}"] = _dmod_B_factory(a)


@register(
    requirements=["dmod_B_dr", "dmod_B_dt", "dmod_B_dz", "grad_rho", "grad_theta", "grad_zeta"],
    attrs=dict(
        long_name="gradient of the modulus of the magnetic field",
        symbol=r"\nabla\left|\mathbf{B}\right|",
    ),
)
def grad_mod_B(ds: xr.Dataset):
    ds["grad_mod_B"] = (
        ds.dmod_B_dr * ds.grad_rho + ds.dmod_B_dt * ds.grad_theta + ds.dmod_B_dz * ds.grad_zeta
    )


@register(
    quantities=("J", "J_contra_r", "J_contra_t", "J_contra_z"),
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


def _mod_factory(v):
    """Factory function for modulus (absolute value) quantities."""

    @register(
        quantities=f"mod_{v}",
        requirements=(v,),
        attrs=dict(
            long_name=f"modulus of the {QUANTITIES[v].attrs[v]['long_name']}",
            symbol=rf"\left|{QUANTITIES[v].attrs[v]['symbol']}\right|",
        ),
    )
    def _mod_v(ds: xr.Dataset):
        ds[f"mod_{v}"] = np.sqrt(xr.dot(ds[v], ds[v], dim="xyz"))

    return _mod_v


# generate functions from factory function
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
    globals()[v] = _mod_factory(v)


# === Straight Field Line Coordinates - Boozer ========================================= #
# dNU_B_dt and dNU_B_dz are overwritten when a boozer transform is performed!


@register(
    requirements=("B", "e_theta", "dLA_dt", "iota", "B_theta_avg", "B_zeta_avg"),
    attrs=dict(
        long_name="poloidal derivative of the Boozer potential computed from the magnetic field",
        symbol=r"\left." + latex_partial(r"\nu_B", "t") + r"\right|_\text{def.}",
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
        symbol=r"\left." + latex_partial(r"\nu_B", "z") + r"\right|_\text{def.}",
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
    Jac_B_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)  # Jac_B / Jac
    ds["e_theta_B"] = (dzB_dz * ds.e_theta - dzB_dt * ds.e_zeta) * Jac_B_Jac


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
    Jac_B_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)  # Jac_B / Jac
    ds["e_zeta_B"] = (dtB_dt * ds.e_zeta - dtB_dz * ds.e_theta) * Jac_B_Jac


def _BJ_boozer_factory(i, BJ):
    """
    Factory function for magnetic field components (if `BJ`="B")
    or current density components (if `BJ`="J") in Boozer coordinates,
    using dot product of cartesian vector of B / J, with tangent basis vector of Boozer coordinates.
    """
    e_i = f"e_{rtz_variables[i]}_B"

    if BJ == "B":
        name = "magnetic field"
    elif BJ == "J":
        name = "current density"

    @register(
        quantities=f"{BJ}_{rtz_variables[i]}_B",
        requirements={e_i, BJ},
        attrs=dict(
            long_name=rf"${rtz_symbols[i]}$ component of the {name} in Boozer coordinates",
            symbol=rf"{BJ}_{{{rtz_symbols[i]}_B}}",
        ),
    )
    def _BJ_boozer(ds: xr.Dataset):
        ds[f"{BJ}_{rtz_variables[i]}_B"] = xr.dot(ds[BJ], ds[e_i], dim="xyz")

    return _BJ_boozer


# generate functions from factory function
for BJ in ["B", "J"]:
    for i in ["r", "t", "z"]:
        globals()[f"{BJ}_{rtz_variables[i]}_B"] = _BJ_boozer_factory(i, BJ)


def _BJ_contra_Boozer_factory(i, BJ):
    """
    Factory function for contravariant components of magnetic field (if `BJ`="B")
    or current density (if `BJ`="J"),in Boozer coordinates.
    Using dot product of cartesian vector of B / J, with gradient of Boozer coordinates.
    """
    grad_i = f"grad_{rtz_variables[i]}_B"
    if BJ == "B":
        name = "magnetic field"
    elif BJ == "J":
        name = "current density"

    @register(
        quantities=f"{BJ}_contra_{i}_B",
        requirements={grad_i, BJ},
        attrs=dict(
            long_name=rf" contravariant ${rtz_symbols[i]}$ component of the {name} in Boozer coordinates",
            symbol=rf"{BJ}^{{{rtz_symbols[i]}_B}}",
        ),
    )
    def _BJ_contra_Boozer(ds: xr.Dataset):
        ds[f"{BJ}_contra_{i}_B"] = xr.dot(ds[BJ], ds[grad_i], dim="xyz")

    return _BJ_contra_Boozer


# generate functions from factory function,
# note that grad_rho = grad_rho_B, so we do not need the rho components
for BJ in ["B", "J"]:
    for i in ["t", "z"]:
        globals()[f"{BJ}_contra_{i}_B"] = _BJ_contra_Boozer_factory(i, BJ)


def _g_ij_B_factory(i, j):
    """Factory function for metric tensor components in Boozer coordinates."""
    e_i = f"e_{rtz_variables[i]}_B"
    e_j = f"e_{rtz_variables[j]}_B"

    @register(
        quantities=f"g_{i}{j}_B",
        requirements={e_i, e_j},
        attrs=dict(
            long_name=f"{i}{j} component of the metric tensor in Boozer coordinates",
            symbol=rf"g_{{{rtz_symbols[i]}_B {rtz_symbols[j]}_B}}",
        ),
    )
    def _g_ij_B(ds: xr.Dataset):
        ds[f"g_{i}{j}_B"] = xr.dot(ds[e_i], ds[e_j], dim="xyz")

    return _g_ij_B


# generate functions from factory function
for i, j in ["rr", "rt", "rz", "tt", "tz", "zz"]:
    globals()[f"g_{i}{j}_B"] = _g_ij_B_factory(i, j)


@register(
    quantities=("k_tt_B", "k_tz_B", "k_zz_B"),
    requirements=["iota", "e_theta", "e_zeta", "k_tt", "k_tz", "k_zz"]
    + sum([[f"dLA_d{ij}", f"dNU_B_d{ij}"] for ij in ("t", "z", "tt", "tz", "zz")], start=[]),
    attrs={
        f"k_{a}{b}_B": dict(
            long_name=f"{a}{b} boozer curvature vector",
            symbol=rf"\mathbf{{k}}_{{{rtz_symbols[a]}_B {rtz_symbols[b]}_B}}",
        )
        for a, b in ["tt", "tz", "zz"]
    },
)
def _k_ij_B(ds: xr.Dataset):
    r"""Factory function for curvature vectors in Boozer coordinates.

    The curvature vector is computed in cartesian space, with

    k_{\alpha\beta} = \frac{\partial}{\partial \alpha} \left(\frac{\partial x}{\partial \beta}\right)

    for the choices of $\alpha,\beta$  being $\vartheta_B,\vartheta_B$, $\vartheta_B,\zeta_B$ and $\zeta_B,\zeta_B$.
    The chain rule is applied to express the quantities in terms of the logical coordinate derivatives,
    together with the derivatives of the Boozer transform.
    """
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz

    dtB_dtt = ds.dLA_dtt + ds.iota * ds.dNU_B_dtt
    dtB_dtz = ds.dLA_dtz + ds.iota * ds.dNU_B_dtz
    dtB_dzz = ds.dLA_dzz + ds.iota * ds.dNU_B_dzz
    dzB_dtt = ds.dNU_B_dtt
    dzB_dtz = ds.dNU_B_dtz
    dzB_dzz = ds.dNU_B_dzz

    JacB_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)
    dJac_JacB_dt = dtB_dtt * dzB_dz + dtB_dt * dzB_dtz - dtB_dtz * dzB_dt - dtB_dz * dzB_dtt
    dJac_JacB_dz = dtB_dtz * dzB_dz + dtB_dt * dzB_dzz - dtB_dzz * dzB_dt - dtB_dz * dzB_dtz
    dJacB_Jac_dtB = JacB_Jac**3 * (dzB_dt * dJac_JacB_dz - dzB_dz * dJac_JacB_dt)
    dJacB_Jac_dzB = JacB_Jac**3 * (dtB_dz * dJac_JacB_dt - dtB_dt * dJac_JacB_dz)

    dt_dtB = JacB_Jac * dzB_dz
    dz_dtB = -JacB_Jac * dzB_dt
    dt_dzB = -JacB_Jac * dtB_dz
    dz_dzB = JacB_Jac * dtB_dt

    dt_dttB = dzB_dz * dJacB_Jac_dtB + JacB_Jac * (dt_dtB * dzB_dtz + dz_dtB * dzB_dzz)
    dt_dtzB = dzB_dz * dJacB_Jac_dzB + JacB_Jac * (dt_dzB * dzB_dtz + dz_dzB * dzB_dzz)
    dt_dzzB = -dtB_dz * dJacB_Jac_dzB - JacB_Jac * (dt_dzB * dtB_dtz + dz_dzB * dtB_dzz)
    dz_dttB = -dzB_dt * dJacB_Jac_dtB - JacB_Jac * (dt_dtB * dzB_dtt + dz_dtB * dzB_dtz)
    dz_dtzB = dtB_dt * dJacB_Jac_dtB + JacB_Jac * (dt_dtB * dtB_dtt + dz_dtB * dtB_dtz)
    dz_dzzB = dtB_dt * dJacB_Jac_dzB + JacB_Jac * (dt_dzB * dtB_dtt + dz_dzB * dtB_dtz)

    ds["k_tt_B"] = dt_dttB * ds.e_theta + dz_dttB * ds.e_zeta
    ds["k_tt_B"] += dt_dtB**2 * ds.k_tt + 2 * dt_dtB * dz_dtB * ds.k_tz + dz_dtB**2 * ds.k_zz

    ds["k_tz_B"] = dt_dtzB * ds.e_theta + dz_dtzB * ds.e_zeta
    ds["k_tz_B"] += dt_dtB * dt_dzB * ds.k_tt + dz_dtB * dz_dzB * ds.k_zz
    ds["k_tz_B"] += (dt_dtB * dz_dzB + dt_dzB * dz_dtB) * ds.k_tz

    ds["k_zz_B"] = dt_dzzB * ds.e_theta + dz_dzzB * ds.e_zeta
    ds["k_zz_B"] += dt_dzB**2 * ds.k_tt + 2 * dt_dzB * dz_dzB * ds.k_tz + dz_dzB**2 * ds.k_zz


@register(
    requirements=("normal", "k_tt_B"),
    attrs=dict(
        long_name="poloidal Boozer component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta_B\theta_B}",
    ),
)
def II_tt_B(ds: xr.Dataset):
    ds["II_tt_B"] = xr.dot(ds.normal, ds.k_tt_B, dim="xyz")


@register(
    requirements=["normal", "k_tz_B"],
    attrs=dict(
        long_name="poloidal-toroidal Boozer component of the second fundamental form",
        symbol=r"\mathrm{II}_{\theta_B\zeta_B}",
    ),
)
def II_tz_B(ds: xr.Dataset):
    ds["II_tz_B"] = xr.dot(ds.normal, ds.k_tz_B, dim="xyz")


@register(
    requirements=["normal", "k_zz_B"],
    attrs=dict(
        long_name="toroidal Boozer component of the second fundamental form",
        symbol=r"\mathrm{II}_{\zeta_B\zeta_B}",
    ),
)
def II_zz_B(ds: xr.Dataset):
    ds["II_zz_B"] = xr.dot(ds.normal, ds.k_zz_B, dim="xyz")


@register(
    requirements=(
        "iota",
        "diota_dr",
        "dLA_dr",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dr",
        "dNU_B_dz",
        "dNU_B_dt",
        "e_rho",
        "e_theta",
        "e_zeta",
    ),
    attrs=dict(
        long_name="radial tangent basis vector in Boozer coordinates",
        symbol=r"\mathbf{e}_{\rho_B}",
    ),
)
def e_rho_B(ds: xr.Dataset):
    dtB_dr = ds.dLA_dr + ds.diota_dr * ds.NU_B + ds.iota * ds.dNU_B_dr
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dr = ds.dNU_B_dr
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    Jac_B_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)  # Jac_B / Jac
    # dr_drB = 1
    dt_drB = Jac_B_Jac * (dtB_dz * dzB_dr - dzB_dz * dtB_dr)
    dz_drB = Jac_B_Jac * (dzB_dt * dtB_dr - dtB_dt * dzB_dr)
    ds["e_rho_B"] = ds.e_rho + dt_drB * ds.e_theta + dz_drB * ds.e_zeta


@register(
    requirements=(
        "iota",
        "diota_dr",
        "dLA_dr",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dr",
        "dNU_B_dz",
        "dNU_B_dt",
        "grad_rho",
        "grad_theta",
        "grad_zeta",
    ),
    attrs=dict(
        long_name="poloidal reciprocal basis vector in Boozer coordinates",
        symbol=r"\nabla\theta_B",
    ),
)
def grad_theta_B(ds: xr.Dataset):
    dtB_dr = ds.dLA_dr + ds.diota_dr * ds.NU_B + ds.iota * ds.dNU_B_dr
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    ds["grad_theta_B"] = dtB_dr * ds.grad_rho + dtB_dt * ds.grad_theta + dtB_dz * ds.grad_zeta


@register(
    requirements=(
        "dNU_B_dr",
        "dNU_B_dt",
        "dNU_B_dz",
        "grad_rho",
        "grad_theta",
        "grad_zeta",
    ),
    attrs=dict(
        long_name="toroidal reciprocal basis vector in Boozer coordinates",
        symbol=r"\nabla\zeta_B",
    ),
)
def grad_zeta_B(ds: xr.Dataset):
    dzB_dr = ds.dNU_B_dr
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    ds["grad_zeta_B"] = dzB_dr * ds.grad_rho + dzB_dt * ds.grad_theta + dzB_dz * ds.grad_zeta


@register(
    requirements=(
        "iota",
        "diota_dr",
        "dLA_dr",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dr",
        "dNU_B_dt",
        "dNU_B_dz",
        "dmod_B_dr",
        "dmod_B_dt",
        "dmod_B_dz",
    ),
    attrs=dict(
        long_name="radial Boozer derivative of the modulus of the magnetic field",
        symbol=r"\frac{\partial\left|\mathbf{B}\right|}{\partial \rho_B}",
    ),
)
def dmod_B_dr_B(ds: xr.Dataset):
    dtB_dr = ds.dLA_dr + ds.diota_dr * ds.NU_B + ds.iota * ds.dNU_B_dr
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dr = ds.dNU_B_dr
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    Jac_B_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)  # Jac_B / Jac
    dt_drB = Jac_B_Jac * (dtB_dz * dzB_dr - dzB_dz * dtB_dr)
    dz_drB = Jac_B_Jac * (dzB_dt * dtB_dr - dtB_dt * dzB_dr)
    ds["dmod_B_dr_B"] = ds.dmod_B_dr + dt_drB * ds.dmod_B_dt + dz_drB * ds.dmod_B_dz


@register(
    requirements=(
        "iota",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dt",
        "dNU_B_dz",
        "dmod_B_dt",
        "dmod_B_dz",
    ),
    attrs=dict(
        long_name="poloidal Boozer derivative of the modulus of the magnetic field",
        symbol=r"\frac{\partial\left|\mathbf{B}\right|}{\partial \theta_B}",
    ),
)
def dmod_B_dt_B(ds: xr.Dataset):
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    Jac_B_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)  # Jac_B / Jac
    dt_dtB = Jac_B_Jac * dzB_dz
    dz_dtB = -Jac_B_Jac * dzB_dt
    ds["dmod_B_dt_B"] = dt_dtB * ds.dmod_B_dt + dz_dtB * ds.dmod_B_dz


@register(
    requirements=(
        "iota",
        "dLA_dt",
        "dLA_dz",
        "dNU_B_dt",
        "dNU_B_dz",
        "dmod_B_dt",
        "dmod_B_dz",
    ),
    attrs=dict(
        long_name="toroidal Boozer derivative of the modulus of the magnetic field",
        symbol=r"\frac{\partial\left|\mathbf{B}\right|}{\partial \zeta_B}",
    ),
)
def dmod_B_dz_B(ds: xr.Dataset):
    dtB_dt = 1 + ds.dLA_dt + ds.iota * ds.dNU_B_dt
    dtB_dz = ds.dLA_dz + ds.iota * ds.dNU_B_dz
    dzB_dt = ds.dNU_B_dt
    dzB_dz = 1 + ds.dNU_B_dz
    Jac_B_Jac = 1 / (dtB_dt * dzB_dz - dtB_dz * dzB_dt)  # Jac_B / Jac
    dt_dzB = -Jac_B_Jac * dtB_dz
    dz_dzB = Jac_B_Jac * dtB_dt
    ds["dmod_B_dz_B"] = dt_dzB * ds.dmod_B_dt + dz_dzB * ds.dmod_B_dz


# === integrals ======================================================================== #
# --- geometry --- #


@register(
    requirements=("Jac",),
    integration=("rho", "theta", "zeta"),
    attrs=dict(long_name="volume", symbol=r"V"),
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
    """
    d/dPhi_n = dr/dPhi_n * d/dr = Phi_0 / dPhi_dr * d/dr
    """
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
    """
    d/dPhi_n = dr/dPhi_n * d/dr = Phi_0 / dPhi_dr * d/dr
    d/dr 1/dPhi_dr = -1/dPhi_dr**2 * dPhi_drr
    """
    ds["dV_dPhi_n2"] = (
        fluxsurface_integral(ds.dJac_dr) * (ds.Phi_edge / ds.dPhi_dr) ** 2
        - fluxsurface_integral(ds.Jac) * ds.Phi_edge**2 / ds.dPhi_dr**3 * ds.dPhi_drr
    )


@register(
    requirements=("e_theta", "e_zeta"),
    attrs=dict(long_name="differential area element", symbol=r"dA"),
)
def dA(ds: xr.Dataset):
    dA = xr.cross(ds.e_theta, ds.e_zeta, dim="xyz")
    ds["dA"] = np.sqrt(xr.dot(dA, dA, dim="xyz"))


@register(
    attrs=dict(long_name="surface area", symbol=r"A_\text{surface}"),
)
def A_surface(ds: xr.Dataset, state: State):
    aux = state.evaluate("dA", rho=1.0, theta="int", zeta="int")
    ds["A_surface"] = fluxsurface_integral(aux.dA).item()


@register(
    attrs=dict(long_name="length of the magnetic axis", symbol=r"L_\text{axis}"),
)
def L_axis(ds: xr.Dataset, state: State):
    aux = state.evaluate("e_zeta", rho=0.0, theta=0.0, zeta="int")
    dL = np.sqrt(xr.dot(aux.e_zeta, aux.e_zeta, dim="xyz"))
    ds["L_axis"] = toroidal_integral(dL).item()


@register(
    requirements=("L_axis",),
    attrs=dict(long_name="effective major radius", symbol=r"r_\text{major,eff}"),
)
def r_major(ds: xr.Dataset):
    ds["r_major"] = ds.L_axis / (2 * np.pi)


@register(
    requirements=("V", "L_axis"),
    attrs=dict(long_name="effective minor radius", symbol=r"r_\text{minor,eff}"),
)
def r_minor(ds: xr.Dataset):
    ds["r_minor"] = np.sqrt(ds.V / (np.pi * ds.L_axis))


@register(
    requirements=("r_major", "r_minor"),
    attrs=dict(long_name="effective aspect ratio", symbol=r"a_\text{eff}"),
)
def aspect_ratio(ds: xr.Dataset):
    ds["aspect_ratio"] = ds.r_major / ds.r_minor


@register(
    requirements=("A_surface", "V", "L_axis"),
    attrs=dict(long_name="effective elongation", symbol=r"E_\text{eff}"),
)
def elongation(ds: xr.Dataset):
    from scipy.optimize import newton
    from gvec.util import ellipse_circumference_factor as ecf

    C = (ds.A_surface / 2 / np.sqrt(np.pi * ds.V * ds.L_axis)).item()
    ds["elongation"] = newton(lambda e: ecf(e) - C, 2)


# --- profiles --- #


@register(
    requirements=("iota",),
    integration=("rho",),
    attrs=dict(long_name="average rotational transform", symbol=r"\overline{\iota}"),
)
def iota_avg(ds: xr.Dataset):
    ds["iota_avg"] = radial_integral(ds.iota)


@register(
    requirements=("iota",),
    integration=("rho",),
    attrs=dict(
        long_name="rotational transform averaged over rho^2", symbol=r"\overline{\iota}_2"
    ),
)
def iota_avg2(ds: xr.Dataset):
    ds["iota_avg2"] = radial_integral(2 * ds.rho * ds.iota)


@register(
    requirements=("mu0", "dPhi_dr", "Jac", "g_tt"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="factor to the toroidal current contribution to the rotational transform",
        description="iota = iota_0 + iota_curr, iota_curr = I_tor * iota_curr_0",
        symbol=r"\iota_{\text{curr},0}",
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
        symbol=r"\iota_\text{curr}",
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
    requirements=("shear",),
    attrs=dict(long_name="average global magnetic shear", symbol=r"\overline{s_g}"),
)
def shear_avg(ds: xr.Dataset):
    ds["shear_avg"] = radial_integral(ds.shear)


@register(
    requirements=("shear",),
    attrs=dict(
        long_name="global magnetic shear averaged over rho^2", symbol=r"\overline{s_g}_2"
    ),
)
def shear_avg2(ds: xr.Dataset):
    ds["shear_avg2"] = radial_integral(2 * ds.rho * ds.shear)


@register(
    requirements=("B", "e_theta"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="flux-surface averaged poloidal magnetic field", symbol=r"\overline{B_\theta}"
    ),
)
def B_theta_avg(ds: xr.Dataset):
    ds["B_theta_avg"] = fluxsurface_integral(xr.dot(ds.B, ds.e_theta, dim="xyz")) / (
        4 * np.pi**2
    )


@register(
    requirements=("B", "dB_dr", "e_theta", "k_rt"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="derivative of the flux-surface averaged poloidal magnetic field",
        symbol=r"\frac{d\overline{B_\theta}}{d\rho}",
    ),
)
def dB_theta_avg_dr(ds: xr.Dataset):
    dB_t_dr = xr.dot(ds.dB_dr, ds.e_theta, dim="xyz") + xr.dot(ds.B, ds.k_rt, dim="xyz")
    ds["dB_theta_avg_dr"] = fluxsurface_integral(dB_t_dr) / (4 * np.pi**2)


@register(
    requirements=("B_theta_avg", "mu0"),
    attrs=dict(long_name="toroidal current enclosed by flux surface", symbol=r"I_\text{tor}"),
)
def I_tor(ds: xr.Dataset):
    ds["I_tor"] = ds.B_theta_avg * 2 * np.pi / ds.mu0


@register(
    requirements=("dB_theta_avg_dr", "mu0"),
    attrs=dict(
        long_name="derivative of the toroidal current enclosed by the flux surface",
        symbol=r"\frac{dI_\text{tor}}{d\rho}",
    ),
)
def dI_tor_dr(ds: xr.Dataset):
    ds["dI_tor_dr"] = ds.dB_theta_avg_dr * 2 * np.pi / ds.mu0


@register(
    requirements=("B", "e_zeta"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="flux-surface averaged toroidal magnetic field", symbol=r"\overline{B_\zeta}"
    ),
)
def B_zeta_avg(ds: xr.Dataset):
    ds["B_zeta_avg"] = fluxsurface_integral(xr.dot(ds.B, ds.e_zeta, dim="xyz")) / (4 * np.pi**2)


@register(
    requirements=("B_zeta_avg", "mu0"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="poloidal current, relative to the magnetic axis", symbol=r"I_\text{pol}"
    ),
)
def I_pol(ds: xr.Dataset):
    ds["I_pol"] = ds.B_zeta_avg * 2 * np.pi / ds.mu0
    ds["I_pol"] = ds.I_pol.sel(rho=0, method="nearest") - ds.I_pol
    if not np.isclose(ds.rho.sel(rho=0, method="nearest"), 0):
        logging.warning(
            f"Computation of `I_pol` uses `rho={ds.rho[0].item():e}` instead of the magnetic axis."
        )


# --- other --- #


@register(
    requirements=("mu0", "gamma", "mod_B", "p", "Jac"),
    integration=("rho", "theta", "zeta"),
    attrs=dict(
        long_name="total MHD energy",
        symbol=r"W_\text{MHD}",
    ),
)
def W_MHD(ds: xr.Dataset):
    ds["W_MHD"] = volume_integral((0.5 * ds.mod_B**2 + (ds.gamma - 1) * ds.mu0 * ds.p) * ds.Jac)


@register(
    requirements=("p", "mod_B", "Jac", "V", "mu0"),
    integration=("rho", "theta", "zeta"),
    attrs=dict(
        long_name="volume-averaged plasma beta",
        symbol=r"\overline{\beta}",
    ),
)
def beta_avg(ds: xr.Dataset):
    beta = 2 * ds.mu0 * ds.p / ds.mod_B**2
    ds["beta_avg"] = volume_integral(beta * ds.Jac) / ds.V


@register(
    attrs=dict(
        long_name="vacuum magnetic well depth",
        symbol=r"d_\text{well}",
    ),
)
def vacuum_magnetic_well_depth(ds: xr.Dataset, state: State):
    aux = state.evaluate("dV_dPhi_n", rho=[1e-4, 1.0], theta="int", zeta="int")
    Vp_edge = aux.dV_dPhi_n.isel(rad=1)
    Vp_axis = aux.dV_dPhi_n.isel(rad=0)
    ds["vacuum_magnetic_well_depth"] = (Vp_axis - Vp_edge) / Vp_axis


@register(
    quantities=("D_Merc", "D_Merc_Shear", "D_Merc_Curr", "D_Merc_Well", "D_Merc_Geod"),
    requirements=(
        "dPhi_dr",
        "dPhi_drr",
        "diota_dr",
        "dp_dr",
        "Phi",
        "chi",
        "Jac",
        "dJac_dr",
        "grad_rho",
        "dB_theta_avg_dr",
        "mu0",
        "J",
        "B",
        "mod_B",
    ),
    attrs={
        "D_Merc": dict(
            long_name="Mercier criterion",
            symbol=r"D_\text{Merc}",
        ),
        "D_Merc_Shear": dict(
            long_name="Shear contribution to the Mercier criterion",
            symbol=r"D_\text{M,Shear}",
        ),
        "D_Merc_Curr": dict(
            long_name="Current contribution to the Mercier criterion",
            symbol=r"D_\text{M,Curr}",
        ),
        "D_Merc_Well": dict(
            long_name="Magnetic well contribution to the Mercier criterion",
            symbol=r"D_\text{M,Well}",
        ),
        "D_Merc_Geod": dict(
            long_name="Geodesic contribution to the Mercier criterion",
            symbol=r"D_\text{M,Geod}",
        ),
    },
)
def D_Merc(ds: xr.Dataset):
    twopi = 2 * np.pi
    diota_dPhi = ds.diota_dr / ds.dPhi_dr
    dp_dPhi = ds.dp_dr / ds.dPhi_dr
    s_chi = np.sign(ds.chi)
    s_Phi = np.sign(ds.Phi)
    dV_dr = fluxsurface_integral(ds.Jac)
    dV_drr = fluxsurface_integral(ds.dJac_dr)
    d2V_dPhi2 = dV_drr / ds.dPhi_dr**2 - dV_dr * ds.dPhi_drr / ds.dPhi_dr**3
    ngradPhi = np.sqrt(xr.dot(ds.grad_rho, ds.grad_rho, dim="xyz")) * ds.dPhi_dr
    dB_theta_avg_dPhi = ds.dB_theta_avg_dr / ds.dPhi_dr
    # dS = xr.cross(ds.e_theta, ds.e_zeta, dim="xyz")
    dS = np.sqrt(xr.dot(ds.grad_rho, ds.grad_rho, dim="xyz")) * ds.Jac
    JBint = fluxsurface_integral(dS * ds.mu0 * xr.dot(ds.J, ds.B, dim="xyz") / ngradPhi**3)
    B2int = fluxsurface_integral(dS * ds.mod_B**2 / ngradPhi**3)
    Bi2int = fluxsurface_integral(dS / ds.mod_B**2 / ngradPhi)
    JB2int = fluxsurface_integral(
        dS * (ds.mu0 * xr.dot(ds.J, ds.B, dim="xyz") / ds.mod_B) ** 2 / ngradPhi**3
    )
    ds["D_Merc_Shear"] = 1 / (16 * np.pi**2) * diota_dPhi**2
    ds["D_Merc_Curr"] = -s_chi / twopi**4 * diota_dPhi * (JBint - dB_theta_avg_dPhi * B2int)
    ds["D_Merc_Well"] = (
        ds.mu0 / twopi**6 * dp_dPhi * (s_Phi * d2V_dPhi2 - ds.mu0 * dp_dPhi * Bi2int) * B2int
    )
    ds["D_Merc_Geod"] = 1 / twopi**6 * (JBint**2 - B2int * JB2int)
    ds["D_Merc"] = (
        ds["D_Merc_Shear"] + ds["D_Merc_Curr"] + ds["D_Merc_Well"] + ds["D_Merc_Geod"]
    )


@register(
    requirements=("mod_B",),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="mirror ratio",
        symbol=r"\Delta_\text{mirror}",
    ),
)
def mirror_ratio(ds: xr.Dataset):
    r"""Compute the mirror ratio of the magnetic field strength on a flux surface.

    The mirror ratio is defined as

    R_mirror = (B_max - B_min) / (B_max + B_min)

    where B_max and B_min are the maximum and minimum values of the magnetic field strength
    on a given flux surface.
    """
    B_max = ds.mod_B.max(dim=("pol", "tor"))
    B_min = ds.mod_B.min(dim=("pol", "tor"))
    ds["mirror_ratio"] = (B_max - B_min) / (B_max + B_min)


@register(
    requirements=("mod_B", "dB_dr", "dB_dt", "dB_dz", "grad_rho", "grad_theta", "grad_zeta"),
    attrs=dict(
        long_name="magnetic gradient scale length",
        symbol=r"L_{\nabla\mathbf{B}}",
    ),
)
def L_gradB(ds: xr.Dataset):
    r"""Compute the magnetic gradient scale length.

    The magnetic gradient scale length is defined as

    L_gradB = sqrt(2) |B| / ||grad B||

    where ||grad B|| is the frobenius norm of the gradient of the magnetic field.
    Details can be found in Kappel et al. PPCF 66 (2024) 025018 DOI:10.1088/1361-6587/ad1a3e.
    """
    gradB = {}  # 3x3 tensor
    for i in "xyz":
        for j in "xyz":
            gradB[i, j] = xr.zeros_like(ds.mod_B)
            for k in ("rho", "theta", "zeta"):
                gradB[i, j] += ds[f"dB_d{k[0]}"].sel(xyz=i) * ds[f"grad_{k}"].sel(xyz=j)
    # frobenius norm
    gradB_normF = np.sqrt(sum(gradB[i, j] ** 2 for i in "xyz" for j in "xyz"))
    ds["L_gradB"] = np.sqrt(2) * ds.mod_B / gradB_normF
