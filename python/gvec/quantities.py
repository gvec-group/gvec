"""
This module defines various quantities and their computation functions for the GVEC package.

The module contains functions that compute different physical quantities such as iota and pressure profiles, coordinate mappings and their derivatives, magnetic field, current density and more.

These functions are registered with the `Evaluations` class using the `@register_compute_func()` decorator.
"""

from .post import Evaluations, State

import logging

import xarray as xr
import numpy as np

register = Evaluations.register_compute_func
rtz_symbols = {"r": r"\rho", "t": r"\theta", "z": r"\zeta"}
rtz_directions = {"r": "radial", "t": "poloidal", "z": "toroidal"}


"""
There are two different kinds of dimensional/coordinate dependencies:
    1) a quantity is defined over some dimensions in the dataset
        * e.g. mu0(), iota(rad), B(vector, rad, pol, tor)
        * sometimes we only care about some of the dimensions, e.g. B(vector, ...)
        * actually all of them could depend on some other dimension, like resolution or time
            * but in the underlying dimension would be `state` and each computation happens only once per state
        * this information is required when assigning the quantity, but not outside the compute function
    2) a quantity makes some assumptions for the dimensions of other quantities
        * e.g. I_tor(rad) needs theta_int, zeta_int, X1 needs rho,theta,zeta or rho,tz_list
        * this information is required when selecting the right compute function, based on the available dimensions
            * but its not just dimensions, e.g. theta_B, zeta_B =/= tz_list, but same compute function

=> let us not overengineer this
* we have two selection criteria: integration points & tensorproduct-tz
-> set flags for them both and change the system later if required
"""


# === helpers ========================================================================== #


def latex_partial(var, deriv):
    return rf"\frac{{\partial {var}}}{{\partial {rtz_symbols[deriv]}}}"


def latex_partial2(var, deriv1, deriv2):
    if deriv1 == deriv2:
        return rf"\frac{{\partial^2 {var}}}{{\partial {rtz_symbols[deriv1]}^2}}"
    return rf"\frac{{\partial^2 {var}}}{{\partial {rtz_symbols[deriv1]}\partial {rtz_symbols[deriv2]}}}"


def latex_partial_smart(var, deriv):
    if len(deriv) == 1:
        return latex_partial(var, deriv[0])
    elif len(deriv) == 2:
        return latex_partial2(var, deriv[0], deriv[1])
    raise TypeError(f"can only handle derivatives up to length 2, got '{deriv}'")


def derivative_name_smart(name, deriv):
    if len(deriv) == 1:
        return f"{rtz_directions[deriv[0]]} derivative of the {name}"
    elif len(deriv) == 2:
        if deriv[0] == deriv[1]:
            return f"second {rtz_directions[deriv[0]]} derivative of the {name}"
        return f"{rtz_directions[deriv[0]]}-{rtz_directions[deriv[1]]} derivative of the {name}"
    raise TypeError(f"can only handle derivatives up to length 2, got '{deriv}'")


# === special ========================================================================== #


@register(
    attrs=dict(long_name="magnetic constant", symbol=r"\mu_0"),
)
def mu0(ds: Evaluations):
    ds["mu0"] = 4 * np.pi * 1e-7


@register(
    attrs=dict(long_name="cartesian vector components", symbol=r"\mathbf{x}"),
)
def vector(ds: Evaluations):
    ds.coords["vector"] = ("vector", ["x", "y", "z"])


# === profiles ========================================================================= #


def _profile(var, evalvar, long_name, symbol):
    """Factory function for profile quantities."""

    @register(
        quantities=var,
        attrs=dict(long_name=long_name, symbol=symbol),
    )
    def profile(ds: Evaluations, state: State):
        ds[var] = ("rad", state.evaluate_profile(evalvar, ds.rho))

    return profile


for var, *args in [
    ("iota", "iota", "rotational transform", r"\iota"),
    (
        "diota_dr",
        "iota_prime",
        "rotational transform gradient",
        r"\frac{d\iota}{d\rho}",
    ),
    ("p", "p", "pressure", r"p"),
    ("dp_dr", "p_prime", "pressure gradient", r"\frac{dp}{d\rho}"),
    ("chi", "chi", "poloidal magnetic flux", r"\chi"),
    ("dchi_dr", "chi_prime", "poloidal magnetic flux gradient", r"\frac{d\chi}{d\rho}"),
    ("Phi", "Phi", "toroidal magnetic flux", r"\Phi"),
    ("dPhi_dr", "Phi_prime", "toroidal magnetic flux gradient", r"\frac{d\Phi}{d\rho}"),
    (
        "dPhi_drr",
        "Phi_2prime",
        "toroidal magnetic flux curvature",
        r"\frac{d^2\Phi}{d\rho^2}",
    ),
    ("Phi_n", "PhiNorm", "normalized toroidal magnetic flux", r"\Phi_n"),
    (
        "dPhi_n_dr",
        "PhiNorm_prime",
        "normalized toroidal magnetic flux gradient",
        r"\frac{d\Phi_n}{d\rho}",
    ),
]:
    globals()[var] = _profile(var, *args)


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
    def base(ds: Evaluations, state: State):
        if set(ds.rho.dims) != {"rad"}:
            raise ValueError(
                f"Expected 'rho' to be of dimension '(rad,)', got {ds.rho.dims}"
            )

        if set(ds.theta.dims) == {"pol"} and set(ds.zeta.dims) == {"tor"}:
            outputs = state.evaluate_base_tens_all(var, ds.rho, ds.theta, ds.zeta)
            for key, value in zip(base.quantities, outputs):
                ds[key] = (("rad", "pol", "tor"), value)

        elif set(ds.theta.dims) | set(ds.zeta.dims) == {"pol", "tor"}:
            stacked = ds[["rho", "theta", "zeta"]].stack(tz=("pol", "tor"))
            thetazeta = np.stack([stacked.theta, stacked.zeta], axis=0)
            outputs = state.evaluate_base_list_tz_all(var, ds.rho, thetazeta)
            for key, value in zip(base.quantities, outputs):
                value = xr.DataArray(value, coords=stacked.coords).unstack("tz")
                ds[key] = value

        elif set(ds.theta.dims) | set(ds.zeta.dims) == {"rad", "pol", "tor"}:
            stacked = ds[["rho", "theta", "zeta"]].stack(tz=("pol", "tor"))
            outputs = []
            for r, rho in enumerate(ds.rho.data):
                surface = stacked.isel(rad=r)
                thetazeta = np.stack([surface.theta, surface.zeta], axis=0)
                outputs.append(state.evaluate_base_list_tz_all(var, [rho], thetazeta))
            for key, value in zip(base.quantities, zip(*outputs)):
                value = xr.DataArray(
                    np.stack(value).squeeze(), coords=stacked.coords
                ).unstack("tz")
                ds[key] = value
        else:
            raise ValueError(
                f"Expected 'theta', 'zeta' to be of dimensions subset of '(rad, pol, tor)', got {ds.theta.dims}, {ds.zeta.dims}"
            )

    return base


for var, long_name, symbol in [
    ("X1", "first reference coordinate", r"X^1"),
    ("X2", "second reference coordinate", r"X^2"),
    ("LA", "straight field line potential", r"\lambda"),
    ("GB", "Boozer potential", r"G_B"),
]:
    globals()[var] = _base(var, long_name, symbol)


@register(
    attrs=dict(long_name="number of field periods", symbol=r"N_{FP}"),
)
def N_FP(ds: Evaluations, state: State):
    ds["N_FP"] = state.nfp


# === mapping ========================================================================== #


@register(
    quantities=("pos", "e_X1", "e_X2", "e_zeta3"),
    requirements=("vector", "X1", "X2", "zeta"),
    attrs=dict(
        pos=dict(long_name="position vector", symbol=r"\mathbf{x}"),
        e_X1=dict(
            long_name="first reference tangent basis vector", symbol=r"\mathbf{e}_{X^1}"
        ),
        e_X2=dict(
            long_name="second reference tangent basis vector",
            symbol=r"\mathbf{e}_{X^2}",
        ),
        e_zeta3=dict(
            long_name="toroidal reference tangent basis vector",
            symbol=r"\mathbf{e}_{\zeta^3}",
        ),
    ),
)
def hmap(ds: Evaluations, state: State):
    outputs = state.evaluate_hmap_only(
        **{
            var: ds[var].broadcast_like(ds.X1).values.flatten()
            for var in hmap.requirements
            if var != "vector"
        }
    )
    for key, value in zip(hmap.quantities, outputs):
        ds[key] = (
            ("vector", "rad", "pol", "tor"),
            value.reshape(3, ds.rho.size, ds.theta.size, ds.zeta.size),
        )


# === metric =========================================================================== #


@register(
    quantities=[
        pattern.format(ij=ij)
        for pattern in ("g_{ij}", "dg_{ij}_dr", "dg_{ij}_dt", "dg_{ij}_dz")
        for ij in ("rr", "rt", "rz", "tt", "tz", "zz")
    ],
    requirements=["X1", "X2", "zeta"]
    + [
        f"d{Xi}_d{j}"
        for j in ("r", "t", "z", "rr", "rt", "rz", "tt", "tz", "zz")
        for Xi in ("X1", "X2")
    ],
    attrs={
        f"g_{ij}": dict(
            long_name=f"{ij} component of the metric tensor", symbol=rf"g_{{{ij}}}"
        )
        for ij in ("rr", "rt", "rz", "tt", "tz", "zz")
    }
    | {
        f"dg_{ij}_d{k}": dict(
            long_name=derivative_name_smart(f"{ij} component of the metric tensor", k),
            symbol=latex_partial(f"g_{{{ij}}}", k),
        )
        for ij in ("rr", "rt", "rz", "tt", "tz", "zz")
        for k in "rtz"
    },
)
def metric(ds: Evaluations, state: State):
    outputs = state.evaluate_metric(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in metric.requirements]
    )
    for key, value in zip(metric.quantities, outputs):
        ds[key] = (
            ds.X1.dims,
            value.reshape(ds.X1.shape),
        )


# === jacobian determinant ============================================================= #


@register(
    quantities=("Jac_h", *(f"dJac_h_d{i}" for i in "rtz")),
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
        "Jac_h": dict(
            long_name="reference Jacobian determinant", symbol=r"\mathcal{J}_h"
        ),
    }
    | {
        f"dJac_h_d{i}": dict(
            long_name=derivative_name_smart("reference Jacobian determinant", i),
            symbol=latex_partial_smart(r"\mathcal{J}_h", i),
        )
        for i in "rtz"
    },
)
def Jac_h(ds: Evaluations, state: State):
    outputs = state.evaluate_jacobian(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in Jac_h.requirements]
    )
    for key, value in zip(Jac_h.quantities, outputs):
        ds[key] = (
            ds.X1.dims,
            value.reshape(ds.X1.shape),
        )


@register(
    quantities=(
        "Jac",
        "Jac_l",
        *(f"dJac{suf}_d{i}" for suf in ["", "_l"] for i in "rtz"),
    ),
    requirements=(
        "Jac_h",
        *(f"dJac_h_d{i}" for i in "rtz"),
        *(
            f"d{Q}_d{i}"
            for Q in "X1 X2".split()
            for i in "r t z rr rt rz tt tz zz".split()
        ),
    ),
    attrs={
        "Jac": dict(long_name="Jacobian determinant", symbol=r"\mathcal{J}"),
        "Jac_l": dict(
            long_name="logical Jacobian determinant", symbol=r"\mathcal{J}_l"
        ),
    }
    | {
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
def Jac(ds: Evaluations):
    ds["Jac_l"] = ds.dX1_dr * ds.dX2_dt - ds.dX1_dt * ds.dX2_dr
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
    ds["Jac"] = ds.Jac_h * ds.Jac_l
    ds["dJac_dr"] = ds.dJac_h_dr * ds.Jac_l + ds.Jac_h * ds.dJac_l_dr
    ds["dJac_dt"] = ds.dJac_h_dt * ds.Jac_l + ds.Jac_h * ds.dJac_l_dt
    ds["dJac_dz"] = ds.dJac_h_dz * ds.Jac_l + ds.Jac_h * ds.dJac_l_dz


# === straight field line coordinates - PEST =========================================== #


@register(
    requirements=("LA",),
    attrs=dict(long_name="poloidal angle in PEST coordinates", symbol=r"\theta_P"),
)
def theta_P(ds: Evaluations):
    ds["theta_P"] = ds.theta + ds.LA


@register(
    requirements=("vector", "theta_sfl", "dLA_dr", "dLA_dt", "dLA_dz"),
    attrs=dict(
        long_name="poloidal reciprocal basis vector in PEST coordinates",
        symbol=r"\nabla \theta_P",
    ),
)
def grad_theta_P(ds: Evaluations):
    ds["grad_theta_P"] = (
        ds.grad_theta * (1 + ds.dLA_dt)
        + ds.grad_rho * ds.dLA_dr
        + ds.grad_zeta * ds.dLA_dz
    )


@register(
    requirements=("Jac", "dLA_dt"),
    attrs=dict(
        long_name="Jacobian determinant in PEST coordinates", symbol=r"\mathcal{J}_P"
    ),
)
def Jac_P(ds: Evaluations):
    ds["Jac_P"] = ds.Jac / (1 + ds.dLA_dt)


# === derived ========================================================================== #


@register(
    requirements=("vector", "e_X1", "e_X2", "dX1_dr", "dX2_dr"),
    attrs=dict(long_name="radial tangent basis vector", symbol=r"\mathbf{e}_\rho"),
)
def e_rho(ds: Evaluations):
    ds["e_rho"] = ds.e_X1 * ds.dX1_dr + ds.e_X2 * ds.dX2_dr


@register(
    requirements=("vector", "e_X1", "e_X2", "dX1_dt", "dX2_dt"),
    attrs=dict(long_name="poloidal tangent basis vector", symbol=r"\mathbf{e}_\theta"),
)
def e_theta(ds: Evaluations):
    ds["e_theta"] = ds.e_X1 * ds.dX1_dt + ds.e_X2 * ds.dX2_dt


@register(
    requirements=("vector", "e_X1", "e_X2", "e_zeta3", "dX1_dz", "dX2_dz"),
    attrs=dict(long_name="toroidal tangent basis vector", symbol=r"\mathbf{e}_\zeta"),
)
def e_zeta(ds: Evaluations):
    ds["e_zeta"] = ds.e_X1 * ds.dX1_dz + ds.e_X2 * ds.dX2_dz + ds.e_zeta3


@register(
    requirements=("vector", "Jac", "e_theta", "e_zeta"),
    attrs=dict(long_name="radial reciprocal basis vector", symbol=r"\nabla\rho"),
)
def grad_rho(ds: Evaluations):
    ds["grad_rho"] = xr.cross(ds.e_theta, ds.e_zeta, dim="vector") / ds.Jac


@register(
    requirements=("vector", "Jac", "e_rho", "e_zeta"),
    attrs=dict(long_name="poloidal reciprocal basis vector", symbol=r"\nabla\theta"),
)
def grad_theta(ds: Evaluations):
    ds["grad_theta"] = xr.cross(ds.e_zeta, ds.e_rho, dim="vector") / ds.Jac


@register(
    requirements=("vector", "Jac", "e_rho", "e_theta"),
    attrs=dict(long_name="toroidal reciprocal basis vector", symbol=r"\nabla\zeta"),
)
def grad_zeta(ds: Evaluations):
    ds["grad_zeta"] = xr.cross(ds.e_rho, ds.e_theta, dim="vector") / ds.Jac


@register(
    quantities=("B", "B_contra_t", "B_contra_z"),
    requirements=(
        "vector",
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
def B(ds: Evaluations):
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
def dB(ds: Evaluations):
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
def J(ds: Evaluations):
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
    ds["J"] = (
        ds.J_contra_r * ds.e_rho
        + ds.J_contra_t * ds.e_theta
        + ds.J_contra_z * ds.e_zeta
    )


def _modulus(v):
    """Factory function for modulus (absolute value) quantities."""

    @register(
        quantities=f"mod_{v}",
        requirements=(v,),
        attrs=dict(
            long_name=f"modulus of the {Evaluations._quantities[v].attrs[v]['long_name']}",
            symbol=rf"\left|{Evaluations._quantities[v].attrs[v]['symbol']}\right|",
        ),
    )
    def mod_v(ds: Evaluations):
        ds[f"mod_{v}"] = np.sqrt(xr.dot(ds[v], ds[v], dim="vector"))

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
]:
    globals()[v] = _modulus(v)

# === Straight Field Line Coordinates - Boozer ========================================= #


@register(
    requirements=("B", "e_theta", "dLA_dt", "iota", "B_theta_avg", "B_zeta_avg"),
    attrs=dict(
        long_name="poloidal derivative of the Boozer potential per definition",
        symbol=r"\left." + latex_partial(r"G_B", "t") + r"\right|_{\text{def}}",
    ),
)
def dGB_dt_def(ds: Evaluations):
    Bt = xr.dot(ds.B, ds.e_theta, dim="vector")
    ds["dGB_dt_def"] = (Bt - ds.B_theta_avg * (1 + ds.dLA_dt)) / (
        ds.iota * ds.B_theta_avg + ds.B_zeta_avg
    )


@register(
    requirements=("B", "e_zeta", "dLA_dz", "iota", "B_theta_avg", "B_zeta_avg"),
    attrs=dict(
        long_name="toroidal derivative of the Boozer potential per definition",
        symbol=r"\left." + latex_partial(r"G_B", "z") + r"\right|_{\text{def}}",
    ),
)
def dGB_dz_def(ds: Evaluations):
    Bz = xr.dot(ds.B, ds.e_zeta, dim="vector")
    ds["dGB_dz_def"] = (Bz - ds.B_theta_avg * ds.dLA_dz - ds.B_zeta_avg) / (
        ds.iota * ds.B_theta_avg + ds.B_zeta_avg
    )


# === integrals ======================================================================== #


@register(
    requirements=("Jac",),
    integration=("rho", "theta", "zeta"),
    attrs=dict(long_name="plasma volume", symbol=r"V"),
)
def V(ds: Evaluations):
    ds["V"] = ds.volume_integral(1, Jac=True)


@register(
    requirements=("Jac", "dPhi_n_dr"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="derivative of the plasma volume w.r.t. normalized toroidal magnetic flux",
        symbol=r"\frac{dV}{d\Phi_n}",
    ),
)
def dV_dPhi_n(ds: Evaluations):
    ds["dV_dPhi_n"] = ds.fluxsurface_integral(1, Jac=True) / ds.dPhi_n_dr


@register(
    requirements=("Jac", "dJac_dr", "dPhi_n_dr"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="second derivative of the plasma volume w.r.t. normalized toroidal magnetic flux",
        symbol=r"\frac{d^2V}{d\Phi_n^2}",
    ),
)
def dV_dPhi_n2(ds: Evaluations):
    ds["dV_dPhi_n2"] = (
        ds.fluxsurface_integral(1, Jac=ds.dJac_dr) / ds.dPhi_n_dr** 2
        # with hardcoded d^2 rho / d Phi_n^2 = -1/(4 rho^3)
        - ds.fluxsurface_integral(1, Jac=ds.Jac) / (4 * ds.rho**3)
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
def minor_major_radius(ds: Evaluations):
    surface_average = ds.volume_integral(1, Jac="Jac_l") / (2 * np.pi)
    ds["minor_radius"] = np.sqrt(surface_average / np.pi)
    ds["major_radius"] = np.sqrt(ds.V / (2 * np.pi * surface_average))


@register(
    requirements=("iota",),
    integration=("rho",),
    attrs=dict(long_name="mean rotational transform", symbol=r"\bar{\iota}"),
)
def iota_mean(ds: Evaluations):
    ds["iota_mean"] = ds.radial_integral("iota", Jac=False)


@register(
    requirements=("mu0", "I_tor", "dPhi_dr", "Jac", "g_tt"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="toroidal current contribution to the rotational transform",
        symbol=r"\iota_{tor}",
    ),
)
def iota_tor(ds: Evaluations):
    Gamma_t = ds.fluxsurface_integral(ds.g_tt / ds.Jac, Jac=False)
    ds["iota_tor"] = 2 * np.pi * ds.mu0 * ds.I_tor / ds.dPhi_dr / Gamma_t


@register(
    requirements=("iota", "diota_dr"),
    attrs=dict(long_name="global magnetic shear", symbol=r"s_g"),
)
def shear(ds: Evaluations):
    ds["shear"] = -ds.rho / ds.iota * ds.diota_dr


@register(
    requirements=("B", "e_theta"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="average poloidal magnetic field", symbol=r"\overline{B_\theta}"
    ),
)
def B_theta_avg(ds: Evaluations):
    ds["B_theta_avg"] = ds.fluxsurface_average(
        xr.dot(ds.B, ds.e_theta, dim="vector"), Jac=False
    )


@register(
    requirements=("B_theta_avg", "mu0"),
    attrs=dict(long_name="toroidal current", symbol=r"I_{tor}"),
)
def I_tor(ds: Evaluations):
    ds["I_tor"] = ds.B_theta_avg * 2 * np.pi / ds.mu0


@register(
    requirements=("B", "e_zeta"),
    integration=("theta", "zeta"),
    attrs=dict(
        long_name="average toroidal magnetic field", symbol=r"\overline{B_\zeta}"
    ),
)
def B_zeta_avg(ds: Evaluations):
    ds["B_zeta_avg"] = ds.fluxsurface_average(
        xr.dot(ds.B, ds.e_zeta, dim="vector"), Jac=False
    )


@register(
    requirements=("B_zeta_avg", "mu0"),
    integration=("theta", "zeta"),
    attrs=dict(long_name="poloidal current", symbol=r"I_{pol}"),
)
def I_pol(ds: Evaluations):
    ds["I_pol"] = ds.B_zeta_avg * 2 * np.pi / ds.mu0
    ds["I_pol"] = ds.I_pol.sel(rho=0, method="nearest") - ds.I_pol
    if not np.isclose(ds.rho.sel(rho=0, method="nearest"), 0):
        logging.warning(
            f"Computation of `I_pol` uses `rho={ds.rho[0].item():e}` instead of the magnetic axis."
        )
