from .post import Evaluations, State

import re

import xarray as xr
import numpy as np

register = Evaluations.register_quantity
rtz_symbols = {"r": r"\rho", "t": r"\theta", "z": r"\zeta"}
rtz_directions = {"r": "radial", "t": "poloidal", "z": "toroidal"}


# === helpers ========================================================================== #


def latex_partial(var, deriv):
    return rf"\frac{{\partial {var}}}{{\partial {rtz_symbols[deriv]}}}"


def latex_partial2(var, deriv1, deriv2):
    if deriv1 == deriv2:
        return rf"\frac{{\partial^2 {var}}}{{\partial {rtz_symbols[deriv1]}^2}}"
    return rf"\frac{{\partial^2 {var}}}{{\partial {rtz_symbols[deriv1]}\partial {rtz_symbols[deriv2]}}}"


# === special ========================================================================== #


@register(coords=())
def mu0(ds: Evaluations):
    ds["mu0"] = 4 * np.pi * 1e-7
    ds.mu0.attrs["long_name"] = "magnetic constant"
    ds.mu0.attrs["symbol"] = r"\mu_0"


# === profiles ========================================================================= #


@register(coords=["rho"])
def iota(ds: Evaluations, state: State):
    ds["iota"] = ("rho", state.evaluate_profile("iota", ds.rho))
    ds.iota.attrs["long_name"] = "rotational transform profile"
    ds.iota.attrs["symbol"] = r"\iota"


@register(coords=["rho"])
def diota_dr(ds: Evaluations, state: State):
    ds["diota_dr"] = ("rho", state.evaluate_profile("iota_prime", ds.rho))
    ds.diota_dr.attrs["long_name"] = "rotational transform gradient profile"
    ds.diota_dr.attrs["symbol"] = r"\frac{d\iota}{d\rho}"


@register(coords=["rho"])
def p(ds: Evaluations, state: State):
    ds["p"] = ("rho", state.evaluate_profile("p", ds.rho))
    ds.p.attrs["long_name"] = "pressure profile"
    ds.p.attrs["symbol"] = r"p"


@register(coords=["rho"])
def dp_dr(ds: Evaluations, state: State):
    ds["dp_dr"] = ("rho", state.evaluate_profile("p_prime", ds.rho))
    ds.dp_dr.attrs["long_name"] = "pressure gradient profile"
    ds.dp_dr.attrs["symbol"] = r"\frac{dp}{d\rho}"


@register(coords=["rho"])
def chi(ds: Evaluations, state: State):
    ds["chi"] = ("rho", state.evaluate_profile("chi", ds.rho))
    ds.chi.attrs["long_name"] = "poloidal magnetic flux profile"
    ds.chi.attrs["symbol"] = r"\chi"


@register(coords=["rho"])
def dchi_dr(ds: Evaluations, state: State):
    ds["dchi_dr"] = ("rho", state.evaluate_profile("chi_prime", ds.rho))
    ds.dchi_dr.attrs["long_name"] = "poloidal magnetic flux gradient profile"
    ds.dchi_dr.attrs["symbol"] = r"\frac{d\chi}{d\rho}"


@register(coords=["rho"])
def Phi(ds: Evaluations, state: State):
    ds["Phi"] = ("rho", state.evaluate_profile("Phi", ds.rho))
    ds.Phi.attrs["long_name"] = "toroidal magnetic flux profile"
    ds.Phi.attrs["symbol"] = r"\Phi"


@register(coords=["rho"])
def dPhi_dr(ds: Evaluations, state: State):
    ds["dPhi_dr"] = ("rho", state.evaluate_profile("Phi_prime", ds.rho))
    ds.dPhi_dr.attrs["long_name"] = "toroidal magnetic flux gradient profile"
    ds.dPhi_dr.attrs["symbol"] = r"\frac{d\Phi}{d\rho}"


@register(coords=["rho"])
def dPhi_drr(ds: Evaluations, state: State):
    ds["dPhi_drr"] = ("rho", state.evaluate_profile("Phi_2prime", ds.rho))
    ds.dPhi_drr.attrs["long_name"] = "toroidal magnetic flux curvature profile"
    ds.dPhi_drr.attrs["symbol"] = r"\frac{d^2\Phi}{d\rho^2}"


@register(coords=["rho"])
def Phi_n(ds: Evaluations, state: State):
    ds["Phi_n"] = ("rho", state.evaluate_profile("PhiNorm", ds.rho))
    ds.Phi_n.attrs["long_name"] = "normalized toroidal magnetic flux profile"
    ds.Phi_n.attrs["symbol"] = r"\Phi_n"


@register(coords=["rho"])
def dPhi_n_dr(ds: Evaluations, state: State):
    ds["dPhi_n_dr"] = ("rho", state.evaluate_profile("PhiNorm_prime", ds.rho))
    ds.dPhi_n_dr.attrs["long_name"] = (
        "normalized toroidal magnetic flux gradient profile"
    )
    ds.dPhi_n_dr.attrs["symbol"] = r"\frac{d\Phi_n}{d\rho}"


# === base ============================================================================= #


@register(
    name=["X1"] + [f"dX1_d{i}" for i in "r t z rr rt rz tt tz zz".split()],
    coords=("rho", "theta", "zeta"),
)
def X1(ds: Evaluations, state: State):
    outputs = state.evaluate_base_tens_all("X1", ds.rho, ds.theta, ds.zeta)
    for key, value in zip(X1.name, outputs):
        ds[key] = (("rho", "theta", "zeta"), value)
    ds.X1.attrs["long_name"] = "first reference coordinate"
    ds.X1.attrs["symbol"] = r"X^1"
    for key in X1.name[1:]:
        deriv = " ".join(rtz_symbols[c] for c in key.split("_")[1][1:])
        ds[key].attrs["long_name"] = "derivative of the first reference coordinate"
        ds[key].attrs["symbol"] = r"\frac{\partial X^1}{\partial" rf"{deriv}}}"


@register(
    name=["X2"] + [f"dX2_d{i}" for i in "r t z rr rt rz tt tz zz".split()],
    coords=("rho", "theta", "zeta"),
)
def X2(ds: Evaluations, state: State):
    outputs = state.evaluate_base_tens_all("X2", ds.rho, ds.theta, ds.zeta)
    for key, value in zip(X2.name, outputs):
        ds[key] = (("rho", "theta", "zeta"), value)
    ds.X1.attrs["long_name"] = "first reference coordinate"
    ds.X1.attrs["symbol"] = r"X^2"
    for key in X1.name[1:]:
        deriv = " ".join(rtz_symbols[c] for c in key.split("_")[1][1:])
        ds[key].attrs["long_name"] = "derivative of the first reference coordinate"
        ds[key].attrs["symbol"] = r"\frac{\partial X^2}{\partial" rf"{deriv}}}"


@register(
    name=["LA"] + [f"dLA_d{i}" for i in "r t z rr rt rz tt tz zz".split()],
    coords=("rho", "theta", "zeta"),
)
def LA(ds: Evaluations, state: State):
    outputs = state.evaluate_base_tens_all("LA", ds.rho, ds.theta, ds.zeta)
    for key, value in zip(LA.name, outputs):
        ds[key] = (("rho", "theta", "zeta"), value)
    ds.LA.attrs["long_name"] = "straight field line potential"
    ds.LA.attrs["symbol"] = r"\lambda"
    for key in LA.name[1:]:
        deriv = " ".join(rtz_symbols[c] for c in key.split("_")[1][1:])
        ds[key].attrs["long_name"] = "derivative of the straight field line potential"
        ds[key].attrs["symbol"] = r"\frac{\partial \lambda}{\partial" rf"{deriv}}}"


# === mapping ========================================================================== #


@register(
    name=("pos", "e_X1", "e_X2", "e_zeta3"),
    coords=("vector", "rho", "theta", "zeta"),
    requirements=("X1", "X2", "zeta"),
)
def hmap(ds: Evaluations, state: State):
    outputs = state.evaluate_hmap_only(
        **{
            var: ds[var].broadcast_like(ds.X1).values.flatten()
            for var in hmap.requirements
        }
    )
    for key, value in zip(hmap.name, outputs):
        ds[key] = (
            hmap.coords,
            value.reshape(3, ds.rho.size, ds.theta.size, ds.zeta.size),
        )

    ds.pos.attrs["long_name"] = "position vector"
    ds.pos.attrs["symbol"] = r"\mathbf{x}"
    ds.e_X1.attrs["long_name"] = "first reference tangent basis vector"
    ds.e_X1.attrs["symbol"] = r"\mathbf{e}_{X^1}"
    ds.e_X2.attrs["long_name"] = "second reference tangent basis vector"
    ds.e_X2.attrs["symbol"] = r"\mathbf{e}_{X^2}"
    ds.e_zeta3.attrs["long_name"] = "toroidal reference tangent basis vector"
    ds.e_zeta3.attrs["symbol"] = r"\mathbf{e}_{\zeta^3}"


# === metric =========================================================================== #


@register(
    name=[
        pattern.format(ij=ij)
        for pattern in ("g_{ij}", "dg_{ij}_dr", "dg_{ij}_dt", "dg_{ij}_dz")
        for ij in ("rr", "rt", "rz", "tt", "tz", "zz")
    ],
    coords=(...,),
    requirements=["X1", "X2", "zeta"]
    + [
        f"d{Xi}_d{j}"
        for j in ("r", "t", "z", "rr", "rt", "rz", "tt", "tz", "zz")
        for Xi in ("X1", "X2")
    ],
)
def metric(ds: Evaluations, state: State):
    idxs = {"r": r"\rho", "t": r"\vartheta", "z": r"\zeta"}
    outputs = state.evaluate_metric(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in metric.requirements]
    )
    for key, value in zip(metric.name, outputs):
        ds[key] = (
            ds.X1.coords,
            value.reshape(ds.X1.shape),
        )
        if m := re.match(r"dg_(\w\w)_d(\w)", key):
            ij, k = m.groups()
            ij = idxs[ij[0]] + idxs[ij[1]]
            k = idxs[k]
            ds[key].attrs["long_name"] = "derivative of a metric coefficient"
            ds[key].attrs["symbol"] = rf"\frac{{\partial g_{{{ij}}}}}{{\partial {k}}}"
        else:
            ij = key.split("_")[1]
            ij = idxs[ij[0]] + idxs[ij[1]]
            ds[key].attrs["long_name"] = "metric coefficient"
            ds[key].attrs["symbol"] = rf"g_{{{ij}}}"


# === jacobian determinant ============================================================= #


@register(
    name=("Jac_h", *(f"dJac_h_d{i}" for i in "r t z".split())),
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
)
def Jac_h(ds: Evaluations, state: State):
    outputs = state.evaluate_jacobian(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in Jac_h.requirements]
    )
    for key, value in zip(Jac_h.name, outputs):
        ds[key] = (
            ds.X1.coords,
            value.reshape(ds.X1.shape),
        )
        if key == "Jac_h":
            ds[key].attrs["long_name"] = "reference Jacobian determinant"
            ds[key].attrs["symbol"] = r"\mathcal{J}_h"
        elif m := re.match(r"dJac_h_d(\w)", key):
            i = m.group(1)
            ds[key].attrs[
                "long_name"
            ] = f"{rtz_directions[i]} derivative of the reference Jacobian determinant"
            ds[key].attrs["symbol"] = latex_partial(r"\mathcal{J}_h", i)
        else:
            raise RuntimeError(f"unexpected key {key}")


@register(
    name=("Jac", "Jac_l", *(f"dJac{suf}_d{i}" for suf in ["", "_l"] for i in "rtz")),
    requirements=(
        "Jac_h",
        *(f"dJac_h_d{i}" for i in "rtz"),
        *(
            f"d{Q}_d{i}"
            for Q in "X1 X2".split()
            for i in "r t z rr rt rz tt tz zz".split()
        ),
    ),
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
    ds.Jac_l.attrs["long_name"] = "logical Jacobian determinant"
    ds.Jac_l.attrs["symbol"] = r"\mathcal{J}_l"
    ds.Jac.attrs["long_name"] = "Jacobian determinant"
    ds.Jac.attrs["symbol"] = r"\mathcal{J}"
    for i in "rtz":
        da = ds[f"dJac_d{i}"]
        da.attrs["long_name"] = (
            f"{rtz_directions[i]} derivative of the Jacobian determinant"
        )
        da.attrs["symbol"] = latex_partial(r"\mathcal{J}", i)
        da = ds[f"dJac_l_d{i}"]
        da.attrs["long_name"] = (
            f"{rtz_directions[i]} derivative of the logical Jacobian determinant"
        )
        da.attrs["symbol"] = latex_partial(r"\mathcal{J}_l", i)


# === derived ========================================================================== #


@register(
    coords=("vector", ...),
    requirements=("e_X1", "e_X2", "dX1_dr", "dX2_dr"),
)
def e_rho(ds: Evaluations):
    da = ds.e_X1 * ds.dX1_dr + ds.e_X2 * ds.dX2_dr
    da.attrs["long_name"] = "radial tangent basis vector"
    da.attrs["symbol"] = r"\mathbf{e}_\rho"
    ds["e_rho"] = da


@register(
    coords=("vector", ...),
    requirements=["e_X1", "e_X2", "dX1_dt", "dX2_dt"],
)
def e_theta(ds: Evaluations):
    da = ds.e_X1 * ds.dX1_dt + ds.e_X2 * ds.dX2_dt
    da.attrs["long_name"] = "poloidal tangent basis vector"
    da.attrs["symbol"] = r"\mathbf{e}_\theta"
    ds["e_theta"] = da


@register(
    coords=("vector", ...),
    requirements=("e_X1", "e_X2", "e_zeta3", "dX1_dz", "dX2_dz"),
)
def e_zeta(ds: Evaluations):
    da = ds.e_X1 * ds.dX1_dz + ds.e_X2 * ds.dX2_dz + ds.e_zeta3
    da.attrs["long_name"] = "toroidal tangent basis vector"
    da.attrs["symbol"] = r"\mathbf{e}_\zeta"
    ds["e_zeta"] = da


@register(
    coords=("vector", ...),
    requirements=("Jac", "e_theta", "e_zeta"),
)
def grad_rho(ds: Evaluations):
    da = xr.cross(ds.e_theta, ds.e_zeta, dim="vector") / ds.Jac
    da.attrs["long_name"] = "radial reciprocal basis vector"
    da.attrs["symbol"] = r"\nabla\rho"
    ds["grad_rho"] = da


@register(
    coords=("vector", ...),
    requirements=("Jac", "e_rho", "e_zeta"),
)
def grad_theta(ds: Evaluations):
    da = xr.cross(ds.e_zeta, ds.e_rho, dim="vector") / ds.Jac
    da.attrs["long_name"] = "poloidal reciprocal basis vector"
    da.attrs["symbol"] = r"\nabla\theta"
    ds["grad_theta"] = da


@register(
    coords=("vector", ...),
    requirements=("Jac", "e_rho", "e_theta"),
)
def grad_zeta(ds: Evaluations):
    da = xr.cross(ds.e_rho, ds.e_theta, dim="vector") / ds.Jac
    da.attrs["long_name"] = "toroidal reciprocal basis vector"
    da.attrs["symbol"] = r"\nabla\zeta"
    ds["grad_zeta"] = da


@register(
    name=("B", "B_contra_t", "B_contra_z"),
    coords=("vector", ...),
    requirements=(
        "iota",
        "dLA_dt",
        "dLA_dz",
        "dPhi_dr",
        "Jac",
        "e_theta",
        "e_zeta",
    ),
)
def B(ds: Evaluations):
    ds["B_contra_t"] = (ds.iota - ds.dLA_dz) * ds.dPhi_dr / ds.Jac
    ds["B_contra_z"] = (1 + ds.dLA_dt) * ds.dPhi_dr / ds.Jac
    ds["B"] = ds.B_contra_t * ds.e_theta + ds.B_contra_z * ds.e_zeta
    ds.B.attrs["long_name"] = "magnetic field"
    ds.B.attrs["symbol"] = r"\mathbf{B}"
    ds.B_contra_t.attrs["long_name"] = "poloidal magnetic field"
    ds.B_contra_t.attrs["symbol"] = r"B^\theta"
    ds.B_contra_z.attrs["long_name"] = "toroidal magnetic field"
    ds.B_contra_z.attrs["symbol"] = r"B^\zeta"


@register(
    name=[f"dB_contra_{i}_d{j}" for i in "tz" for j in "rtz"],
    coords=(...,),
    requirements=[
        "Jac",
        "dPhi_dr",
        "dPhi_drr",
        "iota",
        "diota_dr",
    ]
    + [f"dJac_d{i}" for i in "r t z".split()]
    + [f"dLA_d{i}" for i in "t z rt rz tt tz zz".split()],
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
    for i in "tz":
        for j in "rtz":
            da = ds[f"dB_contra_{i}_d{j}"]
            da.attrs["long_name"] = "derivative of magnetic field"
            da.attrs["symbol"] = latex_partial(f"B^{rtz_symbols[i]}", j)


@register(
    name=["J", "J_contra_r", "J_contra_t", "J_contra_z"],
    coords=(...,),
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

    for i in "rtz":
        da = ds[f"J_contra_{i}"]
        da.attrs["long_name"] = f"contravariant {rtz_directions[i]} current density"
        da.attrs["symbol"] = rf"J^{{{rtz_symbols[i]}}}"
    ds.J.attrs["long_name"] = "current density"
    ds.J.attrs["symbol"] = r"\mathbf{J}"


def _modulus(v):
    @register(
        name=f"mod_{v}",
        coords=(...,),
        requirements=(v,),
    )
    def mod_v(ds: Evaluations):
        da = np.sqrt(xr.dot(ds[v], ds[v], dim="vector"))
        da.attrs["long_name"] = f"modulus of the {ds[v].attrs['long_name']}"
        da.attrs["symbol"] = rf"\left|{ds[v].attrs['symbol']}\right|"
        ds[f"mod_{v}"] = da

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
