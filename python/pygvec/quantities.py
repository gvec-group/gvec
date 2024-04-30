from .post import Evaluations, State
from xarray import Dataset
import xarray as xr

register = Evaluations.register_quantity


# === profiles ================================================================================================================= #


@register(coords=["rho"])
def iota(ds: Dataset, state: State):
    ds["iota"] = ("rho", state.evaluate_profile("iota", ds.rho))
    ds.iota.attrs["long_name"] = "rotational transform profile"
    ds.iota.attrs["symbol"] = r"\iota"


@register(coords=["rho"])
def diota_drho(ds: Dataset, state: State):
    ds["diota_drho"] = ("rho", state.evaluate_profile("iota_prime", ds.rho))
    ds.diota_drho.attrs["long_name"] = "rotational transform gradient profile"
    ds.diota_drho.attrs["symbol"] = r"\frac{d\iota}{d\rho}"


@register(coords=["rho"])
def p(ds: Dataset, state: State):
    ds["p"] = ("rho", state.evaluate_profile("p", ds.rho))
    ds.p.attrs["long_name"] = "pressure profile"
    ds.p.attrs["symbol"] = r"p"


@register(coords=["rho"])
def dp_drho(ds: Dataset, state: State):
    ds["dp_drho"] = ("rho", state.evaluate_profile("p_prime", ds.rho))
    ds.dp_drho.attrs["long_name"] = "pressure gradient profile"
    ds.dp_drho.attrs["symbol"] = r"\frac{dp}{d\rho}"


@register(coords=["rho"])
def chi(ds: Dataset, state: State):
    ds["chi"] = ("rho", state.evaluate_profile("chi", ds.rho))
    ds.chi.attrs["long_name"] = "poloidal magnetic flux profile"
    ds.chi.attrs["symbol"] = r"\chi"


@register(coords=["rho"])
def dchi_drho(ds: Dataset, state: State):
    ds["dchi_drho"] = ("rho", state.evaluate_profile("chi_prime", ds.rho))
    ds.dchi_drho.attrs["long_name"] = "poloidal magnetic flux gradient profile"
    ds.dchi_drho.attrs["symbol"] = r"\frac{d\chi}{d\rho}"


@register(coords=["rho"])
def Phi(ds: Dataset, state: State):
    ds["Phi"] = ("rho", state.evaluate_profile("Phi", ds.rho))
    ds.Phi.attrs["long_name"] = "toroidal magnetic flux profile"
    ds.Phi.attrs["symbol"] = r"\Phi"


@register(coords=["rho"])
def dPhi_drho(ds: Dataset, state: State):
    ds["dPhi_drho"] = ("rho", state.evaluate_profile("Phi_prime", ds.rho))
    ds.dPhi_drho.attrs["long_name"] = "toroidal magnetic flux gradient profile"
    ds.dPhi_drho.attrs["symbol"] = r"\frac{d\Phi}{d\rho}"


@register(coords=["rho"])
def d2Phi_drho2(ds: Dataset, state: State):
    ds["d2Phi_drho2"] = ("rho", state.evaluate_profile("Phi_2prime", ds.rho))
    ds.d2Phi_drho2.attrs["long_name"] = "toroidal magnetic flux curvature profile"
    ds.d2Phi_drho2.attrs["symbol"] = r"\frac{d^2\Phi}{d\rho^2}"


@register(coords=["rho"])
def Phi_n(ds: Dataset, state: State):
    ds["Phi_n"] = ("rho", state.evaluate_profile("PhiNorm", ds.rho))
    ds.Phi_n.attrs["long_name"] = "normalized toroidal magnetic flux profile"
    ds.Phi_n.attrs["symbol"] = r"\Phi_n"


@register(coords=["rho"])
def dPhi_n_drho(ds: Dataset, state: State):
    ds["dPhi_n_drho"] = ("rho", state.evaluate_profile("PhiNorm_prime", ds.rho))
    ds.dPhi_n_drho.attrs["long_name"] = (
        "normalized toroidal magnetic flux gradient profile"
    )
    ds.dPhi_n_drho.attrs["symbol"] = r"\frac{d\Phi_n}{d\rho}"


# === base ===================================================================================================================== #


_X1_out = ["X1", "dX1_drho", "dX1_dtheta", "dX1_dzeta"]


@register(name=_X1_out, coords=["rho", "theta", "zeta"])
def X1(ds: Dataset, state: State):
    outputs = state.evaluate_base_tens_all("X1", ds.rho, ds.theta, ds.zeta)
    for key, value in zip(_X1_out, outputs):
        ds[key] = (("rho", "theta", "zeta"), value)
    ds.X1.attrs["long_name"] = "first reference coordinate"
    ds.X1.attrs["symbol"] = "X^1"
    ds.dX1_drho.attrs["long_name"] = (
        "radial derivative of the first reference coordinate in logical coordinates"
    )
    ds.dX1_drho.attrs["symbol"] = r"\frac{\partial X^1}{\partial \rho}"
    ds.dX1_dtheta.attrs["long_name"] = (
        "poloidal derivative of the first reference coordinate in logical coordinates"
    )
    ds.dX1_dtheta.attrs["symbol"] = r"\frac{\partial X^1}{\partial \theta}"
    ds.dX1_dzeta.attrs["long_name"] = (
        "toroidal derivative of the first reference coordinate in logical coordinates"
    )
    ds.dX1_dzeta.attrs["symbol"] = r"\frac{\partial X^1}{\partial \zeta}"


_X2_out = ["X2", "dX2_drho", "dX2_dtheta", "dX2_dzeta"]


@register(name=_X2_out, coords=["rho", "theta", "zeta"])
def X2(ds: Dataset, state: State):
    outputs = state.evaluate_base_tens_all("X2", ds.rho, ds.theta, ds.zeta)
    for key, value in zip(_X2_out, outputs):
        ds[key] = (("rho", "theta", "zeta"), value)
    ds.X2.attrs["long_name"] = "second reference coordinate"
    ds.X2.attrs["symbol"] = "X^2"
    ds.dX2_drho.attrs["long_name"] = (
        "radial derivative of the second reference coordinate in logical coordinates"
    )
    ds.dX2_drho.attrs["symbol"] = r"\frac{\partial X^2}{\partial \rho}"
    ds.dX2_dtheta.attrs["long_name"] = (
        "poloidal derivative of the second reference coordinate in logical coordinates"
    )
    ds.dX2_dtheta.attrs["symbol"] = r"\frac{\partial X^2}{\partial \theta}"
    ds.dX2_dzeta.attrs["long_name"] = (
        "toroidal derivative of the second reference coordinate in logical coordinates"
    )
    ds.dX2_dzeta.attrs["symbol"] = r"\frac{\partial X^2}{\partial \zeta}"


_LA_out = ["LA", "dLA_drho", "dLA_dtheta", "dLA_dzeta"]


@register(name=_LA_out, coords=["rho", "theta", "zeta"])
def LA(ds: Dataset, state: State):
    outputs = state.evaluate_base_tens_all("LA", ds.rho, ds.theta, ds.zeta)
    for key, value in zip(_LA_out, outputs):
        ds[key] = (("rho", "theta", "zeta"), value)
    ds.LA.attrs["long_name"] = "straight field line potential"
    ds.LA.attrs["symbol"] = r"\lambda"
    ds.dLA_drho.attrs["long_name"] = (
        "radial derivative of the straight field line potential"
    )
    ds.dLA_drho.attrs["symbol"] = r"\frac{\partial \lambda}{\partial \rho}"
    ds.dLA_dtheta.attrs["long_name"] = (
        "poloidal derivative of the straight field line potential"
    )
    ds.dLA_dtheta.attrs["symbol"] = r"\frac{\partial \lambda}{\partial \theta}"
    ds.dLA_dzeta.attrs["long_name"] = (
        "toroidal derivative of the straight field line potential"
    )
    ds.dLA_dzeta.attrs["symbol"] = r"\frac{\partial \lambda}{\partial \zeta}"


# === mapping ================================================================================================================== #


_hmap_out = ["pos", "e_rho", "e_theta", "e_zeta"]
_hmap_reqs = [
    "X1",
    "X2",
    "theta",
    "dX1_drho",
    "dX2_drho",
    "dX1_dtheta",
    "dX2_dtheta",
    "dX1_dzeta",
    "dX2_dzeta",
]


@register(
    name=_hmap_out, coords=["vector", "rho", "theta", "zeta"], requirements=_hmap_reqs
)
def hmap(ds: Dataset, state: State):
    outputs = state.evaluate_hmap(
        *[ds[var].broadcast_like(ds.X1).values.flatten() for var in _hmap_reqs]
    )
    for key, value in zip(_hmap_out, outputs):
        ds[key] = (
            ("vector", "rho", "theta", "zeta"),
            value.reshape(3, ds.rho.size, ds.theta.size, ds.zeta.size),
        )

    ds.pos.attrs["long_name"] = "position vector"
    ds.pos.attrs["symbol"] = r"\mathbf{x}"
    ds.e_rho.attrs["long_name"] = "radial reciprocal basis vector"
    ds.e_rho.attrs["symbol"] = r"\mathbf{x}_\rho"
    ds.e_theta.attrs["long_name"] = "poloidal reciprocal basis vector"
    ds.e_theta.attrs["symbol"] = r"\mathbf{x}_\vartheta"
    ds.e_zeta.attrs["long_name"] = "toroidal reciprocal basis vector"
    ds.e_zeta.attrs["symbol"] = r"\mathbf{x}_\zeta"


# === derived ================================================================================================================== #


@register(
    coords=["rho", "theta", "zeta"],
    requirements=["e_rho", "e_theta", "e_zeta"],
)
def Jac(ds: Dataset):
    ds["Jac"] = xr.dot(
        ds.e_rho, xr.cross(ds.e_theta, ds.e_zeta, dim="vector"), dim="vector"
    )
    ds.Jac.attrs["long_name"] = "Jacobian determinant"
    ds.Jac.attrs["symbol"] = r"\mathcal{J}"


@register(
    coords=["vector", "rho", "theta", "zeta"],
    requirements=[
        "iota",
        "dLA_dtheta",
        "dLA_dzeta",
        "dPhi_drho",
        "Jac",
        "e_theta",
        "e_zeta",
    ],
)
def B(ds: Dataset):
    B_contra_theta = (ds.iota - ds.dLA_dzeta) * ds.dPhi_drho / ds.Jac
    B_contra_zeta = (1 + ds.dLA_dtheta) * ds.dPhi_drho / ds.Jac
    ds["B"] = B_contra_theta * ds.e_theta + B_contra_zeta * ds.e_zeta
    ds.B.attrs["long_name"] = "magnetic field"
    ds.B.attrs["symbol"] = r"\mathbf{B}"
