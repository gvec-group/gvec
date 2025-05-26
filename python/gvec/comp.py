# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""GVEC Postprocessing - Compute Functions"""

from typing import Literal, Iterable, Mapping, MutableMapping, Callable
import logging
import inspect
import re

import numpy as np
import xarray as xr

from .state import State
from . import fourier

# === Globals === #

__all__ = [
    "QUANTITIES",
    "register",
    "compute",
    "table_of_quantities",
    "Evaluations",
    "EvaluationsBoozer",
    "EvaluationsBoozerCustom",
    "radial_integral",
    "fluxsurface_integral",
    "volume_integral",
    "ev2ft",
    "ft_autoremove",
]
QUANTITIES = {}  # dictionary to store the registered quantities (compute functions)


# === helpers ========================================================================== #


rtz_symbols = {"r": r"\rho", "t": r"\theta", "z": r"\zeta"}
rtz_directions = {"r": "radial", "t": "poloidal", "z": "toroidal"}


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


# === Register Compute Functions === #


def register(
    quantities: None | str | Iterable[str] = None,
    requirements: Iterable[str] = (),
    integration: Iterable[str] = (),
    attrs: Mapping = {},
    registry: MutableMapping = QUANTITIES,
):
    """Function decorator to register equilibrium quantities.

    The quantity (compute function) is registered in the QUANTITIES dictionary.
    It contains:
        * a function pointer
        * the name of the computed quantities (used as key in QUANTITIES)
        * the names of required quantities (that should be computed before)
        * the names of the integration axes required for the computation
        * the attributes of the computed quantity (long_name, symbol, etc.)
    """

    def _register(
        func: (
            Callable[[xr.Dataset], xr.Dataset]
            | Callable[[xr.Dataset, State], xr.Dataset]
        ),
    ):
        nonlocal quantities, requirements, integration, attrs
        if quantities is None:
            quantities = [func.__name__]
        if isinstance(quantities, str):
            quantities = [quantities]
        func.quantities = quantities
        func.requirements = requirements
        func.integration = integration
        if len(quantities) == 1 and quantities[0] not in attrs:
            attrs = {quantities[0]: attrs}
        func.attrs = attrs

        for q in quantities:
            if q in registry:
                logging.warning(f"A quantity `{q}` is already registered.")
            registry[q] = func
        return func

    return _register


def table_of_quantities(markdown: bool = False, registry: Mapping = QUANTITIES):
    """
    Generate a table of computable quantities.

    Parameters
    ----------
    markdown : optional
        If True, return the table as a Ipython.Markdown object. Otherwise, return the table as a string.

    Returns
    -------
    str or IPython.display.Markdown
        The table of quantities. If `markdown` is True, the table is returned as an instance of
        IPython.display.Markdown. Otherwise, the table is returned as a string.

    Notes
    -----
    This method generates a table of quantities based on the attributes of the registered quantities.
    The table includes the label, long name, and symbol of each quantity.
    """
    lines = []
    for key, func in sorted(list(registry.items())):
        long_name = func.attrs[key].get("long_name", "")
        symbol = func.attrs[key].get("symbol", "")
        symbol = "$" + symbol.replace("|", r"\|") + "$"
        lines.append((f"`{key}`", long_name, symbol))
    sizes = [max(len(s) for s in col) for col in zip(*lines)]
    txt = f"| {'label':^{sizes[0]}s} | {'long name':^{sizes[1]}s} | {'symbol':^{sizes[2]}s} |\n"
    txt += f"| {'-' * sizes[0]} | {'-' * sizes[1]} | {'-' * sizes[2]} |\n"
    for line in lines:
        txt += f"| {line[0]:^{sizes[0]}s} | {line[1]:^{sizes[1]}s} | {line[2]:^{sizes[2]}s} |\n"
    if markdown:
        from IPython.display import Markdown

        return Markdown(txt)
    else:
        return txt


def compute(
    ev: xr.Dataset,
    *quantities: Iterable[str],
    state: State = None,
    registry: Mapping = QUANTITIES,
) -> xr.Dataset | xr.DataArray:
    """Compute the target equilibrium quantity.

    This method will compute required parameters recursively.
    """
    for quantity in quantities:
        # --- get the compute function --- #
        if quantity in ev:
            continue  # already computed
        if quantity not in registry:
            raise KeyError(f"The quantity `{quantity}` is not registered.")
        func = registry[quantity]
        # --- handle special cases --- #
        # if dLA_dr is requested (and not already present), but LA is already present - this is currently only possible with the boozer transform
        # ToDo: cleanup
        if m := re.match(quantity, r"dLA_d([rtz]{1,2})"):
            if "r" in m.group(1) and "LA" in ev:
                raise ValueError(
                    f"Cannot compute `{quantity}` as it is a radial derivative of lambda, as lambda was recomputed on a surface with a boozer transform."
                )

        # --- handle integration --- #
        # we assume the dimensions are {rad, pol, tor} or {pol, tor}
        # we don't assume which coordinates are associated with which dimensions
        # in particular: (rho, theta, zeta), (rho, theta_B, zeta_B), (rho, alpha, phi_alpha) are all expected
        # some quantities may require integration points in any of {rho, theta, zeta}
        # if the integration points are not present we will create an auxiliary dataset with integration points
        auxcoords = {
            i
            for i in func.integration
            if i not in ev
            or "integration_points" not in ev[i].attrs
            or ev[i].attrs["integration_points"] == "False"
        }
        if auxcoords:
            # --- auxiliary dataset for integration --- #
            logging.info(
                f"Using auxiliary dataset with integration points in {auxcoords} to compute {quantity}."
            )
            if auxcoords > {"rho", "theta", "zeta"}:
                raise ValueError(
                    f"Unsupported integration coordinates for auxiliary dataset: {auxcoords}"
                )
            rho = "int" if "rho" in auxcoords else ev.rho if "rho" in ev else None
            theta = (
                "int" if "theta" in auxcoords else ev.theta if "theta" in ev else None
            )
            zeta = "int" if "zeta" in auxcoords else ev.zeta if "zeta" in ev else None
            obj = Evaluations(rho=rho, theta=theta, zeta=zeta, state=state)
        else:
            obj = ev
        # --- handle requirements --- #
        compute(obj, *func.requirements, state=state, registry=registry)
        # --- compute the quantity --- #
        with xr.set_options(keep_attrs=True):
            if "state" in inspect.signature(func).parameters:
                if state is None:
                    raise ValueError(
                        f"Computation of the quantity `{func.__name__}` requires a state object."
                    )
                func(obj, state)
            else:
                func(obj)
        # --- set attributes --- #
        for q in func.quantities:
            if q in func.attrs:
                obj[q].attrs.update(func.attrs[q])
        # --- handle auxiliary integration dataset --- #
        if auxcoords:
            for q in obj:
                if "weight" in q:
                    continue
                if any([c in auxcoords for c in obj[q].coords]):
                    continue
                ev[q] = (obj[q].dims, obj[q].data, obj[q].attrs)
    if len(quantities) == 1:
        return ev[quantities[0]]
    return ev


# === Create Evaluations Dataset === #


def Evaluations(
    rho: int | Literal["int"] | tuple[float, float, int] | np.ndarray | None = "int",
    theta: (
        int | Literal["int"] | tuple[float, float, int] | np.ndarray | None
    ) = "int",
    zeta: (int | Literal["int"] | tuple[float, float, int] | np.ndarray | None) = "int",
    state: State | None = None,
    nfp: int | None = None,
):
    coords = {}
    # --- get integration points --- #
    if state is not None:
        intp = [state.get_integration_points(q) for q in ["X1", "X2", "LA"]]
        if nfp is not None:
            logging.warning("Both `state` and `nfp` are provided. Disregarding `nfp`.")
        nfp = state.nfp
    # --- parse coordinates --- #
    match rho:
        case xr.DataArray():
            coords["rho"] = rho
        case np.ndarray() | list():
            coords["rho"] = ("rad", rho)
        case "int":
            if state is None:
                raise ValueError("Integration points require a state object.")
            if any(
                [
                    not np.allclose(intp[0][j], intp[i][j])
                    for i in (1, 2)
                    for j in (0, 1)
                ]
            ):
                raise ValueError(
                    "Integration points for rho do not align for X1, X2 and LA."
                )
            coords["rho"] = ("rad", intp[0][0])
            coords["rad_weight"] = ("rad", intp[0][1])
        case int() as num:
            coords["rho"] = ("rad", np.linspace(0, 1, num))
            coords["rho"][1][0] = (
                0.1 * coords["rho"][1][1]
            )  # avoid numerical issues at the magnetic axis
        case (start, stop):
            coords["rho"] = ("rad", np.linspace(start, stop))
        case (start, stop, num):
            coords["rho"] = ("rad", np.linspace(start, stop, num))
        case None:
            pass
        case _:
            raise ValueError(f"Could not parse rho, got {rho}.")
    match theta:
        case xr.DataArray():
            coords["theta"] = theta
        case np.ndarray() | list():
            coords["theta"] = ("pol", theta)
        case "int":
            if state is None:
                raise ValueError("Integration points require a state object.")
            if any(
                [
                    not np.allclose(intp[0][j], intp[i][j])
                    for i in (1, 2)
                    for j in (2, 3)
                ]
            ):
                raise ValueError(
                    "Integration points for theta do not align for X1, X2 and LA."
                )
            coords["theta"] = (
                "pol",
                np.linspace(0, 2 * np.pi, intp[0][2], endpoint=False),
            )
            coords["pol_weight"] = intp[0][3]
        case int() as num:
            coords["theta"] = ("pol", np.linspace(0, 2 * np.pi, num, endpoint=False))
        case (start, stop):
            coords["theta"] = ("pol", np.linspace(start, stop))
        case (start, stop, num):
            coords["theta"] = ("pol", np.linspace(start, stop, num))
        case None:
            pass
        case _:
            raise ValueError(f"Could not parse theta, got {theta}.")
    match zeta:
        case xr.DataArray():
            coords["zeta"] = zeta
        case np.ndarray() | list():
            coords["zeta"] = ("tor", zeta)
        case "int":
            if state is None:
                raise ValueError("Integration points require a state object.")
            if any(
                [
                    not np.allclose(intp[0][j], intp[i][j])
                    for i in (1, 2)
                    for j in (4, 5)
                ]
            ):
                raise ValueError(
                    "Integration points for zeta do not align for X1, X2 and LA."
                )
            coords["zeta"] = (
                "tor",
                np.linspace(0, 2 * np.pi / nfp, intp[0][4], endpoint=False),
            )
            coords["tor_weight"] = intp[0][5]
        case int() as num:
            if nfp is None:
                raise ValueError("Automatic bounds for zeta require `nfp`.")
            coords["zeta"] = (
                "tor",
                np.linspace(0, 2 * np.pi / nfp, num, endpoint=False),
            )
        case (start, stop):
            coords["zeta"] = ("tor", np.linspace(start, stop))
        case (start, stop, num):
            coords["zeta"] = ("tor", np.linspace(start, stop, num))
        case None:
            pass
        case _:
            raise ValueError(f"Could not parse zeta, got {zeta}.")

    # --- init Dataset --- #
    ds = xr.Dataset(coords=coords)

    # --- set attributes & indices --- #
    if "rho" in ds:
        ds.rho.attrs["long_name"] = "Logical radial coordinate"
        ds.rho.attrs["symbol"] = r"\rho"
        ds.rho.attrs["integration_points"] = str(isinstance(rho, str) and rho == "int")
        if ds.rho.dims == ("rad",):
            ds = ds.set_xindex("rho")
    if "theta" in ds:
        ds.theta.attrs["long_name"] = "Logical poloidal angle"
        ds.theta.attrs["symbol"] = r"\theta"
        ds.theta.attrs["integration_points"] = str(
            isinstance(theta, str) and theta == "int"
        )
        if ds.theta.dims == ("pol",):
            ds = ds.set_xindex("theta")
    if "zeta" in ds:
        ds.zeta.attrs["long_name"] = "Logical toroidal angle"
        ds.zeta.attrs["symbol"] = r"\zeta"
        ds.zeta.attrs["integration_points"] = str(
            isinstance(zeta, str) and zeta == "int"
        )
        if ds.zeta.dims == ("tor",):
            ds = ds.set_xindex("zeta")

    if (
        "theta" in ds
        and "zeta" in ds
        and set(ds.theta.dims) >= {"pol", "tor"}
        and set(ds.zeta.dims) >= {"pol", "tor"}
    ):
        ds = ds.set_xindex("theta", "zeta")
    return ds


def radial_integral(quantity: xr.DataArray):
    """Compute the radial integral/average of the given quantity."""
    # --- check for integration points --- #
    if "rad_weight" not in quantity.coords:
        raise ValueError("Radial integral requires integration weights for `rad`.")
    # --- integrate --- #
    return (quantity * quantity.rad_weight).sum("rad")


def fluxsurface_integral(quantity: xr.DataArray):
    """Compute the flux surface integral of the given quantity."""
    # --- check for integration points --- #
    if "pol_weight" not in quantity.coords or "tor_weight" not in quantity.coords:
        raise ValueError(
            "Flux surface average requires integration weights for theta and zeta."
        )
    # --- integrate --- #
    return (quantity * quantity.pol_weight * quantity.tor_weight).sum(("pol", "tor"))


def volume_integral(
    quantity: xr.DataArray,
):
    """Compute the volume integral of the given quantity."""
    # --- check for integration points --- #
    if (
        "rad_weight" not in quantity.coords
        or "pol_weight" not in quantity.coords
        or "tor_weight" not in quantity.coords
    ):
        raise ValueError(
            "Volume integral requires integration weights for rho, theta and zeta."
        )
    # --- integrate --- #
    return (
        quantity * quantity.rad_weight * quantity.pol_weight * quantity.tor_weight
    ).sum(("rad", "pol", "tor"))


def EvaluationsBoozer(
    rho: float | np.ndarray,
    n_theta: int,
    n_zeta: int,
    state: State,
    M: int | None = None,
    N: int | None = None,
    sincos: Literal["sin", "cos", "sincos"] = "sin",
):
    rho = np.asarray(rho)
    theta_B = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    zeta_B = np.linspace(0, 2 * np.pi / state.nfp, n_zeta, endpoint=False)

    squeeze = False
    if rho.ndim == 0:
        rho = np.array([rho])
        squeeze = True
    elif rho.ndim != 1:
        raise ValueError("rho must be 1D")

    ds = xr.Dataset(
        coords=dict(
            rho=("rad", rho),
            theta_B=("pol", theta_B),
            zeta_B=("tor", zeta_B),
        )
    )

    # === Find the logical coordinates of the Boozer grid === #
    stacked = ds[["theta_B", "zeta_B"]].stack(tz=("pol", "tor"))
    tz_B = np.stack([stacked.theta_B, stacked.zeta_B], axis=0)
    sfl_boozer = state.get_boozer(rho, M, N, sincos=sincos)
    tz = state.get_boozer_angles(sfl_boozer, tz_B)
    stacked["theta"] = (("tz", "rad"), tz[0, :, :])
    stacked["zeta"] = (("tz", "rad"), tz[1, :, :])
    ds["theta"] = stacked["theta"].unstack("tz")
    ds["zeta"] = stacked["zeta"].unstack("tz")

    # === Metadata === #
    ds.rho.attrs["long_name"] = "Logical radial coordinate"
    ds.rho.attrs["symbol"] = r"\rho"
    ds.theta_B.attrs["long_name"] = "Boozer straight-fieldline poloidal angle"
    ds.theta_B.attrs["symbol"] = r"\theta_B"
    ds.zeta_B.attrs["long_name"] = "Boozer toroidal angle"
    ds.zeta_B.attrs["symbol"] = r"\zeta_B"
    ds.theta.attrs["long_name"] = "Logical poloidal angle"
    ds.theta.attrs["symbol"] = r"\theta"
    ds.zeta.attrs["long_name"] = "Logical toroidal angle"
    ds.zeta.attrs["symbol"] = r"\zeta"

    # === Indices === #
    # setting them earlier causes issues with the stacking / unstacking
    ds = ds.set_xindex("rho")
    ds = ds.set_xindex("theta_B")
    ds = ds.set_xindex("zeta_B")
    ds = ds.drop_vars("pol")
    ds = ds.drop_vars("tor")

    # === Evaluate LA & NU === #
    # Flatten theta, zeta
    theta = ds.theta.transpose("rad", "pol", "tor").values.reshape(ds.rad.size, -1)
    zeta = ds.zeta.transpose("rad", "pol", "tor").values.reshape(ds.rad.size, -1)

    # Compute base on each radial position
    outputs_la = []
    outputs_nu = []
    for r, rho in enumerate(ds.rho.data):
        thetazeta = np.stack([theta[r, :], zeta[r, :]], axis=0)
        outputs_la.append(
            state.evaluate_boozer_list_tz_all(sfl_boozer, "LA", [r], thetazeta)
        )
        outputs_nu.append(
            state.evaluate_boozer_list_tz_all(sfl_boozer, "NU", [r], thetazeta)
        )

    # Write LA/NU to dataset
    for deriv, value in zip(["", "t", "z", "tt", "tz", "zz"], zip(*outputs_la)):
        if deriv == "":
            var = "LA"
            long_name = "Straight field line potential"
            symbol = r"\lambda"
        else:
            var = f"dLA_d{deriv}"
            long_name = derivative_name_smart("Straight field line potential", deriv)
            symbol = latex_partial_smart(r"\lambda", deriv)
        value = np.stack(value).reshape(ds.rad.size, ds.pol.size, ds.tor.size)
        ds[var] = (
            ("rad", "pol", "tor"),
            value,
            dict(long_name=long_name, symbol=symbol),
        )
    for deriv, value in zip(["", "t", "z", "tt", "tz", "zz"], zip(*outputs_nu)):
        if deriv == "":
            var = "NU_B"
            long_name = "Boozer angular potential"
            symbol = r"\nu_B"
        else:
            var = f"dNU_B_d{deriv}"
            long_name = derivative_name_smart("Boozer angular potential", deriv)
            symbol = latex_partial_smart(r"\nu_B", deriv)
        value = np.stack(value).reshape(ds.rad.size, ds.pol.size, ds.tor.size)
        ds[var] = (
            ("rad", "pol", "tor"),
            value,
            dict(long_name=long_name, symbol=symbol),
        )

    if squeeze:
        ds = ds.squeeze("rad")

    return ds


def EvaluationsBoozerCustom(
    rho: float | np.ndarray,
    theta_B: np.ndarray,
    zeta_B: np.ndarray,
    poloidal: tuple[str, np.ndarray],
    toroidal: tuple[str, np.ndarray],
    state: State,
    M: int | None = None,
    N: int | None = None,
    sincos: Literal["sin", "cos", "sincos"] = "sin",
):
    """Create an Evaluations dataset with a custom Boozer grid.

    This factory function assumes that the target grid still has a poloidal-like and toroidal-like direction
    and that the Boozer coordinates of each grid point are known. They do not have to lie within a single field
    period, nor do they have to be periodic. The Boozer coordinates are used to find the logical coordinates via
    Newton's method.
    """
    rho = np.asarray(rho)
    theta_B = np.asarray(theta_B)
    zeta_B = np.asarray(zeta_B)

    squeeze = False
    if rho.ndim == 0:
        rho = np.array([rho])
        squeeze = True
    if rho.ndim != 1:
        raise ValueError("rho must be 1D")
    if not theta_B.ndim == zeta_B.ndim == 2 or theta_B.shape != zeta_B.shape:
        raise ValueError("theta_B and zeta_B must be 2D of the same shape (pol, tor)")
    if poloidal[1].ndim != 1 or poloidal[1].shape[0] != theta_B.shape[0]:
        raise ValueError(
            "poloidal data must be 1D and have the same length as the first dimension of theta_B/zeta_B"
        )
    if toroidal[1].ndim != 1 or toroidal[1].shape[0] != theta_B.shape[1]:
        raise ValueError(
            "toroidal data must be 1D and have the same length as the second dimension of theta_B/zeta_B"
        )

    ds = xr.Dataset(
        data_vars={
            "theta_B": (("pol", "tor"), theta_B),
            "zeta_B": (("pol", "tor"), zeta_B),
        },
        coords={
            "rho": ("rad", rho),
            poloidal[0]: ("pol", poloidal[1]),
            toroidal[0]: ("tor", toroidal[1]),
        },
    )

    # === Find the logical coordinates of the Boozer grid === #
    stacked = ds[["theta_B", "zeta_B"]].stack(tz=("pol", "tor"))
    tz_B = np.stack([stacked.theta_B, stacked.zeta_B], axis=0)
    sfl_boozer = state.get_boozer(rho, M, N, sincos=sincos)
    tz = state.get_boozer_angles(sfl_boozer, tz_B)
    stacked["theta"] = (("tz", "rad"), tz[0, :, :])
    stacked["zeta"] = (("tz", "rad"), tz[1, :, :])
    ds["theta"] = stacked["theta"].unstack("tz")
    ds["zeta"] = stacked["zeta"].unstack("tz")

    # === Metadata === #
    ds.rho.attrs["long_name"] = "Logical radial coordinate"
    ds.rho.attrs["symbol"] = r"\rho"
    ds.theta_B.attrs["long_name"] = "Boozer straight-fieldline poloidal angle"
    ds.theta_B.attrs["symbol"] = r"\theta_B"
    ds.zeta_B.attrs["long_name"] = "Boozer toroidal angle"
    ds.zeta_B.attrs["symbol"] = r"\zeta_B"
    ds.theta.attrs["long_name"] = "Logical poloidal angle"
    ds.theta.attrs["symbol"] = r"\theta"
    ds.zeta.attrs["long_name"] = "Logical toroidal angle"
    ds.zeta.attrs["symbol"] = r"\zeta"

    # === Indices === #
    # setting them earlier causes issues with the stacking / unstacking
    ds = ds.set_xindex("rho")
    ds = ds.set_xindex(poloidal[0])
    ds = ds.set_xindex(toroidal[0])
    ds = ds.drop_vars("pol")
    ds = ds.drop_vars("tor")

    if squeeze:
        ds = ds.squeeze("rad")

    return ds


# === Fourier Transform === #


def ev2ft(ev, quiet=False):
    m, n = None, None
    data = {}

    if "N_FP" not in ev.data_vars and not quiet:
        logging.warning("recommended quantity 'N_FP' not found in the provided dataset")

    for var in ev.data_vars:
        if ev[var].dims == ():  # scalar
            data[var] = ((), ev[var].data.item(), ev[var].attrs)

        elif ev[var].dims == ("rad",):  # profile
            data[var] = ("rad", ev[var].data, ev[var].attrs)

        elif {"pol", "tor"} <= set(ev[var].dims) <= {"rad", "pol", "tor"}:
            if "rad" in ev[var].dims:
                vft = []
                for r in ev.rad:
                    vft.append(
                        fourier.fft2d(ev[var].sel(rad=r).transpose("pol", "tor").data)
                    )
                vcos, vsin = map(np.array, zip(*vft))
                dims = ("rad", "m", "n")
            else:
                vcos, vsin = fourier.fft2d(ev[var].transpose("pol", "tor").data)
                dims = ("m", "n")

            if m is None:
                m, n = fourier.fft2d_modes(
                    vcos.shape[-2] - 1, vcos.shape[-1] // 2, grid=False
                )

            attrs = {
                k: v
                for k, v in ev[var].attrs.items()
                if k not in {"long_name", "symbol"}
            }
            data[f"{var}_mnc"] = (
                dims,
                vcos,
                dict(
                    long_name=f"{ev[var].attrs['long_name']}, cosine coefficients",
                    symbol=f"{{{ev[var].attrs['symbol']}}}_{{mn}}^c",
                )
                | attrs,
            )
            data[f"{var}_mns"] = (
                dims,
                vsin,
                dict(
                    long_name=f"{ev[var].attrs['long_name']}, sine coefficients",
                    symbol=f"{{{ev[var].attrs['symbol']}}}_{{mn}}^s",
                )
                | attrs,
            )

        elif "xyz" in ev[var].dims and not quiet:
            logging.info(f"skipping quantity '{var}' with cartesian components")

        elif not quiet:
            logging.info(f"skipping quantity '{var}' with dims {ev[var].dims}")

    if "rad" in ev.dims:
        coords = dict(rho=("rad", ev.rho.data, ev.rho.attrs))
    else:
        data["rho"] = ((), ev.rho.item(), ev.rho.attrs)
        coords = {}
    coords |= dict(
        m=(
            "m",
            m if m is not None else [],
            dict(long_name="poloidal mode number", symbol="m"),
        ),
        n=(
            "n",
            n if n is not None else [],
            dict(long_name="toroidal mode number", symbol="n"),
        ),
    )

    ft = xr.Dataset(data, coords=coords)
    if "rad" in ev.dims:
        ft = ft.set_xindex("rho")
    ft.attrs["fourier series"] = (
        "Assumes a fourier series of the form 'v(r, θ, ζ) = Σ v^c_mn(r) cos(m θ - n N_FP ζ) + v^s_mn(r) sin(m θ - n N_FP ζ)'"
    )
    return ft


def ft_autoremove(ft: xr.Dataset, drop=False, **tol_kwargs):
    """autoremove variables which are always close to zero (e.g. due to stellarator symmetry)"""
    selected = []
    for var in ft.data_vars:
        if set(ft[var].dims) >= {"m", "n"} and np.allclose(
            ft[var].data, 0, **tol_kwargs
        ):
            if not drop:
                ft[var] = ((), 0, ft[var].attrs)
            continue
        selected.append(var)
    if drop:
        return ft[selected]
    else:
        return ft
