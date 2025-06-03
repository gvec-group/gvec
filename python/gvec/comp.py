# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""GVEC Postprocessing - Compute Functions"""

from typing import TypeAlias, Literal
from collections.abc import Sequence, Collection, Mapping, MutableMapping, Callable
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
    quantities: None | str | Collection[str] = None,
    requirements: Collection[str] = (),
    integration: Collection[str] = (),
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
    *quantities: Collection[str],
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

CoordinateSpec: TypeAlias = int | float | xr.DataArray | np.ndarray | Sequence


def Evaluations(
    rho: Literal["int"] | CoordinateSpec | None = "int",
    theta: Literal["int"] | CoordinateSpec | None = "int",
    zeta: Literal["int"] | CoordinateSpec | None = "int",
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
        case str() if rho == "int":
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
        case np.ndarray() | Sequence():
            coords["rho"] = ("rad", rho)
        case int() as num:
            coords["rho"] = ("rad", np.linspace(0, 1, num))
            coords["rho"][1][0] = (
                0.1 * coords["rho"][1][1]
            )  # avoid numerical issues at the magnetic axis
        case None:
            pass
        case _:
            raise ValueError(f"Could not parse rho, got {rho}.")
    match theta:
        case xr.DataArray():
            coords["theta"] = theta
        case str() if theta == "int":
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
        case np.ndarray() | Sequence():
            coords["theta"] = ("pol", theta)
        case int() as num:
            coords["theta"] = ("pol", np.linspace(0, 2 * np.pi, num, endpoint=False))
        case None:
            pass
        case _:
            raise ValueError(f"Could not parse theta, got {theta}.")
    match zeta:
        case xr.DataArray():
            coords["zeta"] = zeta
        case str() if zeta == "int":
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
        case np.ndarray() | Sequence():
            coords["zeta"] = ("tor", zeta)
        case int() as num:
            if nfp is None:
                raise ValueError("Automatic bounds for zeta require `nfp`.")
            coords["zeta"] = (
                "tor",
                np.linspace(0, 2 * np.pi / nfp, num, endpoint=False),
            )
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
    rho: CoordinateSpec,
    theta_B: CoordinateSpec,
    zeta_B: CoordinateSpec,
    state: State,
    **boozer_kwargs,
):
    """Create a xarray.Dataset for evaluating GVEC with a grid in Boozer Coordinates.

    Adds LA, NU_B and their angular derivatives to the dataset, as recomputed for the Boozer transform.

    Parameters
    ----------
    rho : int | float | 1D array (DataArray, ndarray, list)
        The specification of the radial, radius-like coordinate.
    theta_B : int | float | 1D array (DataArray, ndarray, list)
        The specification of the poloidal, angle-like Boozer coordinate.
    zeta_B : int | float | 1D array (DataArray, ndarray, list)
        The specification of the poloidal, angle-like Boozer coordinate.
    state : State
        The gvec.State object to create the grid for. Used to perform the Boozer transform.
    """
    match rho:
        case xr.DataArray():
            rho = rho
        case np.ndarray | Sequence():
            rho = ("rad", np.asarray(rho))
        case int():
            rho = ("rad", np.linspace(0, 1, rho + 1)[1:])
        case float():
            rho = ("rad", np.array([rho]))
        case _:
            raise ValueError(f"Could not parse rho, got {rho}.")
    match theta_B:
        case xr.DataArray():
            theta_B = theta_B
        case np.ndarray | Sequence():
            theta_B = ("pol", np.asarray(theta_B))
        case int():
            theta_B = ("pol", np.linspace(0, 2 * np.pi, theta_B, endpoint=False))
        case float():
            theta_B = ("pol", np.array([theta_B]))
        case _:
            raise ValueError(f"Could not parse theta_B, got {theta_B}.")
    match zeta_B:
        case xr.DataArray():
            zeta_B = zeta_B
        case np.ndarray | Sequence():
            zeta_B = ("tor", np.asarray(zeta_B))
        case int():
            zeta_B = (
                "tor",
                np.linspace(0, 2 * np.pi / state.nfp, zeta_B, endpoint=False),
            )
        case float():
            zeta_B = ("tor", np.array([zeta_B]))
        case _:
            raise ValueError(f"Could not parse zeta_B, got {zeta_B}.")

    ds = xr.Dataset(
        coords=dict(
            rho=rho,
            theta_B=theta_B,
            zeta_B=zeta_B,
        )
    )

    # === Find the logical coordinates of the Boozer grid === #
    stacked = ds[["theta_B", "zeta_B"]].stack(tz=("pol", "tor"))
    tz_B = np.stack([stacked.theta_B, stacked.zeta_B], axis=0)
    sfl_boozer = state.get_boozer(ds.rho, **boozer_kwargs)
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
    ds = add_Boozer_LA_NU(ds, state, sfl_boozer)

    return ds


def add_Boozer_LA_NU(ds: xr.Dataset, state: State, sfl_boozer):
    """Add the LA and NU_B variables as computed by the boozer transform to the dataset.

    Helper function for EvaluationsBoozer and related methods.
    """
    # Flatten theta, zeta
    theta = ds.theta.transpose("rad", "pol", "tor").values.reshape(ds.rad.size, -1)
    zeta = ds.zeta.transpose("rad", "pol", "tor").values.reshape(ds.rad.size, -1)

    outputs_la = []
    outputs_nu = []
    # Sequence (list) of sfl_boozer (for each surface)
    if isinstance(sfl_boozer, Sequence):
        for r, rho in enumerate(ds.rho.data):
            thetazeta = np.stack([theta[r, :], zeta[r, :]], axis=0)
            outputs_la.append(
                state.evaluate_boozer_list_tz_all(sfl_boozer[r], "LA", [0], thetazeta)
            )
            outputs_nu.append(
                state.evaluate_boozer_list_tz_all(sfl_boozer[r], "NU", [0], thetazeta)
            )
    # Single sfl_boozer - compute base on each radial position
    else:
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
    return ds


def EvaluationsBoozerCustom(
    rho: CoordinateSpec,
    theta_B: CoordinateSpec,
    zeta_B: CoordinateSpec,
    state: State,
    **boozer_kwargs,
):
    """Create an Evaluations dataset with a custom Boozer grid (e.g. fieldline aligned).

    This factory function assumes that the target grid still has a poloidal-like and toroidal-like direction
    and that the Boozer coordinates of each grid point are known. They do not have to lie within a single field
    period, nor do they have to be periodic. The Boozer coordinates are used to find the logical coordinates via
    Newton's method.

    If a 2D or 3D array for theta_B or zeta_B is passed, the corresponding coordinate for the pol/tor dimension
    needs to be set manually afterwards (e.g. `ev["alpha"] = ("pol", values)` and `ev = ev.set_coords("alpha").set_xindex("alpha")`).

    Parameters
    ----------
    rho : int | float | 1D array (DataArray, ndarray, list)
        The specification of the radial, radius-like coordinate.
    theta_B : int | float | 1D, 2D or 3D array (DataArray, ndarray, list)
        The specification of the poloidal, angle-like Boozer coordinate.
    zeta_B : int | float | 1D, 2D or 3D array (DataArray, ndarray, list)
        The specification of the poloidal, angle-like Boozer coordinate.
    state : State
        The gvec.State object to create the grid for. Used to perform the Boozer transform.
    """
    match rho:
        case xr.DataArray():
            rho = rho
        case np.ndarray() | Sequence():
            rho = np.asarray(rho)
            if rho.ndim != 1:
                raise ValueError(f"rho can only be 1D, but is {rho.ndim}D.")
            rho = ("rad", rho)
        case int():
            rho = ("rad", np.linspace(0, 1, rho + 1)[1:])
        case float():
            rho = ("rad", np.array([rho]))
        case _:
            raise ValueError(f"Could not parse rho, got {rho}.")
    match theta_B:
        case xr.DataArray():
            theta_B = theta_B
        case np.ndarray() | Sequence():
            theta_B = np.asarray(theta_B)
            if theta_B.ndim == 1:
                theta_B = ("pol", theta_B)
            elif theta_B.ndim == 2:
                theta_B = (("pol", "tor"), theta_B)
            elif theta_B.ndim == 3:
                theta_B = (("rad", "pol", "tor"), theta_B)
            else:
                raise ValueError(f"theta_B can only be 1D, 2D, 3D, not {theta_B.ndim}D")
        case int():
            theta_B = ("pol", np.linspace(0, 2 * np.pi, theta_B, endpoint=False))
        case float():
            theta_B = ("pol", np.array([theta_B]))
        case _:
            raise ValueError(f"Could not parse theta_B, got {theta_B}.")
    match zeta_B:
        case xr.DataArray():
            zeta_B = zeta_B
        case np.ndarray() | Sequence():
            zeta_B = np.asarray(zeta_B)
            if zeta_B.ndim == 1:
                zeta_B = ("tor", zeta_B)
            elif zeta_B.ndim == 2:
                zeta_B = (("pol", "tor"), zeta_B)
            elif zeta_B.ndim == 3:
                zeta_B = (("rad", "pol", "tor"), zeta_B)
            else:
                raise ValueError(f"zeta_B can only be 1D, 2D, 3D, not {zeta_B.ndim}D")
        case float():
            zeta_B = ("tor", np.array([zeta_B]))
        case int():
            zeta_B = (
                "tor",
                np.linspace(0, 2 * np.pi / state.nfp, zeta_B, endpoint=False),
            )
        case _:
            raise ValueError(f"Could not parse zeta_B, got {zeta_B}.")

    ds = xr.Dataset(
        coords=dict(
            rho=rho,
        ),
        data_vars=dict(
            theta_B=theta_B,
            zeta_B=zeta_B,
        ),
    )

    # === Find the logical coordinates of the Boozer grid === #
    if "rad" in ds.theta_B.dims or "rad" in ds.zeta_B.dims:  # 3D
        theta = []
        zeta = []
        sfls = []
        for rad, rho in enumerate(ds.rho):
            dsr = ds.isel(rad=rad)
            stacked = dsr[["theta_B", "zeta_B"]]
            stacked = stacked.broadcast_like(stacked).stack(tz=("pol", "tor"))
            tz_B = np.stack([stacked.theta_B, stacked.zeta_B], axis=0)
            sfl_boozer = state.get_boozer([rho], **boozer_kwargs)
            tz = state.get_boozer_angles(sfl_boozer, tz_B)
            sfls.append(sfl_boozer)
            stacked["theta"] = (("tz", "rad"), tz[0, :, :])
            stacked["zeta"] = (("tz", "rad"), tz[1, :, :])
            theta.append(stacked["theta"].unstack("tz"))
            zeta.append(stacked["zeta"].unstack("tz"))
        ds["theta"] = xr.concat(theta, dim="rad")
        ds["zeta"] = xr.concat(zeta, dim="rad")
        ds = add_Boozer_LA_NU(ds, state, sfls)
    else:  # 2D
        stacked = ds[["theta_B", "zeta_B"]].stack(tz=("pol", "tor"))
        tz_B = np.stack([stacked.theta_B, stacked.zeta_B], axis=0)
        sfl_boozer = state.get_boozer(ds.rho, **boozer_kwargs)
        tz = state.get_boozer_angles(sfl_boozer, tz_B)
        stacked["theta"] = (("tz", "rad"), tz[0, :, :])
        stacked["zeta"] = (("tz", "rad"), tz[1, :, :])
        ds["theta"] = stacked["theta"].unstack("tz")
        ds["zeta"] = stacked["zeta"].unstack("tz")
        ds = add_Boozer_LA_NU(ds, state, sfl_boozer)

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
    ds = ds.drop_vars("pol")
    ds = ds.drop_vars("tor")

    if ds.theta_B.dims == ("pol",):
        ds = ds.set_coords("theta_B").set_xindex("theta_B")
    if ds.zeta_B.dims == ("tor",):
        ds = ds.set_coords("zeta_B").set_xindex("zeta_B")

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
