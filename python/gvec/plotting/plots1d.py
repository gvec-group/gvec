# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from typing import Literal
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np

from gvec.core.state import State
from gvec.plotting.utils import (
    _deco_usetex,
    _design_subgrid,
    _extrapolate_axis,
    _subplots,
    _symbol_check,
)


def _plot_line_quantities_from_xarray(
    evaluations, x_axis_values, quantities, subplot_grid, xlabel, plot_kwargs
):
    """
    Plot the quantities from the evaluations object assuming everything is 1D

    Return both the `matplotlib.pyplot.figure` and `Axis` object
    """

    link_xaxis = True
    hide_inner_axis = True

    if subplot_grid is None:
        # If subplot_grid not predetermined we will work it out
        subplot_grid = _design_subgrid(len(quantities))

    f, axs = _subplots(subplot_grid, link_xaxis, False, **plot_kwargs)

    # for i, quantity in enumerate(quantities):
    for i, ax in enumerate(np.asarray(axs).flat):
        ax.plot(
            x_axis_values,
            evaluations[quantities[i]],
        )
        ax.set(
            ylabel=f"${evaluations[quantities[i]].attrs['symbol']}$",
        )

        # If there are multiple plots we need to set the x-axis labels only on the bottom row
        if (hide_inner_axis) & (i - len(axs) + subplot_grid[1] >= 0):
            ax.set_xlabel(f"${xlabel}$")

        ax.gvec_quantity = quantities[i]

    return f, axs


@_deco_usetex
def plot_radial_profile(
    state: State,
    quantities: str | list[str] = ["iota", "p", "I_tor", "I_pol"],
    nrho: int | np.ndarray = 100,
    subplot_grid: list[int] | None = None,
    xaxis: Literal["rho", "rho_squared"] = "rho",
    plot_kwargs: dict = dict(),
):
    """
    Plot the radial profile of given equilibrium quantities.

    Parameters
    ----------
    state : GVEC state object
    quantities : str, list[str], optional
        Default is `["iota","p","I_tor","I_pol"]`.
    nrho : int, numpy.ndarray
        The number of or specific 1D array of radial points to plot at.
        Default is `100`
    subplot_grid : list[int], None, optional
        The grid shape for the subplots. If `None`, grid will be automatically determined.
        Default is `None`.
    xaxis : `"rho"` or `"rho_squared"`, optional
        What quantity to plot on the x axis.
        Default is `"rho"`.

    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
    """

    if isinstance(quantities, str):
        # If plotting a single quantity convert it to a list
        quantities = [quantities]

    evaluations = state.evaluate(*quantities, rho=nrho, theta=0.0, zeta=0.0)

    evaluations = evaluations.sel(theta=0.0).sel(zeta=0.0)

    evaluations = _symbol_check(evaluations, quantities)

    rho = evaluations.rho.data

    if xaxis == "rho_squared":
        xlabel = "$\\rho^2$"
        rho = rho**2
    elif xaxis == "rho":
        xlabel = "$\\rho$"
    else:
        raise ValueError("xaxis must be 'rho' or 'rho_squared'.")

    f, axs = _plot_line_quantities_from_xarray(
        evaluations, rho, quantities, subplot_grid, xlabel, plot_kwargs
    )

    return f, axs


@_deco_usetex
def plot_on_axis(
    state: State,
    quantities: str | list[str] = "mod_B",
    nzeta: int | np.ndarray = 51,
    subplot_grid: list[int] | None = None,
    plot_kwargs: dict = dict(),
):
    """
    Plot a equilibrium quantity (or list of) along the magnetic axis.

    Parameters
    ----------
    state : GVEC state object
    quantities : str, list[str], optional
        Default is "mod_B".
    nzeta : int, numpy.ndarray, optional
        $\zeta$ resolution or array of points to plot at.
        Default is `51`.
    subplot_grid : list[int], None, optional
        The grid shape of `[nrow,ncol]` for the subplots. If `None`, grid will be automatically determined.
        Default is `None`.

    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
    """
    if isinstance(quantities, str):
        # If plotting a single quantity convert it to a list
        quantities = [quantities]

    # Use quadratic extrapolation to obtain values on axis.
    evaluations = _extrapolate_axis(state, quantities, nzeta)

    evaluations = _symbol_check(evaluations, quantities)

    # Since some derived quantities are not defined on axis we will
    # check if there are any NaNs in the dataset
    for quantity in quantities:
        if np.any(np.isnan(evaluations[quantity].data)):
            warn(f"{quantity} has NaNs despite running just off axis.")

    zeta = evaluations.zeta.data

    f, axs = _plot_line_quantities_from_xarray(
        evaluations, zeta, quantities, subplot_grid, "$\zeta$", plot_kwargs
    )

    return f, axs


# === 1D utility functions ===#


def _add_rationals_to_iota_plot(
    state,
    ax,
    limits: tuple[float, float] | None = None,
    n_rationals: int | None = 3,
    n_max: int = 10,
):
    """
    Add high-order rationals as a secondary y-axis to an iota profile plot.
    """
    from math import gcd

    if limits is None:
        limits = ax.get_ylim()

    rationals = []
    for n in range(1, n_max + 1):
        for m in range(1, n * state.nfp + 1):
            if gcd(m, n) != 1:
                continue
            if limits[0] <= m / (n * state.nfp) <= limits[1]:
                rationals.append((m, n * state.nfp))
        if len(rationals) >= n_rationals:
            break
    rationals = sorted(rationals, key=lambda x: x[0] / x[1])

    secax = ax.secondary_yaxis("right")
    secax.set(
        yticks=[m / n for m, n in rationals],
        yticklabels=[f"$\\frac{{{m}}}{{{n}}}$" for m, n in rationals],
    )
    for m, n in rationals:
        ax.axhline(m / n, color="gray", linestyle="--", alpha=0.2)

    return rationals, secax


def add_iota_rationals(
    state,
    figaxs,
    limits: tuple[float, float] | None = None,
    n_rationals: int | None = 3,
    n_max: int = 10,
):
    if isinstance(figaxs, np.ndarray):
        for ax in figaxs.flat:
            if ax.gvec_quantity == "iota":
                _add_rationals_to_iota_plot(state, ax, limits, n_rationals, n_max)
    elif isinstance(figaxs, plt.Axes):
        _add_rationals_to_iota_plot(state, figaxs, limits, n_rationals, n_max)
    elif isinstance(figaxs, tuple):
        for fa in figaxs:
            if isinstance(fa, plt.Axes | np.ndarray):
                add_iota_rationals(state, fa, limits, n_rationals, n_max)
