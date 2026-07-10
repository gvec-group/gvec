# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from typing import Literal
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np

from gvec.plotting.utils import (
    _design_subgrid,
    _subplots,
    _symbol_check,
)


def _plot_line_quantities_from_xarray(
    evaluations, x_axis_values, quantities, subplot_grid, xlabel, plot_kwargs
):
    """
    Plot the quantities from the evaluations object assuming everything is 1D

    Return both the ``matplotlib.pyplot.figure`` and ``Axis`` object
    """

    link_xaxis = True
    hide_inner_axis = True

    if subplot_grid is None:
        # If subplot_grid not predetermined we will work it out
        subplot_grid = _design_subgrid(len(quantities))

    fig, axs = _subplots(subplot_grid, link_xaxis, False, **plot_kwargs)

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
    if hide_inner_axis:
        axs_2d = np.atleast_2d(axs).reshape(subplot_grid)
        for row in range(axs_2d.shape[1]):
            axs_2d[-1, row].set_xlabel(f"${xlabel}$")

    return fig, axs


def plot_radial_profile(
    state,
    quantities: str | list[str] = ["iota", "p", "I_tor", "I_pol"],
    nrho: int | np.ndarray = 101,
    subplot_grid: list[int] | None = None,
    xaxis: Literal["rho", "rho_squared"] = "rho",
    n_rationals: int = 3,
    plot_kwargs: dict = dict(),
):
    """
    Plot the radial profile of given equilibrium quantities.

    Parameters
    ----------
    state : GVEC state object
    quantities : str, list[str], optional
        Default is ``["iota","p","I_tor","I_pol"]``.
    nrho : int, numpy.ndarray
        The number of or specific 1D array of radial points to plot at.
        Default is ``101``
    subplot_grid : list[int], None, optional
        The grid shape for the subplots. If ``None``, grid will be automatically determined.
        Default is ``None``.
    xaxis : ``"rho"`` or ``"rho_squared"``, optional
        What quantity to plot on the x axis.
        Default is ``"rho"``.
    n_rationals : int, optional
        If non-zero, show the largest $n$ rationals on any ``iota`` profile plot.
        Default ``3``
    plot_kwargs: dict, optional
        Any ``**kwargs`` to send to the ``plt.figure()`` function.
        For example ``plot_kwargs={'figsize': (8,8)}``. See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html>`_ for a list of kwargs.

    Returns
    -------
    ``matplotlib.pyplot.figure`` object and ``numpy.ndarray`` of ``matplotlib.axis._axis.Axes`` object(s).
    """

    if isinstance(quantities, str):
        # If plotting a single quantity convert it to a list
        quantities = [quantities]

    evaluations = state.evaluate(*quantities, rho=nrho, theta=0.0, zeta=0.0)

    evaluations = evaluations.sel(theta=0.0).sel(zeta=0.0)

    evaluations = _symbol_check(evaluations, quantities)

    rho = evaluations.rho.data

    if xaxis == "rho_squared":
        xlabel = r"\rho^2"
        rho = rho**2
    elif xaxis == "rho":
        xlabel = r"\rho"
    else:
        raise ValueError("xaxis must be 'rho' or 'rho_squared'.")

    fig, axs = _plot_line_quantities_from_xarray(
        evaluations, rho, quantities, subplot_grid, xlabel, plot_kwargs
    )

    if n_rationals:  # if n_rationals > 0 or not None
        for ax, quantity in zip(np.asarray(axs).flat, quantities):
            if quantity == "iota":
                _add_rationals_to_iota_plot(ax, n_rationals=n_rationals, nfp=state.nfp)

    return fig, axs


def plot_on_axis(
    state,
    quantities: str | list[str] = "mod_B",
    nzeta: int | np.ndarray = 51,
    subplot_grid: list[int] | None = None,
    plot_kwargs: dict = dict(),
):
    """
    Plot a equilibrium quantity (or list of) along the magnetic axis.
    Note that the quantities are always evaluated off-axis (``rho=[1.1e-4,2.2e-4,3.3e-4]``) and extrapolated quadratically to ``rho=0`` (and averaged over ``theta``).

    Parameters
    ----------
    state : GVEC state object
    quantities : str, list[str], optional
        Default is ``"mod_B"``.
    nzeta : int, numpy.ndarray, optional
        ``zeta`` resolution or array of points to plot at.
        Default is ``51``.
    subplot_grid : list[int], None, optional
        The grid shape of ``[nrow,ncol]`` for the subplots. If ``None``, grid will be automatically determined.
        Default is ``None``.
    plot_kwargs: dict, optional
        Any ``**kwargs`` to send to the ``plt.figure()`` function.
        For example ``plot_kwargs={'figsize': (8,8)}``. See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html>`_ for a list of kwargs.

    Returns
    -------
    ``matplotlib.pyplot.figure`` object and ``numpy.ndarray`` of ``matplotlib.axis._axis.Axes`` object(s).
    """
    if isinstance(quantities, str):
        # If plotting a single quantity convert it to a list
        quantities = [quantities]

    # Use quadratic extrapolation to obtain values on axis.
    evaluations = state.evaluate_on_axis(*quantities, theta="int", zeta=nzeta).mean("pol")

    evaluations = _symbol_check(evaluations, quantities)

    # Since some derived quantities are not defined on axis we will
    # check if there are any NaNs in the dataset
    for quantity in quantities:
        if np.any(np.isnan(evaluations[quantity].data)):
            warn(f"{quantity} has NaNs despite extrapolation.")

    zeta = evaluations.zeta.data

    fig, axs = _plot_line_quantities_from_xarray(
        evaluations, zeta, quantities, subplot_grid, r"\zeta", plot_kwargs
    )

    return fig, axs


# === 1D utility functions ===#


def _add_rationals_to_iota_plot(
    ax,
    limits: tuple[float, float] | None = None,
    n_rationals: int = 3,
    nfp: int = 1,
    labels_inside: bool = False,
):
    """
    Add rationals as a secondary y-axis to an iota profile plot. Denominator is restricted to multiples of number of field periods.
    Parameters
    ----------
    state : GVEC state object
    ax : `matplotlib` axis object
    limits : tuple[float, float], optional
        The maximum values of the nominator and denominator. If ``None`` then automatically get the y-limits from the axis.
        Default ``None``
    n_rationals : int, optional
        The maximum number of rationals to add to the plot.
        Default ``3``
    nfp : int, optional
        The number of field periods. Only rationals with numerator that are multiples of this will be plotted.
        Default ``1``
    labels_inside : bool, default: False
        If True, the rational labels will be placed inside the plot area.
    """
    from math import gcd

    if limits is None:
        limits = ax.get_ylim()

    # If we have negative only limits flip them before the computation then flip q/p back afterwards
    pre_mult = 1
    if (limits[0] < 0) & (limits[1] <= 0):
        pre_mult = -1
        limits = limits[::-1]
    limits = [pre_mult * limit for limit in limits]

    rationals = []

    for n in range(1, 1000):
        for m in range(1, 1000):
            if n / m > limits[1] / nfp:
                continue
            if n / m < limits[0] / nfp:
                break
            if gcd(n, m) != 1:
                continue
            rationals.append((n * nfp, m))
            if len(rationals) >= n_rationals:
                break
        if len(rationals) >= n_rationals:
            break

    rationals = sorted(set(rationals), key=lambda x: x[0] / x[1])

    secax = ax.secondary_yaxis("right")
    if pre_mult == 1:
        # If iota is positive this is fine
        secax.set(
            yticks=[n / m for n, m in rationals],
            yticklabels=[f"$\\frac{{{n}}}{{{m}}}$" for n, m in rationals],
        )
        if labels_inside:
            secax.tick_params(direction="in", pad=-10)
    elif pre_mult == -1:
        # if not we need to flip the sign on the ticks, add the sign
        #   and move them inside the plot further
        secax.set(
            yticks=[pre_mult * n / m for n, m in rationals],
            yticklabels=[f"$-\\frac{{{n}}}{{{m}}}$" for n, m in rationals],
        )
        if labels_inside:
            secax.tick_params(direction="in", pad=-17)

    for n, m in rationals:
        ax.axhline(pre_mult * n / m, color="gray", linestyle="--", alpha=0.2)

    return rationals, secax
