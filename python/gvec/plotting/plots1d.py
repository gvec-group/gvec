# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from warnings import warn

# import matplotlib.pyplot as pyplot
import numpy as np
from numpy import ndarray  # type checking

from gvec.core.state import State
from gvec.plotting.utils import (
    _design_subgrid,
    _extrapolate_axis,
    _subplots,
    _symbol_check,
)


def _plot_line_quantities_from_xarray(
    evaluations, x_axis_values, quantities, subplot_grid, xlabel
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

    f, axs = _subplots(subplot_grid, sharex=link_xaxis)

    for i, quantity in enumerate(quantities):
        axs[i].plot(
            x_axis_values,
            evaluations[quantity],
        )
        axs[i].set(
            ylabel=f"${evaluations[quantity].attrs['symbol']}$",
        )

        # If there are multiple plots we need to set the x-axis labels only on the bottom row
        if (hide_inner_axis) & (i - len(axs) + subplot_grid[1] >= 0):
            axs[i].set_xlabel(f"${xlabel}$")

    return f, axs


def plot_radial_profile(
    state: State,
    quantities: str | list = ["iota", "p", "I_tor", "I_pol"],
    nrho: int | ndarray = 100,
    subplot_grid: list | None = None,
    xaxis="rho",
    post_process: dict | None = None,
):
    """
    Plot the radial profile of given equilibrium quantities.

    Parameters
    ----------
    state : GVEC state file
    quantities : str, list, optional
        Default is `["iota","p","I_tor","I_pol"]`.
    nrho : int, numpy.ndarray
        The number of or specific 1D array of radial points to plot at.
        Default is `100`
    subplot_grid : list, None, optional
        The grid shape for the subplots. If `None`, grid will be automatically determined.
        Default is `None`.
    xaxis : `"rho"` or `"rho_squared"`, optional
        What quantity to plot on the x axis. Default is `"rho"`.

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

    f_ax = _plot_line_quantities_from_xarray(evaluations, rho, quantities, subplot_grid, xlabel)

    return f_ax


def plot_on_axis(
    state,
    quantities: str | list = "mod_B",
    nzeta: int | ndarray = 51,
    subplot_grid: list | None = None,
):
    """
    Plot a equilibrium quantity (or list of) along the magnetic axis.

    Parameters
    ----------
    state : GVEC State file
    quantities : str, list, optional
        Default is "mod_B".
    nzeta : int, ndarray, optional
        $\zeta$ resolution or array of points to plot at.
        Default is `51`.
    subplot_grid : list, None, optional
        The grid shape for the subplots. If `None`, grid will be automatically determined.
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
        evaluations, zeta, quantities, subplot_grid, "$\zeta$"
    )

    return f, axs
