from warnings import warn

import matplotlib.pyplot as pyplot
from numpy import any, isnan, ndarray, prod

from gvec.core.state import State

from .utils import (
    _design_subgrid,
    _extrapolate_axis,
    _get_coord_range,
    _subplots,
)

pyplot.rcParams.update({"text.usetex": True})


def _plot_line_quantities_from_dict(plotting_quantities, x_axis_value, subplot_grid, xlabel):
    """
    Plot the quantities from `plotting_quantities`.

    Return both the `matplotlib.pyplot.figure`
    """

    link_xaxis = False
    hide_inner_axis = False

    if (subplot_grid[0] != 1) & (subplot_grid[1] == 1):
        link_xaxis = True
    elif (subplot_grid[0] != 1) & (subplot_grid[1] != 1):
        link_xaxis = True
        hide_inner_axis = True

    # if axis is None:
    f, axs = _subplots(subplot_grid[0], subplot_grid[1], sharex=link_xaxis)
    # else:
    #     axs = axis

    if axs.size == 1:
        # If there is only one axis we plot all quantities on the single axis
        [
            axs[0].plot(x_axis_value, values, label=quantity)
            for (quantity, values) in plotting_quantities.items()
        ]
        axs[0].legend(plotting_quantities.keys())
    else:
        # Plot one value per axis
        for i, quantity in enumerate(plotting_quantities):
            axs[i].plot(x_axis_value, plotting_quantities[quantity])
            # Top right of each subplot will have a text box with the quantity being plotted
            axs[i].annotate(
                quantity,  # TODO: convert this to latex
                xy=(1, 1),
                xycoords="axes fraction",
                xytext=(-0.6, -0.6),
                textcoords="offset fontsize",
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(facecolor="white", edgecolor="black"),
            )

            if (not link_xaxis) | (i == (len(axs) - 1)):
                # If the axis are linked we only need the label on the last plot
                axs[i].set_xlabel(xlabel)
            if (hide_inner_axis) & (i - len(axs) + subplot_grid[1] >= 0):
                # If there are multiple plots we need to set the x-axis labels only on the bottom row
                axs[i].set_xlabel(xlabel)

    # if axis is None:
    return f, axs
    # else:
    # return axs


def _plot_line_quantities_from_xarray(
    evaluations, x_axis_values, quantities, subplot_grid, xlabel
):
    """
    Plot the quantities from the evaluations object assuming everything is 1D

    Return both the `matplotlib.pyplot.figure` and `Axis` object
    """

    link_xaxis = True
    hide_inner_axis = True

    f, axs = _subplots(subplot_grid[0], subplot_grid[1], sharex=link_xaxis)

    for i, quantity in enumerate(quantities):
        axs[i].plot(
            x_axis_values,
            evaluations[quantity],
        )
        # axs[i].annotate(
        #     f"{evaluations[quantity].attrs['long_name']}",
        #     xy=(1, 1),
        #     xycoords="axes fraction",
        #     xytext=(-0.6, -0.6),
        #     textcoords="offset fontsize",
        #     verticalalignment="top",
        #     horizontalalignment="right",
        #     bbox=dict(facecolor="white", edgecolor="black"),
        # )
        axs[i].set(
            ylabel=f"${evaluations[quantity].attrs['symbol']}$",
        )

        # If there are multiple plots we need to set the x-axis labels only on the bottom row
        if (hide_inner_axis) & (i - len(axs) + subplot_grid[1] >= 0):
            axs[i].set_xlabel(xlabel)

    return f, axs


def plot_radial_profile(
    state: State,
    nrho: int | ndarray = 100,
    quantities: str | list = ["iota", "p", "I_tor", "I_pol"],
    subplot_grid: list = [2, 2],
    xaxis="rho",
    # post_process: dict | dict = None, # TODO
):
    """
    Plot the radial profile of given equilibrium quantities.

    Parameters
    ----------
    state: GVEC state file
    nrho: int, numpy.ndarray
        The number of or specific 1D array of radial points to plot at. Default `100`
    quantities: str, list, optional
        Default is `["iota","p","I_tor","I_pol"]`.
    subplot_grid: list, optional
        The grid shape for the subplots. If not provided, the subplot grid will be determined automatically.
        Default is `[2,2]`.
    xaxis: `"rho"` or `"rho_squared"`, optional
        What quantity to plot on the x axis. Default is `"rho"`.
    post_process: dict, optional, --DEPRECIATED--
        `post_process` must be a dict with
            `post_process["quantity to remap"] = [<function>, "quantity name"]`
        such that `post_process["quantity to remap"][0]` is a callable <function> and `post_process["quantity to remap"][1]` is the heading for the subplot.
        The <function> _must_ return a 1D array.
        Default is `None`.

    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
    """

    if isinstance(quantities, str):
        # If plotting a single quantity convert it to a list
        quantities = [quantities]

    rho = _get_coord_range("rho", state.nfp, nrho)
    theta = _get_coord_range("theta", state.nfp, 0.0)
    zeta = _get_coord_range("zeta", state.nfp, 0.0)

    ev = state.evaluate(*quantities, rho=rho, theta=theta, zeta=zeta)

    ev = ev.sel(theta=0.0).sel(zeta=0.0)

    if xaxis == "rho_squared":
        xlabel = "$\\rho^2$"
        rho = rho**2
    elif xaxis == "rho":
        xlabel = "$\\rho$"
    else:
        raise ValueError("xaxis must be 'rho' or 'rho_squared'.")
    # TODO: Plot from xarray like in 020_stellarator.ipynb

    if prod(subplot_grid) != len(quantities):
        subplot_grid = _design_subgrid(len(quantities))
        warn("subplot_grid cannot fit the number of quantities, updating grid to fit.")

    f_ax = _plot_line_quantities_from_xarray(ev, rho, quantities, subplot_grid, xlabel)

    return f_ax


def plot_on_axis(
    state,
    nzeta: int | ndarray = 51,
    quantities: str | list = "mod_B",
    subplot_grid: list = [1, 1],
    post_process: dict = None,
):
    """
    Plot a equilibrium quantity (or list of) along the magnetic axis.

    Parameters
    ----------
    state: GVEC State file
    nzeta: int, ndarray, optional
        $\zeta$ resolution or array of points to plot at. Default 51.
    quantities: str, list, optional
        Default is "mod_B".
    subplot_grid: list, conditionally optional
        The grid shape for the subplots. Default is `[1,1]`. Required if `len(rho)>1`.
    post_process: dict, optional
        `post_process` must be a dict with
            `post_process["quantity to remap"] = [<function>, "quantity name"]`
        such that `post_process["quantity to remap"][0]` is a callable <function> and `post_process["quantity to remap"][1]` is the heading for the subplot.
        The <function> _must_ return a 1D array.
        Default is None

    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
    """
    if isinstance(quantities, str):
        # If plotting a single quantity convert it to a list
        quantities = [quantities]

    zeta = _get_coord_range("zeta", state.nfp, nzeta)

    # Use quadratic extrapolation to obtain values on axis.
    evaluations = _extrapolate_axis(state, quantities, zeta)

    # Since some derived quantities are not defined on axis we will
    # check if there are any NaNs in the dataset
    for quantity in quantities:
        if any(isnan(evaluations[quantity].data)):
            warn(f"{quantity} has NaNs despite running just off axis.")

    f, axs = _plot_line_quantities_from_xarray(
        evaluations, zeta, quantities, subplot_grid, "$\zeta$"
    )

    return f, axs
