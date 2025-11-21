from warnings import warn

import matplotlib.pyplot as pyplot
from numpy import any, array, isnan, linspace, max, min, ndarray, pi

from gvec.core.state import State

from .utils import _get_coord_range, _get_scalars_for_plotting, _subplots

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


def plot_radial_profile(
    state: State,
    nrho: int | ndarray,
    quantities: str | list = "mod_B",
    subplot_grid: list = [1, 1],
    post_process: dict | dict = None,
    xaxis="rho",
):
    """
    Plot the radial profile of given equilibrium quantities.

    Parameters
    ----------
    state: GVEC state file
    nrho: int, numpy.ndarray
        The number of or specific list of radial points to plot at.
    quantities: str, list, optional
        Default is "mod_B".
    subplot_grid: list, conditionally optional
        The grid shape for the subplots. Default is `[1,1]`. Required if `len(rho)>1`.
    share_axis: bool
        Default is True.
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

    rho = _get_coord_range("rho", state.nfp, nrho)
    theta = _get_coord_range("theta", state.nfp, 1)
    zeta = _get_coord_range("zeta", state.nfp, 1)

    ev = state.evaluate(*quantities, rho=rho, theta=theta, zeta=zeta)

    plotting_quantities = _get_scalars_for_plotting(ev, quantities, post_process, "rho")

    if xaxis == "rho_squared":
        xlabel = "$\\rho^2$"
        rho = rho**2
    else:
        xlabel = "$\\rho$"

    f_ax = _plot_line_quantities_from_dict(plotting_quantities, rho, subplot_grid, xlabel)

    return f_ax


def plot_on_axis(
    state,
    nzeta: int | ndarray = 51,
    quantities: str | list = "mod_B",
    subplot_grid: list = [1, 1],
    post_process: dict = None,
    near_axis=False,
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
    share_axis: bool
        Default is True.
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

    if near_axis:
        rho = 1e-8
    else:
        rho = 0.0
    theta = 0.0
    # theta = _get_coord_range("theta", state.nfp, 0.0)
    zeta = _get_coord_range("zeta", state.nfp, nzeta)

    ev = state.evaluate(*quantities, rho=rho, theta=theta, zeta=zeta)
    reevaluate = False

    for quantity in quantities:
        if any(isnan(ev[quantity].data)):
            if near_axis:
                warn(
                    f"{quantity} has NaNs despite running just off axis. Re-evaluating at rho=10^-4."
                )
                reevaluate = True
            else:
                warn(
                    f"NaNs detected in {quantity}. This value is likely zero on axis. Calling with `near_axis=True` may fix this issue by evaluating the quantities at 10^-8."
                )
    if reevaluate:
        ev = state.evaluate(*quantities, rho=1e-4, theta=theta, zeta=zeta)

    plotting_quantities = _get_scalars_for_plotting(ev, quantities, post_process, "zeta")

    f_ax = _plot_line_quantities_from_dict(plotting_quantities, zeta, subplot_grid, "$\\zeta$")

    return f_ax
