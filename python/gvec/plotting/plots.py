from numpy import array, max, min, ndarray, linspace, pi, any, isnan

from warnings import warn

import matplotlib.pyplot as pyplot

from gvec.core.state import State

from .utils import _get_coord_range, _get_scalars_for_plotting


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

    f, axs = pyplot.subplots(*subplot_grid, sharex=link_xaxis)

    # Ensure this is an array so we can iterate over it
    if not isinstance(axs, ndarray):
        axs = array([axs])
    if len(axs) != 1:
        # Otherwise we cannot iterate over this easily
        axs = axs.flatten()

    if axs.size == 1:
        # If there is only one axis we plot all quantities on the single axis
        ax_lines = [
            axs[0].plot(x_axis_value, values, label=quantity)
            for (quantity, values) in plotting_quantities.items()
        ]
        axs[0].legend(plotting_quantities.keys())
    else:
        # Plot one value per axis
        for i, quantity in enumerate(plotting_quantities):
            ax_lines = axs[i].plot(x_axis_value, plotting_quantities[quantity])
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

    return f, axs


def plot_poloidal_plane(
    state: State,
    nrho: int = 51,
    ntheta: int = 51,
    nzeta: int | float | ndarray = 1,
    quantities: str = "mod_B",
    subplot_grid: list = [1, 1],
    share_axis: bool = False,
    show_contours=False,  # TODO: update so a limited number of contours can be specified
    min_rho=0.0,
    # show_origin = False #TODO
):
    """
    plot_poloidal_slice

    Plot a poloidal slice with some equilibrium quantity on it. Defaults to plotting $\|\mathbf{B}\|$

    Parameters
    ----------
    state: GVEC state file
    nrho: int, optional
        The radial resolution of the slices. Default is 51
    ntheta: int, optional
        The poloidal resolution of the slices. Default is 51
    nzeta: int, float, ndarray, optional
        Poloidal slices to plot. Default is 1.
        - int: Number of equally spaced slices.
        - float: The specific slice.
        - ndarray: The specific slices.
    quantities: str
        The quantity to plot. Default is "mod_B"
    subplot_grid: list, conditionally optional
        The grid shape for the subplots. Default is `[1,1]`. Required if `len(nzeta)>1`.
    share_axis: bool
        If true, all subplots will share their `X1` and `X2` axis positions.
        TODO: Fix this so only the minor radius is shared.

    Returns
    -------
    `matplotlib.pyplot.figure` object
    """

    rho = _get_coord_range("rho", state.nfp, nrho, min_val=min_rho)
    theta = _get_coord_range("theta", state.nfp, ntheta)
    zeta = _get_coord_range("zeta", state.nfp, nzeta)

    ev = state.evaluate(*[quantities, "X1", "X2", "LA"], rho=rho, theta=theta, zeta=zeta)

    if any(isnan(ev[quantities])):
        # Sometimes on-axis quantities cannot be computed properly, to avoid erroring out we will just set them to zero
        #   and warn the user that they should adjust if their plots look weird
        ev[quantities].data[isnan(ev[quantities].data)] = 0.0
        warn(
            "NaNs detected in evaluated dataset, these have been set to 0. It is possible that on-axis quantities cannot be evaluated, a minimum rho value can be set with the min_rho input."
        )

    is_scalar = len(ev[quantities].shape) == 1

    # Currently we do not automatically adapt the grid
    # TODO: Try and do this automatically??
    if len(zeta) != 1:
        if subplot_grid == [1, 1]:
            raise ValueError(
                "Since more than one plane is requested please specify the 'subplot_grid' parameter for the subplot grid sizing"
            )

    f, axs = pyplot.subplots(*subplot_grid, sharex=share_axis, sharey=share_axis)

    if not isinstance(axs, ndarray):
        axs = array(axs)  # Makes it easier/possible to loop over
    axs = axs.flatten()

    # All plots will share the same colour scale
    if is_scalar:
        colour_scale = (
            min(ev[quantities]),
            max(ev[quantities]),
        )
    else:
        colour_scale = (
            min(ev[quantities].data[:, :, :]),
            max(ev[quantities].data[:, :, :]),
        )

    # Loop and plot all axes
    for i, ax in enumerate(axs):
        # Check if the quantity is a scalar
        if not is_scalar:  # not a scalar
            plotting_quantity = ev[quantities].data[:, :, i].flatten()
        else:
            plotting_quantity = ev.X1[:, :, i] * 0.0 + ev[quantities]
            plotting_quantity = plotting_quantity.data.flatten()
        f_ax = ax.tricontourf(
            ev.X1.data[:, :, i].flatten(),
            ev.X2.data[:, :, i].flatten(),
            plotting_quantity,
            vmin=colour_scale[0],
            vmax=colour_scale[1],
        )
        if share_axis:
            ax.label_outer()  # Removes any axis labels on subplots on the interior of the grid
        if show_contours:
            # We should plot the rho and \theta^\star contours
            rho_contour = ev.X1[:, :, i] * 0.0 + ev.rho
            rho_levels_vis = linspace(0, 1 - 1e-10, 11)
            ax.contour(
                ev.X1.data[:, :, i],
                ev.X2.data[:, :, i],
                rho_contour,
                rho_levels_vis,
                colors="black",
            )

            theta_star_contour = ev.LA[:, :, i] + ev.theta
            theta_levels_vis = linspace(0, 2 * pi, 10, endpoint=False)
            ax.contour(
                ev.X1.data[:, :, i],
                ev.X2.data[:, :, i],
                theta_star_contour,
                theta_levels_vis,
                colors="red",
            )
            ax.annotate(
                f"$\\zeta={zeta[i]:.2f}$",
                xy=(1, 1),
                xycoords="axes fraction",
                xytext=(-0.6, -0.6),
                textcoords="offset fontsize",
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(facecolor="white", edgecolor="black"),
            )
            if share_axis:
                ax.set(xlabel="X1", ylabel="X2")
                ax.label_outer()

    if share_axis:
        # f.subplots_adjust(wspace=0.05,hspace=0.05)
        f.tight_layout()

    # Adding colourbar
    colourbar_axis = f.add_axes([0.85, 0.15, 0.05, 0.7])
    # Colourbar starts at 0.85
    f.subplots_adjust(right=0.83)

    f.colorbar(f_ax, cax=colourbar_axis, label=quantities)

    return f


def plot_flux_surface(
    state: State,
    nzeta: int,
    ntheta: int,
    rho: float | ndarray = 1.0,
    quantities: str = "mod_B",  # TODO: Allow for multiple quantities instead of multiple surfaces
    subplot_grid: list = [1, 1],
    share_axis: bool = False,
):
    """
    Plot an equilibrium quantity on a toroidal surface at a given `rho` value

    Parameters
    ----------
    state: GVEC state file
    nzeta: int
    ntheta: int
    rho: float, numpy.ndarray, optional
        The flux surface label(s) to plot. Default is 1.0.
    quantities: str, optional
        TODO: Allow for multiple quantities instead of multiple surfaces
    subplot_grid: list, conditionally optional
        The grid shape for the subplots. Default is `[1,1]`. Required if `len(rho)>1`.
    share_axis: bool
        If true, all subplots will share their `x` and `y` axes.

    Returns
    -------
    `matplotlib.pyplot.figure` object
    """

    # Currently we do not automatically adapt the grid
    # TODO: Try and do this automatically
    if isinstance(rho, float):
        if subplot_grid != [1, 1]:
            subplot_grid = [1, 1]
    elif isinstance(rho, (ndarray, list)):
        if subplot_grid == [1, 1]:
            raise ValueError(
                "Since more than one plane is requested please specify the 'subplot_grid' parameter for the subplot grid sizing"
            )

    zeta = _get_coord_range("zeta", state.nfp, nzeta)
    rho = _get_coord_range("rho", state.nfp, rho)
    theta = _get_coord_range("theta", state.nfp, ntheta)

    ev = state.evaluate(quantities, rho=rho, theta=theta, zeta=zeta)

    f, axs = pyplot.subplots(*subplot_grid, sharex=share_axis, sharey=share_axis)

    if not isinstance(axs, ndarray):
        axs = array([axs])  # Makes it easier/possible to loop over

    # Loop over and plot all axes
    for i, ax in enumerate(axs.reshape(-1)):
        ev_i = ev.isel(rad=i)
        f_ax = ax.contour(
            ev_i.theta.data,
            ev_i.zeta.data,
            ev_i[quantities].transpose("pol", "tor").data,
        )
        #   vmin=colour_scale[0], vmax=colour_scale[1])
        if share_axis:
            ax.label_outer()  # Removes any axis labels on subplots on the interior of the grid

    return f


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
    `matplotlib.pyplot.figure` object
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

    f, axs = _plot_line_quantities_from_dict(plotting_quantities, rho, subplot_grid, xlabel)

    return f


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
    `matplotlib.pyplot.figure` object
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

    f, axs = _plot_line_quantities_from_dict(
        plotting_quantities, zeta, subplot_grid, "$\\zeta$"
    )

    return f
