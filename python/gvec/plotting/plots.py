from numpy import array, max, min, ndarray

import matplotlib.pyplot as pyplot

from gvec.core.state import State

from .utils import _get_coord_range, _get_scalars_for_plotting


pyplot.rcParams.update({"text.usetex": True})


def _plot_line_quantities_from_dict(
    plotting_quantities, x_axis_value, subplot_grid, share_axis, xlabel
):
    """
    Plot the quantities from `plotting_quantities`.

    Return both the `matplotlib.pyplot.figure` object and the `array([axes])`.
    """

    f, axs = pyplot.subplots(*subplot_grid, sharex=share_axis)

    # Ensure this is an array so we can iterate over it
    if not isinstance(axs, ndarray):
        axs = array([axs])

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
                quantity,
                xy=(1, 1),
                xycoords="axes fraction",
                xytext=(-0.6, -0.6),
                textcoords="offset fontsize",
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(facecolor="white", edgecolor="black"),
            )
            axs[i].set_xlabel(xlabel)

    return f, axs


def plot_poloidal_slice(
    state: State,
    nrho: int = 51,
    ntheta: int = 51,
    nzeta: int | float | ndarray = 1,
    quantities: str = "mod_B",
    subplot_grid: list = [1, 1],
    share_axis: bool = False,  # TODO: Fix this so that only the minor radius is shared along the x-axis
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

    rho = _get_coord_range("rho", state.nfp, nrho)
    theta = _get_coord_range("theta", state.nfp, ntheta)
    zeta = _get_coord_range("zeta", state.nfp, nzeta)

    ev = state.evaluate(*[quantities, "X1", "X2"], rho=rho, theta=theta, zeta=zeta)

    # Currently we do not automatically adapt the grid
    # TODO: Try and do this automatically??
    if nzeta != 1:
        if subplot_grid == [1, 1]:
            raise ValueError(
                "Since more than one plane is requested please specify the 'subplot_grid' parameter for the subplot grid sizing"
            )

    f, axs = pyplot.subplots(*subplot_grid, sharex=share_axis, sharey=share_axis)

    if not isinstance(axs, ndarray):
        axs = array(axs)  # Makes it easier/possible to loop over

    # All plots will share the same colour scale
    colour_scale = (
        min(ev[quantities].data[:, :, :]),
        max(ev[quantities].data[:, :, :]),
    )

    # Loop and plot all axes
    for i, ax in enumerate(axs.reshape(-1)):
        f_ax = ax.tricontourf(
            ev.X1.data[:, :, i].flatten(),
            ev.X2.data[:, :, i].flatten(),
            ev[quantities].data[:, :, i].flatten(),
            vmin=colour_scale[0],
            vmax=colour_scale[1],
        )
        if share_axis:
            ax.label_outer()  # Removes any axis labels on subplots on the interior of the grid

    if share_axis:
        # f.subplots_adjust(wspace=0.05,hspace=0.05)
        f.tight_layout()

    # Adding colourbar
    colourbar_axis = f.add_axes([0.85, 0.15, 0.05, 0.7])
    # Colourbar starts at 0.85
    f.subplots_adjust(right=0.83)

    f.colorbar(f_ax, cax=colourbar_axis)

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
    share_axis: bool = True,
    post_process: dict | dict = None,
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

    f, axs = _plot_line_quantities_from_dict(
        plotting_quantities, rho, subplot_grid, share_axis, "$\\rho$"
    )

    return f


def plot_on_axis(
    state,
    nzeta: int | ndarray = 51,
    quantities: str | list = "mod_B",
    subplot_grid: list = [1, 1],
    share_axis: bool = True,
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

    rho = 0.0
    theta = 0.0
    # theta = _get_coord_range("theta", state.nfp, 0.0)
    zeta = _get_coord_range("zeta", state.nfp, nzeta)

    ev = state.evaluate(*quantities, rho=rho, theta=theta, zeta=zeta)

    plotting_quantities = _get_scalars_for_plotting(ev, quantities, post_process, "zeta")

    f, axs = _plot_line_quantities_from_dict(
        plotting_quantities, zeta, subplot_grid, share_axis, "$\\zeta$"
    )

    return f
