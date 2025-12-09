from warnings import warn

import matplotlib.pyplot as pyplot
from numpy import amax, amin, any, isnan, linspace, max, min, ndarray, pi, prod

from gvec.core.state import State

from gvec.plotting.utils import _design_subgrid, _get_coord_range, _subplots

pyplot.rcParams.update({"text.usetex": True})


def plot_poloidal_plane(
    state: State,
    nrho: int = 21,
    ntheta: int = 51,
    zeta: int | float | ndarray = 9,
    quantities: str = "mod_B",
    subplot_grid: list = [3, 3],
    share_axis: bool = False,
    st_contours: list = [0, 0],
    theta_contour_style: str = "theta",
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
    zeta: int, float, ndarray, optional
        Poloidal slices to plot. Default is 9.
        - int: Number of equally spaced slices.
        - float: The specific slice.
        - ndarray: The specific slices.
    quantities: str
        The quantity to plot. Default is "mod_B"
    subplot_grid: list, conditionally optional
        The grid shape for the subplots. Default is `[3,3]`. Required if `zeta>1`.
    share_axis: bool
        If true, all subplots will share their `X1` and `X2` axis positions.
    st_contours: list
        Number of ``$\theta$`` and ``$\zeta$`` contours to plot.
    theta_contour_style: str
        Either `"theta"` or `"theta_star"`.

    Returns
    -------
    `matplotlib.pyplot.figure` object
    """

    if theta_contour_style not in ["theta", "theta_star"]:
        raise ValueError("theta_contour_type must be 'theta' or 'theta_star'")

    rho = _get_coord_range("rho", state.nfp, nrho)
    theta = _get_coord_range("theta", state.nfp, ntheta)
    zeta_eval = _get_coord_range("zeta", state.nfp, zeta)

    ev = state.evaluate(*[quantities, "X1", "X2", "LA"], rho=rho, theta=theta, zeta=zeta_eval)

    if any(isnan(ev[quantities])):
        # Sometimes on-axis quantities cannot be computed properly, to avoid erroring out we will just set them to zero
        #   and warn the user that they should adjust if their plots look weird
        ev[quantities].data[isnan(ev[quantities].data)] = 0.0
        warn(
            "NaNs detected in evaluated dataset, these have been set to 0. It is possible that on-axis quantities cannot be evaluated, a minimum rho value can be set with the min_rho input."
        )

    if prod(subplot_grid) != len(zeta_eval):
        warn(
            f"subplot_grid is {subplot_grid} but number of slices is {len(zeta_eval)}, updating subplot grid to fit."
        )
        subplot_grid = _design_subgrid(len(zeta_eval))
        warn(f"new grid is {subplot_grid}.")

    f, axs = _subplots(*subplot_grid, sharex=share_axis, sharey=share_axis)

    is_scalar = len(ev[quantities].shape) == 1  # TODO: Why did I do this?
    # All plots will share the same colour scale
    if is_scalar:
        colour_scale = (
            min(ev[quantities]),
            max(ev[quantities]),
        )
    else:
        colour_scale = (
            amin(ev[quantities].data),
            amax(ev[quantities].data),
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
        if st_contours[0] > 0:
            # We should plot the rho
            rho_contour = ev.X1[:, :, i] * 0.0 + ev.rho  # TODO: expand dims?
            rho_levels_vis = linspace(0, 1 - 1e-10, st_contours[0])
            f_ax = ax.contour(
                ev.X1.data[:, :, i],
                ev.X2.data[:, :, i],
                rho_contour,
                rho_levels_vis,
                colors="black",
            )

        if st_contours[1] > 0:
            # and/or \theta or \theta^\star contours
            if theta_contour_style == "theta":
                prefactor = 0.0
            else:
                prefactor = 1.0
            theta_contour = prefactor * ev.LA[:, :, i] + ev.theta
            theta_levels_vis = linspace(0, 2 * pi, st_contours[1], endpoint=False)
            ax.contour(
                ev.X1.data[:, :, i],
                ev.X2.data[:, :, i],
                theta_contour,
                theta_levels_vis,
                colors="red",
            )

        ax.annotate(
            f"$\\zeta={zeta_eval[i]:.2f}$",
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

    return f, axs


def plot_on_flux_surface(
    state: State,
    nzeta: int,
    ntheta: int,
    rho: float | ndarray = 1.0,
    quantities: str
    | list = "mod_B",  # TODO: Allow for multiple quantities instead of multiple surfaces
    subplot_grid: list = [1, 1],
    share_axis: bool = True,
    levels: int | ndarray = 10,
    boozer: bool = True,
    filled_contours=False,
):
    """
    Plot an equilibrium quantity over the two angles $(\\vartheta, \\zeta)$ of a flux surface at a given `rho` value

    TODO: Add boozer option

    Parameters
    ----------
    state: GVEC state file
    nzeta: int
    ntheta: int
    rho: float, numpy.ndarray, optional
        The flux surface label(s) to plot. Default is 1.0.
    quantities: str, optional
        TODO: Allow for multiple quantities instead of multiple surfaces -- DONE: NEEDS TEST
    subplot_grid: list, conditionally optional
        The grid shape for the subplots. Default is `[1,1]`. Required if `len(rho)>1`.
    share_axis: bool, optional
        If true, all subplots will share their `x` and `y` axes.
        Default `True`.
    levels: int, numpy.ndarray, optional
        If `int` then chooses number of levels in the contour plot.
        If an `numpy.ndarray` then plots contours at given values.
        Default is `10`
    boozer: bool
        Plot boozer surfaces or regular ``$\theta-\zeta$`` coordinates.
        Default is `True`
    filled_contours: bool
        Use `contour` (False) or `contourf` (True) in plotting.
        Default is `False` (filled contour)


    # TODO: maybe add kwargs

    Returns
    -------
    `matplotlib.pyplot.figure` object
    """

    zeta = _get_coord_range("zeta", state.nfp, nzeta)
    rho_eval = _get_coord_range("rho", state.nfp, rho)
    theta = _get_coord_range("theta", state.nfp, ntheta)

    if isinstance(quantities, list) and (len(rho_eval) > 1):
        raise ValueError(
            "You can either plot multiple quantities on a single surface or a single quantity on multiple surfaces. Either quantities is a list or the number of rho positions is > 1."
        )
    if isinstance(quantities, str):
        quantities_eval = [quantities]
        nplots = len(rho_eval)
    elif isinstance(quantities, list):
        quantities_eval = quantities
        nplots = len(quantities)

    if prod(subplot_grid) != nplots:
        warn(
            f"subplot_grid is {subplot_grid} but number of slices is {nplots}, updating subplot grid to fit."
        )
        subplot_grid = _design_subgrid(nplots)
        warn(f"new grid is {subplot_grid}.")

    # Plotting boozer needs different evaluate function
    if boozer:
        ev = state.evaluate_sfl(
            *quantities_eval, rho=rho_eval, theta=theta, zeta=zeta, sfl="boozer"
        )
    else:
        ev = state.evaluate(*quantities_eval, rho=rho_eval, theta=theta, zeta=zeta)

    f, axs = _subplots(*subplot_grid, sharex=share_axis, sharey=share_axis)

    for i, ax in enumerate(axs):
        # We should not change the slice if we are plotting multiple quantities
        # and we should select the quantity to plot
        if isinstance(quantities, list):
            ev_i = ev.isel(rad=0)
            quantity = quantities[i]
        else:
            ev_i = ev.isel(rad=i)
            quantity = quantities

        # Method switching
        if filled_contours:
            contour_method = ax.contourf
        else:
            contour_method = ax.contour

        contour_method(
            ev_i.theta.data,
            ev_i.zeta.data,
            ev_i[quantity].transpose("pol", "tor").data,
            levels=levels,
        )
        if isinstance(quantities, list):
            ax.annotate(
                f"${ev_i[quantity].attrs['symbol']}$",
                xy=(1, 1),
                xycoords="axes fraction",
                xytext=(-0.6, -0.6),
                textcoords="offset fontsize",
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(facecolor="white", edgecolor="black"),
            )

        if share_axis:
            # Removes any axis labels on subplots on the interior of the grid
            ax.label_outer()

    return f, axs
