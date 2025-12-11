# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from typing import Literal
from warnings import warn

import numpy as np
from numpy import ndarray  # type checking

from gvec.core.state import State
from gvec.plotting.utils import _design_subgrid, _subplots, _symbol_check


def plot_poloidal_plane(
    state: State,
    quantity: str = "mod_B",
    /,
    nrho: int = 21,
    ntheta: int = 51,
    zeta: int | float | ndarray = 9,
    subplot_grid: list | None = None,
    share_axis: bool = False,
    st_contours: list = [0, 0],
    theta_contour_style: str = "theta",
):
    """
    plot_poloidal_plane

    Plot a poloidal plane with some equilibrium quantity on it. Defaults to plotting $\|\mathbf{B}\|$

    Parameters
    ----------
    state : GVEC state file
    quantity : str, optional
        The quantity to plot. Default is "mod_B"
    nrho : int, optional
        The radial resolution of the slices. Default is 51
    ntheta : int, optional
        The poloidal resolution of the slices. Default is 51
    zeta : int, float, ndarray, optional
        Poloidal slices to plot. Default is 9.
        - int: Number of equally spaced slices.
        - float: The specific slice.
        - ndarray: The specific slices.
    subplot_grid : list, None, optional
        The grid shape for the subplots. If `None`, grid will be automatically determined. Default is `None`.
    share_axis : bool
        If true, all subplots will share their `X1` and `X2` axis positions.
    st_contours : list
        Number of ``$\theta$`` and ``$\zeta$`` contours to plot.
    theta_contour_style : str
        Either `"theta"` or `"theta_star"`.

    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
    """

    if theta_contour_style not in ["theta", "theta_star"]:
        raise ValueError("theta_contour_type must be 'theta' or 'theta_star'")

    evaluations = state.evaluate(quantity, "X1", "X2", "LA", rho=nrho, theta=ntheta, zeta=zeta)

    zeta_eval = evaluations.zeta.data

    if np.any(np.isnan(evaluations[quantity])):
        # Sometimes on-axis quantity cannot be computed properly, to avoid erroring out we will just set them to zero
        #   and warn the user that they should adjust if their plots look weird
        evaluations[quantity].data[np.isnan(evaluations[quantity].data)] = 0.0
        warn(
            "NaNs detected in evaluated dataset, these have been set to 0. It is possible that on-axis quantity cannot be evaluated, a minimum rho value can be set with the min_rho input."
        )

    if not subplot_grid:
        subplot_grid = _design_subgrid(len(zeta_eval))

    f, axs = _subplots(subplot_grid, sharex=share_axis, sharey=share_axis)

    is_scalar = len(evaluations[quantity].shape) == 1
    # All plots will share the same colour scale
    # Need to make sure we use the correct function for getting the min and max
    if is_scalar:
        colour_scale = (
            np.min(evaluations[quantity]),
            np.max(evaluations[quantity]),
        )
    else:
        colour_scale = (
            np.amin(evaluations[quantity].data),
            np.amax(evaluations[quantity].data),
        )

    # Loop and plot all axes
    for i, ax in enumerate(axs):
        # Check if the quantity is a scalar
        if not is_scalar:  # not a scalar
            plotting_quantity = evaluations[quantity].data[:, :, i].flatten()
        else:
            plotting_quantity = evaluations.X1[:, :, i] * 0.0 + evaluations[quantity]
            plotting_quantity = plotting_quantity.data.flatten()

        f_ax = ax.tricontourf(
            evaluations.X1.data[:, :, i].flatten(),
            evaluations.X2.data[:, :, i].flatten(),
            plotting_quantity,
            vmin=colour_scale[0],
            vmax=colour_scale[1],
        )
        if share_axis:
            ax.label_outer()  # Removes any axis labels on subplots on the interior of the grid
        if st_contours[0] > 0:
            # We should plot the rho
            rho_contour = evaluations.X1[:, :, i] * 0.0 + evaluations.rho  # TODO: expand dims?
            rho_levels_vis = np.linspace(0, 1 - 1e-10, st_contours[0])
            f_ax = ax.contour(
                evaluations.X1.data[:, :, i],
                evaluations.X2.data[:, :, i],
                rho_contour,
                rho_levels_vis,
                colors="black",
            )

        if st_contours[1] > 0:
            # \theta or \theta^\star contours
            if theta_contour_style == "theta":
                prefactor = 0.0
            else:
                prefactor = 1.0
            theta_contour = prefactor * evaluations.LA[:, :, i] + evaluations.theta
            theta_levels_vis = np.linspace(0, 2 * np.pi, st_contours[1], endpoint=False)
            ax.contour(
                evaluations.X1.data[:, :, i],
                evaluations.X2.data[:, :, i],
                theta_contour,
                theta_levels_vis,
                colors="red",
            )

        # The slice label will be added as an annotation to the top right of the subplot
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

        # Remove interior axis labels and add axis labels to the boundary plots
        if share_axis:
            ax.set(xlabel="X1", ylabel="X2")
            ax.label_outer()

    if share_axis:
        f.tight_layout()

    # Adding colourbar
    colourbar_axis = f.add_axes([0.85, 0.15, 0.05, 0.7])
    # Colourbar starts at 0.85
    f.subplots_adjust(right=0.83)

    evaluations = _symbol_check(evaluations, [quantity])

    f.colorbar(f_ax, cax=colourbar_axis, label=f"${evaluations[quantity].symbol}$")

    return f, axs


def plot_on_flux_surface(
    state: State,
    quantities: str | list[str] = "mod_B",
    rho: float | ndarray | list = 1.0,
    ntheta: int = 11,
    nzeta: int = 11,
    subplot_grid: list | None = None,
    share_axis: bool = True,
    levels: int | ndarray = 10,
    sfl: Literal["pest", "boozer"] | None = "boozer",
    filled_contours=False,
):
    """
    Plot an equilibrium quantity over the two angles $(\\vartheta, \\zeta)$ of a flux surface at (a) given `rho` value(s). Alternatively, plot
    multiple `quantities` on a single `rho` surface.

    Parameters
    ----------
    state : GVEC state file
    quantities: str, list[str], optional
        Plot either a single quantitiy on a number of `rho` surfaces or a number of `quantities` on a single `rho` surface.
        Default `mod_B`
    rho : float, numpy.ndarray, optional
        The flux surface label(s) to plot. Default is 1.0.
    ntheta : int
        Default is 11
    nzeta : int
        Default is 11
    subplot_grid : list, optional
        The grid shape for the subplots. If `None`, grid will be automatically determined. Default is `None`.
    share_axis : bool, optional
        If true, all subplots will share their `x` and `y` axes.
        Default `True`.
    levels : int, numpy.ndarray, optional
        If `int` then chooses number of levels in the contour plot.
        If an `numpy.ndarray` then plots contours at given values.
        Default is `10`
    sfl : str
        Plot surfaces in `"boozer"`, `"pest"` or regular ``$\theta-\zeta$`` coordinates.
        Default is `"boozer"`
    filled_contours : bool
        Use `contour` (False) or `contourf` (True) in plotting.
        Default is `False` (filled contour)


    Returns
    -------
    `matplotlib.pyplot.figure` object
    """

    if isinstance(rho, int) or isinstance(rho, float):
        rho_len = rho
    else:
        rho_len = len(rho)

    if isinstance(quantities, list) and (rho_len > 1):
        raise ValueError(
            "You can either plot multiple quantities on a single surface or a single quantity on multiple surfaces. Either quantities is a list or the number of rho positions is > 1."
        )

    if isinstance(quantities, str):
        quantities_eval = [quantities]
        nplots = rho_len
    elif isinstance(quantities, list):
        quantities_eval = quantities
        nplots = len(quantities)

    if subplot_grid is None:
        subplot_grid = _design_subgrid(nplots)

    # Plotting boozer needs different evaluate function
    if sfl:
        evaluations = state.evaluate_sfl(
            *quantities_eval, rho=rho, theta=ntheta, zeta=nzeta, sfl=sfl
        )
    else:
        evaluations = state.evaluate(*quantities_eval, rho=rho, theta=ntheta, zeta=nzeta)

    f, axs = _subplots(subplot_grid, sharex=share_axis, sharey=share_axis)

    evaluations = _symbol_check(evaluations, quantities_eval)

    for i, ax in enumerate(axs):
        # We should not change the slice if we are plotting multiple quantities
        # and we should select the quantity to plot
        if isinstance(quantities, list):
            evaluations_i = evaluations.isel(rad=0)
            quantity = quantities[i]
        else:
            evaluations_i = evaluations.isel(rad=i)
            quantity = quantities

        # Method switching
        if filled_contours:
            contour_method = ax.contourf
        else:
            contour_method = ax.contour

        f_ax = contour_method(
            evaluations_i.theta.data,
            evaluations_i.zeta.data,
            evaluations_i[quantity].transpose("pol", "tor").data,
            levels=levels,
        )
        if isinstance(quantities, list):
            # Specify the quantity type in the top right corner if multiple quantities were requested
            ax.annotate(
                f"${evaluations_i[quantity].attrs['symbol']}$",
                xy=(1, 1),
                xycoords="axes fraction",
                xytext=(-0.6, -0.6),
                textcoords="offset fontsize",
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(facecolor="white", edgecolor="black"),
            )

        if isinstance(quantities, list):
            f.colorbar(f_ax, ax=ax)  # , label=f"${evaluations[quantity].symbol}$")

        if share_axis:
            # Removes any axis labels on subplots on the interior of the grid
            ax.label_outer()

    if isinstance(quantities, str):
        # Adding colourbar if single quantity was requested
        colourbar_axis = f.add_axes([0.85, 0.15, 0.05, 0.7])
        # Colourbar starts at 0.85
        f.subplots_adjust(right=0.83)

        f.colorbar(f_ax, cax=colourbar_axis, label=f"${evaluations[quantity].symbol}$")

    return f, axs
