# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from typing import Literal
from warnings import warn

import numpy as np

from gvec.core.state import State
from gvec.plotting.utils import _design_subgrid, _subplots, _symbol_check


def plot_poloidal_plane(
    state: State,
    quantity: str = "mod_B",
    /,
    nrho: int = 21,
    ntheta: int = 51,
    zeta: int | float | np.ndarray = 9,
    subplot_grid: list[int] | None = None,
    share_axis: bool = False,
    rho_contours: int = 4,
    theta_contours: int = 8,
    sfl: Literal["pest"] | None = "pest",
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
        The radial resolution of the slices.
        Default is `51`
    ntheta : int, optional
        The poloidal resolution of the slices.
        Default is `51`
    zeta : int, float, ndarray, optional
        The number of equally spaced slices (int), the specific `zeta` value (float) or values (np.ndarray).
        Default is `9`.
    subplot_grid : list[int], None, optional
        The grid shape for the subplots. If `None`, grid will be automatically determined.
        Default is `None`.
    share_axis : bool
        If `True`, all subplots will share their `X1` and `X2` axis positions.
        Default `False`
    rho_contours : int, optional
        The number of ``$\rho$`` contours to plot.
        Default `4`
    theta_contours : int, optional
        The number of ``$\theta$`` contours to plot.
        Default `8`
    sfl : `"pest"` or `None`, optional
        Plot the `theta` contours or `pest` contours (``$\theta^\star$``).
        Default `"pest"`.

    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
    """

    evaluations = state.evaluate(quantity, "X1", "X2", "LA", rho=nrho, theta=ntheta, zeta=zeta)

    if rho_contours:
        ev_rho_contours = state.evaluate(
            "X1", "X2", rho=rho_contours, theta=np.linspace(0, 2 * np.pi, ntheta), zeta=zeta
        )
    if theta_contours:
        if sfl == "pest":
            ev_theta_contours = state.evaluate_sfl(
                "X1",
                "X2",
                rho=np.linspace(0, 1, nrho),
                theta=theta_contours,
                zeta=zeta,
                sfl="pest",
            )
        else:
            ev_theta_contours = state.evaluate(
                "X1", "X2", rho=np.linspace(0, 1, nrho), theta=theta_contours, zeta=zeta
            )

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
        if rho_contours:
            # We should plot the rho
            ax.plot(
                ev_rho_contours.X1[:, :, i].T, ev_rho_contours.X2[:, :, i].T, "w", linewidth=1.0
            )

        if theta_contours:
            # \theta or \theta^\star contours
            ax.plot(
                ev_theta_contours.X1[:, :, i], ev_theta_contours.X2[:, :, i], "w", linewidth=1.0
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

    # Adding colourbar
    evaluations = _symbol_check(evaluations, [quantity])
    f.colorbar(f_ax, ax=axs.ravel().tolist(), label=f"${evaluations[quantity].symbol}$")

    return f, axs


def plot_on_flux_surface(
    state: State,
    quantities: str | list[str] = "mod_B",
    rho: float | np.ndarray | list = 1.0,
    ntheta: int = 51,
    nzeta: int = 51,
    subplot_grid: list[int] | None = None,
    share_axis: bool = True,
    levels: int | np.ndarray | list = 10,
    sfl: Literal["pest", "boozer"] | None = "boozer",
    style: Literal["contour", "filled-contour"] = "contour",
    **boozer_kwargs,
):
    """
    Plot an equilibrium quantity over the two angles $(\\vartheta, \\zeta)$ of a flux surface at (a) given `rho` value(s). Alternatively, plot
    multiple `quantities` on a single `rho` surface.

    Parameters
    ----------
    state : GVEC state file
    quantities : str, list[str], optional
        Plot either a single quantitiy on a number of `rho` surfaces or a number of `quantities` on a single `rho` surface.
        Default `mod_B`
    rho : float, numpy.ndarray, list, optional
        The flux surface label(s) to plot.
        Default is `1.0`.
    ntheta : int
        Resolution in `theta`.
        Default is `11`
    nzeta : int
        Resolution in `zeta`.
        Default is `11`
    subplot_grid : list[int], optional
        The grid shape for the subplots. If `None`, grid will be automatically determined.
        Default is `None`.
    share_axis : bool, optional
        If `True`, all subplots will share their `x` and `y` axes.
        Default `True`.
    levels : int, numpy.ndarray, optional
        If `int` then chooses number of levels in the contour plot. If an `numpy.ndarray` or `list` then plots contours at given values.
        Default is `10`
    sfl : str, optional
        Plot surfaces in `"boozer"`, `"pest"` or regular ``$\theta-\zeta$`` coordinates.
        Default is `"boozer"`
    style : str, optional
        Use `"contour"` (False) or `"filled-contour"` (True) in plotting.
        Default is `"filled-contour"` (filled contour).
    boozer_kwargs : optional
        Keyword arguments for the case where `boozer` is used.


    Returns
    -------
    `matplotlib.pyplot.figure` object and `numpy.ndarray` of `matplotlib.axis._axis.Axes` object(s).
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
        share_axis = False
    elif isinstance(quantities, list):
        quantities_eval = quantities
        nplots = len(quantities)

    if subplot_grid is None:
        subplot_grid = _design_subgrid(nplots)

    # Plotting boozer needs different evaluate function
    if sfl:
        evaluations = state.evaluate_sfl(
            *quantities_eval, rho=rho, theta=ntheta, zeta=nzeta, sfl=sfl, **boozer_kwargs
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
        if style == "filled-contour":
            contour_method = ax.contourf
        elif style == "contour":
            contour_method = ax.contour
        else:
            raise ValueError("style must be 'filled-contour' or 'contour'.")

        if sfl == "boozer":
            theta_vals = evaluations_i.theta_B
            zeta_vals = evaluations_i.zeta_B
        elif sfl == "pest":
            theta_vals = evaluations_i.theta_P
            zeta_vals = evaluations_i.zeta
        else:
            theta_vals = evaluations_i.theta
            zeta_vals = evaluations_i.zeta

        f_ax = contour_method(
            zeta_vals.data,
            theta_vals.data,
            evaluations_i[quantity].transpose("pol", "tor").data,
            levels=levels,
        )
        if isinstance(quantities, list):
            plot_label = f"${evaluations_i[quantity].attrs['symbol']}$"
        else:
            plot_label = f"$\\rho={evaluations.rho.data[i]}$"
            # Specify the quantity type in the top right corner if multiple quantities were requested
        ax.annotate(
            plot_label,
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

        ax.set(
            xlabel=f"${zeta_vals.attrs['symbol']}$", ylabel=f"${theta_vals.attrs['symbol']}$"
        )
        if share_axis:
            # Removes any axis labels on subplots on the interior of the grid
            ax.label_outer()

    if isinstance(quantities, str):
        # Adding colourbar if single quantity was requested
        f.colorbar(f_ax, ax=axs.ravel().tolist(), label=f"${evaluations[quantity].symbol}$")

    return f, axs
