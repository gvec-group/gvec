# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from typing import Literal
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np

from gvec.core.state import State
from gvec.plotting.utils import _deco_usetex, _design_subgrid, _subplots, _symbol_check


@_deco_usetex
def plot_poloidal_plane(
    state: State,
    quantity: None | str = "mod_B",
    nrho: int = 21,
    ntheta: int = 51,
    zeta: int | float | np.ndarray | list[float] = 9,
    subplot_grid: list[int] | None = None,
    share_axis: bool = False,
    rho_contours: int = 4,
    rho_contours_color: str | None = None,
    theta_contours: int = 8,
    theta_contours_color: str | None = None,
    sfl: Literal["pest"] | None = "pest",
    plot_kwargs: dict = dict(),
):
    """
    Plot a poloidal plane with some equilibrium quantity on it. Defaults to plotting $|B|$

    Parameters
    ----------
    state : GVEC state object
    quantity : str, optional
        The quantity to plot. Default is ``"mod_B"``. If ``None``, no contours are plotted.
    nrho : int, optional
        The radial resolution of the slices.
        Default is ``51``
    ntheta : int, optional
        The poloidal resolution of the slices.
        Default is ``51``
    zeta : int, float, ndarray, optional
        The number of equally spaced slices (``int``), the specific ``zeta`` value (``float``) or values (``np.ndarray``).
        Default is ``9``.
    subplot_grid : list[int], None, optional
        The grid shape for the subplots. If ``None``, grid will be automatically determined.
        Default is ``None``.
    share_axis : bool
        If ``True``, all subplots will share their ``X1`` and ``X2`` axis positions.
        Default ``False``
    rho_contours : int, optional
        The number of ``rho`` contours to plot. If ``0``, no contours are plotted.
        Default ``4``
    rho_contours_color : str, optional
        The color of the ``rho`` contours.
        Default ``"white"``
    theta_contours : int, optional
        The number of ``theta`` contours to plot. If ``0``, no contours are plotted.
        Default ``8``
    theta_contours_color : str, optional
        The color of the ``theta`` contours.
        Default ``"white"``
    sfl : ``"pest"`` or ``None``, optional
        Plot the ``theta`` contours or ``pest`` contours ($\\theta^\\star$).
        Default ``"pest"``.
    plot_kwargs: dict, optional
        Any ``**kwargs`` to send to the ``plt.figure()`` function.
        For example ``plot_kwargs={'figsize': (8,8)}``. See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html>`_ for a list of kwargs.

    Returns
    -------
    ``matplotlib.pyplot.figure`` object and ``numpy.ndarray`` of ``matplotlib.axis._axis.Axes`` object(s).
    """

    if rho_contours_color is None:
        rho_contours_color = "black" if quantity is None else "white"
    if theta_contours_color is None:
        theta_contours_color = "black" if quantity is None else "white"

    quantities = ["X1", "X2", "LA"]
    if quantity is not None:
        quantities.append(quantity)

    ev_contour = state.evaluate(
        *quantities, rho=nrho, theta=np.linspace(0, 2 * np.pi, ntheta), zeta=zeta
    )

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

    zeta_eval = ev_contour.zeta.data
    if quantity is not None:
        if np.any(np.isnan(ev_contour[quantity])):
            # Sometimes on-axis quantity cannot be computed properly, to avoid erroring out we will just set them to zero
            #   and warn the user that they should adjust if their plots look weird
            ev_contour[quantity].data[np.isnan(ev_contour[quantity].data)] = 0.0
            warn(
                "NaNs detected in evaluated dataset, these have been set to 0. It is possible that on-axis quantity cannot be evaluated, a minimum rho value can be set with the min_rho input."
            )

    if not subplot_grid:
        subplot_grid = _design_subgrid(len(zeta_eval))

    fig, axs = _subplots(subplot_grid, share_axis, share_axis, **plot_kwargs)

    if quantity is not None:
        # All plots will share the same colour scale
        plotting_quantity = ev_contour[quantity].broadcast_like(ev_contour.X1)
        colour_scale = (plotting_quantity.min().item(), plotting_quantity.max().item())

    # if plotting_quantity.ndim
    # print(plotting_quantity.ndim)
    if "xyz" in plotting_quantity.coords:
        raise ValueError(
            f"Plotting quantity must be a scalar field but {quantity} is a vector field."
        )

    # Loop and plot all axes
    for i, ax in enumerate(np.asarray(axs).flat):
        if quantity is not None:
            # Check if the quantity is a scalar
            # if not is_scalar:  # not a scalar
            # plotting_quantity = ev_contour[quantity].data[:, :, i]  # .flatten()
            # else:
            # plotting_quantity = ev_contour.X1[:, :, i] * 0.0 + ev_contour[quantity]
            # plotting_quantity = plotting_quantity  # .flatten()

            f_ax = ax.contourf(
                ev_contour.X1.isel(tor=i),  # .flatten(),
                ev_contour.X2.isel(tor=i),  # .flatten(),
                plotting_quantity.isel(tor=i),
                vmin=colour_scale[0],
                vmax=colour_scale[1],
            )
        if share_axis:
            ax.label_outer()  # Removes any axis labels on subplots on the interior of the grid
        if rho_contours:
            # We should plot the rho
            ax.plot(
                ev_rho_contours.X1[:, :, i].T,
                ev_rho_contours.X2[:, :, i].T,
                color=rho_contours_color,
                linewidth=1.0,
            )

        if theta_contours:
            # \theta or \theta^\star contours
            ax.plot(
                ev_theta_contours.X1[:, :, i],
                ev_theta_contours.X2[:, :, i],
                color=theta_contours_color,
                linewidth=1.0,
            )

        # The slice label will be added as an annotation to the top right of the subplot
        zeta_angle = zeta_eval[i] / np.pi
        # ax.annotate(
        #     f"$\\zeta={zeta_angle:.2f} \\pi$",
        #     xy=(1, 1),
        #     xycoords="axes fraction",
        #     xytext=(-0.6, -0.6),
        #     textcoords="offset fontsize",
        #     verticalalignment="top",
        #     horizontalalignment="right",
        #     bbox=dict(facecolor="white", edgecolor="black"),
        # )
        ax.set_title(f"$\\zeta={zeta_angle:.2f} \\pi$")

        # Remove interior axis labels and add axis labels to the boundary plots
        if share_axis:
            ax.set(xlabel="X1", ylabel="X2")
            ax.label_outer()
            ax.set_aspect("equal")
        else:
            # Bufer the axis limits and make sure all things are square
            xlims = ax.get_xlim()
            ylims = ax.get_ylim()

            margin = 1.05

            if np.diff(ylims) > np.diff(xlims):
                ax.set_xlim((np.sum(xlims) + (-1, 1) * np.diff(ylims) * margin) / 2)
                ax.set_ylim((np.sum(ylims) + (-1, 1) * np.diff(ylims) * margin) / 2)
            else:
                ax.set_xlim((np.sum(xlims) + (-1, 1) * np.diff(xlims) * margin) / 2)
                ax.set_ylim((np.sum(ylims) + (-1, 1) * np.diff(xlims) * margin) / 2)

            ax.set_box_aspect(1)

    if quantity is not None:
        # Adding colourbar
        ev_contour = _symbol_check(ev_contour, [quantity])
        fig.colorbar(
            f_ax, ax=np.asarray(axs).ravel().tolist(), label=f"${ev_contour[quantity].symbol}$"
        )

    return fig, axs


@_deco_usetex
def plot_on_flux_surface(
    state: State,
    quantities: str | list[str] = "mod_B",
    rho: float | np.ndarray | list = 1.0,
    ntheta: int = 51,
    nzeta: int = 51,
    subplot_grid: list[int] | None = None,
    share_axis: bool = True,
    share_contours: bool = False,
    levels: int | np.ndarray | list = 10,
    sfl: Literal["pest", "boozer"] | None = "boozer",
    style: Literal["contour", "filled-contour"] = "contour",
    plot_kwargs: dict = dict(),
    **boozer_kwargs,
):
    """
    Plot an equilibrium quantity over the two angles $(\\vartheta, \\zeta)$ of a flux surface at (a) given ``rho`` value(s). Alternatively, plot
    multiple ``quantities`` on a single ``rho`` surface.

    Parameters
    ----------
    state : GVEC state object
    quantities : str, list[str], optional
        Plot either a single quantitiy on a number of `rho` surfaces or a number of `quantities` on a single `rho` surface.
        Default ``mod_B``
    rho : float, numpy.ndarray, list, optional
        The flux surface label(s) to plot.
        Default is ``1.0``.
    ntheta : int
        Resolution in ``theta``.
        Default is ``11``
    nzeta : int
        Resolution in ``zeta``.
        Default is `11`
    subplot_grid : list[int], optional
        The grid shape for the subplots. If ``None``, grid will be automatically determined.
        Default is ``None``.
    share_axis : bool, optional
        If ``True``, all subplots will share their ``x`` and ``y`` axes.
        Default ``True``.
    levels : int, numpy.ndarray, optional
        If ``int`` then chooses number of levels in the contour plot. If an ``numpy.ndarray`` or ``list`` then plots contours at given values.
        Default is ``10``
    sfl : str, optional
        Plot surfaces in ``"boozer"``, ``"pest"`` or regular $(\\theta,\\zeta)$ (``"None"``) coordinates.
        Default is ``"boozer"``
    style : str, optional
        Use ``"contour"`` (``False``) or ``"filled-contour"`` (``True``) in plotting.
        Default is ``"filled-contour"`` (filled contour).
    plot_kwargs: dict, optional
        Any ``**kwargs`` to send to the ``plt.figure()`` function.
        For example ``plot_kwargs={'figsize': (8,8)}``. See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html>`_ for a list of kwargs.
    boozer_kwargs : optional
        Keyword arguments for the case where ``boozer`` is used.


    Returns
    -------
    ``matplotlib.pyplot.figure`` object and ``numpy.ndarray`` of ``matplotlib.axis._axis.Axes`` object(s).
    """

    if isinstance(rho, int):
        rho_len = rho
    elif isinstance(rho, float):
        rho_len = 1
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
    theta = np.linspace(0.0, 2 * np.pi, ntheta)
    zeta = np.linspace(0.0, 2 * np.pi / state.nfp, nzeta)
    if sfl:
        evaluations = state.evaluate_sfl(
            *quantities_eval, rho=rho, theta=theta, zeta=zeta, sfl=sfl, **boozer_kwargs
        )
    else:
        evaluations = state.evaluate(*quantities_eval, rho=rho, theta=theta, zeta=zeta)

    fig, axs = _subplots(subplot_grid, share_axis, share_axis, **plot_kwargs)

    evaluations = _symbol_check(evaluations, quantities_eval)

    # The actual plotting bit
    for i, ax in enumerate(np.asarray(axs).flat):
        # We should not change the slice if we are plotting multiple quantities
        # and we should select the quantity to plot
        if isinstance(quantities, list):
            evaluations_i = evaluations.isel(rad=0)
            quantity = quantities_eval[i]
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
            zeta_vals.data / (2 * np.pi),
            theta_vals.data / (2 * np.pi),
            evaluations_i[quantity].transpose("pol", "tor").data,
            levels=levels,
        )

        if isinstance(quantities, list):
            # Specify the quantity type in the top right corner if multiple quantities were requested
            plot_label = f"${evaluations_i[quantity].attrs['symbol']}$"
            fig.colorbar(f_ax, ax=ax)  # , label=f"${evaluations[quantity].symbol}$")
        else:
            plot_label = f"$\\rho={evaluations.rho.data[i]}$"
            if not share_contours:
                fig.colorbar(f_ax, ax=ax)

        # ax.annotate(
        #     plot_label,
        #     xy=(1, 1),
        #     xycoords="axes fraction",
        #     xytext=(-0.6, -0.6),
        #     textcoords="offset fontsize",
        #     verticalalignment="top",
        #     horizontalalignment="right",
        #     bbox=dict(facecolor="white", edgecolor="black"),
        # )
        ax.set_title(plot_label)

        # zeta_angle = zeta_eval[i] / np.pi

        ax.set(
            xlabel=rf"${zeta_vals.attrs['symbol']} / (2\pi)$",
            ylabel=rf"${theta_vals.attrs['symbol']} / (2\pi)$",
        )

        # Removes any axis labels on subplots on the interior of the grid
        ax.label_outer()

    if isinstance(quantities, str) and share_contours:
        # Adding colourbar if single quantity was requested
        fig.colorbar(f_ax, ax=axs.ravel().tolist(), label=f"${evaluations[quantity].symbol}$")

    return fig, axs


@_deco_usetex
def plot_fourier_on_surface(
    state: State,
    quantity: str = "mod_B",
    rho: float = 1.0,
    ntheta: int = 101,
    nzeta: int = 101,
    sfl: Literal["pest", "boozer"] | None = None,
    limit: float | None = 1e-15,
):
    """
    Diagnostic plot for plotting the Fourier modes of a given quantitity on a flux surface.

    Parameters
    ----------
    state : GVEC state object
    quantities : str, optional
        Quantitiy to plot.
        Default ``mod_B``
    rho : float, optional
        The flux surface label(s) to plot.
        Default is ``1.0``.
    ntheta : int, optional
        Resolution in ``theta``.
        Default is ``101``
    nzeta : int, optional
        Resolution in ``zeta``.
        Default is ``101``
    sfl : str, optional
        Plot surfaces in ``"boozer"``, ``"pest"`` or regular $\\theta-\\zeta$ (``None``) coordinates.
        Default is ``"boozer"``
    limit: float, optional
        Cut-off value for the Fourier amplitudes to plot.
        Default is ``1e-15``


    Returns
    -------
    ``matplotlib.pyplot.figure`` object and ``numpy.ndarray`` of ``matplotlib.axis._axis.Axes`` object(s).
    """
    from gvec.core.compute import ev2ft

    quantities = [quantity, "N_FP"]
    if sfl is None:
        evaluations = state.evaluate(*quantities, rho=rho, theta=ntheta, zeta=nzeta)
    else:
        evaluations = state.evaluate_sfl(
            *quantities, rho=rho, theta=ntheta, zeta=nzeta, sfl=sfl
        )
    evaluations = evaluations[quantities]

    evaluations = _symbol_check(evaluations, quantities)

    evft = ev2ft(evaluations).squeeze().sortby("n")

    levels = np.linspace(-14, 0, 8)
    symbol = evaluations[quantity].attrs.get("symbol", f"\\mathrm{{{quantity}}}")

    fig, axs = _subplots([1, 2], True, True)
    for ax, suffix in zip(axs, ["mnc", "mns"]):
        c = ax.contourf(
            evft.n,
            evft.m,
            np.log10(np.maximum(np.abs(evft[f"{quantity}_{suffix}"]), 1e-16)),
            levels=levels,
            extend="both",
        )
    fig.colorbar(c, ax=axs, label=rf"$\log_{{10}}|{symbol}|$")
    coords = {
        None: "logical coordinates",
        "pest": "PEST coordinates",
        "boozer": "Boozer coordinates",
    }
    fig.suptitle(
        f"Fourier spectrum of ${symbol}$ on flux surface $\\rho={rho}$ in {coords[sfl]}"
    )
    for ax in axs:
        ax.set(xlabel="$n$")
    axs[0].set(
        ylabel="$m$",
        title="cosine components",
    )
    axs[1].set(
        title="sine components",
    )

    if limit is not None:
        power = np.sqrt(evft[f"{quantity}_mnc"] ** 2 + evft[f"{quantity}_mns"] ** 2)
        limit_m = power.m.where((power > limit).sum(dim="n") == 0).min().item()
        limit_m = np.nanmax([limit_m, 5])
        limit_n1 = (
            power.n.where((power > limit).sum(dim="m") == 0).where(power.n > 0).min().item()
        )
        limit_n2 = (
            power.n.where((power > limit).sum(dim="m") == 0).where(power.n < 0).max().item()
        )
        limit_n = np.nanmax([abs(limit_n1), abs(limit_n2), 5])
        for ax in axs:
            ax.set(
                xlim=(-limit_n - 1, limit_n + 1),
                ylim=(0, limit_m + 1),
            )

    return fig, axs
