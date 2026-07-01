# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from typing import Literal
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors

from gvec.core.state import State
from gvec.plotting.utils import _design_subgrid, _subplots, _symbol_check


def plot_poloidal_plane(
    state: State,
    quantity: None | str = "mod_B",
    rho: int | np.ndarray | list[float] = 21,
    ntheta: int = 51,
    zeta: int | float | np.ndarray | list[float] = None,
    subplot_grid: list[int] | None = None,
    share_axis: bool = False,
    rho_contours: int = 4,
    rho_contours_color: str | None = None,
    theta_contours: int = 8,
    theta_contours_color: str | None = None,
    sfl: Literal["pest"] | None = "pest",
    levels: int = 10,
    colorbar_bounds: tuple = (None, None),
    colorbar_scale: Literal["linear", "log"] = "linear",
    plot_kwargs: dict = dict(),
):
    """
    Plot a poloidal plane with some equilibrium quantity on it. Defaults to plotting $|B|$

    Parameters
    ----------
    state : GVEC state object
    quantity : str, optional
        The quantity to plot. Default is ``"mod_B"``. If ``None``, no contours are plotted.
    rho : int, optional
        The radial resolution of the slices. Either by number, equally space, (``int``), or values (``np.ndarray``).
        Default is ``21``
    ntheta : int, optional
        The poloidal resolution of the slices.
        Default is ``51``
    zeta : int, float, ndarray, optional
        The number of equally spaced slices (``int``), the specific ``zeta`` value (``float``) or values (``np.ndarray``).
        Default is ``None``. Note: If None and ``n=0`` then ``zeta=1``, otherwise ``zeta=9``.
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
    levels : int, optional
        Choose number of levels in the contour plot.
        Default is ``10``
    colorbar_bounds: tuple, optional
        The bounds for the color scale of the plotted quantity, as a tuple of (min, max). Default is (None,None), then the bound is determined from the minimum and maximum of the plotted quantity. If only one is set to None, then that bound is determined from the plotted quantity.
    colorbar_scale: str, optional
        Scale of the colorbar, either ``"linear"`` or ``"log"`` for logarithmic scaling.
        Default is ``"linear"``.
    plot_kwargs: dict, optional
        Any ``**kwargs`` to send to the ``plt.figure()`` function.
        For example ``plot_kwargs={'figsize': (8,8)}``. See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html>`_ for a list of kwargs.

    Returns
    -------
    ``matplotlib.pyplot.figure`` object and ``numpy.ndarray`` of ``matplotlib.axis._axis.Axes`` object(s).

    Remark
    ------
    If you want to plot a new quantity that is not in the gvec list of quantities: Before calling the plotting, define a function that adds the new quantity to the dataset ``ds``, and use a decorator to register it. The decorator contains a list of required existing quantities and attributes of the new quantity (like the symbol), and the function defines how to compute it from the required quantities. Once the quantity is registered, you can use it for plotting.

    Here an example:

    >>> @gvec.core.compute.register(
    ...    requirements=("mod_B","mu0"),
    ...    attrs=dict(
    ...        long_name="magnetic pressure",
    ...        symbol=r"|B|^2/(2\mu_0)",
    ...    ),
    >>> )
    >>> def magnetic_pressure(ds):
    ...    ds["magnetic_pressure"] = ds.mod_B**2/(2*ds.mu0)

    You can nest this and register a second quantity that uses the previously registered in its requirements.
    """

    if rho_contours_color is None:
        rho_contours_color = "black" if quantity is None else "white"
    if theta_contours_color is None:
        theta_contours_color = "black" if quantity is None else "white"

    quantities = ["X1", "X2", "LA"]
    if quantity is not None:
        quantities.append(quantity)

    if zeta is None:
        # If zeta is not set, determine it based on the number of toroidal modes
        x1_n_max = state.parameters["X1_mn_max"][1]
        x2_n_max = state.parameters["X2_mn_max"][1]
        n_max = max(x1_n_max, x2_n_max)
        if n_max == 0:
            zeta = 1
        else:
            zeta = 9

    ev_contour = state.evaluate(
        *quantities, rho=rho, theta=np.linspace(0, 2 * np.pi, ntheta), zeta=zeta
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
                rho=rho,
                theta=theta_contours,
                zeta=zeta,
                sfl="pest",
            )
        else:
            ev_theta_contours = state.evaluate(
                "X1", "X2", rho=rho, theta=theta_contours, zeta=zeta
            )

    zeta_eval = ev_contour.zeta.data
    if quantity is not None:
        if np.any(np.isnan(ev_contour[quantity])):
            # Sometimes on-axis quantity cannot be computed properly, to avoid erroring out we will just set them to zero
            #   and warn the user that they should adjust if their plots look weird
            ev_contour[quantity].data[np.isnan(ev_contour[quantity].data)] = 0.0
            warn("NaNs detected in evaluated dataset, these have been set to 0.")

    if not subplot_grid:
        subplot_grid = _design_subgrid(len(zeta_eval))

    fig, axs = _subplots(subplot_grid, share_axis, share_axis, **plot_kwargs)

    if quantity is not None:
        # All plots will share the same color scale
        plotting_quantity = ev_contour[quantity].broadcast_like(ev_contour.X1)
        # Make sure we are not trying to plot a vector field
        if "xyz" in plotting_quantity.coords:
            raise ValueError(
                f"Plotting quantity must be a scalar field but {quantity} is a vector field."
            )

        data_bounds = (plotting_quantity.min().item(), plotting_quantity.max().item())
        color_bounds = data_bounds
        # replace bounds if user-defined `colorbar_bounds` is set.
        if colorbar_bounds[0] is not None:
            color_bounds = (colorbar_bounds[0], color_bounds[1])
        if colorbar_bounds[1] is not None:
            color_bounds = (color_bounds[0], colorbar_bounds[1])
        if colorbar_scale == "log":
            if color_bounds[0] <= 0:
                raise ValueError(
                    f"Log colorbar scale is not possible with non-positive bounds, bounds used: {color_bounds}. Make sure bounds are positive and non-zero."
                )
            color_norm = colors.LogNorm(vmin=color_bounds[0], vmax=color_bounds[1])
        else:
            color_norm = colors.Normalize(vmin=color_bounds[0], vmax=color_bounds[1])
        # explicitly define a discrete colormap, used for both the contour plots and the colorbar.
        discrete_cmap = plt.get_cmap("viridis").resampled(levels)
        cbar_args = {}
        if data_bounds[1] > color_bounds[1]:
            discrete_cmap.set_over("red")
            cbar_args["extend"] = "max"

    # Loop and plot all axes
    for i, ax in enumerate(np.asarray(axs).flat):
        if quantity is not None:
            ax.contourf(
                ev_contour.X1.isel(tor=i),
                ev_contour.X2.isel(tor=i),
                plotting_quantity.isel(tor=i),
                levels=levels,
                norm=color_norm,
                cmap=discrete_cmap,
            )
        if rho_contours:
            # We should plot the rho contours
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

        ax.set_title(f"$\\zeta={zeta_angle:.3f} \\pi$")

        if share_axis:
            # Remove interior axis labels (numbers, not ticks)
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

    # always add X1,X2 axis labels to the left/lower boundary plots
    axs_2d = np.atleast_2d(axs).reshape(subplot_grid)
    for row in range(axs_2d.shape[1]):
        axs_2d[-1, row].set_xlabel("$X^1$")
    for col in range(axs_2d.shape[0]):
        axs_2d[col, 0].set_ylabel("$X^2$")

    if quantity is not None:
        # Adding colorbar
        ev_contour = _symbol_check(ev_contour, [quantity])
        fig.colorbar(
            cm.ScalarMappable(norm=color_norm, cmap=discrete_cmap),
            ax=np.asarray(axs).ravel().tolist(),
            label=f"${ev_contour[quantity].symbol}$",
            **cbar_args,  # extend with red if data>vmax
        )

    return fig, axs


def plot_on_flux_surface(
    state: State,
    quantities: str | list[str] = "mod_B",
    rho: float | np.ndarray | list = 1.0,
    ntheta: int = 51,
    nzeta: int = 51,
    subplot_grid: list[int] | None = None,
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
        Additional keyword arguments to pass to the ``get_boozer`` method of the ``state`` object.
        These can be used to specify the Boozer transform parameters.
        For example the maximum mode number factor ``boozer_kwargs={'MNfactor': 3}``.


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

    # Work out the number of plots we want so we can design the grid if `subplot_grid` is not set
    if isinstance(quantities, str):
        quantities_eval = [quantities]
        nplots = rho_len
    elif isinstance(quantities, list):
        quantities_eval = quantities
        nplots = len(quantities)
    else:
        raise ValueError("quantities must be a string or a list of strings to evaluate.")

    if subplot_grid is None:
        subplot_grid = _design_subgrid(nplots)

    # Plotting boozer needs different evaluate function
    theta = np.linspace(0.0, 2 * np.pi, ntheta)
    zeta = np.linspace(0.0, 2 * np.pi / state.nfp, nzeta)
    if sfl:
        # Reduce some cost of the boozer computation
        if sfl == "boozer":
            if "MNfactor" not in boozer_kwargs:
                boozer_kwargs["MNfactor"] = 3
            if "radial_derivative" not in boozer_kwargs:
                boozer_kwargs["radial_derivative"] = False

        evaluations = state.evaluate_sfl(
            *quantities_eval, rho=rho, theta=theta, zeta=zeta, sfl=sfl, **boozer_kwargs
        )
    else:
        evaluations = state.evaluate(*quantities_eval, rho=rho, theta=theta, zeta=zeta)

    fig, axs = _subplots(subplot_grid, True, True, **plot_kwargs)

    evaluations = _symbol_check(evaluations, quantities_eval)

    # use a discrete cmap for the contours and the colorbar
    discrete_cmap = plt.get_cmap("viridis").resampled(levels)
    if isinstance(quantities, str) and share_contours:
        # Ensure the colormap is consistent across plots
        plotting_quantity = evaluations[quantities]
        # Make sure we are not trying to plot a vector field
        if "xyz" in plotting_quantity.coords:
            raise ValueError(
                f"Plotting quantity must be a scalar field but {quantities} is a vector field."
            )
        color_scale = (plotting_quantity.min().item(), plotting_quantity.max().item())
        color_norm = colors.Normalize(vmin=color_scale[0], vmax=color_scale[1])
    else:
        color_norm = None

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
            norm=color_norm,
            cmap=discrete_cmap,
        )

        if isinstance(quantities, list):
            # Specify the quantity type in the top right corner if multiple quantities were requested
            plot_label = f"${evaluations_i[quantity].attrs['symbol']}$"
            fig.colorbar(f_ax, ax=ax)
        else:
            plot_label = f"$\\rho={evaluations.rho.data[i]}$"
            if not share_contours:
                fig.colorbar(f_ax, ax=ax, label=f"${evaluations[quantity].symbol}$")

        ax.set(
            title=plot_label,
            xlabel=rf"${zeta_vals.attrs['symbol']} / (2\pi)$",
            ylabel=rf"${theta_vals.attrs['symbol']} / (2\pi)$",
        )

        # Removes any axis labels on subplots on the interior of the grid
        ax.label_outer()

    if isinstance(quantities, str) and share_contours:
        # Adding colorbar if single quantity was requested
        fig.colorbar(
            cm.ScalarMappable(color_norm, cmap=discrete_cmap),
            ax=np.asarray(axs).ravel().tolist(),
            label=f"${evaluations[quantity].symbol}$",
        )

    return fig, axs


def plot_fourier_on_surface(
    state: State,
    quantity: str = "mod_B",
    rho: float = 1.0,
    ntheta: int = 101,
    nzeta: int = 101,
    sfl: Literal["pest", "boozer"] | None = None,
    limit: float | None = 1e-15,
    contour: bool = False,
    plot_kwargs: dict[str] = {},
    **boozer_kwargs,
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
    contour: bool, default: False
        If ``True``, plot filled contours, otherwise plot on a rectangular grid.
    plot_kwargs: dict, optional
        Any ``**kwargs`` to send to the ``plt.figure()`` function.
        For example ``plot_kwargs={'figsize': (8,8)}``. See the `matplotlib documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html>`_ for a list of kwargs.
    boozer_kwargs : optional
        Additional keyword arguments to pass to the ``get_boozer`` method of the ``state`` object.
        These can be used to specify the Boozer transform parameters.
        For example the maximum mode number factor ``boozer_kwargs={'MNfactor': 3}``.


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
            *quantities, rho=rho, theta=ntheta, zeta=nzeta, sfl=sfl, **boozer_kwargs
        )
    evaluations = evaluations[quantities]

    evaluations = _symbol_check(evaluations, quantities)

    evft = ev2ft(evaluations).squeeze().sortby("n")

    levels = np.linspace(-14, 0, 8)
    symbol = evaluations[quantity].attrs.get("symbol", f"\\mathrm{{{quantity}}}")
    if not contour:
        cmap = cm.get_cmap("plasma", 7)
        norm = colors.LogNorm(vmin=1e-14, vmax=1)

    fig, axs = _subplots([1, 2], sharex=True, sharey=True, **plot_kwargs)
    for ax, suffix in zip(axs, ["mnc", "mns"]):
        if contour:
            c = ax.contourf(
                evft.n,
                evft.m,
                np.log10(np.maximum(np.abs(evft[f"{quantity}_{suffix}"]), 1e-16)),
                levels=levels,
                extend="both",
            )
        else:
            c = ax.pcolormesh(
                evft.n,
                evft.m,
                np.abs(evft[f"{quantity}_{suffix}"]),
                shading="nearest",
                cmap=cmap,
                norm=norm,
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
        limit_n1 = (
            power.n.where((power > limit).sum(dim="m") == 0).where(power.n > 0).min().item()
        )
        limit_n2 = (
            power.n.where((power > limit).sum(dim="m") == 0).where(power.n < 0).max().item()
        )
        if not np.isnan(limit_m):
            limit_m = np.max([limit_m, 5])
            for ax in axs:
                ax.set(ylim=(0, limit_m + 1))
        if not np.isnan(limit_n1) and not np.isnan(limit_n2):
            limit_n = np.max([abs(limit_n1), abs(limit_n2), 5])
            for ax in axs:
                ax.set(xlim=(-limit_n - 1, limit_n + 1))

    return fig, axs
