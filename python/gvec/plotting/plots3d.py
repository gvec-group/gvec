# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

import numpy as np

from gvec.core.state import State
from typing import Literal


def plot_3d_surface(
    state: State,
    quantity: str = "mod_B",
    rho: float | np.ndarray | list = 1.0,
    ntheta: int = 41,
    nzeta: int = 51,
    period: Literal["single", "half", "full"] = "single",
    zeta_contours: int = 0,
    to_file: str | None = None,
    surface_kwargs: dict = dict(),
    fig_in=None,
):
    """
    Generate a 3D surface plot with the quantity provided by ``quantity`` on it at a given ``rho`` position and
    poloidal and toroidal resolution ``ntheta`` and ``nzeta`` per field period.

    Parameters
    ----------
    state : GVEC state object
    quantity : str, optional
        Equilibrium quantity to plot on the surface.
        Default is ``"mod_B"``.
    rho : float, optional
        Radial position of the surface.
        Default ``1.0``.
    ntheta : int
        Poloidal resolution.
        Default ``41``.
    nzeta : int
        Toroidal resolution, per field period
        Default ``51``.
    period : str
        Plot the ``"full"`` surface, a ``"single"`` field period or ``"half"`` a field period.
        Default ``"single"``
    zeta_contours : int
        Number of contour lines of zeta coordinate on each surface, per field period. Default is ``0`` (no contours). Not counting the last contour ``max(zeta)``, which is always added.   Note: Adapts ``nzeta`` to ensure equally spaced contours.
    to_file : str
        If a string, will automatically save the plot to a file with the given input in the current working directory. Recommended to use this if the plots don't display.
        Default is ``None``
    surface_kwargs : dict, optional
        Keyword arguments for the ``Surface`` function of ``plotly``, i.e. ``dict(opacity=0.5)``
    fig_in : ``plotly.graph_objects.Figure`` object, optional
        If provided, the plot adds to this figure instead of creating a new one. Default is ``None``.

    Returns
    -------
    ``plotly.plot`` object
    """
    import plotly.offline as plotly_offline
    from plotly import graph_objects as go

    if period == "full":
        zeta_end = 2.0 * np.pi
        nzeta = nzeta * state.nfp
        nzeta_contours = zeta_contours * state.nfp
    elif period == "single":
        zeta_end = 2.0 * np.pi / state.nfp
    elif period == "half":
        zeta_end = np.pi / state.nfp
        nzeta = nzeta // 2
        nzeta_contours = zeta_contours // 2
    else:
        raise ValueError("period must be 'full', 'single' or 'half'.")

    if zeta_contours > 0:
        nzeta_skip = max(round((nzeta - 1) / nzeta_contours), 1)
        zeta = np.linspace(0.0, zeta_end, nzeta_skip * nzeta_contours + 1)
    else:
        zeta = np.linspace(0.0, zeta_end, nzeta)
    theta = np.linspace(0.0, 2 * np.pi, ntheta)

    evaluations = state.evaluate(
        "pos", quantity, rho=np.asarray(rho).tolist(), theta=theta, zeta=zeta
    )
    if fig_in is not None:
        assert isinstance(fig_in, go.Figure), (
            "fig_in must be a plotly.graph_objects.Figure object."
        )
        fig = fig_in
    else:
        fig = go.Figure()
    min_value, max_value = (
        np.amin(evaluations[quantity].data),
        np.amax(evaluations[quantity].data),
    )

    for irho, rhoval in enumerate(evaluations.rho.data):
        evaluation = evaluations.isel(rad=irho)
        fig.add_trace(
            go.Surface(
                x=evaluation.pos.sel(xyz="x"),
                y=evaluation.pos.sel(xyz="y"),
                z=evaluation.pos.sel(xyz="z"),
                surfacecolor=evaluation[quantity].broadcast_like(evaluation.pos.sel(xyz="x")),
                # colorbar_title_text=f"${evaluation[quantity].attrs['symbol']}$",
                # colorbar_title_text=f"${evaluation[quantity].attrs.get('symbol', quantity)}$",
                colorbar_title_text=quantity,
                showscale=(irho == 0),
                cmin=min_value,
                cmax=max_value,
                **surface_kwargs,
            )
        )
        if zeta_contours > 0:
            for z in zeta[::nzeta_skip]:
                fig.add_trace(
                    go.Scatter3d(
                        x=evaluation.pos.sel(xyz="x").sel(zeta=z),
                        y=evaluation.pos.sel(xyz="y").sel(zeta=z),
                        z=evaluation.pos.sel(xyz="z").sel(zeta=z),
                        mode="lines",
                        line=dict(color="black", width=2),
                        showlegend=False,
                    )
                )
    # Ensure the figure correctly scales the axis to the data
    fig["layout"]["scene"]["aspectmode"] = "data"

    if isinstance(to_file, str):  # See why plotly is not showing
        plotly_offline.plot(fig, filename=to_file)

    return fig
