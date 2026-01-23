# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

import numpy as np

from gvec.core.state import State


def plot_3d_surface(
    state: State,
    quantity: str = "mod_B",
    rho: float | np.ndarray | list = 1.0,
    ntheta: int = 41,
    nzeta: int = 51,
    period: str = "single",
    to_file: str | None = None,
    surface_kwargs: dict = dict(),
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
    surface_kwargs : dict, optional
        Keyword arguments for the surface plot. i.e. `dict(opacity=0.5)`
    to_file : str
        If a string, will automatically save the plot to a file with the given input in the current working directory. Recommended to use this if the plots don't display.
        Default is ``None``

    Returns
    -------
    ``plotly.plot`` object
    """
    import plotly.offline as plotly_offline
    from plotly import graph_objects

    if period == "full":
        zeta_end = 2.0 * np.pi
        nzeta = nzeta * state.nfp
    elif period == "single":
        zeta_end = 2.0 * np.pi / state.nfp
    elif period == "half":
        zeta_end = np.pi / state.nfp
        nzeta = nzeta // 2
    else:
        raise ValueError("period must be 'full', 'single' or 'half'.")

    zeta = np.linspace(0.0, zeta_end, nzeta)
    theta = np.linspace(0.0, 2 * np.pi, ntheta)

    evaluations = state.evaluate(
        "pos", quantity, rho=np.asarray(rho).tolist(), theta=theta, zeta=zeta
    )

    fig = graph_objects.Figure()
    min_value, max_value = (
        np.amin(evaluations[quantity].data),
        np.amax(evaluations[quantity].data),
    )

    for irho, rhoval in enumerate(evaluations.rho.data):
        evaluation = evaluations.isel(rad=irho)
        fig.add_trace(
            graph_objects.Surface(
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
    # Ensure the figure correctly scales the axis to the data
    fig["layout"]["scene"]["aspectmode"] = "data"

    if isinstance(to_file, str):  # See why plotly is not showing
        plotly_offline.plot(fig, filename=to_file)

    return fig
