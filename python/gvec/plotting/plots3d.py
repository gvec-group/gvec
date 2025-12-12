# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

import numpy as np

from gvec.core.state import State


def plot_3d_surface(
    state: State,
    quantity: str = "mod_B",
    rho: float = 1.0,
    ntheta: int = 41,
    nzeta: int = 51,
    full_surface: bool = False,
    to_file: str | None = None,
):
    """
    Generate a 3D surface plot with the quantity provided by `quantity` on it at a given `rho` position and
    poloidal and toroidal resolution `ntheta` and `nzeta`.

    Parameters
    ----------
    state : GVEC state file
    quantity : str, optional
        Equilibrium quantity to plot on the surface.
        Default is "mod_B".
    rho : float, optional
        Radial position of the surface.
        Default `1.0`.
    ntheta : int
        Poloidal resolution.
        Default `41`.
    nzeta : int
        Toroidal resolution.
        Default `51`.
    to_file : str
        If a string, will automatically save the plot to a file with the given input in the current working directory. Recommended to use this if the plots don't display.
        Default is `None`.

    Returns
    -------
    `plotly.plot` object
    """
    import plotly.offline as plotly_offline
    from plotly import graph_objects

    nfp = state.nfp
    if full_surface:
        nfp = 1.0

    zeta = np.linspace(0.0, 2 * np.pi / nfp, nzeta)
    theta = np.linspace(0.0, 2 * np.pi, ntheta)

    evaluation = state.evaluate(quantity, rho=[rho], theta=theta, zeta=zeta).sel(rho=rho)

    plt = graph_objects.Figure()

    plt.add_trace(
        graph_objects.Surface(
            x=evaluation.pos.sel(xyz="x"),
            y=evaluation.pos.sel(xyz="y"),
            z=evaluation.pos.sel(xyz="z"),
            surfacecolor=evaluation[quantity],
            # colorbar_title_text=f"${evaluation[quantity].attrs['symbol']}$",
            # colorbar_title_text=f"${evaluation[quantity].attrs.get('symbol', quantity)}$",
            colorbar_title_text=quantity,
        )
    )

    if isinstance(to_file, str):  # See why plotly is not showing
        plotly_offline.plot(plt, filename=to_file)

    return plt
