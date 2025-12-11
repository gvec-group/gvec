# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

import numpy as np

from gvec.core.state import State


def plot_3d_surface(
    state: State,
    rho: float,
    ntheta: int,
    nzeta: int,
    quantity: str = "mod_B",
    full_surface=False,
    to_file: str | None = None,
):
    """
    Generate a 3D surface plot with the quantity provided by `equilibrium_quantity` on it at a given `rho` position.

    Parameters
    ----------
    state : GVEC state file
    rho : float
        Radial position of the surface
    ntheta : int
        Poloidal resolution
    nzeta : int
        Toroidal resolution
    equilibrium_quantity : String, optional
        default is "mod_B"
    to_file : str
        If a string, will automatically save the plot to a file with the given input in the current working directory. Recommended to use this if the plots don't display.
        default is None

    Returns
    -------
    `plotly.plot` object


    plot a 3D surface at a given rho. By default plots the exterior boundary
    """
    import plotly.offline as plotly_offline
    from plotly import graph_objects

    nfp = state.nfp
    if full_surface:
        nfp = 1.0

    zeta = np.linspace(0.0, 2 * np.pi / nfp, nzeta)

    ev = state.evaluate(quantity, rho=[rho], theta=ntheta, zeta=zeta).sel(rho=rho)

    # if existing_plot is None:
    plt = graph_objects.Figure()
    # else:
    # plt = existing_plot
    plt.add_trace(
        graph_objects.Surface(
            x=ev.pos.sel(xyz="x"),
            y=ev.pos.sel(xyz="y"),
            z=ev.pos.sel(xyz="z"),
            surfacecolor=ev[quantity],
            colorbar_title_text=quantity,
        )
    )

    if isinstance(to_file, str):  # TODO See why plotly is not showing
        plotly_offline.plot(plt, filename=to_file)

    return plt


def plot_boundary(state, nzeta, ntheta, quantity="mod_B", full_surface=True, to_file=None):
    """
    Plot the boundary of a GVEC solve.

    Wrapper around `plot_3d_surface` with `rho=1.0`. See `help(plot_3d_surface)` for information on inputs.
    """
    return plot_3d_surface(
        state,
        1.0,
        ntheta,
        nzeta,
        quantity,
        full_surface=full_surface,
        to_file=to_file,
    )
