from plotly import graph_objects
import plotly.offline as plotly_offline

from gvec.core.state import State

from numpy import meshgrid

from .utils import _get_coord_range


def plot_3d_surface(
    state: State,
    nzeta: int,
    ntheta: int,
    rho: float,
    equilibrium_quantity: str = "mod_B",
    full_surface=False,
    offline: bool = False,
):
    """
    Generate a 3D surface plot with the quantity provided by `equilibrium_quantity` on it at a given `rho` position.

    Parameters
    ----------
    state: GVEC state file
    nzeta: int
        Toroidal resolution
    ntheta: int
        Poloidal resolution
    rho: float
        Radial position of the surface
    equilibrium_quantity: String, optional
        default is "mod_B"
    offline: bool
        If true, will automatically save the plot to a file in the current working directory.
        default is False

    Returns
    -------
    `plotly.plot` object


    plot a 3D surface at a given rho. By default plots the exterior boundary
    """

    nfp = state.nfp
    if full_surface:
        nfp = 1.0

    theta = _get_coord_range("theta", nfp, ntheta)
    zeta = _get_coord_range("zeta", nfp, nzeta)

    ev = state.evaluate(equilibrium_quantity, rho=[rho], theta=theta, zeta=zeta).sel(rho=rho)

    # if existing_plot is None:
    plt = graph_objects.Figure()
    # else:
    # plt = existing_plot
    plt.add_trace(
        graph_objects.Surface(
            x=ev.pos.sel(xyz="x"),
            y=ev.pos.sel(xyz="y"),
            z=ev.pos.sel(xyz="z"),
            surfacecolor=ev[equilibrium_quantity],
            colorbar_title_text=equilibrium_quantity,
        )
    )

    if offline:  # TODO See why plotly is not showing
        plotly_offline.plot(plt)
        print(
            "Plotly figure saved to working directory. To avoid saving call with `offline=False`."
        )

    return plt


def plot_boundary(state, nzeta, ntheta, equilibrium_quantity="mod_B", full_surface=True):
    """
    Plot the boundary of a GVEC solve.

    Parameters
    ----------
    state
    nzeta
    ntheta
    equilibrium_property

    TODO: Docstring

    Wrapper around `plot_3d_surface` with `rho=1.0`. See `help(plot_3d_surface)` for information on inputs.
    """
    return plot_3d_surface(
        state, nzeta, ntheta, 1.0, equilibrium_quantity, full_surface=full_surface
    )
