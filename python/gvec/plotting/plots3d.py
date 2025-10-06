from plotly import graph_objects
import plotly.offline as plotly_offline

from gvec.core.state import State

from gvec_plotting.utils import _get_coord_range


def plot_3d_surface(
    state: State,
    nzeta: int,
    ntheta: int,
    equilibrium_quantity: str = "mod_B",
    rho: float = 1.0,
    offline: bool = True,
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
    equilibrium_quantity: String, optional
        default is "mod_B"
    rho: float
        default is 1.0
    offline: bool
        If true, will automatically save the plot to a file in the current working directory.
        default is True

    Returns
    -------
    `plotly.plot` object


    plot a 3D surface at a given rho. By default plots the exterior boundary
    """

    theta = _get_coord_range("theta", state.nfp, ntheta)
    zeta = _get_coord_range("zeta", state.nfp, nzeta)

    ev = state.evaluate(equilibrium_quantity, rho=[rho], theta=theta, zeta=zeta).sel(
        rho=rho
    )

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
        )
    )

    if offline:  # TODO See why plotly is not showing
        plotly_offline.plot(plt)
        print(
            "Plotly figure saved to working directory. To avoid saving call with `offline=False`."
        )

    return plt


def plot_boundary(state, nzeta, ntheta, equilibrium_quantity="mod_B"):
    """
    Plot the boundary of a GVEC solve.

    Parameters
    ----------
    state
    nzeta
    ntheta
    equilibrium_property

    Wrapper around `plot_3d_surface` with `rho=1.0`. See `help(plot_3d_surface)` for information on inputs.
    """
    return plot_3d_surface(state, nzeta, ntheta, equilibrium_quantity)
