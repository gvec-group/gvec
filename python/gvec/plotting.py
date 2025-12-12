# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
from gvec.core.state import State
import logging
import numpy as np
from numpy.typing import ArrayLike

try:
    import matplotlib.pyplot as plt

    has_matplotlib = True
except ImportError:
    has_matplotlib = False

logger = logging.getLogger(__name__)


def plot_zeta_cuts(
    state: State,
    rho: ArrayLike = [0.01, 0.1, 0.5, 1.0],
    theta: ArrayLike | None = None,
    zeta: ArrayLike | None = None,
    axs: ArrayLike | None = None,
    figsize: tuple = (8, 8),
    **kwargs,
):
    """Plot polodial cross-sections with const. rho contours.

    Parameters
    ----------
    state : State
        GVEC state object used for the visualization.
    rho : ArrayLike, optional
        Radial contour levels, by default [0.01,0.1,0.5,1.0]
    theta : ArrayLike | None, optional
        Poloidal visualization points, by default 64.
    zeta : ArrayLike | None, optional
        Toroidal positions where cross-sections are drawn, by default three equally spaced.
    axs : ArrayLike | None, optional
        Axes to draw the contours into, by default None. If None, a new figure with axes is returned.
    figsize : tuple, optional
        Figure size as used by matplotlib.pyplot, by default (8,8).

    Returns
    -------
    tuple(fig,axs)
        Figure and axes of the plot, if axs was none, else just axs.
    """

    if not has_matplotlib:
        logger.warning("Failed to import matplotlib.pyplot! Is it installed?")
        return

    if theta is None:
        theta = np.linspace(0, 2 * np.pi, 64, endpoint=True)
    if zeta is None:
        zeta = np.linspace(0, 2 * np.pi / state.nfp, 3, endpoint=False)

    ev = state.evaluate(
        "X1",
        "X2",
        rho=rho,
        theta=theta,
        zeta=zeta,
    )
    return_ax = False
    if axs is None:
        return_ax = True
        fig, axs = plt.subplots(len(ev.zeta), tight_layout=True, figsize=figsize, dpi=300)
    for i, ax_i in enumerate(axs):
        for rho in ev.rho:
            ax_i.plot(
                ev.X1.sel(rho=rho, zeta=ev.zeta[i]),
                ev.X2.sel(rho=rho, zeta=ev.zeta[i]),
                **kwargs,
            )
        ax_i.axis("equal")
        ax_i.set(
            title=f"$\\zeta$={ev.zeta[i].data:.2f}",
            xlabel="X1",
            ylabel="X2",
        )
    if return_ax:
        return (fig, axs)
    else:
        return axs
