# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
try:
    import matplotlib.pyplot as plt

    plt.rcParams["text.usetex"] = True
except ImportError:
    pass

from gvec.plotting.plots1d import plot_on_axis, plot_radial_profile
from gvec.plotting.plots2d import plot_on_flux_surface, plot_poloidal_plane
from gvec.plotting.plots3d import plot_3d_surface
