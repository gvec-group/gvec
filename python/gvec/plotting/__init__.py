import matplotlib.pyplot as plt

from .plots1d import plot_on_axis, plot_radial_profile
from .plots2d import plot_on_flux_surface, plot_poloidal_plane
from .plots3d import plot_3d_surface, plot_boundary

plt.rcParams.update({"text.usetex": True})
