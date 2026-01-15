# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT

from gvec.core.state import State
from gvec.plotting.plots1d import plot_on_axis, plot_radial_profile
from gvec.plotting.plots2d import plot_on_flux_surface, plot_poloidal_plane
from gvec.plotting.plots3d import plot_3d_surface

State.plot_on_axis = plot_on_axis
State.plot_radial_profile = plot_radial_profile
State.plot_on_flux_surface = plot_on_flux_surface
State.plot_poloidal_plane = plot_poloidal_plane
State.plot_3d_surface = plot_3d_surface
