# %% [markdown]
#
# We present some examples for how to use the basic plotting functions in the python interface for GVEC. Note that
# these are not intended for advanced plotting, but for quickly visualising some basic properties.
# We will use the same elliptic tokamak GVEC case for the plotting examples.

# %%

import gvec
import gvec.plotting as gp
import numpy as np

# %%


params = {
    "ProjectName": "test_plotting",
    "which_hmap": 1,
    "PhiEdge": 1.0,
    "iota": {"type": "polynomial", "coefs": [0.625, 0.35]},
    "pres": {"type": "polynomial", "coefs": [1.0, -1.0], "scale": 1000.0},
    "nfp": 3,
    "X1_b_cos": {(0, 0): 3.0, (1, 0): 1.0, (1, 1): 0.4},
    "X2_b_sin": {(1, 0): 1.0, (1, 1): -0.4, (0, 1): -0.25},
    "init_average_axis": True,
    "sgrid_nElems": 2,
    "X1_mn_max": [3, 3],
    "X2_mn_max": [3, 3],
    "LA_mn_max": [3, 3],
    "X1X2_deg": 5,
    "LA_deg": 5,
    "totalIter": 1000,
    "minimize_tol": 1.0e-6,
}


try:
    state = gvec.find_state()
except ValueError:
    run = gvec.run(params)
    state = run.state


# %% [markdown]
#
# Note that all 1D and 2D plotting functions return the `matplotlib.figure.Figure` object, so you may edit the output if you wish.
#
# # 1D plots
# Two are two convience functions for plotting 1D properties.
# The first is to plot radial profiles, i.e. scalar functions of $\rho$. Note that only scalar quantities can be evaluated directly.

# %%

p_plot_radial_profile, ax_plot_radial_profile = gp.plot_radial_profile(state, 11)
p_plot_radial_profile.show()

# %% [markdown]
# Properties along the magnetic axis can be plotted with `plot_on_axis`. Since some derived quantities cannot be evaluated on axis,
# this will use quadratic extrapolation to generate the values for plotting. Also note that although we have specified the subplot grid here, this is not required unless
# you want a specific shape.

p_on_axis, ax_on_axis = gp.plot_on_axis(state, quantities=["mod_B", "X1"], subplot_grid=[2, 1])
p_on_axis.show()

# %% [markdown]
# # 2D plotting
#
# ## Poloidal slice plots
#
# The second and third inputs is the number of discrete $\theta$ and $\zeta$ points.
# The ``$X^1$`` and ``$X^2$`` values on the x and y axes can be shared by specifying `share_axis=True`.
# Note that the inputs for the 1D and 2D plots are very similar.

# %%

p_poloidal_slice_mod_B, ax_poloidal_slice_mod_B = gp.plot_poloidal_plane(
    state, 11, 11, share_axis=True, zeta=4
)
p_poloidal_slice_mod_B.show()


# %% [markdown]
# ## Flux surface plots
# As in the case with the 1D plots, the poloidal slice and flux surface plots share the same inputs.
# Note that multiple quantities can be evaluated and at different positions. By default the outermost flux surface,
# i.e. the boundary of the equilibrium, will be plotted.

# %%

p_plot_flux_surface, ax_plot_flux_surface = gp.plot_on_flux_surface(state, 11, 11)
p_plot_flux_surface.show()

# %% [markdown]
# Note that we can also evaluate different values on different flux surfaces. by specifying the radial positions of the flux surfaces with the
# fourth input and, and the quanitites with a keyword argument. Note that the quantities should be scalar values, as with the 1D plots.

# %%

p_plot_flux_surface_grid, ax_plot_flux_surface = gp.plot_on_flux_surface(
    state,
    11,
    11,
    np.array([0.3, 0.6]),
    subplot_grid=[2, 1],
)
p_plot_flux_surface_grid.show()

# %%
# By default, `plot_on_flux_surface` plots in Boozer coordinates. This can be switched to ``$\theta-\zeta$`` coordinates by setting `boozer=False`.
# We can also switch to filled contours

p_plot_flux_surface_grid, ax_plot_flux_surface = gp.plot_on_flux_surface(
    state, 11, 11, boozer=False, filled_contours=True
)
p_plot_flux_surface_grid.show()

# %%
# Instead of plotting multiple slices we can instead plot various quantities on the same surface.

# %%

p_plot_flux_surface_grid, ax_plot_flux_surface = gp.plot_on_flux_surface(
    state, 11, 11, quantities=["mod_B", "Jac"]
)
p_plot_flux_surface_grid.show()

# %% [markdown]
# # 3D plotting
#
# For 3D plotting we use [plotly](https://plotly.com/python/) as a backend.
# To plot the boundary we only need the state file and the resolution of the plot. Note that we can also specify the quantity we want to plot with
# the `quantities` keyword, by default $\|B\|$ will be plotted on the boundary.
#
# 3D plots return the `plotly.graph_objs._figure.Figure` object.
#

# %%

p_3d_boundary = gp.plot_boundary(state, 11, 11)
p_3d_boundary.show()


# %% [markdown]
# The `plot_boundary` function is implemented as a convience function sets `rho=1.0` as input to the function `plot_3d_surface`.
# If for some reason _plotly_ does not `show`, there is an optional input, `offline` (default False),
# which can be set to `True` to save the plot to your current working directory.

# %%

p_3d_surface = gp.plot_3d_surface(state, 11, 11, 0.5)
p_3d_surface.show()


# %% [markdown]
# # Advanced options
#
# As a final note there are some hidden inputs to the 1D and 2D plotting functions, specifically the keyword `post_process` option.
# If the plotting routines `plot_radial_profile` or `plot_on_axis` recieve `quantities` which are _not_ scalars,
# the `post_process` option can be used to adapt them for plotting.
#
# In the example below we will write a `lambda` function which takes the `ev.B` value in a GVEC `state.evaluate` object
# and spits out ``B_y``. The input to `post_process` should be a dictionary with the format,
# `{"name of field to perform operation on": {"function": <function to process>, "name": "name of new value"]}`
# the function should return a 1D numpy `ndarray`.
# This is limited, since it only takes a single value to do something to it, but it is not intended for
# complicated plotting.
#
# Note that we the "symbol" key is not required. The "name" will be used if "symbol" is not present.

# %%
post_process = {"B": {"function": lambda val: val.sel(xyz="y"), "name": "B_y", "symbol": "B_y"}}

p_line, ax_line = gp.plot_radial_profile(
    state,
    51,
    quantities=["mod_B", "B"],
    subplot_grid=[1, 2],
    post_process=post_process,
)
p_line.show()
