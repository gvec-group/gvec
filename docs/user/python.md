# Postprocessing with python
GVEC has *Python* bindings (referred to as *pyGVEC*) to run gvec and evaluate gvec equilibria from *Python*, while relying on the compiled *Fortran* functions to perform the actual computations.
This ensures post-processing is consistent with computing the equilibrium and also improves performance.


## Installation

Please follow the instructions for installing [**gvec and its python bindings**](install.md)

and confirm successful installation with
```python
import gvec
```

## The `gvec.State` class

The central object for evaluating a GVEC equilibrium is the `gvec.State` class.

A state can be loaded from a given *parameter-* and *statefile* with
```python
state = gvec.load_state("parameter.ini", "EXAMPLE_State_0001_00001000.dat")
```
or from a run directory with:
```python
state = gvec.find_state("rundir")
```

If a run directory contains multiple states, all can be loaded with
```python
states = gvec.find_states("rundir")
```

A state can also be created directly from a given *parameter-* and *statefile* with
```python
state = gvec.State("parameter.ini", "EXAMPLE_State_0001_00001000.dat")
```
or initialised with a parameter dictionary using:
```python
state = gvec.State.new(parameters, "rundir")
```

## Evaluating a state

A variety of quantities can be evaluated for a given state using the `state.evaluate` function.
This function creates an evaluations dataset (the grid on which the state is evaluated)
and recursively computes all desired quantities and any quantities they depend on:
```python
ev = state.evaluate("pos", "B", rho=[0.1, 0.5, 0.9], theta=20, zeta=0.0)
```

A list of all computable quantities is given [below](#available-quantities-for-evaluation).

The grid can be specified explicitly (with an array or float) or automatically with either an integer
for linear spacing (within one field period) or `"int"` for the integration points required by GVEC.

An existing dataset can be extended using `gvec.compute`:
```python
gvec.compute(ev, "J", "V", state=state)
```

The evaluations dataset `ev` is an [xarray](https://docs.xarray.dev/) `Dataset`, which groups several variables,
their coordinates and metadata, similar to the *netCDF* format or *pandas dataframes*.

The individual quantities can be accessed using `ev.B` or `ev["B"]` and are enriched with coordinate information
which can be used to a variety of operations, e.g.:
```python
import xarray as xr

mean_B = ev.B.mean(dim=("rad","pol","tor"))
B2 = xr.dot(ev.B, ev.B, dim="xyz")
JxB = xr.cross(ev.J, ev.B, dim="xyz")
```
Indexing is best done using `ev.B.sel(rho=0.5)` (by value) or `ev.B.isel(rad=0)` (by position)
and the raw `numpy` array can be extracted with `ev.B.values`.
To ensure a particular order, use `ev.B.transpose("rad", "pol", "tor", "xyz").values`.
To convert a scalar (e.g. the volume `V`) into a python scalar use the `ev.V.item()` method.

## Computing integrals

Specifying `"int"` for `rho`, `theta` and `zeta` in the `Evaluations` function, chooses the grid to be the integration points used by GVEC internally.
The functions `radial_integral`, `fluxsurface_integral` and `volume_integral` can then be used to perform integration with the apropriate weights.
E.g. to compute the volume averaged plasma beta, one could use:

```python
ev = state.evaluate("mod_B", "mu0", "p", "Jac", "V", rho="int", theta="int", zeta="int")
beta = ev.p / (ev.mod_B**2 / (2 * ev.mu0))
beta_avg = gvec.volume_integral(beta * ev.Jac) / ev.V
```

## Low-level access

The `gvec.State` class is automatically *bound* to the fortran library when needed and allows evaluating a variety of quantities.
E.g. the number of field periods can be accessed with `state.nfp` and the maximum fourier modes with `state.get_mn_max()`.
For a full list of availble methods see the [API](<dev/api/core-state>).

The debugging output for the postprocessing can be accessed with `state.stdout`.

## Boozer transform

To evaluate the equilibrium in Boozer angles, you can use `gvec.evaluate_sfl(..., sfl="boozer")`, which performs a Boozer transform to obtain a set of $\vartheta,\zeta$ points which correspond to your desired grid in $\vartheta_B,\zeta_B$.
The evaluations with this new dataset work the same as above, note however that the suffixes `t` and `z` still refer to components/derivatives with respect to $\vartheta,\zeta$.
Some additional quantities, like `B_contra_t_B` or `e_zeta_B` are now also available.
```python
ev = state.evaluate_sfl("B_contra_t_B", "B_contra_z_B", rho=[0.1, 0.5, 0.9], theta=20, zeta=40, MNfactor=5)
```

Similar to `Evaluations`, the grid of coordinates can be specified with an integer for equidistant spacing or explicitly with a list, array or DataArray.
The optional `MNfactor` (default `5`) sets the maximum fourier modes for the boozer transform `M`, `N` to the specified multiple of the highest modenumber of $X^1,X^2,\lambda$.
`M`, `N` can also be specified directly.

```{note}
The Boozer transform recomputes $\lambda$ with a higher resolution (to satisfy the integrability condition for $\nu_B$)!
Therefore some quantities will differ between the equilibrium evaluation and Boozer evaluation.

In particular $\langle B_\vartheta \rangle, \langle B_\zeta \rangle$ will differ from $B_{\vartheta_B},B_{\zeta_B}$ by an offset.

Currently the Boozer transform is performed for each surface individually and radial derivatives are therfore not available.
This means that that $\frac{\partial \mathbf{B}}{\partial \rho}$ and $\mathbf{J}$ are not available!
```

### Field-aligned grid

Some applications require a fieldline-aligned grid, which can be generated using `EvaluationsBoozerCustom`:
```python
import numpy as np
import gvec

state = gvec.State("parameter.ini", "EXAMPLE_State_0001_00001000.dat")
rho = [0.5, 1.0]  # radial positions
alpha = np.linspace(0, 2 * np.pi, 20, endpoint=False)  # fieldline label
phi = np.linspace(0, 2 * np.pi / state.nfp, 101)  # angle along the fieldline

# evaluate the rotational transform (fieldline angle) on the desired surfaces
iota = state.evaluate("iota", rho=rho, theta=None, zeta=None).iota

# 3D toroidal and poloidal arrays that correspond to fieldline coordinates for each surface
theta_B = alpha[None, :, None] + iota.data[:, None, None] * phi[None, None, :]

# create the grid
ev = gvec.EvaluationsBoozerCustom(rho=rho, theta_B=theta_B, zeta_B=phi, state=state, MNfactor=5)

# set the fiedline label as poloidal coordinate & index
ev["alpha"] = ("pol", alpha)
ev["alpha"].attrs = dict(symbol=r"\alpha", long_name="fieldline label")
ev = ev.set_coords("alpha").set_xindex("alpha")

state.compute(ev, "B")
```

## Available Quantities for Evaluation
The following table contains the quantities that can be evaluated with the python bindings.

```{include} ../generators/quantities.md
```
