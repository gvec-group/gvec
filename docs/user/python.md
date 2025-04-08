# Python Bindings
GVEC has *Python* bindings (referred to as *pyGVEC*) to run gvec and evaluate gvec equilibria from *Python*, while relying on the compiled *Fortran* functions to perform the actual computations.
This ensures post-processing is consistent with computing the equilibrium and also improves performance.


## Installation

Please follow the instructions for installing [**gvec and its python bindings**](install.md).


## Python interface to gvec state

The low-level interface is provided with a `State` class, which can be instantiated from a given *parameter-* and *statefile*:
```python
from gvec import State

with State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    ...
```

A high-level interface for evaluations is provided in form of a [xarray](https://docs.xarray.dev/) `Dataset` that can be generated with the `Evaluations` factory:
```python
from gvec import State, Evaluations

with State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    ev = Evaluations(rho=[0.1, 0.5, 0.9], theta="int", zeta="int", state=state)
    state.compute(ev, "B")
```
Here the additional arguments configure the points in the radial, poloidal and toroidal direction respectively, with `"int"` selecting the integration points that are used internally by GVEC.
The `ev` object is an instance of the `xarray.Dataset` and the individual `xarray.DataArray`s can then be accessed using `ev.B` or `ev["B"]`.
A `xarray.Dataset` closely mirrors the structure of netCDF, grouping several variables with named dimensions and coordinates as well as metadata.
The `state.compute` method can be used to compute various quantities that are then added to the `Dataset`.
Here *pyGVEC* takes care of computing all the required intermediate quantities, which are also added to the `Dataset`.

### Boozer transform

To evaluate the equilibrium in Boozer angles, you can use `BoozerEvaluations`, which performs a Boozer transform to obtain a set of $\vartheta,\zeta$ points which correspond to your desired grid in $\vartheta_B,\zeta_B$.
The evaluations with this new dataset works the same as above, not however that the suffixes `t` and `z` still refer to components/derivatives with respect to $\vartheta,\zeta$.
```python
from gvec import State, Evaluations

with State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    ev = gvec.EvaluationsBoozer(rho=[0.1, 0.5, 0.9], n_theta=51, n_zeta=51, state=state)
    state.compute(ev, "B")
```

## Available Quantities for Evaluation
The following table contains the quantities that can be evaluated with the python bindings.

```{include} ../generators/quantities.md
```
