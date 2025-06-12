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

To evaluate the equilibrium in Boozer angles, you can use `EvaluationsBoozer`, which performs a Boozer transform to obtain a set of $\vartheta,\zeta$ points which correspond to your desired grid in $\vartheta_B,\zeta_B$.
The evaluations with this new dataset work the same as above, note however that the suffixes `t` and `z` still refer to components/derivatives with respect to $\vartheta,\zeta$.
Some additional quantities, like `B_contra_t_B` or `e_zeta_B` are now also available.
```python
from gvec import State, Evaluations

with State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    ev = gvec.EvaluationsBoozer(rho=[0.1, 0.5, 0.9], theta_B=20, zeta_B=40, state=state, MNfactor=5)
    state.compute(ev, "B")
```

Similar to `Evaluations`, the grid of coordinates can be specified with an integer for equidistant spacing or explicitly with a list, array or DataArray.
The optional `MNfactor` sets the maximum fourier modes for the boozer transform `M`, `N` to the specified multiple of the highest modenumber of $X^1,X^2,\lambda$.
`M`, `N` can also be specified directly.

```{note}
The Boozer transform recomputes $\lambda$ with a higher resolution (to satisfy the integrability condition for $\nu_B$)!
Therefore some quantities will differ between the equilibrium evaluation and Boozer evaluation.

In particular $\langle B_\theta \rangle, \langle B_\zeta \rangle$ will differ from $B_{\theta_B},B_{\zeta_B}$ by an offset.
```

### Fieldline aligned grid

Some applications require a fieldline-aligned grid, which can be generated using `EvaluationsBoozerCustom`:
```python
import numpy as np
import gvec

with gvec.State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    rho = [0.5, 1.0]  # radial positions
    alpha = np.linspace(0, 2 * np.pi, 20, endpoint=False)  # fieldline label
    phi = np.linspace(0, 2 * np.pi / state.nfp, 101)  # angle along the fieldline

    # evaluate the rotational transform (fieldline angle) on the desired surfaces
    ev = gvec.Evaluations(rho=rho, theta=None, zeta=None, state=state)
    state.compute(ev, "iota")

    # 3D toroidal and poloidal arrays that correspond to fieldline coordinates for each surface
    theta_B = alpha[None, :, None] + ev.iota.data[:, None, None] * phi[None, None, :]

    # create the grid
    ev = gvec.EvaluationsBoozerCustom(rho=rho, theta_B=theta_B, zeta_B=phi, state=state, MNfactor=5)

    # set the fiedline label as coordinate & index
    ev["alpha"] = ("pol", alpha)
    ev["alpha"].attrs = dict(symbol=r"\alpha", long_name="fieldline label")
    ev = ev.set_coords("alpha").set_xindex("alpha")

    state.compute(ev, "B")
```

## Available Quantities for Evaluation
The following table contains the quantities that can be evaluated with the python bindings.

```{include} ../generators/quantities.md
```
