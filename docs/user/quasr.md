# QUASR and G-Frame

## GVEC-QUASR interface
:::{note}
The QUASR interface requires [`simsopt`](https://github.com/hiddenSymmetries/simsopt) to be installed.
:::

The *QUAsi-symmetric Stellarator Repository* [QUASR](https://quasr.flatironinstitute.org/) [^QUASR1] [^QUASR2] [^QUASR3] is a database of curl-free stellarators optimized for volume quasi-symmetry.

A QUASR configuration can be loaded with
```{code} bash
pygvec load-quasr ID -v
```
where `ID` is replaced with the desired configuration and `-v` is for verbose. Alternatively
```{code} bash
pygvec load-quasr -s FILE -v
```
can be used instead, to load a boundary from a `simsopt` compatible JSON file (e.g. manually downloaded from QUASR).
With
```{code} bash
pygvec load-quasr -f FILE -v
```
the cartesian boundary data is read directly from the supplied `netCDF` file (in this case `simsopt` is also not required).

In the script, depending on the type of input, the boundary surface data is used to construct a *G-Frame*, which is described in [^GFrame]. The script writes a `netCDF` file containing the *G-Frame* and the boundary cross-sections, and along with it a first [GVEC parameter file](./gvec-parameter-list.md).

* The `-v` is for verbose mode, which prints the steps executed by the script to screen.
* The `--clean` parameter sets the desired tolerance to remove modes from the input surface (cleanup), and therefore changes the boundary surface!
* The `--stellsym` flag imposes stellarator symmetry for the input surface. It is recommended to first check without this flag and see if necessary and applicable. Changes the boundary surface!
* The `--cutoff` parameter sets a desired cutoff toroidal mode number for the G-Frame construction (does not change the boundary surface).
* The `--tol` parameter sets the desired tolerance for determining minimal necessary Fourier modes for the output boundary X1,X2.
* The `--nt` and `--nz` parameters set the number of points in $\vartheta$ and $\zeta$ respectively for one field period, from which a *G-Frame* as well as the boundary cross-sections are computed. The points exclude the periodic endpoint and should be chosen to be odd.
* With `--save-xyz` the cartesian boundary data can be saved as a `netCDF` file.
* All parameters can be seen with `pygvec load-quasr --help`.

## Construction from other surface data

One can also use other surfaces from which to construct a *G-Frame* and its corresponding cross-sections, providing the surface data only.

The surface cartesian position $(x,y,z)(\vartheta_i,\zeta_j)$ must be evaluated
at a meshgrid on the full torus, excluding the periodic endpoint:

$\vartheta_i=2\pi \frac{i}{n_\vartheta},i=0\dots,n_\vartheta-1,\quad \zeta_j=2\pi\frac{j}{n_\zeta},j=0,\dots,n_\zeta-1$

$n_\vartheta,n_\zeta$ should be chosen to be odd.

Below an example, with an `eval_surface` function that evaluates the surface in cartesian coordinates,
provided the number of points `ntheta` and `nzeta` over one field period and the number of field periods `nfp`.

```python
import gvec
import numpy as np
theta=np.linspace(0,2*np.pi,ntheta,endpoint=False)
zeta=np.linspace(0,2*np.pi,nzeta*nfp,endpoint=False)
xyz=np.zeros((nzeta,ntheta,3))
for j in range(nzeta):
        for i in range(ntheta):
            xyz[j,i,:] = eval_surface(theta[i],zeta[j])
gvec.scipts.quasr.save_xyz(xyz,nfp,'surface.nc')
```
Once the surface is written to file, the above script can be used:
```bash
pygvec load-quasr -f surface.nc -v
```
The construction of the G-Frame uses the toroidal angle $\zeta$ from the data, a recommended choice is to use the Boozer angle for $\zeta$.
:::{note}
If the surface data is computed from cylinder coordinates, the construction will not change this, and the resulting cross-sections will not change shape!
:::

## Visualization of the G-Frame file

The script above produces a `netCDF` file containing the *G-Frame* and the corresponding boundary. One can visualize both via VTK.

```python
import gvec
gvec.vtk.gframe_to_vtk("netcdf-file-gframe.nc")
```
This produces a `.vts` file for the axis of the G-Frame and of the boundary surface, which can be viewed e.g. in Paraview. More visualization options are available, see docstring of `gframe_to_vtk`.


<!--- References -->

[^QUASR1]: A. Giuliani, F. Wechsung, A. Cerfon, G. Stadler, M. Landreman (2022). Single-stage gradient-based stellarator coil design: Optimization for near-axis quasi-symmetry. Journal of Computational Physics, 459, 111147, [DOI: 10.1016/j.jcp.2022.111147](https://doi.org/10.1016/j.jcp.2022.111147).

[^QUASR2]: A. Giuliani, F. Wechsung, G. Stadler, A. Cerfon, M. Landreman (2022). Direct computation of magnetic surfaces in Boozer coordinates and coil optimization for quasisymmetry. Journal of Plasma Physics, 88 (4), 905880401, [DOI: 10.1017/S0022377822000563](https://doi.org/10.1017/S0022377822000563).

[^QUASR3]: A. Giuliani, F. Wechsung, A. Cerfon, M. Landreman, G. Stadler (2023). Direct stellarator coil optimization for nested magnetic surfaces with precise quasi-symmetry. Phys. Plasmas, 30 (4), 042511, [DOI: 10.1063/5.0129716](https://doi.org/10.1063/5.0129716).

[^GFrame]: F. Hindenlang, G. Plunk, O. Maj (2025). Computing MHD equilibria of stellarators with a flexible coordinate frame. Plasma Physics and Controlled Fusion, 67, 045002, [DOI: 10.1088/1361-6587/adba11](https://doi.org/10.1088/1361-6587/adba11).
