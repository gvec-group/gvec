# QUASR

The *QUAsi-symmetric Stellarator Repository* [QUASR](https://quasr.flatironinstitute.org/) [^QUASR1] [^QUASR2] [^QUASR3] is a database of curl-free stellarators optimized for volume quasi-symmetry.

A QUASR configuration can be loaded with
```{code} bash
pygvec load-quasr ID
```
where `ID` is replaced with the desired configuration. Alternatively
```{code} bash
pygvec load-quasr -f FILE
```
can be used to instead load the configuration from a simsopt compatible JSON file manually downloaded from QUASR.

This script will download the requested QUASR configuration, generate a *G-Frame* to be used as $h$-map and the boundary representation in that *G-Frame*, as well as a [GVEC parameter file](./gvec-parameter-list.md).

The `--tol` parameter sets the desired tolerance of the boundary representation which directly impacts the necessary degrees of freedom and therefore computational speed.
The `--nt` and `--nz` parameters set the number of points in $\vartheta$ and $\zeta$ respectively which are used to construct the *G-Frame*.
Other parameters can be seen with `pygvec load-quasr --help`.

<!--- References -->

[^QUASR1]: A. Giuliani, F. Wechsung, A. Cerfon, G. Stadler, M. Landreman (2022). Single-stage gradient-based stellarator coil design: Optimization for near-axis quasi-symmetry. Journal of Computational Physics, 459, 111147, [DOI: 10.1016/j.jcp.2022.111147](https://doi.org/10.1016/j.jcp.2022.111147).

[^QUASR2]: A. Giuliani, F. Wechsung, G. Stadler, A. Cerfon, & M. Landreman (2022). Direct computation of magnetic surfaces in Boozer coordinates and coil optimization for quasisymmetry. Journal of Plasma Physics, 88 (4), 905880401, [DOI: 10.1017/S0022377822000563](https://doi.org/10.1017/S0022377822000563).

[^QUASR3]: A. Giuliani, F. Wechsung, A. Cerfon, M. Landreman, G. Stadler (2023). Direct stellarator coil optimization for nested magnetic surfaces with precise quasi-symmetry. Phys. Plasmas, 30 (4), 042511, [DOI: 10.1063/5.0129716](https://doi.org/10.1063/5.0129716).
