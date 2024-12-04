# GVEC

## Overview

GVEC (Galerkin Variational Equilibrium Code) is an open-source software for
the generation of three-dimensional ideal MHD equilibria.
The ideas are strongly based on on [VMEC](https://princetonuniversity.github.io/STELLOPT/VMEC) (Hirshman & Whitson, 1983).

The main features of GVEC are

* Use of modern **object-oriented FORTRAN**
* **Radial High Order B-spline** discretization
* **Number of Fourier modes** for each variable $X^1,X^2,\lambda$ can be chosen separately
* **Flexible choice of the mapping** between the space $\left(X^1,X^2,\zeta\right) \mapsto \left(x,y,z\right)$ (in VMEC fixed to $\left(R,Z,\phi\right)\mapsto\left(x,y,z\right)$)
  to find equilibria in complex-shaped domains (magnetic islands, knotted domains...)
* A VMEC generated netcdf file can be used for initialization of GVEC.
* GVEC can convert its final solution to several formats (VMEC netcdf, specific converters)

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany.

The list of contributors is found in [CONTRIBUTORS.md](CONTRIBUTORS.md).
Outside contributions are always welcome!

## License

GVEC is released under the terms of the GNU General Public License v3.0.
For the full license terms see the included [license](LICENSE) file.

Copyright 2017-2024 GVEC Contributors | Max Planck Institute for Plasma Physics

## Documentation

 * a pdf that documents the [theory and implementation details](https://gitlab.mpcdf.mpg.de/gvec-group/GVEC_doc/blob/master/GVEC_prototype/GVEC_prototype.pdf) of GVEC
 * [user and developer documentation](https://gvec-group.pages.mpcdf.de/gvec) built with *sphinx*
 * auto-generated [fortran code documentation](https://gvec-group.pages.mpcdf.de/gvec/ford/index.html) built with [FORD](https://forddocs.readthedocs.io/en/latest/)

## Installation & Getting started

See the documentation on [Installation](https://gvec-group.pages.mpcdf.de/gvec/user/install.html) and [Getting Started](https://gvec-group.pages.mpcdf.de/gvec/user/getting-started.html).

## Reporting Bugs & Contributing to GVEC

Please report any bugs you find on the [Issue tracker](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/issues).

Contributions are always welcome, best get into contact directly with the maintainers.
Also see the relevant [documentation](https://gvec-group.pages.mpcdf.de/gvec/user/index.html).
