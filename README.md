# GVEC - 3D MHD Equilibrium Solver

[![License](https://img.shields.io/badge/license-MIT-red)](./LICENSE.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15026781.svg)](https://doi.org/10.5281/zenodo.15026781)

## Overview

GVEC (Galerkin Variational Equilibrium Code) is an open-source software for
the generation of three-dimensional ideal magnetohydrodynamic (MHD) equilibria.

The main features of GVEC are:

* The equilibrium is found by **minimizing the MHD energy** under the assumption of closed nested flux surfaces. This approach is based on [VMEC](https://princetonuniversity.github.io/STELLOPT/VMEC) (Hirshman & Whitson, 1983).
* The equilibrium is found with a **fixed plasma boundary shape** and given radial profiles (rotational transform and pressure).
* **High Order B-spline** discretization in the radial direction with **a smooth representation** at the magnetic axis.
* **Double-angle Fourier** representation in the poloidal and toroidal direction of the flux surfaces, with different resolutions for each solution variable $X^1,X^2,\lambda$. Stellarator symmetry may be explicitly imposed.
* **Flexible choice of the mapping** between the logical and cartesian space $\left(X^1,X^2,\zeta\right) \mapsto \left(x,y,z\right)$
  to find equilibria in complex-shaped domains (magnetic islands, knotted domains...).

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany.

The list of contributors is found in [CONTRIBUTORS.md](CONTRIBUTORS.md).
Outside contributions are always welcome!

## Documentation

 * [user and developer documentation](https://gvec-group.pages.mpcdf.de/gvec) built with *sphinx*
   * [Installation](https://gvec-group.pages.mpcdf.de/gvec/user/install.html)
   * [Getting Started](https://gvec-group.pages.mpcdf.de/gvec/user/getting-started.html)
   * [Theoretical considerations](https://gvec-group.pages.mpcdf.de/gvec/user/theory.html)
 * auto-generated [fortran code documentation](https://gvec-group.pages.mpcdf.de/gvec/ford/index.html) built with [FORD](https://forddocs.readthedocs.io/en/latest/)

## Installation & Getting started

See the documentation on [Installation](https://gvec-group.pages.mpcdf.de/gvec/user/install.html) and [Getting Started](https://gvec-group.pages.mpcdf.de/gvec/user/getting-started.html).

## Reporting Bugs & Contributing to GVEC

Please report any bugs you find on the [Issue tracker](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/issues).

Contributions are always welcome, best get into contact directly with the maintainers.
Also see the relevant [documentation](https://gvec-group.pages.mpcdf.de/gvec/user/index.html).

## License

GVEC is released under the terms of the [MIT License](https://spdx.org/licenses/MIT.html).
For the full license terms see the included [LICENSE.txt](LICENSE.txt) file.

Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics

Parts of this software are licensed differently:
* `src/base/bsplines/` is part of [SeLaLib](https://github.com/selalib/selalib/) and licensed with `CECILL-B`.
* `src/mod_timings.f90` & `src/perf2timings.f90` are wrappers for the [ftimings](https://gitlab.mpcdf.mpg.de/loh/ftimings) library, licensed with `LGPL-3.0-only`.
* `src/globals/cla.f90` is [CLAF90](https://ingria.ceoas.oregonstate.edu/fossil/CLAF90) licensed with a modified `MIT` license.

## Citing GVEC

If you use GVEC in your work, please be sure to cite the following Zenodo record:

https://zenodo.org/records/15026781

with the DOI: [`10.5281/zenodo.15026781`](https://doi.org/10.5281/zenodo.15026781)

A bibtex entry is found in [CITATION.bib](CITATION.bib).


## References

This is a list of references in which the GVEC equilibrium solver was utilized:

|   |        |
|----- |------- |
| [HPM25] | Florian Hindenlang, Gabriel G Plunk, and Omar Maj. *Computing MHD equilibria of stellarators with a flexible coordinate frame*. Plasma Physics and Controlled Fusion, 67(4):045002, mar 2025. doi:10.1088/1361-6587/adba11.|
| [PDR+25] | Gabriel G Plunk, Michael Drevlak, Eduardo Rodríguez, Robert Babin, Alan Goodman, and Florian Hindenlang. *Back to the figure-8 stellarator*. Plasma Physics and Controlled Fusion, 67(3):035025, feb 2025. doi:10.1088/1361-6587/adb64b.|
| [PDS+23] | Jonas Puchmayr, Mike G Dunne, Erika Strumberger, Matthias Willensdorfer, Hartmut Zohm, and Florian Hindenlang. *Helical mode localization and mode locking of ideal MHD instabilities in magnetically perturbed tokamak plasmas*. Nuclear Fusion, 2023. |
| [MND+20] | *Maurice Maurer, A Banon Navarro, Tilman Dannert, Marco Restelli, Florian Hindenlang, Tobias Goerler, Daniel Told, Denis Jarema, Gabriele Merlo, and Frank Jenko*. GENE-3D: a global gyrokinetic turbulence code for stellarators. Journal of Computational Physics, 420:109694, 2020.|
| [NRH+22] | Nikita Nikulsin, Rohan Ramasamy, Matthias Hoelzl, Florian Hindenlang, Erika Strumberger, Karl Lackner, Sibylle Guenter, JOREK Team, and others. *JOREK3D: an extension of the JOREK nonlinear MHD code to stellarators*. Physics of Plasmas, 2022.|
| [NMP+20] | A Banon Navarro, G Merlo, G G Plunk, P Xanthopoulos, A Von Stechow, A Di Siena, M Maurer, F Hindenlang, F Wilms, and F Jenko. *Global gyrokinetic simulations of ITG turbulence in the magnetic configuration space of the Wendelstein 7-X stellarator*. Plasma Physics and Controlled Fusion, 62(10):105005, 2020.|
| [WNM+21] | Felix Wilms, Alejandro Bañón Navarro, Gabriele Merlo, Leonhard Leppin, Tobias Görler, Tilman Dannert, Florian Hindenlang, and Frank Jenko. *Global electromagnetic turbulence simulations of W7-X-like plasmas with GENE-3D*. Journal of Plasma Physics, 87(6):905870604, 2021. |
