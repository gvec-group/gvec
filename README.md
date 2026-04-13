# GVEC - 3D MHD Equilibrium Solver

[![readthedocs](https://img.shields.io/badge/Read%20the%20Docs-8CA1AF?logo=readthedocs&labelColor=gray)](https://gvec.readthedocs.io)
[![gitlab](https://img.shields.io/badge/GitLab-FC6D26?logo=gitlab&labelColor=gray)](https://gitlab.mpcdf.mpg.de/gvec-group/gvec)
[![github-mirror](https://img.shields.io/badge/GitHub%20mirror-gray?logo=github&labelColor=gray)](https://github.com/gvec-group/gvec)
[![matrix-chat](https://img.shields.io/badge/chat-gvec%3Ampg.de-0DBD8B?logo=matrix)](https://matrix.to/#/#gvec:mpg.de)

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.09670/status.svg)](https://doi.org/10.21105/joss.09670)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.15026780.svg)](https://doi.org/10.5281/zenodo.15026780)
[![GitLab Release](https://img.shields.io/gitlab/v/release/gvec-group%2Fgvec?gitlab_url=https%3A%2F%2Fgitlab.mpcdf.mpg.de&logo=gitlab)](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/releases)
[![PyPI](https://img.shields.io/pypi/v/gvec?logo=PyPI&label=PyPI)](https://pypi.org/project/gvec/)

[![License](https://img.shields.io/badge/license-MIT-green?logo=open-source-initiative)](./LICENSE.txt)
[![python-version](https://img.shields.io/pypi/pyversions/gvec?logo=python)](https://pypi.org/project/gvec)
[![pre-commit](https://img.shields.io/badge/pre--commit-gray?logo=pre-commit)](https://github.com/pre-commit/pre-commit)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

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

 * [user and developer documentation](https://gvec.readthedocs.io/latest) built with *sphinx* and hosted on *readthedocs*.
   * [Installation](https://gvec.readthedocs.io/latest/user/install.html)
   * [Getting Started](https://gvec.readthedocs.io/latest/tutorials/getting-started.html)
   * [Theoretical considerations](https://gvec.readthedocs.io/latest/user/theory.html)
 * auto-generated [fortran code documentation](https://gvec.readthedocs.io/latest/ford/index.html) built with [FORD](https://forddocs.readthedocs.io/en/latest/)

## Statement of need

MHD equilibrium solutions are the basis for a number of high fidelity plasma physics models and associated codes.
For example, they provide the initial conditions for linear and nonlinear MHD solvers (e.g. CASTOR3D ([Puchmayr et al., 2023](https://doi.org/10.1088/1741-4326/acdd12)), CAS3D ([Schwab, 1993](https://doi.org/10.1063/1.860656)), Jorek3D ([Nikulsin et al., 2022](https://doi.org/10.1063/5.0087104)), Struphy ([Holderied et al., 2022](https://doi.org/10.1016/j.jcp.2021.110143)), M3D-C1 ([Jardin, 2004](https://doi.org/10.1016/j.jcp.2004.04.004))), or the magnetic field for particle orbit tracing (e.g. SIMPLE ([Albert et al., 2020](https://doi.org/10.1016/j.jcp.2019.109065))) and turbulence simulations (e.g. BOUT++ ([Shanahan et al., 2024]()), GENE ([Bañón Navarro et al., 2020](https://doi.org/10.1088/1361-6587/aba858))).
3D MHD equilibria are directly used to analyse potential stellarator configurations in optimisation frameworks, such as SIMSOPT ([Landreman et al., 2021](https://doi.org/10.21105/joss.03525)) or STELLOPT ([Lazerson et al., 2020](https://doi.org/10.11578/dc.20180627.6)).

GVEC has a flexible coordinate frame, allowing it to represent boundary shapes beyond those possible with the standard cylindrical coordinates used by many equilibrium codes. This can be used for example to optimise the boundary shape of the figure-8 stellarator in [Plunk et al. (2025)](https://doi.org/10.1088/1361-6587/adb64b).

## Installation & Getting started

GVEC is available on [PyPI](https://pypi.org/project/gvec/):
```bash
pip install gvec
```

For required libraries, other installation methods and more details see the documentation on [Installation](https://gvec.readthedocs.io/latest/user/install.html) and [Getting Started](https://gvec.readthedocs.io/latest/tutorials/getting-started.html).

## Reporting Bugs & Contributing to GVEC

The GVEC development repository is hosted at the [MPCDF-GitLab](https://gitlab.mpcdf.mpg.de/gvec-group/gvec).
There you can also find the [issue tracker](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/issues), [merge requests](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/merge_requests) and [release notes](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/releases).

A public mirror of GVEC is hosted on [GitHub](https://github.com/gvec-group/gvec).
There everyone can open [issues](https://github.com/gvec-group/gvec/issues) and [pull requests](https://github.com/gvec-group/gvec/pulls).

**Contributions are always welcome!**
More information on how to contribute can be found in the [developer documentation](https://gvec.readthedocs.io/latest/dev/index.html).

If you don't know where to start, don't hesitate to contact [the maintainers](https://www.ipp.mpg.de/person/139885/5497858) directly or join the [matrix-chat channel](https://matrix.to/#/#gvec:mpg.de) for questions and discussions.

## License

GVEC is released under the terms of the [MIT License](https://spdx.org/licenses/MIT.html).
For the full license terms see the included [LICENSE.txt](LICENSE.txt) file.

Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics

Parts of this software are licensed differently:
* `src/base/bsplines/` is part of [SeLaLib](https://github.com/selalib/selalib/) and licensed with `CECILL-B`.
* `src/mod_timings.f90` & `src/perf2timings.f90` are wrappers for the [ftimings](https://gitlab.mpcdf.mpg.de/loh/ftimings) library, licensed with `LGPL-3.0-only`.
* `src/globals/cla.f90` is [CLAF90](https://ingria.ceoas.oregonstate.edu/fossil/CLAF90) licensed with a modified `MIT` license.

## Citing GVEC

If you use GVEC in your work, please cite:
> Hindenlang et al., (2026). GVEC: A flexible 3D MHD equilibrium solver. Journal of Open Source Software, 11(120), 9670, [doi:10.21105/joss.09670](https://doi.org/10.21105/joss.09670)

A corresponding bibtex entry is found in [CITATION.bib](CITATION.bib).

To cite a specific version of GVEC, you can use the corresponding archive on [Zenodo](https://doi.org/10.5281/zenodo.15026780).

To reference the *G-Frame*, please cite:
> Hindenlang et al., (2025). Computing MHD equilibria of stellarators with a flexible coordinate frame. Plasma Physics and Controlled Fusion, 67(4), 045002, [doi:10.1088/1361-6587/adba11](https://doi.org/10.1088/1361-6587/adba11)

## References

This is a list of selected references in which the GVEC equilibrium solver was utilized:

| | |
|-|-|
| [HMB+26] | Florian Hindenlang, Omar Maj, Robert Babin, Robert Köberl, Dean Muir, Tiago Tamissa Ribeiro, Markus Rampp, and Eric Sonnendrücker. *GVEC: A flexible 3D MHD equilibrium solver*. Journal of Open Source Software, 11(120):9670, apr 2026. doi:10.21105/joss.09670.|
| [HPM+25] | Florian Hindenlang, Gabriel G Plunk, and Omar Maj. *Computing MHD equilibria of stellarators with a flexible coordinate frame*. Plasma Physics and Controlled Fusion, 67(4):045002, mar 2025. doi:10.1088/1361-6587/adba11.|
| [PDR+25] | Gabriel G Plunk, Michael Drevlak, Eduardo Rodríguez, Robert Babin, Alan Goodman, and Florian Hindenlang. *Back to the figure-8 stellarator*. Plasma Physics and Controlled Fusion, 67(3):035025, feb 2025. doi:10.1088/1361-6587/adb64b.|
| [PDS+23] | Jonas Puchmayr, Mike G Dunne, Erika Strumberger, Matthias Willensdorfer, Hartmut Zohm, and Florian Hindenlang. *Helical mode localization and mode locking of ideal MHD instabilities in magnetically perturbed tokamak plasmas*. Nuclear Fusion, 2023. |
| [MND+20] | Maurice Maurer, A Banon Navarro, Tilman Dannert, Marco Restelli, Florian Hindenlang, Tobias Goerler, Daniel Told, Denis Jarema, Gabriele Merlo, and Frank Jenko. *GENE-3D: a global gyrokinetic turbulence code for stellarators.* Journal of Computational Physics, 420:109694, 2020.|
| [NRH+22] | Nikita Nikulsin, Rohan Ramasamy, Matthias Hoelzl, Florian Hindenlang, Erika Strumberger, Karl Lackner, Sibylle Guenter, JOREK Team, and others. *JOREK3D: an extension of the JOREK nonlinear MHD code to stellarators*. Physics of Plasmas, 2022.|
| [NMP+20] | A Banon Navarro, G Merlo, G G Plunk, P Xanthopoulos, A Von Stechow, A Di Siena, M Maurer, F Hindenlang, F Wilms, and F Jenko. *Global gyrokinetic simulations of ITG turbulence in the magnetic configuration space of the Wendelstein 7-X stellarator*. Plasma Physics and Controlled Fusion, 62(10):105005, 2020.|
| [WNM+21] | Felix Wilms, Alejandro Bañón Navarro, Gabriele Merlo, Leonhard Leppin, Tobias Görler, Tilman Dannert, Florian Hindenlang, and Frank Jenko. *Global electromagnetic turbulence simulations of W7-X-like plasmas with GENE-3D*. Journal of Plasma Physics, 87(6):905870604, 2021. |
