---
myst:
  html_meta:
    "description lang=en": |
      Top-level documentation for GVEC, with links to the rest of the site..
html_theme.primary_secondary.remove: true
html_theme.sidebar_secondary.remove: true
---


<div style="text-align: center;">

# The Galerkin Variational Equilibrium Code
## A flexible 3D MHD equilibrium solver
</div>

[GVEC](https://gitlab.mpcdf.mpg.de/gvec-group/gvec) is an open-source software for the generation of three-dimensional ideal MHD equilibria.

```{grid} 2
:gutter: 2
:class-container: gallery-directive

:::{grid-item-card} Inspired by VMEC
Ideas are strongly based on [VMEC](https://princetonuniversity.github.io/STELLOPT/VMEC) (Hirshman & Whitson, 1983).
:::
:::{grid-item-card} Python bindings
Installable with `pip`. Python bindings for running, postprocessing and integration with other tools.
:::
:::{grid-item-card} Radial B-Splines
Radial discretization using B-Splines of arbitrary polynomial degree. Fourier series in poloidal and toroidal direction with different maximum modenumber for each variable.
:::
:::{grid-item-card} Multiple Interfaces
Initialize with a VMEC netCDF output or interface with other codes: JOREK, CASTOR3D, GENE...
:::
:::{grid-item-card} Flexible Mapping
Choice of the mapping $(X^1,X^2,\zeta) \mapsto (x,y,z)$, not restricted to $(R,Z,\phi)$, but e.g. a [generalized Frenet frame](#g-frame).
:::
:::{grid-item-card} Modern Fortran
Use of modern object-oriented Fortran
:::
```

```{figure} static/frenet_n2-12_bfield.png
:width: 70 %
:align: center

The magnetic field of a two-fieldperiod QI-stellarator configuration (configuration taken from [Hindenlang et al. (2025)](https://doi.org/10.1088/1361-6587/adba11)).
```

```{grid} 2 2 2 4
:gutter: 2
:class-container: gallery-directive

:::{grid-item-card} Installation
:link: user/install
:link-type: doc
:class-card: intro-card
:img-top: static/fa-wrench-solid-full.svg
:::
:::{grid-item-card} Getting Started
:link: tutorials/getting-started
:link-type: doc
:class-card: intro-card
:img-top: static/fa-person-running-solid-full.svg
:::
:::{grid-item-card} User Guide
:link: user/index
:link-type: doc
:class-card: intro-card
:img-top: static/fa-book-open-solid-full.svg
:::
:::{grid-item-card} Developer Guide
:link: dev/index
:link-type: doc
:class-card: intro-card
:img-top: static/fa-terminal-solid-full.svg
:::
```
## Citing GVEC

If you use GVEC in your work, please cite:
> Hindenlang et al., (2026). GVEC: A flexible 3D MHD equilibrium solver. Journal of Open Source Software, 11(120), 9670, [doi:10.21105/joss.09670](https://doi.org/10.21105/joss.09670)

To cite a specific version of GVEC, you can use the corresponding archive on [Zenodo](https://doi.org/10.5281/zenodo.15026780).

## Statement of need

MHD equilibrium solutions are the basis for a number of high fidelity plasma physics models and associated codes.
For example, they provide the initial conditions for linear and nonlinear MHD solvers (e.g. CASTOR3D ([Puchmayr et al., 2023](https://doi.org/10.1088/1741-4326/acdd12)), CAS3D ([Schwab, 1993](https://doi.org/10.1063/1.860656)), Jorek3D ([Nikulsin et al., 2022](https://doi.org/10.1063/5.0087104)), Struphy ([Holderied et al., 2022](https://doi.org/10.1016/j.jcp.2021.110143)), M3D-C1 ([Jardin, 2004](https://doi.org/10.1016/j.jcp.2004.04.004))), or the magnetic field for particle orbit tracing (e.g. SIMPLE ([Albert et al., 2020](https://doi.org/10.1016/j.jcp.2019.109065))) and turbulence simulations (e.g. BOUT++ ([Shanahan et al., 2024]()), GENE ([Bañón Navarro et al., 2020](https://doi.org/10.1088/1361-6587/aba858))).
3D MHD equilibria are directly used to analyse potential stellarator configurations in optimisation frameworks, such as SIMSOPT ([Landreman et al., 2021](https://doi.org/10.21105/joss.03525)) or STELLOPT ([Lazerson et al., 2020](https://doi.org/10.11578/dc.20180627.6)).

GVEC has a flexible coordinate frame, allowing it to represent boundary shapes beyond those possible with the standard cylindrical coordinates used by many equilibrium codes. This can be used for example to optimise the boundary shape of the figure-8 stellarator in [Plunk et al. (2025)](https://doi.org/10.1088/1361-6587/adb64b).

## Contact

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany.

The list of contributors is found in <project:dev/CONTRIBUTORS.md>.


```{toctree}
:maxdepth: 2
:hidden:

User Guide <user/index>
Developer Guide <dev/index>
```
