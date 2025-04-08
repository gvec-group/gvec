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

::::{grid} 2
:::{grid-item-card}  Inspired by VMEC
Ideas are strongly based on [VMEC](https://princetonuniversity.github.io/STELLOPT/VMEC) (Hirshman & Whitson, 1983).
:::
:::{grid-item-card}  Python bindings
Installable with `pip`. Python bindings for running, postprocessing and integration with other tools.
:::{grid-item-card}  Radial B-Splines
Radial discretization using B-Splines of arbitrary polynomial degree. Fourier series in poloidal and toroidal direction with different maximum modenumber for each variable.
:::
:::{grid-item-card}  Multiple Interfaces
Initialize with a VMEC netCDF output and convert to VMEC netCDF or interface with other specific converters.
:::
:::{grid-item-card}  Flexible Mapping
Choice of the mapping $(X^1,X^2,\zeta) \mapsto (x,y,z)$, not restricted to $(R,Z,\phi)$, but e.g. a generalized Frenet frame.
:::
:::{grid-item-card}  Modern Fortran
Use of modern object-oriented Fortran
:::
::::

```{figure} static/frenet_n2-12_bfield.png
:width: 70 %
:align: center

The magnetic field of a two-fieldperiod QI-stellarator configuration (configuration taken from [[HPM25]](https://doi.org/10.1088/1361-6587/adba11)).
```

## User Guide

```{toctree}
:maxdepth: 2
:titlesonly:

user/index
```

## Developer Guide

```{toctree}
:maxdepth: 2
:titlesonly:

dev/index
```

## Fortran API

Automatic Fortran code documentation generated with [ford](https://forddocs.readthedocs.io).

[Fortran Code Documentation](ford/index.html){.external}

## Contact

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany.

The list of contributors is found in <project:dev/CONTRIBUTORS.md>.
