---
myst:
  html_meta:
    "description lang=en": |
      Top-level documentation for GVEC, with links to the rest of the site..
html_theme.sidebar_secondary.remove: true
---

# GVEC

The Galerkin Variational Equilibrium Code.
An open-source software for the generation of three-dimensional ideal MHD equilibria.

::::{grid} 2
:::{grid-item-card}  Inspired by VMEC
Ideas are strongly based on [VMEC](https://princetonuniversity.github.io/STELLOPT/VMEC) (Hirshman & Whitson, 1983).
:::
:::{grid-item-card}  Radial B-Splines
Radial discretization using B-Splines of arbitrary polynomial degree. Fourier series in poloidal and toroidal direction with different maximum modenumber for each variable.
:::
:::{grid-item-card}  Flexible Mapping
Choice of the mapping $(X^1,X^2,\zeta) \mapsto (x,y,z)$, not restricted to $(R,Z,\phi)$.
:::
:::{grid-item-card}  Multiple Interfaces
Initialize with a VMEC netCDF output and convert to VMEC netCDF or interface with other specific converters.
:::
:::{grid-item-card}  Modern Fortran
Use of modern object-oriented Fortran
:::
:::{grid-item-card}  Python Postprocessing
*Python bindings for easy postprocessing and integration with other tools.* (in development)
:::
::::


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

[Fortran Code Documentation](/_static/ford/index.html){.external}

## Contact

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics 
in Garching, Germany.

The list of contributors is found in <project:dev/CONTRIBUTORS.md>.
