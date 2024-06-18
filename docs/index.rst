.. GVEC documentation master file, created by
   sphinx-quickstart on Thu Jun 13 14:10:57 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GVEC: Galerkin Variational Equilibrium Code
===========================================

GVEC is an open-source software for the generation of three-dimensional ideal MHD equilibria.
The ideas are strongly based on `VMEC <https://princetonuniversity.github.io/STELLOPT/VMEC>`_ (Hirshman & Whitson, 1983).

The main features of GVEC are

* **Radial discretization using B-Splines** of arbitrary polynomial degree
* **Number of Fourier modes** for each variable :math:`X^1,X^2,\lambda` can be chosen separately
* **Flexible choice of the mapping** between the space :math:`\left(X^1,X^2,\zeta\right) \mapsto \left(x,y,z\right)` (in VMEC fixed to `(R,Z,phi)-->(x,y,z)` ) 
  to find equilibria in complex-shaped domains (magnetic islands, knotted domains...)
* Use of modern **object-oriented FORTRAN**
* **Python postprocessing** tools (in development)
* A VMEC generated netcdf file can be used for initialization of GVEC.
* GVEC can convert its final solution to several formats (VMEC netcdf, specific converters) 

Contact
-------

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics 
in Garching, Germany.

The list of contributors is found in :doc:`CONTRIBUTORS`.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 1

   Getting Started <README>
   Installation <INSTALL>
   dev
   Contributors <CONTRIBUTORS>
