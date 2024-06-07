project: GVEC &mdash; 3D MHD equilibrium code 
summary: Galerkin Variational Equilibrium Code
author: Florian Hindenlang et. al.
email: florian.hindenlang@ipp.mpg.de
website: https://gitlab.mpcdf.mpg.de/gvec-group/gvec 
license: by
src_dir: ../src/
output_dir: ./_build
exclude_dir: 
page_dir: ./static
mathjax_config: ./MathJax-latex-macros.js
predocmark: >
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.


[GVEC](https://gitlab.mpcdf.mpg.de/gvec-group/gvec) (Galerkin Variational Equilibrium Code) is an open-source software for
the generation of three-dimensional ideal MHD equilibria.
The ideas are strongly based on on the VMEC code.

The main features of GVEC are

* Use of modern **object-oriented FORTRAN**
* **Radial High Order Finite Element** discretization: Splines with continuity \(C^{deg-1}\) or discontinuous polynomials
* **Number of Fourier modes** for each variable \(X^1,X^2,\lambda\) can be chosen separately
* **Flexible choice of the mapping** between the space \(\left(X^1,X^2,\zeta\right) \rightarrow \left(x,y,z\right)\) (in VMEC fixed to \(\left(R,Z,\phi\right) \rightarrow \left(x,y,z\right)\)) 
  to find equilibria in complex-shaped domains (magnetic islands, knotted domain...)
* Also a VMEC generated netcdf outputfile can be used for initialization of GVEC.


GVEC has been developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
lead by Prof. Eric Sonnendruecker at the Max-Planck Institute for Plasma Physics 
in Garching, Germany.

