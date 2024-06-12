project: GVEC
summary: GVEC (Galerkin Variational Equilibrium Code) is an open-source software for the generation of three-dimensional ideal MHD equilibria.
author: Florian Hindenlang et. al.
email: florian.hindenlang@ipp.mpg.de
project_gitlab: https://gitlab.mpcdf.mpg.de/gvec-group/gvec 
doc_license: by
src_dir: ../src/
include: ../src/
output_dir: ./_build
page_dir: ./static
mathjax_config: ./MathJax-latex-macros.js
predocmark: >
display: public
         protected
         private
source: true
graph: true
search: true
macro: TEST
       LOGIC=.true.

The main features of GVEC are

* Solves for the inverse coordinate mapping using the minimization of an energy functional, inspired by [VMEC](https://princetonuniversity.github.io/STELLOPT/VMEC) (Hirshman & Whitson, 1983)
* Use of modern **object-oriented FORTRAN**
* **Radial High Order Finite Element** discretization: Splines with continuity \(C^{deg-1}\) or discontinuous polynomials
* **Number of Fourier modes** for each variable \(X^1,X^2,\lambda\) can be chosen separately
* **Flexible choice of the mapping** between the space \(\left(X^1,X^2,\zeta\right) \rightarrow \left(x,y,z\right)\) (in VMEC fixed to \(\left(R,Z,\phi\right) \rightarrow \left(x,y,z\right)\)) 
  to find equilibria in complex-shaped domains (magnetic islands, knotted domain...)
* Also a VMEC generated netcdf outputfile can be used for initialization of GVEC.

GVEC is being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics 
in Garching, Germany.

## Getting Started
* [README](page/README.html)
* [Installation](page/INSTALL.html)
