project: GVEC &mdash; 3D MHD equilibrium code 
summary: Galerkin Variational Equilibrium Code
author: Omar Maj, Florian Hindenlang
email: hindenlang@gmail.com
website: https://gitlab.mpcdf.mpg.de/gvec-group/gvec 
git: https://gitlab.mpcdf.mpg.de/gvec-group/gvec  
license: by
src_dir: ./../src/
output_dir: ./../../ford-gvec-doc
exclude_dir: ./../../ford-gvec-doc
             ./../../ford-gvec-doc/src
             ./../src/lbfgsb
mathjax_config: ./../ford-config/MathJax-latex-macros.js
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
The ideas are strongly based on on the VMEC code, 
see [VMEC wiki pages](https://bitbucket.org/lazerson_princeton/stellopt/wiki/VMEC).

The main features of GVEC are

* Use of modern **object-oriented FORTRAN**
* **Radial High Order Finite Element** discretization: Splines with continuity `C^(deg-1)` or discontinuous polynomials
* **Number of Fourier modes** for each variable `X^1,X^2,lambda` can be chosen separately
* **Flexible choice of the mapping** between the space `(X^1,X^2,zeta)--> (x,y,z)` (in VMEC fixed to `(R,Z,phi)-->(x,y,z)` ) 
  to find equilibria in complex-shaped domains (magnetic islands, knotted domain...)
* Also a VMEC generated netcdf outputfile can be used for initialization of GVEC.

Some parts of GVEC are reused from the open-source code [HOPR](https://github.com/fhindenlang/hopr).

GVEC has been developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
lead by Prof. Eric Sonnendruecker at the Max-Planck Institute for Plasma Physics 
in Garching, Germany.

