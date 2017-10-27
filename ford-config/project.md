project: GVEC &mdash; 3D MHD equilibrium code 
summary: Galerkin Variational Equilibrium Code
author: Omar Maj, Florian Hindenlang
email: hindenlang@gmail.com
website: https://gitlab.mpcdf.mpg.de/ipphinde/gvec 
bitbucket: https://gitlab.mpcdf.mpg.de/ipphinde/gvec  
project_bitbucket: 
license: by
src_dir: ./../src/
output_dir: ./../doc
exclude_dir: ./../doc
             ./../doc/src
mathjax_config: ./../ford-config/MathJax-latex-macros.js
predocmark: >
display: public
         protected
         private
source: false
graph: false
search: false
macro: TEST
       LOGIC=.true.


GVEC (Galerkin Variational Equilibrium Code) is an open-source software for
the generation of three-dimensional ideal MHD equilibria.
The ideas are strongly based on on the VMEC code,
see [VMEC wiki pages](http://vmecwiki.pppl.wikispaces.net/VMEC).
Also a VMEC generated netcdf outputfile can be used for initialization of GVEC.

GVEC reuses parts of the open-source code [HOPR](https://github.com/fhindenlang/hopr).

GVEC has been developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
lead by Prof. Eric Sonnendruecker at the Max-Planck Institute for Plasma Physics
in Garching, Germany.

