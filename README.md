### GVEC 

#### Overview

GVEC (Galerkin Variational Equilibrium Code) is an open-source software for
the generation of three-dimensional ideal MHD equilibria.
The ideas are strongly based on on the VMEC code, 
see [VMEC wiki pages](https://princetonuniversity.github.io/STELLOPT/VMEC).

The main features of GVEC are

* Use of modern **object-oriented FORTRAN**
* **Radial High Order Finite Element** discretization: Splines with continuity `C^(deg-1)` or discontinuous polynomials
* **Number of Fourier modes** for each variable `X^1,X^2,lambda` can be chosen separately
* **Flexible choice of the mapping** between the space `(X^1,X^2,zeta)--> (x,y,z)` (in VMEC fixed to `(R,Z,phi)-->(x,y,z)` ) 
  to find equilibria in complex-shaped domains (magnetic islands, knotted domain...)
* Also a VMEC generated netcdf outputfile can be used for initialization of GVEC.

GVEC has been developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
lead by Prof. Eric Sonnendruecker at the Max-Planck Institute for Plasma Physics 
in Garching, Germany.

The list of contributors is found in [CONTRIBUTORS.md](CONTRIBUTORS.md)

### Documentation

A pdf that documents the theory and implementation details of GVEC  **[ can be found here](https://gitlab.mpcdf.mpg.de/gvec-group/GVEC_doc/blob/master/GVEC_prototype/GVEC_prototype.pdf)**.

Direct code documentation can be found **[here.](http://gvec-group.pages.mpcdf.de/ford-gvec-doc)**
It is generated using [FORD](https://github.com/cmacmackin/ford) 
and the input file for running ford is in `ford-config/project.md`.

### License

GVEC is Copyright (C) 2017, F. Hindenlang, O. Maj, E. Sonnendruecker 

and is released under the terms of the GNU General Public License v3.0. 
For the full license terms see the included license file [LICENSE.md](LICENSE.md).

### Installation

For [installation instructions see INSTALL.md](INSTALL.md).


### Run a GVEC example

After [installation](INSTALL.md), the binary executable should be found in `build/bin/gvec`. 
There are several example input files named `paramter.ini`, which are found the folders `ini/w7x` `ini/jet_single_null` `ini/test` .
For execution, go into one of these folders and execute for example the following commands
``` 
  cd ini/w7x
  ../../build/bin/gvec parameter.ini |tee log
``` 
which pipes the screen output also into the file `log`.


### Visualization

For line plots, csv datafiles are generated (we like to use [veusz](https://veusz.github.io/) for plotting). 

For 3D data, we write `.vtu` files 
that can be visualized in [paraview](https://www.paraview.org).


### Object-Oriented Programming in FORTRAN

Here is a recommendation for a tutorial on how to program in an object-oriented way
with [polymorphism in fortran](https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd).

