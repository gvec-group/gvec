# GVEC 

## Overview

GVEC (Galerkin Variational Equilibrium Code) is an open-source software for
the generation of three-dimensional ideal MHD equilibria.
The ideas are strongly based on on the VMEC code, 
see [VMEC wiki pages](https://princetonuniversity.github.io/STELLOPT/VMEC).

The main features of GVEC are

* Use of modern **object-oriented FORTRAN**
* **Radial High Order B-spline** discretization: Splines with continuity `C^(deg-1)` or discontinuous polynomials
* **Number of Fourier modes** for each variable `X^1,X^2,lambda` can be chosen separately
* **Flexible choice of the mapping** between the space `(X^1,X^2,zeta)--> (x,y,z)` (in VMEC fixed to `(R,Z,phi)-->(x,y,z)` ) 
  to find equilibria in complex-shaped domains (magnetic islands, knotted domains...)
* A VMEC generated netcdf file can be used for initialization of GVEC.
* GVEC can convert its final solution to several formats (VMEC netcdf, specific converters) 

GVEC has been developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
lead by Prof. Eric Sonnendruecker at the Max-Planck Institute for Plasma Physics 
in Garching, Germany.

The list of contributors is found in [CONTRIBUTORS.md](CONTRIBUTORS.md)

## Installation

After cloning this repository, follow the installation process: [For instructions see INSTALL.md](INSTALL.md).

## Documentation

A pdf that documents the theory and implementation details of GVEC  **[ can be found here](https://gitlab.mpcdf.mpg.de/gvec-group/GVEC_doc/blob/master/GVEC_prototype/GVEC_prototype.pdf)**.

Direct code documentation can be found **[here.](http://gvec-group.pages.mpcdf.de/ford-gvec-doc)**
It is generated using [FORD](https://forddocs.readthedocs.io/en/latest/) 
and the input file for running ford is in `ford-config/project.md`.

### License

GVEC is Copyright (C) 2017, F. Hindenlang, O. Maj, E. Sonnendruecker 

and is released under the terms of the GNU General Public License v3.0. 
For the full license terms see the included license file [LICENSE.md](LICENSE.md).

### Testing

GVEC has been equipped with automatic CI testing on the gitlab.mpcdf.mpg.de, using shared MPCDF the gitlab runners to execute the tests. 
More details on the CI setup are found in [test-CI/README-CI.md](test-CI/README-CI.md).

The CI controls builds of the code, then calls pytest for running it and checking the results (requires `python >3.10` to be installed!).
The pytest feature allows to locally run the same tests. More details and examples on running the tests with pytest are found in [test-CI/README.md](test-CI/README.md).

Some pytest commands can also be executed directly after the [cmake install process](INSTALL.md). Simply change to the build directory, and execute:
```bash
ctest -T test --output-on-failure -R
```

### Run a GVEC example

After [installation](INSTALL.md), the binary executable should be found in `build/bin/gvec`. 
There are several example input files named `paramter.ini`, which are found in a subfolder of `test-CI/examples`.
For execution, go into one of these folders and execute for example the following commands
```bash
  cd test-CI/examples/ellipstell_lowres
  ../../../build/bin/gvec parameter.ini |tee log
``` 
(pipes the screen output also into the file `log`)

You can also restart a simulation by using on of the restart files (`*_State_*.dat`). 
Before the restart, resolution parameters in the `.ini` file can be changed, so that the new iterations will be on a finer grid, for example, or with more modes. The restart is triggered by simply adding the restart filename as an argument to the execution command, for example:
```bash 
  ../../build/bin/gvec parameter.ini TOKSY_State_0000_00000200.dat |tee log
``` 
Then the first integer (`_0000_`) will be incremented for the newly written restart files. 


### Visualization

For line plots, csv datafiles are generated (we like to use [veusz](https://veusz.github.io/) for plotting). 

For 3D data, we write `.vtu` files 
that can be visualized in [paraview](https://www.paraview.org).


### Object-Oriented Programming in FORTRAN

Here is a recommendation for a tutorial on how to program in an object-oriented way
with [polymorphism in fortran](https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd).

