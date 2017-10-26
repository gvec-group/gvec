# GVEC INSTALLATION PROCEDURE


## Prerequisites

GVEC supports Linux-based systems only requires a x86_64
compliant platform and has been tested on the following platforms:

- Ubuntu 16.04 or newer


### Compilers

GVEC requires a C and a Fortran 2003 compliant compiler,
compilers tested with GVEC include

- GNU Compiler Collection 4.6 or newer
- Intel C/Fortran Compiler 12 or newer (recommended)

GVEC furthermore requires CMake 3.0+ as a build system.

### Libraries

The following libraries are required, if not mentioned
otherwise, including their development headers:

- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- Fortran netcdf library


## Compiling GVEC

GVEC supports CMake as a build system, which should be
available on most systems. The previously available
custom Makefile suport has been removed.
For compiling GVEC, create a new sub-directory,
e.g. "build" . Inside that directory execute
 
   CC=<C-Compiler> FC=<Fortran-Compiler>  ccmake ../

Here you can specify library paths and options. If no
preinstallied libraries for netcdf are found, an error occurs
Press <c> to configure and <g> to create the Makefiles.
Finally compile GVEC in the build folder by typing `make`.

### Libraries

Under ubuntu, the following packages should be installed:

- gcc
- g++
- gfortran
- liblapack-dev
- lib1g-dev
- libnetcdff-dev

### Visualization

For line plots, csv datafiles are generated. For plotting 2D and 3D data, we use [paraview](https://www.paraview.org).
