# GVEC INSTALLATION PROCEDURE


## Prerequisites

GVEC supports Linux-based systems only requires a x86\_64
compliant platform and has been tested on the following platforms:

- Ubuntu 16.04 or newer


### Compilers

GVEC requires a C and a Fortran 2003 compliant compiler,
compilers tested with GVEC include

- GNU Compiler Collection 4.6 or newer
- Intel C/Fortran Compiler 12 or newer (recommended)
- CMake 3.0+ as a build system

### Libraries

The following libraries are required, if not mentioned
otherwise, including their development headers:

- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- Fortran netcdf library


Under ubuntu, the following packages should be installed:

- `cmake` and `cmake-curses-gui`
- `gcc`,`g++` and `gfortran`
- `liblapack3` and `liblapack-dev`
- `zlib1g-dev`
- `libnetcdf-dev` and `libnetcdff-dev` (for VMEC netcdf datafile readin only)

## Compiling GVEC

GVEC supports CMake as a build system, which should be
available on most systems. The previously available
custom Makefile suport has been removed.
For compiling GVEC, create a new sub-directory,
e.g. `build` 
``` 
   mkdir build ; cd build
```
Inside that directory execute
``` 
   ccmake ../
``` 
Here you can specify library paths and options. Press "enter" to change options.
Press "c" to configure and "g" to create the Makefiles.
If `BUILD_CDF=ON` and no preinstallied libraries for netcdf are found, an error occurs...

Finally compile GVEC in the build folder by typing 
```
  make
```
   
 




