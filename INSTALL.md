# GVEC INSTALLATION PROCEDURE


## Prerequisites

GVEC supports Linux-based systems only requires a x86\_64
compliant platform and has been tested on the following platforms:

- Ubuntu 16.04 or newer
- MPCDF cobra


### Compilers

GVEC requires a C and a Fortran 2003 compliant compiler,
compilers tested with GVEC include

- GNU Compiler Collection 5.4 or newer
- Intel C/Fortran Compiler 17 or newer (recommended)
- CMake 3.5+ as a build system

### Libraries

The following libraries are required, if not mentioned
otherwise, including their development headers:

- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- netcdf library (Fortran & serial!)



#### Under Linux (Ubuntu)

Under ubuntu, the following packages should be installed:

- `cmake` and `cmake-curses-gui`
- `gcc`,`g++` and `gfortran`
- `liblapack3` and `liblapack-dev`
- `zlib1g-dev`
- `libnetcdf-dev` and `libnetcdff-dev` (for VMEC netcdf datafile readin only)

## Compiling GVEC

GVEC supports CMake as a build system, which should be
available on most systems. 

**IMPORTANT:**
Before executing cmake, be sure that you have all modules (netcdf must be compiled in serial)
and be sure to export an environment variable `FC` to point to the compiler. For intel for example:
``` 
   export FC=`which ifort`
```
or for GNU compiler
``` 
   export FC=`which gfortran`
```

To compile GVEC, create a new sub-directory that can have any name, e.g. `build` 
``` 
   mkdir build ; cd build
```

Inside that directory execute
``` 
   ccmake ../
``` 
Here you can specify options. Press "enter" to change options.
Press "c" to configure and "g" to create the Makefiles. 
If `BUILD_NETCDF=ON` and no preinstallied libraries for netcdf are found, an error occurs...

In the main `CMakeList.txt` file, some pre-defined setups (library paths) for different architectures are controlled 
by setting the  `CMAKE_HOSTNAME` to `cobra/hydra/...` .

Finally compile GVEC in the build folder by typing (`-j` compiles in parallel)
```
  make -j
```
   
#### Compiling on Cobra cluster (Dec. 2019)

Load the modules and export the fortran compiler : 

```
   module purge 
   module load git cmake intel mkl hdf5-serial
   module load netcdf-serial
   export FC=`which ifort`
```

Follow the steps above and be sure to set in ccmake the `CMAKE_HOSTNAME` to `cobra`.

If you want to run gvec with the OpenMP parallelization, be sure to set the desired number of threads:
```
   export OMP_NUM_THREADS=??
```

#### Compiling on Cobra cluster with MPI (Oct 2022)

Load the modules and export the fortran compiler : 

```
    module purge
    module load git cmake 
    module load intel/19.1.3 mkl hdf5-serial
    module load impi/2019.9
    export FC=`which mpiifort`
```
Example for an interactive run:
```
    srun  --nodes=1 --ntasks-per-node=#MPItasks --cpus-per-task=$OMP_NUM_THREADS #cores -p interactive --time=TIME_LESS_THAN_2HOURS --mem=MEMORY_LESS_THAN_32G build/bin/gvec ini/performance/parameter_tiagoMPI.ini
```

