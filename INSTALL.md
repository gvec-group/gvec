# GVEC INSTALLATION PROCEDURE

Quick-links:
*  [Install on a MPCDF cluster](#install-on-a-mpcdf-cluster) 
*  [Install on a Linux](#install-on-linux) 
*  [Install on a Mac](#install-on-a-mac) 
*  [Compile GVEC](#compile-gvec)

## Prerequisites
The following libraries are required, if not mentioned
otherwise, including their development headers:
- cmake 
- git lfs
- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- netcdf library (Fortran & serial!)

**Note**: `git lfs` has been recently introduced, once you have gvec cloned and the branch checked out, be sure to run in the gvec directory 
```bash
git lfs install
git lfs pull
```


## Install on a MPCDF cluster

*tested on raven and cobra updated 2024/02/26*

GVEC has been installed and tested on the MPCDF cobra and raven clusters

### Compiler and Libraries

To setup the compiler and load the modules, we use a prepared shell script in `CI_setup`, that depends on the compiler:
*   For Intel compiler on raven or cobra, run for example:
    ```bash
    . CI_setup/mpcdf_setup_intel
    ```
*   For Gnu compiler on raven, run for example:
    ```bash
    . CI_setup/mpcdf_setup_gnu
    ```

**Note**: `git lfs` has been recently introduced, once you have gvec cloned and the branch checked out, be sure to run in the gvec directory 
```bash
git lfs install
git lfs pull
```

### Compilation

...continue with section [Compile GVEC](#compile-gvec)

## Install on a Mac

*tested on a mac with homebrew, updated 2024/02/26*

### Compilers

GVEC requires a C and a Fortran 2003 compliant compiler,
compilers tested with GVEC include

- GNU Compiler Collection 9 or newer
- Intel C/Fortran Compiler 17 or newer (recommended)
- CMake 3.5+ as a build system

### Libraries

Using homebrew (`brew install`) on a Mac, the following packages should be installed:
- `cmake` 
- `git-lfs`
- `netcdf-fortran`
- `gcc` (possibly no need to install explicitly)
- `lapack` (possibly no need to install explicitly)

**Note**: `git lfs` has been recently introduced, once you have gvec cloned and the branch checked out, be sure to run in the gvec directory 
```bash
git lfs install
git lfs pull
```

### Compilation

GVEC supports CMake as a build system, which should be available on most systems. 

...continue with section [Compile GVEC](#compile-gvec)

## Install on linux

*tested on ubuntu, updated 2024/02/26*

### Compilers

GVEC requires a C and a Fortran 2003 compliant compiler,
compilers tested with GVEC include

- GNU Compiler Collection 9 or newer
- Intel C/Fortran Compiler 17 or newer (recommended)
- CMake 3.5+ as a build system

### Libraries

Under ubuntu, the following packages should be installed:

- `cmake` and `cmake-curses-gui`
- `git-lfs`
- `gcc`,`g++` and `gfortran`
- `liblapack3` and `liblapack-dev`
- `zlib1g-dev`
- `libnetcdf-dev` and `libnetcdff-dev`

**Note**: `git lfs` has been recently introduced, once you have gvec cloned and the branch checked out, be sure to run in the gvec directory 
```bash
git lfs install
git lfs pull
```

### Compilation

**IMPORTANT:** Before executing cmake, be sure that you have all libraries (netcdf must be compiled in serial). It might also be necessary to export an environment variable `FC` to point to the compiler. 
*  For gnu compiler:
   ```bash 
   export FC=`which gfortran`
   ```
*  For Intel compiler
   ```bash
   export FC=`which ifort`
   ```
...continue with section [Compile GVEC](#compile-gvec)

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


## Compile GVEC

* The standard way of compiling GVEC is using cmake presets, see [here](#configure-and build-with-cmake-presets)
* There is also an interactive way with ccmake , see [here](#configure-and-build-interactively)

### Configure and build with cmake presets

With Cmake version > 3.22, the CMakePresets feature can be used to configure and then build the code. 

1.  Start from the GVEC directory with
    ```bash
    cmake --list-presets
    ```
    to show a list of presets (defined `CMakePresets.json` and `CMakeUserPresets.json`).
1.  Select a preset and specify the `build` directory (the build directory can have any name). 
    *   **on Linux**: 
        ```bash
        cmake --preset gvec_config_release -B build
        ```
    *   **on mac (homebrew)**
        ```bash
        cmake --preset gvec_config_release_mac_brew -B build
        ```
1.  Then compile with  (`-j` compiles in parallel)
    ```bash
    cmake --build build -j
    ```

Further, the user can also create own presets by creating his own preset file `CMakeUserPresets.json` in the GVEC directory. Be careful to only add new entries with new names, as they must be different from those in `CMakePresets.json`. For example compiling on a mac in debug mode:
```json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 22,
    "patch": 0
    },
  "configurePresets": [
      {
          "name": "gvec_config_debug_mac",
          "displayName": "GVEC configure: default debug build on a MAC",
          "hidden": false,
          "cacheVariables": {
              "CMAKE_BUILD_TYPE": "Debug",
              "COMPILE_GVEC": "On",
              "CMAKE_HOSTNAME": "mac_brew",
              "LINK_GVEC_TO_NETCDF": "On",
              "USE_OPENMP": "On",
              "COMPILE_GVEC_AS_STATIC_LIB": "On"
          }
      }
  ]
}
```
The user presets then appear also on the list of presets. 

**Note:** The preset files allow building the code in **VScode** with "CMake" and "CMake Tools" extensions.

### Configure and build interactively

To compile GVEC interactively (needs `ccmake` command):

1.  create a new subdirectory that can have any name, e.g. `build` 
    ```bash 
    mkdir build ; cd build
    ```
1.  Inside that directory execute
    ```bash
    ccmake ../
    ``` 
    `ccmake` gives you a visual setup on the terminal. 
    *  Press "enter" to change options, and press "enter" again to fix the change
    *  Press "c" to configure and "g" to create the Makefiles. 
    *  If `BUILD_NETCDF=ON` and no preinstalled libraries for netcdf are found, an error occurs...
    * On a Mac, be sure to activate `COMPILE_GVEC_AS_STATIC_LIB=ON` (in ccmake, toggle to view all variables by typing `t`)
    *  In the main `CMakeList.txt` file, some pre-defined setups (library paths) for different architectures are controlled 
       by setting the  `CMAKE_HOSTNAME` to `cobra`/`raven`/`mac_brew`/`mac_ports`/`tokp`/.. .
1.  Finally, compile GVEC in the build directory by typing (`-j` compiles in parallel)
    ```bash
    make -j
    ```

## Run GVEC with OpenMP

If you run gvec with the OpenMP parallelization, be sure to set the desired number of threads as an environment variable:
   ```bash
   #replace ??? by the number of threads you want to use
   export OMP_NUM_THREADS=???
   ```

## Running tests

After compilation, you can quickly run some tests via `ctest`, that then calls the `pytest` environment of GVEC (requires `python >3.10` to be installed!). 

Change to the build directory, and execute:
```bash
ctest -T test --output-on-failure -R
```
