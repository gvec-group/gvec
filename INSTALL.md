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

- cmake 
- git lfs
- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- netcdf library (Fortran & serial!)



#### Under Linux (Ubuntu, Feb.2024)

Under ubuntu, the following packages should be installed:

- `cmake` and `cmake-curses-gui`
- `git-lfs`
- `gcc`,`g++` and `gfortran`
- `liblapack3` and `liblapack-dev`
- `zlib1g-dev`
- `libnetcdf-dev` and `libnetcdff-dev`

#### Under Mac-OS (homebrew, Feb.2024)

Using homebrew (`brew install`) on a Mac, the following packages should be installed:

- `cmake` 
- `git-lfs`
- `gcc`
- `lapack`
- `netcdf` and `netcdf-fortran`

**Note**: `git lfs` has been recently introduced, if you have gvec checked out, be sure to run in the gvec folder 
```bash
  git lfs install
```
## Compiling GVEC

GVEC supports CMake as a build system, which should be
available on most systems. 

**IMPORTANT:**
Before executing cmake, be sure that you have all modules (netcdf must be compiled in serial)
and be sure to export an environment variable `FC` to point to the compiler. For gnu compiler for example:
```bash 
   export FC=`which gfortran`
```
or for Intel compiler
```bash
   export FC=`which ifort`
```

### Configure and build with "CMakePresets"

With Cmake version > 3.22, the CMakePresets feature can be used to configure and then build the code. 

Start from the GVEC directory with
```bash
  cmake --list-presets
```
to show a list of presets (defined `CMakePresets.json`).
 
Then a preset can be chosen and configured for a specified `build` directory (the build directory can have any name). 
```bash
  ### ON LINUX
  cmake --preset gvec_config_release -B build
  ### ON MAC (HOMEBREW)
  cmake --preset gvec_config_release_mac -B build
```
and then is compiled with
```bash
  cmake --build build
```

Further, the user can also create own presets by creating his own preset file `CMakeUserPresets.json` in the GVEC directory. The names must be different from those in `CMakePresets.json`. 
For example compiling on a mac in debug mode:
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
}
```
The user presets then appear also on the list of presets. 

**Note:** The preset files allow building the code in **VScode** with "CMake" and "CMake Tools" extensions.
### Compile with interactive cmake variables

To compile GVEC, create a new sub-directory that can have any name, e.g. `build` 
```bash 
   mkdir build ; cd build
```

Inside that directory execute
```bash
   ccmake ../
``` 

`ccmake` gives you a visual setup on the terminal. Press "enter" to change options.
Press "c" to configure and "g" to create the Makefiles. 
If `BUILD_NETCDF=ON` and no preinstallied libraries for netcdf are found, an error occurs...

In the main `CMakeList.txt` file, some pre-defined setups (library paths) for different architectures are controlled 
by setting the  `CMAKE_HOSTNAME` to `cobra/raven/draco/...` .

Finally compile GVEC in the build folder by typing (`-j` compiles in parallel)
```bash
  make -j
```

### Compiling on MPCDF cluster (raven, Feb. 2024)

Load the modules from the prepared shell script `CI_setup/MPCDF_setup_intel`

```bash
./CI_setup/raven_setup_intel
```
Follow the steps above and be sure to set in ccmake the `CMAKE_HOSTNAME` is set correctly (`raven`).

If you run gvec with the OpenMP parallelization, be sure to set the desired number of threads as an environment variable:
```bash
   export OMP_NUM_THREADS=???
```

## Running tests with ctests

You can run now several tests via ctest, that then calls the pytest environment, (requires `python >3.10` to be installed!). 

Change to the build directory, and execute:
```bash
ctest -T test --output-on-failure -R
```
