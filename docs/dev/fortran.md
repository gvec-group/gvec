# Fortran-only

Before the python bindings were introduced, GVEC was Fortran-only. Here, installation and use of the Fortran-only version is described, which might be helpful from a developers' perspective.

## Install Fortran executable of GVEC with `cmake`

```{warning}
This will only install the Fortran executables of GVEC, which are mostly used for testing. The python bindings are not installed with this option.
```

The standard way of compiling GVEC is using cmake presets, but there is also an interactive way with ccmake.

:::::{note}
Before executing cmake, be sure that you have all libraries (netcdf must be compiled in serial). It might also be necessary to export an environment variable `FC` to point to the compiler.

::::{tab-set}
:sync-group: compiler


:::{tab-item} GNU
:sync: gnu

```bash
export FC=`which gfortran`
```
:::

:::{tab-item} Intel ifx
:sync: intel

```bash
export FC=`which ifx`
```

:::
:::{tab-item} Intel ifort
:sync: intel

```bash
export FC=`which ifort`
```

:::
::::
:::::

::::::{tab-set}
:::::{tab-item} CMake Presets
<!-- ### Configure and build with cmake presets -->

With Cmake version > 3.22, the CMakePresets feature can be used to configure and then build the code.

1.  Start from the GVEC directory with
    ```bash
    cmake --list-presets
    ```
    to show a list of presets (defined `CMakePresets.json` and `CMakeUserPresets.json`).
1.  Select a preset and specify the `build` directory (the build directory can have any name).

    ::::{tab-set}
    :sync-group: os

    :::{tab-item} Linux
    :sync: linux

    ```bash
    cmake --preset gvec_config_release -B build
    ```

    :::
    :::{tab-item} MacOS
    :sync: mac

    ```bash
    cmake --preset gvec_config_release_mac_brew -B build
    ```

    :::
    ::::

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

:::{note}
The preset files allow building the code in **VScode** with "CMake" and "CMake Tools" extensions.
:::

:::::
:::::{tab-item} CCMake Interactive
<!-- ### Configure and build interactively -->

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

:::::
::::::

Now GVEC should be installed!
You should find the `gvec` and `gvec_post` binaries in `build/bin/` and can continue with [](/tutorials/getting-started).



## Run the GVEC Fortran executable
```{warning}
For users, we suggest using pygvec instead of the GVEC Fortran executable!
It has  has a simpler parameterfile (ending with `.ini`) and thus less features than using `pygvec` with TOML/YAML parameterfiles.
It is mainly used for testing purposes.
```

1) To install the Fortran executable of GVEC, follow the [installation instructions](/user/install).
2) The binary executables `gvec` and `gvec_post` should now be found in `build/bin/`.
3) GVEC is configured with a custom parameter file, typically called `parameter.ini`.
Example parameter files are found in `ini/` or `test-CI/examples/`

### Running GVEC

There are several test example input files named `parameter.ini`, which are found in a subfolder of [`test-CI/examples` {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/test-CI/examples/).

*   For execution, go into one of these folders and execute for example the following commands
    ```bash
    cd test-CI/examples/ellipstell_lowres
    ../../../build/bin/gvec parameter.ini |tee log
    # (|tee pipes the screen output also into the file `log`)
    ```
*   You can also restart a simulation by using one of the restart files (`*_State_*.dat`).
    Before the restart, resolution parameters in the `.ini` file can be changed, so that the new iterations will be on a finer grid, for example, or with more modes. The restart is triggered by simply adding the restart filename as an argument to the execution command, for example:
    ```bash
    ../../build/bin/gvec parameter.ini ELLIPSTELL_State_0000_00000200.dat |tee log
    ```
    Then the first integer (`_0000_`) will be incremented for the newly written restart files.

#### Run GVEC with OpenMP

If you run gvec with the OpenMP parallelization, be sure to set the desired number of threads as an environment variable:
   ```bash
   #replace ??? by the number of threads you want to use
   export OMP_NUM_THREADS=???
   ```

### Running tests

After compilation, you can quickly run some tests via `ctest`, that then calls the `pytest` environment of GVEC (requires `python >3.10` to be installed!).

Change to the build directory, and execute:
```bash
ctest -T test --output-on-failure -R
```

### Visualization

Using the python interface, any statefile can be loaded and visualized using the `ipython` notebook [`visu.ipynb`](<path:../../python/examples/visu.ipynb>) (view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/visu.ipynb)).

For line plots, csv datafiles are generated.

For 3D visualization data, it is possible to write `*visu*.vtu` files, that can be visualized in [paraview](https://www.paraview.org). There is an option to write visualization data in netcdf, `*visu*.nc`, which can be read for example in python.
