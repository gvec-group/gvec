# Installation

<!-- ::::{tab-set}
:sync-group: host

:::{tab-item} General
:sync: general

:::
:::{tab-item} MPCDF
:sync: mpcdf

:::
:::{tab-item} Ubuntu
:sync: ubuntu

:::
:::{tab-item} MacOS
:sync: mac

:::
:::: -->
<!-- ::::{tab-set}
:sync-group: os

:::{tab-item} Linux
:sync: linux

:::
:::{tab-item} MacOS
:sync: mac

:::
:::: -->
<!-- ::::{tab-set}
:sync-group: compiler

:::{tab-item} Intel
:sync: intel

:::
:::{tab-item} GNU
:sync: gnu

:::
:::: -->

:::{note}
These instructions have last been updated on 2025-02-04 and tested with ubuntu/raven/cobra/macOS.
:::

## Prerequisites

GVEC requires a C and a Fortran 2003 compliant compiler.
Compilers tested with GVEC include:

- GNU Compiler Collection 9 or newer
- Intel C/Fortran Compiler 17 or newer (recommended)
- CMake 3.5+ as a build system

::::::{tab-set}
:sync-group: host

:::::{tab-item} General
:sync: general

Additionally GVEC requires:

- git
- cmake
- libc6
- zlib
- BLAS/LAPACK (or compatible, e.g. ATLAS, MKL)
- netcdf library (Fortran & serial!)

:::::
:::::{tab-item} cluster
:sync: mpcdf

On MPCDF clusters (raven, viper) we can use the `module` system to manage dependencies.
There are prepared scripts in `CI_setup` that need to be sourced to load the modules, depending on the machine and the compiler.
For example on raven with the latest intel compiler, use
```bash
. CI_setup/raven_setup_ifx-2025
```
or on viper with the latest gfortran compiler, use
```bash
. CI_setup/viper_setup_gfortran-13
```

:::::
:::::{tab-item} Ubuntu
:sync: ubuntu

Install the following packages using `apt`:

- `cmake` and `cmake-curses-gui`
- `gcc`,`g++` and `gfortran`
- `liblapack3` and `liblapack-dev`
- `zlib1g-dev`
- `libnetcdf-dev` and `libnetcdff-dev`

:::::
:::::{tab-item} MacOS
:sync: mac

Install the following packages using homebrew (`brew install`)

- `cmake`
- `netcdf-fortran`
- `gcc` (possibly no need to install explicitly)
- `lapack` (possibly no need to install explicitly)

:::::
::::::

## Install GVEC with `pip`

The latest release of GVEC is available on [PyPI](https://pypi.org/project/gvec/) and can be installed with:
```bash
pip install gvec
```

Currently GVEC is only available as an `sdist` (source-distribution) and not as a `wheel` (pre-compiled). This means `pip` will download GVEC and compile it with CMake.

Installing GVEC with `pip` will also install the python bindings and a number of scripts, see:
```bash
pygvec --help
```

### Installing pre-releases

The latest version of the `develop` branch of GVEC is available in the GitLab package registry and can be installed with:
```bash
pip install gvec --index-url https://gitlab.mpcdf.mpg.de/api/v4/projects/1395/packages/pypi/simple
```

You can also tell `pip` to install directly from the repository:
```bash
pip install git+ssh://git@gitlab.mpcdf.mpg.de/gvec-group/gvec.git
```

And even specify a different branch, using `@branch_name`, e.g.:
```bash
pip install git+ssh://git@gitlab.mpcdf.mpg.de/gvec-group/gvec.git@main
```

#### Troubleshooting

* `Cannot open include file 'netcdf.inc'` during installation
  * `gvec/src/vmec/vmec_readin.f90(494): error #5102: Cannot open include file 'netcdf.inc'`
  * gvec likely did not recognize your hostname and cannot find netCDF. Possible fixes:
    * explicitly tell gvec which preset to use, e.g. `pip install git+ssh://git@gitlab.mpcdf.mpg.de/gvec-group/gvec.git --config-settings=cmake.define.CMAKE_HOSTNAME=tokp`
    * disable netCDF: `pip install git+ssh://git@gitlab.mpcdf.mpg.de/gvec-group/gvec.git --config-settings=cmake.define.LINK_GVEC_TO_NETCDF=Off`
* `INTEL_MKL_ERROR` when trying to use pyGVEC
  * `INTEL MKL ERROR: /usr/lib/x86_64-linux-gnu/libmkl_avx2.so: undefined symbol: mkl_sparse_optimize_bsr_trsm_i8.`
  * `Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so or libmkl_def.so.`
  * you can try `export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libmkl_def.so:/usr/lib/x86_64-linux-gnu/libmkl_avx2.so:/usr/lib/x86_64-linux-gnu/libmkl_core.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_thread.so:/usr/lib/x86_64-linux-gnu/libiomp5.so`
* `undefined reference to EVP_KDF_CTX` when trying to import pyGVEC and using `conda`
  * `/lib64/libk5crypto.so.3: undefined reference to EVP_KDF_CTX_new_id@OPENSSL_1_1_1b`
  * this can be caused by the `conda` environment conflicting with system libraries. You can try: `export LD_PRELOAD="/usr/lib64/libcrypto.so /usr/lib64/libssl.so"`
* on Fedora linux, it seems that the gfortran compiler needs an additional flag, in `CMakeList.txt`:
  ```
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/lib64/gfortran/modules")
  ```

## Getting GVEC with `git` & `CMake`

1.  Clone the repository to a local directory:
    ::::{tab-set}

    :::{tab-item} HTTPS

    Either you use the `https` address to clone, which always requests username and password when connecting to the mpcdf gitlab:
    ```bash
    git clone https://gitlab.mpcdf.mpg.de/gvec-group/gvec.git
    ```

    :::
    :::{tab-item} SSH

    Or you use the `git@` address to clone, in which case you have to copy your **public ssh key** of your machine to your profile on `gitlab.mpcdf.mpg.de` (typically its a `.pub` file in `~/.ssh`, generated by `ssh-keygen`). Then you do not have to specify any username or password.
    ```bash
    git clone git@gitlab.mpcdf.mpg.de:gvec-group/gvec.git
    ```

    :::
    ::::

1.  Now enter the created folder of the clone ( `cd gvec` )
1.  To change to the **latest release**, simply checkout `main` branch:
    ```bash
    git checkout main
    ```
    A list of branches on the remote is shown with `git branch -r`. You can change to a different branch of the GVEC repository (for a specific feature not yet available in the release) with:
    ```bash
    git checkout NAME_OF_BRANCH
    ```

### Installing the python bindings manually
You can install the gvec python package **manually**, from a cloned repository.
We strongly recommend to **always** use a clean virtual environment for the installation, e.g.
```bash
cd gvec
python -m venv .venv
source .venv/bin/activate
```

Then you can install the gvec python package manually with
```bash
pip install .[dev,examples] -v
```

## Install GVEC with `cmake`

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
You should find the `gvec` and `gvec_post` binaries in `build/bin/` and can continue with [](getting-started).
