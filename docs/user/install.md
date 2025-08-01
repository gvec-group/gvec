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

The python bindings require Python 3.10 or newer.

::::::{tab-set}
:sync-group: host

:::::{tab-item} General
:sync: general

Additionally GVEC requires:

- git
- cmake
- pkg-config
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
- `python3.10` (or newer)

:::::
:::::{tab-item} MacOS
:sync: mac

Install the following packages using homebrew (`brew install`)

- `cmake`
- `netcdf-fortran`
- `pkgconf`
- `gcc` (possibly no need to install explicitly)
- `lapack` (possibly no need to install explicitly)
- `python@3.10` (or newer)

:::::
::::::

## Install GVEC with `pip`

The latest release of GVEC is available on [PyPI](https://pypi.org/project/gvec/) and can be installed with:
```bash
pip install gvec
```

Currently GVEC is only available as an `sdist` (source-distribution) and not as a `wheel` (pre-compiled).
This means `pip` will download GVEC and compile it with CMake.
To see more details of the compilation process you can specify the `--verbose` flag.
For the compilation to work properly you need to have the prerequisite system packages installed/loaded.

Installing GVEC with `pip` will also install the python bindings and a number of scripts, see:
```bash
pygvec --help
```

### Installing pre-releases

You can also tell `pip` to install directly from the repository:
```bash
pip install git+https://gitlab.mpcdf.mpg.de/gvec-group/gvec.git
```

And even specify a different branch, using `@branch_name`, e.g.:
```bash
pip install git+https://gitlab.mpcdf.mpg.de/gvec-group/gvec.git@main
```

#### Troubleshooting

* `no such option --config-settings` during installation
  * upgrade pip with `pip install pip --upgrade`
* `Cannot open include file 'netcdf.inc'` during installation
  * `gvec/src/vmec/vmec_readin.f90(494): error #5102: Cannot open include file 'netcdf.inc'`
  * gvec cannot find netCDF
* `INTEL_MKL_ERROR` when trying to use pyGVEC
  * `INTEL MKL ERROR: /usr/lib/x86_64-linux-gnu/libmkl_avx2.so: undefined symbol: mkl_sparse_optimize_bsr_trsm_i8.`
  * `Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so or libmkl_def.so.`
  * you can try `export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libmkl_def.so:/usr/lib/x86_64-linux-gnu/libmkl_avx2.so:/usr/lib/x86_64-linux-gnu/libmkl_core.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_thread.so:/usr/lib/x86_64-linux-gnu/libiomp5.so`
* `undefined reference to EVP_KDF_CTX` when trying to `import gvec` and using `conda`
  * `/lib64/libk5crypto.so.3: undefined reference to EVP_KDF_CTX_new_id@OPENSSL_1_1_1b`
  * this can be caused by the `conda` environment conflicting with system libraries. You can try: `export LD_PRELOAD="/usr/lib64/libcrypto.so /usr/lib64/libssl.so"`
* on Fedora linux, it seems that the gfortran compiler needs an additional flag, in `CMakeList.txt`:
  ```
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/lib64/gfortran/modules")
  ```

## Getting GVEC with `git`

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
python3 -m venv .venv
source .venv/bin/activate
```

Then you can install the gvec python package manually with
```bash
pip install .[dev,examples] -v
```
