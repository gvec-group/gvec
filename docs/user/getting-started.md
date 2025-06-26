# Getting Started

## Run GVEC via its python bindings

First, please follow the installation instructions for installing [gvec with python bindings](install.md). Details on the python bindings are given [here](python.md).

We have prepared a simple elliptic tokamak example as a `ipython` notebook: [`tokamak_gvecrun_and_visualize.ipynb`](<path:../../python/examples/tokamak_gvecrun/tokamak_gvecrun_and_visualize.ipynb>).
All files for this example are also part of the repository at `python/examples/tokamak_gvecrun/` (view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/tokamak_gvecrun)).

Note that the kernel for the ipython notebook should be chosen as the virtual environment where the gvec python package is installed.

Here, we mention the main steps from the notebook to run gvec and post-process the result.
1.  Load the package with
    ```python
    os.environ["OMP_NUM_THREADS"]="2"
    import gvec
    ```
    Note that the number of openMP threads must be set before the import.

1.  We define the parameters in a dictionary
    ```python
    params = {
        "ProjectName": "GVEC_elliptok",
        ...
    }
    ```
1.  To run the simulation, we use `gvec.run` and specify a subdirectory for storing the results.
    ```python
    runpath = Path(f"run_{1:02d}")
    run = gvec.run(params, runpath=runpath)
    ```
1.  The final equilibrium solution is kept in the `run` object as a `run.state`, which can be loaded and evaluated using `gvec.State`
    ```python
    state = run.state
    rho = np.linspace(0, 1, 20)
    theta = np.linspace(0, 2 * np.pi, 50)
    ev = state.evaluate("X1", "X2", "LA", "iota", "p", rho=rho, theta=theta, zeta=[0.0])
    ```
    Here, the visualization grid in the logical coordinates `rho,theta,zeta` has to be provided. The `ev` contains all computed variables as an `xarray` dataset, which are then plotted. A list of the available output variables is printed with
    ```python
    gvec.table_of_quantities(markdown=True)
    ```

Instead of running GVEC within the notebook, one can also export the parameter dictionary to a TOML parameterfile
```python
gvec.util.write_parameters(params, "parameter.toml")
```
and then run GVEC from the command line:

```bash
pygvec run parameter.toml
```

1. More visualization examples are provided in the `ipython` notebook [`visu.ipynb`](<path:../../python/examples/visu.ipynb>) (view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/visu.ipynb)).


## Run GVEC with prescribed toroidal current

To run GVEC with a prescribed toroidal current profile one has to utilize the python-bindings and switch to the more flexible TOML or YAML style parameter files. The additional parameters that need to be set are `I_tor` and `picard_current`. To run GVEC with such a TOML/YAML input file, one can just calls from the command line:
```bash
pygvec run parameter.toml
```
Example parameter files are given below.
Using `picard_current = "auto"` selects the default algorithm for the current optimization with the specified `totalIter` and `minimize_tol`.
For more information and more detailed control over the current optimization (including the possibility to do resolution refinement) see the [stages](stages.md) section.

::::{tab-set}
:::{tab-item} TOML

```{code-block} toml
:caption: `parameter.toml`
# GVEC parameter file for W7X
ProjectName = "W7X"
minimize_tol = 1.0e-06

...

totalIter = 5000
picard_current = "auto"

[I_tor]
type = "polynomial"
coefs = [0.0]

[X1_b_cos]
"(0, 0)" = 5.5
"(0, 1)" = 0.2354

...
```

Full example: [`parameter.toml`](<path:../../python/examples/current_profile/parameter.toml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_profile/parameter.toml))

:::

:::{tab-item} YAML
```{code-block} yaml
:caption: `parameter.yaml`
# GVEC parameter file for W7X
ProjectName: W7X
minimize_tol: 1.0e-06

...

totalIter: 5000
picard_current: auto

I_tor:
  type: polynomial
  coefs: [0.0]

X1_b_cos:
  (0, 0): 5.5
  (0, 1): 0.2354

...
```

Full example: [`parameter.yaml`](<path:../../python/examples/current_profile/parameter.yaml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_profile/parameter.yaml))

:::
::::



## Run the GVEC Fortran executable
```{warning}
For users, we suggest using pygvec instead of the GVEC Fortran executable!
It has  has a simpler parameterfile (ending with `.ini`) and thus less features than using `pygvec` with TOML/YAML parameterfiles.
It is mainly used for testing purposes.
```

1) To install the Fortran executable of GVEC, follow the [installation instructions](install).
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
