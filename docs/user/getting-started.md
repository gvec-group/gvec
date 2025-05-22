# Getting Started

## Run GVEC via its python bindings

First, please follow the installation instructions for installing [gvec with python bindings](install.md). Details on the python bindings are given [here](python.md).

We have prepared a circular tokamak example as a `ipython` notebook: [`run_and_visualize_gvec.ipynb`](<path:../../python/examples/gvecrun_tokamak/run_and_visualize_gvec.ipynb>).
All files for this example are also part of the repository at `python/examples/gvecrun_tokamak/` (view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/gvecrun_tokamak)).

Note that the kernel for the ipython notebook should be chosen as the virtual environment where the gvec python package is installed.

Here, we mention the main steps from the notebook to run gvec and post-process the result.
1.  Load the package with
    ```python
    os.environ["OMP_NUM_THREADS"]="4"
    import gvec
    ```
    Note that the number of openMP threads must be set before the import.

1.  To run gvec, a template [`parameter.ini`](<path:../../python/examples/gvecrun_tokamak/parameter.ini>) file is necessary.
    First, a sub-directory `runpath="run_01"` for the gvec run is created.

1.  From the template, one can modify its parameters using `gvec.util.adapt_parameter_file` and providing a dictionary with (key,value) pairs.
    Here, the template parameter file is modified and written to the run directory:
    ```python
    gvec.util.adapt_parameter_file(template, runpath / "parameter.ini", **params)
    ```
1.  To run the simulation, one changes to the run directory and executes `gvec.run`
    ```python
    with chdir(runpath):
        gvec.run("parameter.ini", stdout_path="stdout.txt")
    ```
1.  The final equilibrium solution is written to a state file, which can be loaded and evaluated using `gvec.State`
    ```python
    statefile = sorted(runpath.glob("*State*.dat"))[-1]
    with gvec.State(runpath / "parameter.ini", statefile) as state:
        rho = np.linspace(0, 1, 20)  # radial visualization points
        theta = np.linspace(0, 2 * np.pi, 50)  # poloidal visualization points
        ev = gvec.Evaluations(rho=rho, theta=theta, zeta=1, state=state)
        state.compute(ev, "X1", "X2", "LA","iota","p")
    ```
    Here, the visualization grid in the logical coordinates `rho,theta,zeta` has to be provided. The `ev` contains all computed variables, which are then plotted. A list of the available output variables is printed with
    ```python
    gvec.comp.table_of_quantities(markdown=True)
    ```

1. More visualization examples are provided in the `ipython` notebook [`visu.ipynb`](<path:../../python/examples/visu.ipynb>) (view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/visu.ipynb)).

## Run GVEC via the command line

1) To install GVEC, follow the [installation instructions](install).
2) The binary executables `gvec` and `gvec_post` should now be found in `build/bin/`.
3) GVEC is configured with a custom parameter file, typically called `parameter.ini`.
Example parameter files are found in `ini/` or `test-CI/examples/`

### Running GVEC

There are several example input files named `parameter.ini`, which are found in a subfolder of [`test-CI/examples` {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/test-CI/examples/).

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

For 3D visualization data, we write `*visu*.vtu` files, that can be visualized in [paraview](https://www.paraview.org). There is an option to write visualization data in netcdf, `*visu*.nc`, which can be read for example in python.

### Current constraint & increasing resolution

:::{warning}
This feature is experimental!
:::

The python bindings for GVEC provide a simple wrapper to run GVEC:
```bash
pygvec run parameter.ini
```

They also allow running GVEC with a *current constraint* (strictly speaking a current-optimization using picard iterations) and refinement.

For this the parameters need to be specified in YAML or TOML files, which are more flexible than the classic GVEC-INI files.
The parameters use the same keys, with a different syntax for specifying the boundary and axis coefficients.
In addition the parameters for `iota`, `pres` and `sgrid` are grouped together
and two new groups of parameters, `Itor` and `stages` are available.
If `Itor` and `stages` are not present in the parameterfile, a TOML or YAML parameterfile is equivalent to a `.ini` file.

The `stages` parameter is a list of stages, which are executed in order, for which each parameters can be selected that replace the default values.
The `Itor` parameter defines the target toroidal current profile (with the same syntax as the other profiles, but currently only with `type: polynomial`).
GVEC will now use picard iterations to optimize the input `iota` profile, such that the resulting current converges towards the target profile.
The `iota` parameter is now only used for the initial profile in the first run of the first stage.
The `runs` parameter within `stages` sets the number of picard iterations (GVEC restarts) to be done within each stage.

::::{tab-set}
:::{tab-item} TOML

```{code-block} toml
:caption: `parameter.toml`
# GVEC parameter file for W7X
ProjectName = "W7X"
whichInitEquilibrium = 0
minimize_tol = 1.0e-07

...

stages = [
    {runs = 10, maxIter = 10},
    {runs = 10, maxIter = 1000},
]

[Itor]
type = "polynomial"
scale = 1000
coefs = [0.0]

[iota]
type = "polynomial"
coefs = [0.9]

[X1_b_cos]
"(0, 0)" = 1.0
"(0, 1)" = 0.1

...
```

Full example: [`parameter.toml`](<path:../../python/examples/current_constraint/parameter.toml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_constraint/parameter.toml))

:::

:::{tab-item} YAML
```{code-block} yaml
:caption: `parameter.yaml`
# GVEC parameter file for W7X
ProjectName: W7X
whichInitEquilibrium: 0
minimize_tol: 1.0e-07

...

stages:
-   runs: 20
    maxIter: 10
-   runs: 10
    maxIter: 1000

Itor:
  type: polynomial
  scale: 1000
  coefs: [0.0]
iota:
  type: polynomial
  coefs: [0.9]

X1_b_cos:
  (0, 0): 1.0
  (0, 1): 0.1

...
```

Full example: [`parameter.yaml`](<path:../../python/examples/current_constraint/parameter.yaml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_constraint/parameter.yaml))

:::
::::
