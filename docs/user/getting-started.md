# Getting Started

1) To install GVEC, follow the [installation instructions](install).
2) The binary executables `gvec` and `gvec_post` should now be found in `build/bin/`
3) GVEC is configured with a custom parameter file, typically called `parameter.ini`.
Example parameter files are found in `ini/` or `test-CI/examples/`

## Running GVEC

There are several example input files named `parameter.ini`, which are found in a subfolder of `test-CI/examples`.

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

## Visualization

For line plots, csv datafiles are generated (we like to use [veusz](https://veusz.github.io/) for plotting). 

For 3D data, we write `.vtu` files that can be visualized in [paraview](https://www.paraview.org).
