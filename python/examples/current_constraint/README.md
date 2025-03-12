# Run GVEC with zero-current-optimization and several stages

The `parameter.toml` file in this example specifies a W7X equilibrium solve which is optimized for zero toroidal current using picard iterations
This mode of operation is also called *current constraint*.
The radial resoultion and maximum iterations of each GVEC run (picard iteration) is set in two stages.

To run this example, install GVEC with python bindings and then execute:
```bash
pygvec run parameter.toml -p
```
where the `-p` option creates two plots for the convergence of the profiles and MHD energy.

To see more command line options execute:
```bash
pygvec run -h
```
