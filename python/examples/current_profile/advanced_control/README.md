# Run GVEC with automatic zero-current-optimization

The `parameter.toml` file in this example specifies a W7X equilibrium solve which is optimized for zero toroidal current using picard iterations.
This mode of operation is also called *current optimization*. In the example the required stages are automatically generated.

To run this example, install GVEC with python bindings and then execute:
```bash
pygvec run parameter.toml -p
```
where the `-p` option creates plots for the convergence of the profiles MHD energy and forces.

To see more command line options execute:
```bash
pygvec run -h
```
