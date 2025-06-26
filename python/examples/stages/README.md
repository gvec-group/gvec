# Run GVEC with several stages

The `parameter.toml` file in this example specifies a W7X equilibrium solve which uses stages to increase the radial resolution while simultaneously reducing the force tolerance (`minimize_tol`).

To run this example, install GVEC with python bindings and then execute:
```bash
pygvec run parameter.toml -p
```
or
```bash
pygvec run parameter.yaml -p
```
where the `-p` option creates two plots for the convergence of the profiles and MHD energy.

To see more command line options execute:
```bash
pygvec run -h
```
