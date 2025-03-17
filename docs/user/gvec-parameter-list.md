# GVEC parameter list

This is a list of all parameters that can be set in the GVEC parameter file.
We group the parameters by the following (sub-)categories:

1. [Initialization](#initialization-parameters)
1. [Discretization](#discretization-parameters)
1. [Boundary and axis](#boundary-and-axis-parameters)
1. [Profile specification](#profile-parameters)
1. [Minimizer](#minimizer-parameters)
1. [hmap](#hmap-parameters)
1. [Visualization/Output](#visualization-and-output-parameters)


## Initialization parameters

GVEC has basically two modes of initialization, see `whichInitequilibrium` parameter.
One choice is starting from scratch, by specifying the boundary, initial axis guess and profiles in the GVEC parameter file.
The other choice allows starting from an existing VMEC solution (`wout` file), which is read in GVEC, the boundary and magnetic axis are taken, as well as the profiles. This also sets some defaults for the GVEC discretization parameters, to match what has been used in VMEC, like the Fourier series of the variables. Many things, like discretization parameters, the boundary, axis and profiles, can then still be changed by setting them explicitly in the GVEC parameterfile.


```{include} ../generators/parameters-initialization.md
```

## Discretization parameters

There are three solution variables in GVEC, $X^1$, $X^2$ and $\lambda$ (in the parameter names, they are simply `X1`, `X2` and `LA`). They are scalar unknowns, described by a B-Spline discretization in the radial direction $\rho$ and a double-periodic Fourier series in the poloidal angle $\vartheta$ and the toroidal angle $\zeta$, see [coordinate conventions](#coordinate-conventions).

The B-Spline is defined by a polynomial degree $p$ and a radial grid, which is for now provided by the number of elements and a grid type. The Fourier series is defined by maximum mode numbers $m_\text{max},n_\text{max}$ in poloidal and toroidal direction, respectively. Also, they can be specified as sine, cosine or full sine/cosine series.

Note that each variable has its own discretization parameters, which can be set separately.


```{include} ../generators/parameters-discretization.md
```

## Boundary and axis parameters
For the variables $X^1$ and $X^2$, the boundary and initial guess for the axis must be specified.

One possibility is to provide the coefficients of the sine and cosine Fourier modes with poloidal mode number $m$ and toroidal mode number $n$. The initial guess for the axis is also given as cosine and sine Fourier mode coefficients, but obviously only toroidal modes. Alternatively, for simple enough geometries, the axis can also be guessed by setting `init_axis_average=True`.

Another possibility is to provide a specific dataset via a netcdf file, that contains the positions of the boundary on a regular grid in $\vartheta$ and $\zeta$, so $X^1(\vartheta_i,\zeta_j), X^2(\vartheta_i,\zeta_j)$ values.


```{include} ../generators/parameters-bcs.md
```

## Profile parameters
The profiles for the rotational transform (`iota`) and pressure (`pres`) are defined over a normalized radial coordinate $s=\rho^2 \in [0,1]$, where the radial coordinate $\rho\sim\sqrt{\Phi}$ is proportional to the square root of the torodial magnetic flux, so $\rho=0$ is the magnetic axis. The profile can be defined either by a polynomial, by B-spline coefficients and a knot sequence, or by a set of point-value pairs which are then interpolated.

```{include} ../generators/parameters-profiles.md
```

## Minimizer parameters

In GVEC, the total MHD energy $W_{MHD}$ is minimized using a gradient based method. The iteration stops if the norms of the gradients is below a given threshold, or at the maximum number of iterations. One can write intermediate results to disk in a given interval of iterations, and also log information on major quantities during the minimization in an iteration interval to a log-file.


```{include} ../generators/parameters-minimizer.md
```

## hmap parameters
The hmap is an exchangeable function that specifies how the $(q^1,q^2,q^3)=(X^1,X^2,\zeta)$ variables are mapped to Cartesian coordinates $(x,y,z)$. The simplest map is a straight periodic cylinder, or the cylinder coordinates. There is a more academic knot-map, and then two types of axis-following frames, the Frenet frame given only by a closed 3D curve, provided as a Fourier series in $R,Z$. And the generalized frame, or G-Frame, where a curve and the normal and binormal vector are given as at point positions/vectors in cartesian coordinates, which are provided via a dataset in a netcdf file.

```{include} ../generators/parameters-hmap.md
```

## Visualization and output parameters

Even though the visualization can now be done via the python bindings, the fortran executable also allows for a number of visualization output files during the minimization run, with  `*visu*.csv` for 1D data and as paraview files (`*visu*.vtu`) for 3D data. The 3D data can also be written as a netcdf file (`*visu*.nc`).

```{include} ../generators/parameters-visualization.md
```
