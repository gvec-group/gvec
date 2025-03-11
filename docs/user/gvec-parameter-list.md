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

```{include} ../generators/parameters-initialization.md
```

## Discretization parameters

```{include} ../generators/parameters-discretization.md
```

## Boundary and axis parameters
```{include} ../generators/parameters-bcs.md
```

## Profile parameters
```{include} ../generators/parameters-profiles.md
```

## Minimizer parameters

In GVEC, the total MHD energy $W_{MHD}$ is minimized using a gradient based method.

```{include} ../generators/parameters-minimizer.md
```

## hmap parameters

```{include} ../generators/parameters-hmap.md
```

## Visualization and output parameters

```{include} ../generators/parameters-visualization.md
```
