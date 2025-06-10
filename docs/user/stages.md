# Stages

:::{warning}
This feature is experimental!
:::

The python bindings for GVEC provide a simple wrapper to run GVEC:
```bash
pygvec run parameter.ini
```

They also allow running GVEC with a *current constraint* (strictly speaking a current-optimization using picard iterations) and refinement.

For this, the parameters need to be specified in YAML or TOML files, which are more flexible than the classic GVEC-INI files.
The parameters use the same keys, with a different syntax for specifying the boundary and axis coefficients.
In addition the parameters for `iota`, `pres` and `sgrid` are grouped together
and three new groups of parameters, `I_tor`, `picard_current` and `stages` are available.
If `I_tor`, `picard_current` and `stages` are not present in the parameterfile, a TOML or YAML parameterfile is equivalent to an `.ini` file.

The `stages` parameter is a list of stages, which are executed in order, for which each parameter can be selected to replace the globally defined or default values.
Note that each stage inherits the base parameters, but does not take into account any previous stages (that is, after the stage the parameters will revert back to the global/default parameters).
The exception to this is `totalIter` (as described below), `iota` and `init_LA` (which is always `False` when restarting).

 Furthermore, `picard_current` defines the parameters for the algorithm when running GVEC with a fixed toroidal current profile. Both, `picard_current` and `I_tor` are required to run GVEC with fixed toroidal current.

A `pygvec` run with stages will produce as an output a directory `{ProjectName}_gvec_stages` with subdirectories containing the numbered individual GVEC runs of each stage, as well as the `parameter_{ProjectName}_final.ini` and `{ProjectName}_State_final.dat` files, which are the ini and last state file of the last restart in the last stage. These latter files can then be used for further analysis or restarts.
Note that `{ProjectName}` is the project name set in the parameter file.

When running GVEC with stages one abort criterion is again the number of iterations. The limit on the total iterations over all stages and restarts is set trough the parameter `totaliter`. Note that `totaliter` is different from the usual `maxIter`. When using `stages`, `maxIter` limits the maximum number of iterations per restart. Therefore, `maxIter` can be changed during each stage, however, `totaliter` will be kept fixed to its initial value.

### Increasing resolution

To demonstrate one intended use of `stages`, we will increase the radial resolution while simultaneously decreasing `minimize_tol` (i.e. improving the equilibrium solution) with a fixed `iota` profile. Note that the stages are independent from one another, except for `iota`. That is, if a parameter is not set during a stage, that parameter will fall back to its value outside of the stage. For example, in the parameter files below the global value for `minimize_tol` is $10^{-7}$, but we specify it during each stage and therefore it is replaced during each stage.
The example below showcases how the corresponding input files would look like when using `.toml` or `.yaml`:

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
    {minimize_tol = 1e-3, sgrid.nElems = 3},
    {minimize_tol = 1e-5, sgrid.nElems = 10},
    {minimize_tol = 1e-6, sgrid.nElems = 20},
]

[iota]
type = "polynomial"
coefs = [
    -0.8625290502868942, 0.08116648327976568, -0.3057372847655277,
    0.4672872124759052, -0.23677929291598848, -3.126329344369636,
    10.14720008596784, -14.253993484428593, 9.742801872387513,
    -2.657588003523321
    ]

[X1_b_cos]
"(0, 0)" = 5.5
"(0, 1)" = 0.2354

...
```

Full example: [`parameter.toml`](<path:../../python/examples/stages/parameter.toml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/stages/parameter.toml))

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
- minimize_tol: 0.001
  sgrid:
    nElems: 3
- minimize_tol: 1.0e-05
  sgrid:
    nElems: 10
- minimize_tol: 1.0e-06
  sgrid:
    nElems: 20

iota:
  type: polynomial
  coefs: [
    -0.8625290502868942, 0.08116648327976568, -0.3057372847655277,
    0.4672872124759052, -0.23677929291598848, -3.126329344369636,
    10.14720008596784, -14.253993484428593, 9.742801872387513,
    -2.657588003523321
    ]

X1_b_cos:
  (0, 0): 5.5
  (0, 1): 0.2354

...
```

Full example: [`parameter.yaml`](<path:../../python/examples/stages/parameter.yaml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/stages/parameter.yaml))

:::
::::
### Current constraint basics
TL;DR:
  - Set final `minimize_tol`, e.g. $10^{-5}$
  - Set `I_tor` with the same syntax as `pres` or `iota`
  - Set `picard_current="auto"`

Given `I_tor` and `picard_current`, GVEC will use picard iterations to optimize the `iota` profile, such that the resulting toroidal current converges towards the prescribed `I_tor` profile. A `iota` profile is not required but can still be provided. In this case the `iota` parameters will act as an initial guess. As mentioned above, `I_tor` has the same parameters as the other two profiles, `iota` and `pres`. Via `picard_current` we can control the behavior of the picard iterations. Per default `picard_current` is set to `off`, which corresponds to a fixed `iota` run. For a fixed current profile run, we can set `picard_current="auto"`. Via this mode, a set of stages will be automatically generated to converge both `I_tor` as well as the force tolerance specified by `minimize_tol`. Note, however, that with this `"auto"` mode the `stages` parameter must not be set in the parameter file. The generated stages can be found in the `{ProjectName}_gvec_stages/parameter_{ProjectName}.stages.toml` file.

::::{tab-set}
:::{tab-item} TOML

```{code-block} toml
:caption: `parameter.toml`
# GVEC parameter file for W7X
ProjectName = "W7X"
whichInitEquilibrium = 0
minimize_tol = 1.0e-06

...

picard_current = "auto"

[I_tor]
type = "polynomial"
coefs = [0.0]

[X1_b_cos]
"(0, 0)" = 5.5
"(0, 1)" = 0.2354

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
minimize_tol: 1.0e-06

...
picard_current: auto

I_tor:
  type: polynomial
  coefs: [0.0]

X1_b_cos:
  (0, 0): 5.5
  (0, 1): 0.2354

...
```

Full example: [`parameter.yaml`](<path:../../python/examples/current_constraint/parameter.yaml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_constraint/parameter.yaml))

:::
::::

### Current constraint advanced control

Instead of `picard_current="auto"`, we can also set the two parameters that influence the picard iterations manually. The two options are `target` an `iota_tol`. The former specifies which optimization targets are considered during a picard iteration. The latter specifies how well we want the current-constraint to be fulfilled.

Since we are technically optimizing for the prescribed `I_tor`, it can be useful to allow small deviations from `I_tor`, similarly to allowing deviations from the force balance via `minimize_tol`. As the underlying algorithm utilizes the current contribution to the rotational transform profile $\iota$, and $\iota$ is without units, we specify this deviation in terms of a tolerance on the (targeted) $\iota$: `iota_tol`. The `picard_current="auto"` mode will always try to get this tolerance below $10^{-10}$.

Given `iota_tol`, we can now choose to either aggressively optimize for `iota_tol` or choose to optimize for `minimize_tol` first and then try to also fulfill `iota_tol`. The former would require one to set `target="iota"` whereas the latter corresponds to `target="iota_and_force"`. Generally, the `target="iota"` option is intended to be used with low `maxIter` and low `iota_tol` in an initial stage, if no prior knowledge on $\iota$ is present. If the initial guess for $\iota$ is reasonable, using `target="iota_and_force"` with a low value for `iota_tol` is recommended. Such a stage is typically a follow up to a `target="iota"` stage.

The example below demonstrates the use of `picard_current` with stages. It mimics the behavior of `picard_current="auto"` but also performs refinement during the stages:

::::{tab-set}
:::{tab-item} TOML

```{code-block} toml
:caption: `parameter.toml`
# GVEC parameter file for W7X
ProjectName = "W7X"
whichInitEquilibrium = 0
minimize_tol = 1.0e-06

maxIter = 1000
totaliter = 5000
...

stages = [
    {minimize_tol = 1e-3, sgrid.nElems = 3, picard_current={target="iota",iota_tol=1e-3}, maxIter = 10},
    {minimize_tol = 1e-5, sgrid.nElems = 10, picard_current={iota_tol=1e-6}},
    {minimize_tol = 1e-6, sgrid.nElems = 20},
]

[I_tor]
type = "polynomial"
coefs = [0.0]

[picard_current]
target = "iota_and_force"
iota_tol = 1e-10


[X1_b_cos]
"(0, 0)" = 5.5
"(0, 1)" = 0.2354

...
```

Full example: [`parameter.toml`](<path:../../python/examples/current_constraint/advanced_control/parameter.toml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_constraint/advanced_control/parameter.toml))

:::

:::{tab-item} YAML
```{code-block} yaml
:caption: `parameter.yaml`
# GVEC parameter file for W7X
ProjectName: W7X
whichInitEquilibrium: 0
minimize_tol: 1.0e-06

maxIter: 1000
totaliter: 5000

...

stages:
- minimize_tol: 0.001
  sgrid:
    nElems: 3
  picard_current:
    target: iota
    iota_tol: 0.001
  maxIter: 10
- minimize_tol: 1.0e-05
  sgrid:
    nElems: 10
  picard_current:
    iota_tol: 1.0e-06
- minimize_tol: 1.0e-06
  sgrid:
    nElems: 20

picard_current:
  target: iota_and_force
  iota_tol: 1.0e-10

I_tor:
  type: polynomial
  coefs: [0.0]

X1_b_cos:
  (0, 0): 5.5
  (0, 1): 0.2354

...
```

Full example: [`parameter.yaml`](<path:../../python/examples/current_constraint/advanced_control/parameter.yaml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_constraint/advanced_control/parameter.yaml))

:::
::::
