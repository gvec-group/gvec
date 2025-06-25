# Stages

The python bindings for GVEC provide a simple wrapper to run GVEC:
```bash
pygvec run parameter.ini
```

They also allow running GVEC with a prescribed current profile and refinement.

For this, the parameters **need to be specified in YAML or TOML files**, which are more flexible than the classic GVEC-INI files. The parameters use the same keys, with a different syntax for specifying the boundary and axis coefficients.

Additionally, the parameters for `iota`, `pres` and `sgrid` are grouped together
and three **new groups of parameters**,`I_tor` `picard_current` and `stages` are available.
If `I_tor`, `picard_current` and `stages` are not present in the parameterfile, a TOML or YAML parameterfile is equivalent to an `.ini` file.

The `stages` parameter is a list of stages, which are executed in order, for which each parameter can be selected to replace the globally defined or default values.
Note that each stage inherits the base parameters, but does not take into account any previous stages (that is, after the stage the parameters will revert back to the global/default parameters).
The exception to this is `totalIter` (as described below), `iota` and `init_LA` (which is always `False` when restarting).

 Furthermore, `picard_current` defines the parameters for the algorithm when running GVEC with a fixed toroidal current profile. Both, `picard_current` and `I_tor` are required to run GVEC with fixed toroidal current.

A `pygvec` run with stages will produce as an output a directory `{ProjectName}_gvec_stages` with subdirectories containing the numbered individual GVEC runs of each stage, as well as the `parameter_{ProjectName}_final.ini` and `{ProjectName}_State_final.dat` files, which are the `.ini` and last state file of the last restart in the last stage. These latter files can then be used for further analysis or restarts.
Note that `{ProjectName}` is the project name set in the parameter file.

When running GVEC with stages, one abort criterion is again the number of iterations. The limit on the total iterations over all stages and restarts is set through the parameter `totalIter`. Note that `totalIter` is different from the usual `maxIter`. When using `stages`, `maxIter` limits the maximum number of iterations per restart. Therefore, `maxIter` can be changed during each stage, however, `totalIter` will be kept fixed to its initial value.

### Increasing resolution

To demonstrate one intended use of `stages`, we will increase the radial resolution while simultaneously decreasing `minimize_tol` (i.e. improving the equilibrium solution) with a fixed `iota` profile. Note that the stages are independent of one another, except for `iota`. That is, if a parameter is not set within a stage, that parameter will fall back to its value specified outside. For example, in the parameter files below the global value for `minimize_tol` is $10^{-7}$, but we specify it during each stage, thus it is replaced during each stage.
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

## Prescribing a toroidal current profile
### Theory

Instead of fixing the $\iota$ profile, leading to a fixed poloidal flux $\chi'(\rho) =  \iota(\rho)\Phi'(\rho)$, we could also fix the toroidal current profile $I_{\text{tor}}(\rho)$. The toroidal current profile is defined as

$$
  I_{\text{tor}}(\hat{\rho}) := \int_0^{\hat{\rho}}\int_0^{2\pi} \mathcal{J} J^{\zeta}\, \mathrm{d}\rho \mathrm{d}\vartheta,
$$

with $J^\zeta$ denoting the contravariant toroidal current density and $\mathcal{J}$ the usual Jacobian determinant. Expressing $J^\zeta$ in terms of the covariant magnetic field components yields

$$
I_{\text{tor}}(\hat{\rho}) = \frac{1}{\mu_0} \int_0^{\hat{\rho}} \int_0^{2\pi} \mathcal{J}\frac{1}{\mathcal{J}}\left(\frac{\partial B_\vartheta}{\partial\rho}-\frac{\partial B_\rho}{\partial\vartheta}\right)\, \mathrm{d}\rho \,\mathrm{d}\vartheta.
$$

 Due to the integration with respect to $\vartheta$ and periodicity the term with $\frac{\partial B_\rho}{\partial\vartheta}$ vanishes. Hence, we can express the toroidal current profile as

$$
  I_{\text{tor}}(\rho) = \frac{1}{\mu_0} \int_0^{2\pi} \left . B_\vartheta\right |_{\rho} \, \mathrm{d}\vartheta \frac{1}{2\pi\mu_0} \int_0^{2\pi}  \int_0^{2\pi}  \left . B_\vartheta\right |_{\rho} \, \mathrm{d}\vartheta \, \mathrm{d}\zeta= \frac{1}{2\pi\mu_0} {\left\langle B_\vartheta |_{\rho} \right\rangle}.
$$

$B_\vartheta$ can be obtained via (see also magnetic field in  [theory](theory.md)):

$$
\begin{align*}
B_\vartheta &= g_{\vartheta\vartheta}B^\vartheta + g_{\vartheta\zeta}B^\zeta\\
& = \frac{g_{\vartheta\vartheta}}{\mathcal{J}}\left(\chi'-\Phi'\frac{\partial\lambda}{\partial\zeta}\right)+\frac{g_{\vartheta\zeta}}{\mathcal{J}}\Phi'\left(1+\frac{\partial\lambda}{\partial\vartheta}\right).
\end{align*}
$$

Inserting this into the equation for $I_{\text{tor}}$ lets us solve for $\chi'$. As $\Phi'$ is known, we can also solve for $\iota = \frac{\chi'}{\Phi'}$. By introducing the quantities $\Gamma_\vartheta=\frac{g_{\vartheta\vartheta}}{\mathcal{J}}$ and $\Gamma_\zeta=\frac{g_{\vartheta\zeta}}{\mathcal{J}}$ the final expression for $\iota$ reads:

$$
  \iota(\rho) = \frac{2\pi\mu_0}{\Phi'(\rho)\left\langle\Gamma_\vartheta|_{\rho}\right\rangle}I_{\text{tor}}(\rho)+\frac{\Phi'(\rho)}{\left\langle\Gamma_\vartheta|_{\rho}\right\rangle}\left\langle\left.\left(\Gamma_\vartheta\frac{\partial\lambda}{\partial\zeta}-\Gamma_\zeta\left(1+\frac{\partial\lambda}{\partial\vartheta}\right)\right)\right|_{\rho}\right\rangle
$$

Here, we can identify two contributions to the rotational transform; one independent of $I_{\text{tor}}$ which we refer to as $\iota_0$ and one depending on $I_{\text{tor}}$ referred to here as $\iota_{\text{curr}}$:

$$
\begin{align*}
  \iota_{\text{curr}}(\rho) &:= \frac{2\pi\mu_0}{\Phi'(\rho)\left\langle\Gamma_\vartheta|_{\rho}\right\rangle}I_{\text{tor}}(\rho),\\
  \iota_{0}(\rho)&:= \frac{\Phi'(\rho)}{\left\langle\Gamma_\vartheta|_{\rho}\right\rangle}\left\langle\left.\left(\Gamma_\vartheta\frac{\partial\lambda}{\partial\zeta}-\Gamma_\zeta\left(1+\frac{\partial\lambda}{\partial\vartheta}\right)\right)\right|_{\rho}\right\rangle.
\end{align*}
$$

When prescribing a toroidal current profile $I_{\text{tor}}(\rho)$ in GVEC, the code will start with an initial guess for the $\iota$ profile and do a Picard iteration to find the $\iota$ profile that yields the desired toroidal current profile. Each Picard iteration corresponds to an energy minimization with a fixed $\iota(\rho)$, starting from the solution of the previous Picard iteration. We prescribe the $\iota$ profile as the target $\iota_T=\iota_0+\iota_{\text{curr}}$, using the given current profile and the previous GVEC solution for the terms $\Gamma_\vartheta,\Gamma_\zeta,\lambda$. After each Picard iteration, we can check for the difference $\Delta\iota$ between the target and $\iota(\rho)$ computed from the GVEC solution only. Finally, for convergence, one has to check that both $\Delta\iota$ and the MHD equilibrium forces are below the desired tolerance. $\Delta\iota$ is calculated from a set of $N$ evaluations at different radial positions $\rho_i \in [0,1]$:

$$
\Delta\iota := \sqrt{\frac{1}{N}\sum_i^N\left(\iota_T\left(\rho_i\right)-\iota\left(\rho_i\right)\right)²}.
$$

The general approach is then:
1. Prescribe $I_{\text{tor}}(\rho)$ and an initial guess for $\iota$ (the default initial guess is $\iota(\rho)=0$)
2. Minimize energy with fixed $\iota(\rho)$ until the force tolerance is reached or the maximum number of iterations for this minimization is reached.
3. Calculate $\Delta\iota$
4. if $\Delta\iota$ is below the tolerance for $\iota$, finish; else, update $\iota(\rho) = \iota_0(\rho)+\iota_{\text{curr}}(\rho)$ using $I_{\text{tor}}(\rho)$ and continue.
5. Repeat from 2.

#### The default approach
In general $\iota$ could be updated in each step of the energy minimization. However, we found that with the given implementation it is typically better for convergence to delay this update until the DOF have been sufficiently modified by the energy minimization, for a reasonable guess for $\iota$. Thus, the default algorithm for GVEC with prescribed current profile will automatically generate stages that can be separated into two phases. In the initial phase the energy minimization with fixed $\iota$ is run for only a maximum of $10$ iterations (`maxIter=10`) until $\iota$ is updated. This ensures that $\iota$ is quickly adapted if the initial guess for $\iota$ is poor. That is, during this phase $\iota_T$ is the prime target for the optimization. If $\Delta\iota$ falls below a tolerance of $10^{-3}$ the second phase is started. In this phase we aggressively target the prescribed force tolerance (`minimize_tol`) while keeping $\iota$ fixed. Only when the force tolerance is reached, $\iota$ is updated. With the updated $\iota$, energy minimization with fixed $\iota$ is performed until the force tolerance is reached again. This process is repeated until both $\Delta\iota$ and the forces are below the prescribed tolerance or the computational budget is exhausted. This procedure could be achieved simply within two separate stages. However, numerical experiments have shown that ramping up the tolerances for $\Delta\iota$ and the forces is advantageous. Therefore, the second phase is automatically split into several stages where the force and $\Delta\iota$ tolerance is increased until the prescribed tolerances are reached.

### Basic control

The following parameters can be set to control the current profile optimization:
  - Set `picard_current="auto"`
  - Set final `minimize_tol`, e.g. $10^{-5}$
  - Set `I_tor` with the same syntax as `pres` or `iota`


Given `I_tor` and `picard_current`, GVEC will use Picard iterations as described above to find the `iota` profile where the resulting toroidal current profile matches the prescribed `I_tor` profile. A `iota` profile is not required, but can still be provided. In this case the `iota` parameters will act as an initial guess. As mentioned above, `I_tor` has the same profile definition parameters as the other two profiles, `iota` and `pres`.

Via `picard_current` we can control the behavior of the Picard iterations. Per default `picard_current` is set to `off`, which corresponds to a fixed `iota` run. For a fixed current profile run, we can set `picard_current="auto"`. Via this mode, a set of stages will be automatically generated to converge both `I_tor` as well as the force tolerance specified by `minimize_tol`. Note, however, that with this `"auto"` mode, the `stages` parameter must not be set in the parameter file. The generated stages can be found in the `{ProjectName}_gvec_stages/parameter_{ProjectName}.stages.toml` file.

::::{tab-set}
:::{tab-item} TOML

```{code-block} toml
:caption: `parameter.toml`
# GVEC parameter file for W7X
ProjectName = "W7X"
whichInitEquilibrium = 0
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
whichInitEquilibrium: 0
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

### Advanced control

Instead of `picard_current="auto"`, we can also set the two parameters that influence the Picard iterations manually. The two key options are `target` and `iota_tol`. The former specifies which targets are considered during a Picard iteration. The latter specifies how well the current-constraint has to be fulfilled.

Since we are technically optimizing for the prescribed `I_tor`, it can be useful to allow small deviations from `I_tor`, similarly to allowing deviations from the force balance via `minimize_tol`. As the underlying algorithm utilizes $\iota_{\text{curr}}$, and $\iota$ is without units, we specify this deviation in terms of a tolerance on the (targeted) rotational transform $\iota_T$: `iota_tol`. The `picard_current="auto"` mode will always try to get this tolerance below $10^{-10}$.

Given `iota_tol`, we can now choose `target="iota"` to either aggressively optimize for `iota_tol`, or we can choose `target="iota_and_force"` to optimize for `minimize_tol` first and then try to also fulfill `iota_tol`. Generally, the `target="iota"` option is intended to be used with low `maxIter` and low `iota_tol` in an initial stage, if no prior knowledge on $\iota$ is present. If the initial guess for $\iota$ is reasonable, using `target="iota_and_force"` with a low value for `iota_tol` is recommended. Such a stage is typically follows a `target="iota"` stage.

The final parameter that might be set for `picard_current` is `maxRestarts`. This parameter limits the maximum number of restarts that can be performed during a stage. Per default `maxRestarts=30`. Its intended use is to have an abort criterion if the targeted `iota_tol` can not be reached during a stage, similar to `maxIter` and `totalIter`.

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
totalIter = 5000
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
maxRestarts = 30

[X1_b_cos]
"(0, 0)" = 5.5
"(0, 1)" = 0.2354

...
```

Full example: [`parameter.toml`](<path:../../python/examples/current_profile/advanced_control/parameter.toml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_profile/advanced_control/parameter.toml))

:::

:::{tab-item} YAML
```{code-block} yaml
:caption: `parameter.yaml`
# GVEC parameter file for W7X
ProjectName: W7X
whichInitEquilibrium: 0
minimize_tol: 1.0e-06

maxIter: 1000
totalIter: 5000

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
  maxRestarts: 30

I_tor:
  type: polynomial
  coefs: [0.0]

X1_b_cos:
  (0, 0): 5.5
  (0, 1): 0.2354

...
```

Full example: [`parameter.yaml`](<path:../../python/examples/current_profile/advanced_control/parameter.yaml>)
(view [online {fab}`square-gitlab`](https://gitlab.mpcdf.mpg.de/gvec-group/gvec/-/blob/develop/python/examples/current_profile/advanced_control/parameter.yaml))

:::
::::
