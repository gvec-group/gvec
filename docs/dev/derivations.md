# Derivations for computable Quantities

This page collects some derivations for computable quantities:

Throughout this document we use $\alpha,\beta$ as generic indices for the coordinates $\rho,\theta,\zeta$ unless otherwise specified.

## Flux aligned coordinates

$$
\begin{align}
\vec{e}_\alpha &= \pdv{\vec{x}}{\alpha} \\
\vec{e}^\alpha &= \grad \alpha \\
\Jac &= \vec{e}_\rho \cdot (\vec{e}_\theta \times \vec{e}_\zeta) \\
\Jac &= \qty(\grad \rho \cdot \grad \theta \times \grad \zeta)^{-1} \\
\end{align}
$$

$$
\begin{align}
\vec{e}_\rho &= \Jac \grad\theta \times \grad\zeta \\
\vec{e}_\theta &= \Jac \grad\zeta \times \grad\rho \\
\vec{e}_\zeta &= \Jac \grad\rho \times \grad\theta \\
\end{align}
$$

$$
\begin{align}
\vec{k}_{\alpha\beta} &:= \pdv[2]{\vec{x}}{\alpha}{\beta} \\
&= \pdv{\vec{e}_\alpha}{\beta} = \pdv{\vec{e}_\beta}{\alpha} \\
\end{align}
$$

## Magnetic field $\vec{B}$


$$
\begin{align}
B^\rho &= 0 \\
B^\thet &= \frac{1}{\Jac} \dPhidr \left(\iota-\pdv{\lambda}{\zeta} \right) \\
B^\zeta &= \frac{1}{\Jac} \dPhidr \left(1+\pdv{\lambda}{\thet} \right) \\
\end{align}
$$


$\Rightarrow$ `B_contra_t`, `B_contra_z`

## Derivatives of $\modB$


$$
\begin{align}
\modB^2 &= B_\thet B^\thet + B_\zeta B^\zeta\\
\modB^2 &= B^\thet B^\thet g_{\thet\thet} + 2 B^\thet B^\zeta g_{\thet\zeta} + B^\zeta B^\zeta g_{\zeta\zeta} \\
\pdv{\modB^2}{\alpha} &= 2 B^\thet \pdv{B^\thet}{\alpha} g_{\thet\thet} + B^\thet B^\thet \pdv{g_{\thet\thet}}{\alpha} \\
                                        &+ 2 \left(\pdv{B^\thet}{\alpha} B^\zeta + B^\thet \pdv{B^\zeta}{\alpha}\right) g_{\thet\zeta}
                                        + B^\thet B^\zeta \pdv{g_{\thet\zeta}}{\alpha} \\
                                        &+ 2 B^\zeta \pdv{B^\zeta}{\alpha} g_{\zeta\zeta} + B^\zeta B^\zeta \pdv{g_{\zeta\zeta}}{\alpha} \\
\pdv{\modB}{\alpha} &= \frac{1}{2\modB} \pdv{\modB^2}{\alpha} = \frac{\vec{B}}{\modB} \cdot \pdv{\vec{B}}{\alpha} \\
\nabla \modB &= \pdv{\modB}{\rho} \nabla\rho + \pdv{\modB}{\thet} \nabla\thet + \pdv{\modB}{\zeta} \nabla\zeta \\
\end{align}
$$


for $\alpha\in\left\{\rho,\thet,\zeta\right\}$.

$\Rightarrow$ `dmodB_dr`, `dmodB_dt`, `dmodB_dz`, `grad_modB`

## Derivatives of $\vec{B}$

As $\vec{B}$ is a vector, its derivative $\grad\vec{B}$ is a matrix, which complicates things a bit, as the matrix product is not commutative. We denote the *outer product* with $\otimes$.

$$
\begin{align}
\grad \vec{B} &:= \grad \otimes \vec{B} \\
&= \sum_\alpha \grad\alpha \otimes \pdv{\vec{B}}{\alpha} \\
&= \grad\rho \otimes \pdv{\vec{B}}{\rho} + \grad\thet \otimes \pdv{\vec{B}}{\thet} + \grad\zeta \otimes \pdv{\vec{B}}{\zeta} \\
\pdv{\vec{B}}{\alpha} &:= \sum_\beta \pdv{\qty(B^\beta \vec{e}_\beta)}{\alpha} \\
&= \sum_\beta \pdv{B^\beta}{\alpha} \vec{e}_\beta + B^\beta \vec{k}_{\beta\alpha} \\
&= \cancel{\pdv{B^\rho}{\alpha}} \erho + \cancel{B^\rho} \vec{k}_{\rho\alpha}
+ \pdv{B^\thet}{\alpha} \ethet + B^\thet \vec{k}_{\thet\alpha}
+ \pdv{B^\zeta}{\alpha} \ezeta + B^\zeta \vec{k}_{\zeta\alpha} \\
\nabla \vec{B} &= \sum_{\alpha\beta} \grad\alpha \otimes \left(\pdv{B^\beta}{\alpha} \vec{e}_\beta + B^\beta \vec{k}_{\beta\alpha} \right)
\end{align}
$$

for $\alpha,\beta\in\left\{\rho,\thet,\zeta\right\}$.

Due to `xarray`'s limitations with multiple dimensions of the same name, $\grad\vec{B}$ is not directly available as a computable quantity, and its components should be used instead.

$\Rightarrow$ `dB_dr`, `dB_dt`, `dB_dz`

---

From

$$
\begin{align}
B^\thet &= \frac{1}{\Jac} \dPhidr \left(\iota-\pdv{\lambda}{\zeta} \right) \\
B^\zeta &= \frac{1}{\Jac} \dPhidr \left(1+\pdv{\lambda}{\thet} \right) \\
\end{align}
$$

it follows

$$
\begin{align}
\pdv{B^\thet}{\rho} &=
    \frac{1}{\Jac} \frac{d^2\Phi}{d\rho^2} \left(\iota-\pdv{\lambda}{\zeta} \right)
    + \frac{1}{\Jac} \dPhidr \left(\frac{d\iota}{d\rho}-\pdv{^2\lambda}{\rho\partial\zeta} \right)
    - \pdv{\Jac}{\rho} \frac{1}{\Jac^2} \dPhidr \left(\iota-\pdv{\lambda}{\zeta} \right) \\
\pdv{B^\thet}{\thet} &=
    - \frac{1}{\Jac} \dPhidr \pdv{^2\lambda}{\thet\partial\zeta}
    - \pdv{\Jac}{\thet} \frac{1}{\Jac^2} \dPhidr \left(\iota-\pdv{\lambda}{\zeta} \right) \\
\pdv{B^\thet}{\zeta} &=
    - \frac{1}{\Jac} \dPhidr \pdv{^2\lambda}{\zeta^2}
    - \pdv{\Jac}{\zeta} \frac{1}{\Jac^2} \dPhidr \left(\iota-\pdv{\lambda}{\zeta} \right) \\
\pdv{B^\zeta}{\rho} &=
    \frac{1}{\Jac} \frac{d^2\Phi}{d\rho^2} \left(1+\pdv{\lambda}{\thet} \right)
    + \frac{1}{\Jac} \dPhidr \pdv{^2\lambda}{\rho\partial\thet}
    - \pdv{\Jac}{\rho} \frac{1}{\Jac^2} \dPhidr \left(1+\pdv{\lambda}{\thet} \right) \\
\pdv{B^\zeta}{\thet} &=
    \frac{1}{\Jac} \dPhidr \pdv{^2\lambda}{^2\thet}
    - \pdv{\Jac}{\thet} \frac{1}{\Jac^2} \dPhidr \left(1+\pdv{\lambda}{\thet} \right) \\
\pdv{B^\zeta}{\zeta} &=
    \frac{1}{\Jac} \dPhidr \pdv{^2\lambda}{\thet\partial\zeta}
    - \pdv{\Jac}{\zeta} \frac{1}{\Jac^2} \dPhidr \left(1+\pdv{\lambda}{\thet} \right) \\
\end{align}
$$

$\Rightarrow$ `dB_contra_t_dr`, `dB_contra_t_dt`, `dB_contra_t_dz`, `dB_contra_z_dr`, `dB_contra_z_dt`, `dB_contra_z_dz`

---

The gradient $\nabla\vec{B}$ can be used to compute the *magnetic gradient length scale* $L_{\nabla\vec{B}}$ (`L_gradB`). Details are found in *John Kappel et al 2024 PPCF 66 025018* [DOI:10.1088/1361-6587/ad1a3e](https://www.doi.org/10.1088/1361-6587/ad1a3e).


## Normalized magnetic field

The normalized magnetic field (i.e. the unit vector along the magnetic field) is defined as

$$
\begin{align}
\vec{b} &:= \frac{\vec{B}}{\modB} \\
\end{align}
$$

Its gradient (a matrix) is then:

$$
\begin{align}
\grad\vec{b} &= \sum_\alpha \grad\alpha \otimes \pdv{\vec{b}}{\alpha} \\
\pdv{\vec{b}}{\alpha} &= \frac{1}{\modB}\pdv{\vec{B}}{\alpha} - \frac{\vec{B}}{\modB^2}\pdv{\modB}{\alpha} \\
&= \frac{1}{\modB}\qty(\pdv{\vec{B}}{\alpha} - \vec{b} \qty(\vec{b} \cdot \pdv{\vec{B}}{\alpha})) \\
\end{align}
$$

$\Rightarrow$ `db_dr`, `db_dt`, `db_dz`

## Geodesic curvature

The fieldline curvature vector $\vec{\kappa}_B$ is defined as

$$
\begin{align}
\vec{\kappa}_B &:= \vec{b} \cdot \grad\vec{b} \\
&= \sum_\alpha \qty(\frac{\vec{B}}{\modB} \cdot \grad\alpha) \pdv{\vec{b}}{\alpha} \\
\end{align}
$$

The geodesic curvature $\kappa_G$ is defined as

$$
\kappa_G &:= \vec{\kappa}_B \cdot \qty(\grad\rho \times \vec{b}) \\
$$

Which could be rewritten as

$$
\begin{align}
\kappa_G &= \vec{b} \cdot \qty(\grad \vec{b}) \cdot \qty(\grad\rho \times \vec{b}) \\
&= \sum_\alpha \qty(\frac{\vec{B}}{\modB} \cdot \grad\alpha) \pdv{\vec{b}}{\alpha} \cdot \qty(\grad\rho \times \vec{b}) \\
\end{align}
$$

$\Rightarrow$ `kappa_B`, `kappa_G`

## Effective geometric quantities

The plasma volume $V$, surface are $A_\text{surface}$ and length of the magnetic axis $L_\text{axis}$ are defined as

$$
\begin{align}
V &= \int_0^1 d\rho \int_0^{2\pi} d\theta \int_0^{2\pi} d\zeta \Jac \\
A_\text{surface} &= \left. \int_0^{2\pi} d\theta \int_0^{2\pi} d\zeta \left|\vec{e}_\theta \times \vec{e}_\zeta\right| \right|_{\rho=1}\\
L_\text{axis} &= \left. \int_0^{2\pi} d\zeta \left|\vec{e}_\zeta\right| \right|_{\rho=0,\theta=0}\\
\end{align}
$$

All three are independent of the coordinate frame and parametrization, and $V,A_\text{surface}$ only depend on the boundary shape and not the full equilibrium.

From these we can compute the *effective minor radius*, *effective major radius*, *effective aspect ratio* and *effective elongation*, by relating them to an equivalent torus with elliptical cross section.
In particular we define the equivalent torus as a torus with circular axis and constant elliptical cross-section, which has the same volume, surface area and axis length as the configuration of interest.

$$
\begin{align}
V &= 2\pi^2 r_\text{minor,eff}^2 r_\text{major,eff} \\
A_\text{surface} &= 4\pi^2 r_\text{minor,eff} r_\text{major,eff} \tilde{C}(E_\text{eff}) \\
L_\text{axis} &= 2\pi r_\text{major,eff}
\end{align}
$$

with the effective elongation $E_\text{eff} := \frac{a}{b}$ defined as the ratio of the cross-sections semi-major axis and semi-minor axis ($a \geq b$).
The circumference of the ellipse does not have a closed form.
We use Ramanujan's approximation

$$
\begin{align}
C &= 2\pi r_\text{minor,eff} \tilde{C} \\
\tilde{C} &= \frac{E_\text{eff} + 1}{2\sqrt{E_\text{eff}}} \left[ 1 + 3 \frac{h}{10 + \sqrt{4 - 3h}} \right] \\
h &:= \frac{(E_\text{eff} - 1)^2}{(E_\text{eff} + 1)^2}.
\end{align}
$$

We invert these formulas to obtain

$$
\begin{align}
r_\text{major,eff} &= \frac{L_\text{axis}}{2\pi} \\
r_\text{minor,eff} &= \sqrt{\frac{V}{\pi L_\text{axis}}} \\
a_\text{eff} &= \frac{r_\text{major,eff}}{r_\text{minor,eff}} \\
\tilde{C}(E_\text{eff}) &= \frac{A_\text{surface}}{2\sqrt{\pi V L_\text{axis}}},
\end{align}
$$

where we find $E_\text{eff}$ from the value of $\tilde{C}$ using the Newton-Raphson method.

$\Rightarrow$ `V`, `A_surface`, `L_axis`, `r_minor`, `r_major`, `aspect_ratio`, `elongation`

## Vacuum magnetic well depth

We define the *vacuum magnetic well depth* as

$$
\begin{align}
d_\text{well} &= \frac{\frac{dV}{d\Phi_n}(\rho=0) - \frac{dV}{d\Phi_n}(\rho=1)}{\frac{dV}{d\Phi_n}(\rho=0)}
\end{align}
$$

positive values of $d_\text{well}$ indicate $\frac{d^2V}{d\Phi_n^2} < 0$ which is favorable for stability.

$\Rightarrow$ `vacuum_magnetic_well_depth`

## Mercier criterion

We follow the formulas reported in *Landreman & Jorge* (2020), given in *Bauer et al.* (1984); *Ichiguchi et al.* (1993).

$$
D_\text{Merc} = D_\text{M,Shear} + D_\text{M,Curr} + D_\text{M,Well} + D_\text{M,Geod},
$$

where a positive value of $D_\text{Merc}$ indicates stability and

$$
\begin{align}
D_\text{M,Shear} &= \frac{1}{16\pi^2} \pqty{\dv{\iota}{\Phi}}^2 \\
D_\text{M,Curr} &= -\frac{s_\chi}{\pqty{2\pi}^4} \dv{\iota}{\Phi} \int\dd S \frac{\vec{\Xi} \cdot \vec{B}}{\abs{\grad\Phi}^3} \\
D_\text{M,Well} &= \frac{\mu_0}{\pqty{2\pi}^6} \dv{p}{\Phi} \pqty{s_\Phi \dv[2]{V}{\Phi} - \mu_0 \dv{p}{\Phi} \int dS \frac{1}{\modB^2 \abs{\grad\Phi}}}\int\dd S \frac{\modB^2}{\abs{\grad\Phi}^3} \\
D_\text{M,Geod} &= \frac{1}{\pqty{2\pi}^6} \bqty{
    \pqty{\int\dd S \frac{\mu_0 \vec{J}\cdot\vec{B}}{\abs{\grad\Phi}^3}}^2
   -\pqty{\int\dd S \frac{\modB^2}{\abs{\grad\Phi}^3}}
    \pqty{\int\dd S \frac{\pqty{\mu_0 \vec{J}\cdot\vec{B}}^2}{\modB^2 \abs{\grad\Phi}^3}}
} \\
\end{align}
$$

with

$$
\begin{align}
s_\chi &= \text{sgn}{\chi} \\
s_\Phi &= \text{sgn}{\Phi} \\
\dd S &= \abs{\grad\rho} \abs{\Jac} \dd\theta \dd\zeta = \abs{\vec{e}_\theta \times \vec{e}_\zeta} \dd\theta \dd\zeta\\
\vec{\Xi} &= \mu_0 \vec{J} - \dv{\aqty{B_\theta}}{\Phi} \vec{B} \\
\frac{\vec{\Xi} \cdot \vec{B}}{\abs{\grad\Phi}^3} &= \frac{\mu_0 \vec{J}\cdot\vec{B}}{\abs{\grad\Phi}^3} - \dv{\aqty{B_\theta}}{\Phi} \frac{\modB^2}{\abs{\grad\Phi}^3} \\
\grad\Phi &= \dv{\Phi}{\rho} \grad\rho \\
\dv{\iota}{\Phi} &= \dv{\iota}{\rho} \pqty{\dv{\Phi}{\rho}}^{-1} \\
\dv{p}{\Phi} &= \dv{p}{\rho} \pqty{\dv{\Phi}{\rho}}^{-1} \\
\dv{\aqty{B_\theta}}{\Phi} &= \dv{\aqty{B_\theta}}{\rho} \pqty{\dv{\Phi}{\rho}}^{-1} \\
\dv[2]{V}{\Phi} &= \dv[2]{V}{\rho} \pqty{\dv{\Phi}{\rho}}^{-2} - \dv{V}{\rho} \dv[2]{\Phi}{\rho} \pqty{\dv{\Phi}{\rho}}^{-3} \\
\end{align}
$$

Note that all four terms scale with $\pqty{\dv{\Phi}{\rho}}^{-2}$.

$\Rightarrow$ `D_Merc`, `D_Merc_Shear`, `D_Merc_Curr`, `D_Merc_Well`, `D_Merc_Geod`

## Transformation to Boozer coordinates

In terms of the GVEC coordinates $(\rho,\thet,\zeta)$ the Boozer transform is given as

$$
\begin{align}
\rho_B &= \rho \\
\thet_B &= \thet +\lambda(\rho,\thet,\zeta) +\iota(\rho)\nu_B(\rho,\thet,\zeta) \\
\zeta_B &= \zeta \nu_B(\rho,\thet,\zeta) \\
\end{align}
$$

The derivatives are

$$
\begin{align}
\pdv{\rho_B}{\rho} &= 1, &
\pdv{\thet_B}{\rho} &= \pdv{\lambda}{\rho}+\pdv{\iota}{\rho}\nu_B + \iota\pdv{\nu_B}{\rho} ,&
\pdv{\zeta_B}{\rho} &= \pdv{\nu_B}{\rho} \\
\pdv{\rho_B}{\thet} &=  0 ,&
\pdv{\thet_B}{\thet} &= 1 + \pdv{\lambda}{\thet} + \iota\pdv{\nu_B}{\thet} ,&
\pdv{\zeta_B}{\thet} &= \pdv{\nu_B}{\thet} \\
\pdv{\rho_B}{\zeta} &= 0  ,&
\pdv{\thet_B}{\zeta} &= \pdv{\lambda}{\zeta} + \iota\pdv{\nu_B}{\zeta} ,&
\pdv{\zeta_B}{\zeta} &= 1 + \pdv{\nu_B}{\zeta} \\
\end{align}
$$

We can compute the ratio of the Jacobian determinants

$$\begin{align}
\frac{\Jac}{\Jac_B} &= \frac{\erho \cdot(\ethet \times \ezeta)}{\vec{e}_{\rho_B} \cdot(\vec{e}_{\thet_B} \times \vec{e}_{\zeta_B})}  \\
&= \pdv{\thet_B}{\thet}\pdv{\zeta_B}{\zeta} - \pdv{\thet_B}{\zeta}\pdv{\zeta_B}{\thet} \\
\end{align}
$$

with that we can express the inverse derivatives

$$
\begin{align}
\pdv{\rho}{\rho_B} &= 1, &
\pdv{\rho}{\thet_B} &= 0  ,&
\pdv{\rho}{\zeta_B} &= 0 \\
 \pdv{\thet}{\rho_B} &= \frac{\Jac_B}{\Jac}\left(\pdv{\thet_B}{\zeta} \pdv{\zeta_B}{\rho} - \pdv{\zeta_B}{\zeta} \pdv{\thet_B}{\rho}\right ),&
\pdv{\thet}{\thet_B} &= \frac{\Jac_B}{\Jac} \pdv{\zeta_B}{\zeta},&
\pdv{\thet}{\zeta_B} &= -\frac{\Jac_B}{\Jac} \pdv{\thet_B}{\zeta} \\
\pdv{\zeta}{\rho_B} &=  \frac{\Jac_B}{\Jac}\left(\pdv{\zeta_B}{\thet} \pdv{\thet_B}{\rho} - \pdv{\thet_B}{\thet} \pdv{\zeta_B}{\rho} \right ),&
\pdv{\zeta}{\thet_B} &=  -\frac{\Jac_B}{\Jac}\pdv{\zeta_B}{\thet},&
\pdv{\zeta}{\zeta_B} &=  \frac{\Jac_B}{\Jac}  \pdv{\thet_B}{\thet}\\
\end{align}
$$

The basis vectors are then computed as:

$$
\begin{align}
\vec{e}_{\rho_B} &=&  \erho + &\pdv{\thet}{\rho_B}\ethet + \pdv{\zeta}{\rho_B}\ezeta \\
\vec{e}_{\thet_B} &=&  &\pdv{\thet}{\thet_B}\ethet + \pdv{\zeta}{\thet_B}\ezeta \\
\vec{e}_{\zeta_B} &=&  &\pdv{\thet}{\zeta_B}\ethet + \pdv{\zeta}{\zeta_B}\ezeta \\
\end{align}
$$

and equivalently

$$
\begin{align}
\nabla \rho_B &= \nabla \rho \\
\nabla \thet_B &= \pdv{\thet_B}{\rho}\nabla \rho + \pdv{\thet_B}{\thet}\nabla \thet + \pdv{\thet_B}{\zeta}\nabla \zeta \\
\nabla \zeta_B &= \pdv{\zeta_B}{\rho}\nabla \rho + \pdv{\zeta_B}{\thet}\nabla \thet + \pdv{\zeta_B}{\zeta}\nabla \zeta \\
\end{align}
$$

## Transformation to PEST coordinates

In terms of the GVEC coordinates $(\rho,\thet,\zeta)$ the PEST transform is given as

$$
\begin{align}
\rho_P &= \rho \\
\thet_B &= \thet + \lambda(\rho,\thet,\zeta) \\
\zeta_P &= \zeta \\
\end{align}
$$

The derivatives are

$$
\begin{align}
\pdv{\rho_P}{\rho} &= 1, &
\pdv{\thet_P}{\rho} &= \pdv{\lambda}{\rho},&
\pdv{\zeta_P}{\rho} &= 0,\\
\pdv{\rho_P}{\thet} &=  0,&
\pdv{\thet_P}{\thet} &= 1 + \pdv{\lambda}{\thet},&
\pdv{\zeta_P}{\thet} &= 0,\\
\pdv{\rho_P}{\zeta} &= 0,&
\pdv{\thet_P}{\zeta} &= \pdv{\lambda}{\zeta},&
\pdv{\zeta_P}{\zeta} &= 1.\\
\end{align}
$$

The reciprocal basis vectors are then

$$
\begin{align}
\grad \rho_P &= \grad \rho \\
\grad \thet_P &= \pdv{\thet_P}{\rho}\grad \rho + \pdv{\thet_P}{\thet}\grad \thet + \pdv{\thet_P}{\zeta}\grad \zeta \\
\grad \zeta_P &= \grad \zeta \\
\end{align}
$$

The Jacobian determinant is therefore

$$
\begin{align}
\Jac_P &= \frac{1}{\grad\rho_P \cdot \grad\theta_P \times \grad\zeta_P} \\
&= \frac{1}{\grad\rho \cdot \pdv{\thet_P}{\thet}\grad\thet \times \grad\zeta} \\
&= \frac{\Jac}{\pdv{\thet_P}{\thet}} \\
&= \frac{\Jac}{1 + \pdv{\lambda}{\thet}}. \\
\end{align}
$$

with that we can express the inverse derivatives

$$
\begin{align}
\pdv{\rho}{\rho_P} &= 1, &
\pdv{\rho}{\thet_P} &= 0  ,&
\pdv{\rho}{\zeta_P} &= 0 \\
 \pdv{\thet}{\rho_P} &= -\pdv{\thet_P}{\rho}\left(\pdv{\thet_P}{\thet}\right)^{-1},&
\pdv{\thet}{\thet_P} &= \left(\pdv{\thet_P}{\thet}\right)^{-1}, &
\pdv{\thet}{\zeta_P} &= -\pdv{\thet_P}{\zeta}\left(\pdv{\thet_P}{\thet}\right)^{-1} \\
\pdv{\zeta}{\rho_P} &=  0,&
\pdv{\zeta}{\thet_P} &=  0,&
\pdv{\zeta}{\zeta_P} &=  1\\
\end{align}
$$

The basis vectors are then computed as:

$$
\begin{align}
\vec{e}_{\rho_P} &=&  \erho + &\pdv{\thet}{\rho_P}\ethet  \\
\vec{e}_{\thet_P} &=&  &\pdv{\thet}{\thet_P}\ethet  \\
\vec{e}_{\zeta_P} &=&  &\pdv{\thet}{\zeta_P}\ethet + \ezeta \\
\end{align}
$$
