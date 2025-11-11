# Derivations for computable Quantities

This page collects some derivations for computable quantities:

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
\pdv{\modB}{\alpha} &= \frac{1}{2\modB} \pdv{\modB^2}{\alpha} \\
\nabla \modB &= \pdv{\modB}{\rho} \nabla\rho + \pdv{\modB}{\thet} \nabla\thet + \pdv{\modB}{\zeta} \nabla\zeta \\
\end{align}
$$


for $\alpha\in\left\{\rho,\thet,\zeta\right\}$.

$\Rightarrow$ `dmodB_dr`, `dmodB_dt`, `dmodB_dz`, `grad_modB`

## Derivatives of $\vec{B}$

$$
\begin{align}
\nabla \vec{B} &:= \pdv{\vec{B}}{\vec{x}}
= \sum_\alpha \pdv{\vec{B}}{\alpha} \nabla\alpha
= \pdv{\vec{B}}{\rho} \nabla\rho + \pdv{\vec{B}}{\thet} \nabla\thet + \pdv{\vec{B}}{\zeta} \nabla\zeta \\
\pdv{\vec{B}}{\alpha} &:= \sum_\beta \pdv{B^\beta}{\alpha} \vec{e}_\beta + B^\beta \vec{k}_{\beta\alpha} \\
&= \cancel{\pdv{B^\rho}{\alpha}} \erho + \cancel{B^\rho} \vec{k}_{\rho\alpha}
+ \pdv{B^\thet}{\alpha} \ethet + B^\thet \vec{k}_{\thet\alpha}
+ \pdv{B^\zeta}{\alpha} \ezeta + B^\zeta \vec{k}_{\zeta\alpha} \\
\nabla \vec{B} &= \sum_{\alpha\beta} \left(\pdv{B^\beta}{\alpha} \vec{e}_\beta + B^\beta \vec{k}_{\beta\alpha} \right) \nabla\alpha
\end{align}
$$

for $\alpha,\beta\in\left\{\rho,\thet,\zeta\right\}$.

Note that $\vec{e}_\beta \nabla\alpha$ and $\vec{k}_{\beta\alpha} \nabla\alpha$ are *outer products* (?) and $\nabla\vec{B} \in  \R^3\times\R^3$ (?).

Due to `xarray`'s limitations with multiple dimensions of the same name, $\nabla\vec{B}$ is not directly available as a computable quantitiy, and its components should be used instead.

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
