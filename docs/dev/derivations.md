# Derivations for computable Quantities

This page collects some derivations for computable quantities:

## Magnetic field $\vec{B}$

$$
\begin{align}
B^\rho &= 0 \\
B^\thet &= \frac{1}{\Jac} \dPhidr \left(\iota-\frac{\partial\lambda}{\partial\zeta} \right) \\
B^\zeta &= \frac{1}{\Jac} \dPhidr \left(1+\frac{\partial\lambda}{\partial\thet} \right) \\
\end{align}
$$

$\Rightarrow$ `B_contra_t`, `B_contra_z`

## Derivatives of $\modB$

$$
\begin{align}
\modB^2 &= B_\thet B^\thet + B_\zeta B^\zeta\\
\modB^2 &= B^\thet B^\thet g_{\thet\thet} + 2 B^\thet B^\zeta g_{\thet\zeta} + B^\zeta B^\zeta g_{\zeta\zeta} \\
\frac{\partial\modB^2}{\partial \alpha} &= 2 B^\thet \frac{\partial B^\thet}{\partial \alpha} g_{\thet\thet} + B^\thet B^\thet \frac{\partial g_{\thet\thet}}{\partial \alpha} \\
                                        &+ 2 \left(\frac{\partial B^\thet}{\partial \alpha} B^\zeta + B^\thet \frac{\partial B^\zeta}{\partial \alpha}\right) g_{\thet\zeta}
                                        + 2 B^\thet B^\zeta \frac{\partial g_{\thet\zeta}}{\partial \alpha} \\
                                        &+ 2 B^\zeta \frac{\partial B^\zeta}{\partial \alpha} g_{\zeta\zeta} + B^\zeta B^\zeta \frac{\partial g_{\zeta\zeta}}{\partial \alpha} \\
\frac{\partial\modB}{\partial \alpha} &= \frac{1}{2\modB} \frac{\partial\modB^2}{\partial \alpha} \\
\nabla \modB &= \frac{\partial\modB}{\partial\rho} \nabla\rho + \frac{\partial\modB}{\partial\thet} \nabla\thet + \frac{\partial\modB}{\partial\zeta} \nabla\zeta \\
\end{align}
$$

for $\alpha\in\left\{\rho,\thet,\zeta\right\}$.

$\Rightarrow$ `dmodB_dr`, `dmodB_dt`, `dmodB_dz`, `grad_modB`

## Derivatives of $\vec{B}$

$$
\begin{align}
\nabla \vec{B} &:= \frac{\partial \vec{B}}{\partial \vec{x}}
= \sum_\alpha \frac{\partial \vec{B}}{\partial \alpha} \nabla\alpha
= \frac{\partial \vec{B}}{\partial \rho} \nabla\rho + \frac{\partial \vec{B}}{\partial \thet} \nabla\thet + \frac{\partial \vec{B}}{\partial \zeta} \nabla\zeta \\
\frac{\partial \vec{B}}{\partial \alpha} &:= \sum_\beta \frac{\partial B^\beta}{\partial \alpha} \vec{e}_\beta + B^\beta \vec{k}_{\beta\alpha} \\
&= \cancel{\frac{\partial B^\rho}{\partial \alpha}} \erho + \cancel{B^\rho} \vec{k}_{\rho\alpha}
+ \frac{\partial B^\thet}{\partial \alpha} \ethet + B^\thet \vec{k}_{\thet\alpha}
+ \frac{\partial B^\zeta}{\partial \alpha} \ezeta + B^\zeta \vec{k}_{\zeta\alpha} \\
\nabla \vec{B} &= \sum_{\alpha\beta} \left(\frac{\partial B^\beta}{\partial \alpha} \vec{e}_\beta + B^\beta \vec{k}_{\beta\alpha} \right) \nabla\alpha \\
\end{align}
$$

for $\alpha,\beta\in\left\{\rho,\thet,\zeta\right\}$.

Note that $\vec{e}_\beta \nabla\alpha$ and $\vec{k}_{\beta\alpha} \nabla\alpha$ are *outer products* (?) and $\nabla\vec{B} \in  \R^3\times\R^3$ (?).

Due to `xarray`'s limitations with multiple dimensions of the same name, $\nabla\vec{B}$ is not directly available as a computable quantitiy, and its components should be used instead.

$\Rightarrow$ `dB_dr`, `dB_dt`, `dB_dz`

---

$$
\begin{align}
\lVert\nabla\vec{B}\rVert^2 &:= \nabla\vec{B} : \nabla\vec{B} \\
&= \sum_{ij} \sum_{\alpha\beta} \left(\frac{\partial \vec{B}}{\partial \alpha}\right)_i \left(\nabla\alpha\right)_j \left(\frac{\partial \vec{B}}{\partial \beta}\right)_i \left(\nabla\beta\right)_j
\end{align}
$$
???

for $\alpha,\beta\in\left\{\rho,\thet,\zeta\right\}$.

$\Rightarrow$ `grad_B_norm2`

---

From

$$
\begin{align}
B^\thet &= \frac{1}{\Jac} \dPhidr \left(\iota-\frac{\partial\lambda}{\partial\zeta} \right) \\
B^\zeta &= \frac{1}{\Jac} \dPhidr \left(1+\frac{\partial\lambda}{\partial\thet} \right) \\
\end{align}
$$

it follows

$$
\begin{align}
\frac{\partial B^\thet}{\partial \rho} &=
    \frac{1}{\Jac} \frac{d^2\Phi}{d\rho^2} \left(\iota-\frac{\partial\lambda}{\partial\zeta} \right)
    + \frac{1}{\Jac} \dPhidr \left(\frac{d\iota}{d\rho}-\frac{\partial^2\lambda}{\partial\rho\partial\zeta} \right)
    - \frac{\partial\Jac}{\partial\rho} \frac{1}{\Jac^2} \dPhidr \left(\iota-\frac{\partial\lambda}{\partial\zeta} \right) \\
\frac{\partial B^\thet}{\partial \thet} &=
    - \frac{1}{\Jac} \dPhidr \frac{\partial^2\lambda}{\partial\thet\partial\zeta}
    - \frac{\partial\Jac}{\partial\thet} \frac{1}{\Jac^2} \dPhidr \left(\iota-\frac{\partial\lambda}{\partial\zeta} \right) \\
\frac{\partial B^\thet}{\partial \zeta} &=
    - \frac{1}{\Jac} \dPhidr \frac{\partial^2\lambda}{\partial\zeta^2}
    - \frac{\partial\Jac}{\partial\zeta} \frac{1}{\Jac^2} \dPhidr \left(\iota-\frac{\partial\lambda}{\partial\zeta} \right) \\
\frac{\partial B^\zeta}{\partial \rho} &=
    \frac{1}{\Jac} \frac{d^2\Phi}{d\rho^2} \left(1+\frac{\partial\lambda}{\partial\thet} \right)
    + \frac{1}{\Jac} \dPhidr \frac{\partial^2\lambda}{\partial\rho\partial\thet}
    - \frac{\partial\Jac}{\partial\rho} \frac{1}{\Jac^2} \dPhidr \left(1+\frac{\partial\lambda}{\partial\thet} \right) \\
\frac{\partial B^\zeta}{\partial \thet} &=
    \frac{1}{\Jac} \dPhidr \frac{\partial^2\lambda}{\partial^2\thet}
    - \frac{\partial\Jac}{\partial\thet} \frac{1}{\Jac^2} \dPhidr \left(1+\frac{\partial\lambda}{\partial\thet} \right) \\
\frac{\partial B^\zeta}{\partial \zeta} &=
    \frac{1}{\Jac} \dPhidr \frac{\partial^2\lambda}{\partial\thet\partial\zeta}
    - \frac{\partial\Jac}{\partial\zeta} \frac{1}{\Jac^2} \dPhidr \left(1+\frac{\partial\lambda}{\partial\thet} \right) \\
\end{align}
$$

$\Rightarrow$ `dB_contra_t_dr`, `dB_contra_t_dt`, `dB_contra_t_dz`, `dB_contra_z_dr`, `dB_contra_z_dt`, `dB_contra_z_dz`


## export GIST files for GENE

?-handed coordinate systems $\left(s,\alpha,\theta\right)$ and $\left(x^1, x^2, x^3\right) = \left(a \sqrt{s}, a\frac{\sqrt{s_0}}{q_0}\alpha, a\theta\right)$.

With respect to a ?-handed $\left(\rho, \vartheta, \zeta\right)$, the fieldline coordinates are defined as

$$
\begin{align}
s &= \rho^2 &
\rho &= \sqrt{s} \\
\alpha &= \sigma\vartheta/\iota - \sigma\zeta &
\zeta &= \theta/\iota - \sigma\alpha \\
\theta &= \sigma\vartheta &
\vartheta &= \sigma\theta \\
\end{align}
$$

<!-- which can be obtained from the right-handed GVEC angles by flipping the toroidal or poloidal direction. -->

$$
\begin{align}
\frac{\partial s}{\partial \rho} &= 2\rho &
\frac{\partial \rho}{\partial s} &= \frac{1}{2\sqrt{s}} \\
\frac{\partial \alpha}{\partial \rho} &= -\sigma \frac{d\iota}{d\rho} \frac{\vartheta}{\iota^2} &
\frac{\partial \zeta}{\partial s} &= -\frac{d\iota}{ds} \frac{\theta}{\iota^2} \\
\frac{\partial \alpha}{\partial \vartheta} &= \sigma\frac{1}{\iota} &
\frac{\partial \zeta}{\partial \theta} &= \frac{1}{\iota} \\
\frac{\partial \alpha}{\partial \zeta} &= -\sigma &
\frac{\partial \zeta}{\partial \alpha} &= -\sigma \\
\frac{\partial \theta}{\partial \vartheta} &= \sigma &
\frac{\partial \vartheta}{\partial \theta} &= \sigma \\
\end{align}
$$

$$
\begin{align}
\nabla s &= 2\rho \nabla\rho \\
\nabla\alpha &= \sigma\left( -\frac{d\iota}{d\rho} \frac{\vartheta}{\iota^2} \nabla\rho + \frac{1}{\iota} \nabla\vartheta - \nabla\zeta \right)\\
\nabla\theta &= \sigma\nabla\vartheta \\

\nabla x^1 &= \frac{a}{2\sqrt{s}} \nabla s \\
\nabla x^2 &= \frac{a\sqrt{s_0}}{q_0} \nabla\alpha \\
\nabla x^3 &= a\nabla\theta \\
\end{align}
$$

$$
\begin{align}
\Jac^{-1}_{s\alpha\theta} &:= \nabla s \cdot \nabla\alpha \times \nabla\theta \\
&= 2\rho \nabla\rho \cdot \sigma \left(-\frac{d\iota}{d\rho} \frac{\vartheta}{\iota^2} \nabla\rho + \frac{1}{\iota} \nabla\vartheta - \nabla\zeta\right) \times \sigma \nabla\vartheta \\
&= -2\rho \nabla\rho \cdot \nabla\zeta \times \nabla\vartheta \\
&= 2\rho \Jac^{-1}_{\rho\vartheta\zeta} \\

\Jac^{-1}_{123} &:= \nabla x^1 \cdot \nabla x^2 \times \nabla x^3 \\
&= \frac{a^3 \sqrt{s_0}}{2\sqrt{s}q_0} \Jac^{-1}_{s\alpha\theta} \\
&= \frac{a^3}{2q_0} \Jac^{-1}_{s\alpha\theta} \\
\end{align}
$$

$$
\begin{align}
\frac{\partial Q}{\partial s} &= \frac{1}{2\sqrt{s}} \frac{\partial Q}{\partial \rho} - \frac{d\iota}{ds}\frac{\theta}{\iota^2} \frac{\partial Q}{\partial \zeta} \\
\frac{\partial Q}{\partial \alpha} &= -\sigma \frac{\partial Q}{\partial \zeta} \\
\frac{\partial Q}{\partial \theta} &= \sigma \frac{\partial Q}{\partial \vartheta} + \frac{1}{\iota} \frac{\partial Q}{\partial \zeta}\\
\end{align}
$$

---

GVEC Boozer coordinates $\left(\rho,\theta_B,\zeta_B\right)$ are right-handed, straight-fieldline coordinates.
To transform into fieldline coordinates $\left(\rho_\alpha,\alpha,\phi_\alpha\right)$:

$$
\begin{align}
\rho_\alpha &= \rho_B \\
\alpha &= \theta_B - \iota \zeta_B \\
\phi_\alpha &= \theta_B
\end{align}
$$

$$
\begin{align}
\rho_B &= \rho_\alpha \\
\theta_B &= \phi_\alpha \\
\zeta_B &= \left(\phi_\alpha - \alpha\right) / \iota
\end{align}
$$

$$
\begin{align}
\pdv{\rho_\alpha}{\rho_B} &= 1 \\
\pdv{\alpha}{\rho_B} &= - \frac{d\iota}{d\rho}\zeta \\
\pdv{\alpha}{\theta_B} &= 1 \\
\pdv{\alpha}{\zeta_B} &= \iota \\
\pdv{\phi_\alpha}{\theta_B} &= 1
\end{align}
$$

$$
\begin{align}
\pdv{\rho_B}{\rho_\alpha} &= 1 \\
\pdv{\theta_B}{\phi_\alpha} &= 1 \\
\pdv{\zeta_B}{\rho_\alpha} &= \left(\alpha - \phi_\alpha\right) \frac{d\iota}{d\rho} / \iota^2
\pdv{\zeta_B}{\alpha} &= -1 / \iota \\
\pdv{\zeta_B}{\phi_\alpha} &= 1 / \iota
\end{align}
$$

$$
\begin{align}
\nabla\rho_\alpha &= \nabla\rho_B \\
\nabla\alpha &= - \frac{d\iota}{d\rho}\zeta \nabla\rho_B + \nabla\theta_B + \iota \nabla\zeta_B \\
\nabla\phi_\alpha &= \nabla\theta_B
\end{align}
$$

$$
\begin{align}
\vec{e}_{\rho_\alpha} &= \vec{e}_{\rho_B} + \left(\alpha - \phi_\alpha\right) \frac{d\iota}{d\rho} / \iota^2 \vec{e}_{\zeta_B} \\
\vec{e}_{\alpha} &= - 1 / \iota \vec{e}_{\zeta_B} \\
\vec{e}_{\phi_\alpha} &= \vec{e}_{\theta_B} + 1 / \iota \vec{e}_{\zeta_B}
\end{align}
$$

$$
\begin{align}
B^{\rho_B} &= 0 \\
B^{\theta_B} &= \frac{1}{\Jac_{B}}\frac{d\Phi}{d\rho}\iota \\
B^{\zeta_B} &= \frac{1}{\Jac_{B}}\frac{d\Phi}{d\rho}
\end{align}
$$

$$
\begin{align}
\vec{B} &= B^{\theta_B} \vec{e}_{\theta_B} + B^{\zeta_B} \vec{e}_{\zeta_B}
&= \frac{1}{\Jac_{B}} \frac{d\Phi}{d\rho} \left( \iota\vec{e}_{\theta_B} + \vec{e}_{\zeta_B} \right)
B^{\rho_\alpha} &= 0 \\
B^{\alpha} &= \\
B^{\phi_\alpha} &=
\end{align}
$$

---

## 2nd try

### Boozer coordinates

Boozer coordinates $\left(\rho,\theta_B,\zeta_B\right)$ are a particular kind of straight-fieldline coordinates, where the following relations hold for the magnetic field $\vec{B}$:

$$
\begin{align}
\vec{B}\cdot\vec{e}^{\rho_B} = B^{\rho_B} &= 0 \\
\vec{B}\cdot\vec{e}^{\theta_B} = B^{\theta_B} &= \frac{1}{\Jac_{B}}\frac{d\Phi}{d\rho}\iota \\
\vec{B}\cdot\vec{e}^{\zeta_B} = B^{\zeta_B} &= \frac{1}{\Jac_{B}}\frac{d\Phi}{d\rho} \\
\vec{B}\cdot\nabla\theta_B = B_{\theta_B}(\rho_B) \\
\vec{B}\cdot\nabla\zeta_B = B_{\zeta_B}(\rho_B) \\
\end{align}
$$

Here $\rho_B\in[0,1]$ is a radial, radius-like coordinate, typically chosen such that $\rho^2$ is the normalized toroidal flux.
$\theta_B, \zeta_B \in [0,2\pi)$ are the poloidal and toroidal periodic, angle-like coordinates.

### changing coordinate directions

In GVEC Boozer coordinates are right-handed and follow the directions of the logical (flux aligned) coordinates $\left(\rho,\theta,\zeta\right)$.
The directions of the two angular coordinates are arbitrary and need to be flipped in certain cases.

Flipping toroidally $\zeta \to -\zeta$ changes $\Phi \to -\Phi, \iota \to -\iota, \Jac \to -\Jac$.

Flipping poloidally $\theta \to -\theta$ changes $\chi \to -\chi, \iota \to -\iota, \Jac \to -\Jac$.

The geometry also needs to be flipped accordingly and all derived quantities change, but the formulas stay the same.

### Fieldline coordinates

The transformation to the fieldline coordinates $\left(\rho_\alpha,\alpha,\phi_\alpha\right)$ can be defined as:

$$
\begin{align}
\rho_\alpha &= \rho_B \\
\alpha &= \theta_B / \iota - \zeta_B \\
\phi_\alpha &= \theta_B
\end{align}
$$

$$
\begin{align}
\pdv{\rho_\alpha}{\rho_B} &= 1 &
\pdv{\rho_\alpha}{\theta_B} &= 0 &
\pdv{\rho_\alpha}{\zeta_B} &= 0 \\
\pdv{\alpha}{\rho_B} &= -\frac{\theta_B}{\iota^2} \frac{d\iota}{d\rho} &
\pdv{\alpha}{\theta_B} &= \frac{1}{\iota} &
\pdv{\alpha}{\zeta_B} &= -1 \\
\pdv{\phi_\alpha}{\rho_B} &= 0 &
\pdv{\phi_\alpha}{\theta_B} &= 1 &
\pdv{\phi_\alpha}{\zeta_B} &= 0 \\
\end{align}
$$

The inverse transformation is then:

$$
\begin{align}
\rho_B &= \rho_\alpha \\
\theta_B &= \phi_\alpha \\
\zeta_B &= \phi_\alpha / \iota - \alpha
\end{align}
$$

$$
\begin{align}
\pdv{\rho_B}{\rho_\alpha} &= 1 &
\pdv{\rho_B}{\alpha} &= 0 &
\pdv{\rho_B}{\phi_\alpha} &= 0 \\
\pdv{\theta_B}{\rho_\alpha} &= 0 &
\pdv{\theta_B}{\alpha} &= 0 &
\pdv{\theta_B}{\phi_\alpha} &= 1 \\
\pdv{\zeta_B}{\rho_\alpha} &= -\frac{\phi_\alpha}{\iota^2} \frac{d\iota}{d\rho} &
\pdv{\zeta_B}{\alpha} &= -1 &
\pdv{\zeta_B}{\phi_\alpha} &= \frac{1}{\iota} \\
\end{align}
$$

The reciprocal and tangential basis vectors are then:

$$
\begin{align}
\nabla\rho_\alpha &= \nabla\rho_B \\
\nabla\alpha &= -\frac{\theta_B}{\iota^2} \frac{d\iota}{d\rho} \nabla\rho_B + \frac{1}{\iota} \nabla\theta_B - \nabla\zeta_B \\
\nabla\phi_\alpha &= \nabla\theta_B
\end{align}
$$

$$
\begin{align}
\vec{e}_{\rho_\alpha} &= \vec{e}_{\rho_B} -\frac{\phi_\alpha}{\iota^2} \frac{d\iota}{d\rho} \vec{e}_{\zeta_B} \\
\vec{e}_{\alpha} &= -\vec{e}_{\zeta_B} \\
\vec{e}_{\phi_\alpha} &= \vec{e}_{\theta_B} + \frac{1}{\iota} \vec{e}_{\zeta_B}
\end{align}
$$

And the Jacobian determinant is

$$
\begin{align}
\Jac_\alpha &= \vec{e}_{\rho_\alpha} \cdot \vec{e}_\alpha \times \vec{e}_{\phi_\alpha} \\
&= \left(\vec{e}_{\rho_B} -\frac{\phi_\alpha}{\iota^2} \frac{d\iota}{d\rho} \vec{e}_{\zeta_B}\right)
\cdot \left(-\vec{e}_{\zeta_B}\right)
\times \left( \vec{e}_{\theta_B} + \frac{1}{\iota} \vec{e}_{\zeta_B} \right) \\
&= \Jac_B.
\end{align}
$$

The magnetic field takes the form

$$
\begin{align}
\vec{B} &= B^{\theta_B} \vec{e}_{\theta_B} + B^{\zeta_B} \vec{e}_{\zeta_B} \\
&= \frac{1}{\Jac_{B}} \frac{d\Phi}{d\rho} \left( \iota\vec{e}_{\theta_B} + \vec{e}_{\zeta_B} \right) \\
&= \frac{\iota}{\Jac_{B}} \frac{d\Phi}{d\rho} \vec{e}_{\phi_\alpha} \\
B^{\rho_\alpha} &= 0 \\
B^{\alpha} &= 0 \\
B^{\phi_\alpha} &= \frac{\iota}{\Jac_{B}} \frac{d\Phi}{d\rho} = \frac{1}{\Jac_B} \frac{d\chi}{d\rho}, \\
\end{align}
$$

which shows that the magnetic field is tangential to the fieldline defined by constant $\rho,\alpha$.

The derivatives of some quantity $Q$ in fieldline coordinates are then:

$$
\begin{align}
\pdv{Q}{\rho_\alpha} &= \pdv{Q}{\rho_B} -\frac{\phi_\alpha}{\iota^2} \frac{d\iota}{d\rho} \pdv{Q}{\zeta_B} \\
\pdv{Q}{\alpha} &= -\pdv{Q}{\zeta_B} \\
\pdv{Q}{\phi_\alpha} &= \pdv{Q}{\theta_B} + \frac{1}{\iota} \pdv{Q}{\zeta_B}
\end{align}
$$

### GIST quantities

The GIST coordinate system is defined as

$$\left(x^1,x^2,x^3\right) = \left(\rho, \rho_0 \iota_0 \alpha, \phi_\alpha\right),$$

where

$$\iota_0 = \iota(\rho_0)$$

and $\rho_0$ is the position of the flux surface of interest.

The normalized GIST coordinate system is defined as

$$\left(\hat{x}^1,\hat{x}^2,\hat{x}^3\right) = \left(ax^1, ax^2, ax^3\right).$$

The reciprocal basis vectors in the normalized GIST coordinate system are then

$$
\begin{align}
\nabla\hat{x}^1 &= \hat\nabla x^1 = a \nabla\rho \\
\nabla\hat{x}^2 &= \hat\nabla x^2 = a \rho_0\iota_0 \nabla\alpha \\
\nabla\hat{x}^3 &= \hat\nabla x^3 = a \nabla\phi_\alpha,
\end{align}
$$

with the normalized gradient $\hat\nabla = a\nabla$.
The associated components of the metric tensor, defined as

$$g^{ij} = \nabla i \cdot \nabla j$$

are

$$
\begin{align}
\hat{g}^{11} &= a^2 g^{\rho\rho} &
\hat{g}^{12} &= a^2 \rho_0 \iota_0 g^{\rho\alpha} \\
& &
\hat{g}^{22} &= a^2 \rho_0^2 \iota_0^2 g^{\alpha\alpha}.
\end{align}
$$

The Jacobian determinant is

$$\hat\Jac = \frac{\Jac_\alpha}{\rho_0\iota_0 a^3}$$

and the magnetic field is defined by:

$$
\begin{align}
\vec{B} &:= B_{ref} \hat\nabla x^1 \times \hat\nabla x^2 \\
&= B_{ref} a^2\rho_0\iota_0 \nabla\rho \times \nabla\alpha \\
&= B_{ref} \frac{a^2\rho_0\iota_0}{\Jac_\alpha} \vec{e}_{\phi_\alpha} \\
\vec{B} &= \frac{\iota}{\Jac_{B}} \frac{d\Phi}{d\rho} \vec{e}_{\phi_\alpha} \\
B_{ref} &= \frac{d\Phi}{d\rho} \frac{1}{a^2\rho_0} \\
\end{align}
$$

Note that with the usual convention of $\Phi = \Phi_{edge}\rho^2$, it follows

$$\frac{d\Phi}{d\rho} = 2\rho\Phi_{edge} \quad\Rightarrow\quad B_{ref} = \frac{2\Phi_{edge}}{a^2}.$$

Some quantities of interest are (see relevant literature on GENE for definitions and motivation):

$$
\begin{align}
L_1 &= - \frac{1}{B_{ref}\iota_0\rho_0} \left[ \pdv{\modB}{\alpha}
+ \frac{g^{\rho\rho}g^{\alpha\phi}-g^{\rho\alpha}g^{\rho\phi}}{g^{\rho\rho}g^{\alpha\alpha}-g^{\rho\alpha}g^{\rho\alpha}}
\pdv{\modB}{\phi_\alpha} \right] \\
L_2 &= \frac{1}{B_{ref}} \left[ \pdv{\modB}{\rho_\alpha}
+ \frac{g^{\alpha\alpha}g^{\rho\phi}-g^{\rho\alpha}g^{\alpha\phi}}{g^{\rho\rho}g^{\alpha\alpha}-g^{\rho\alpha}g^{\rho\alpha}}
\pdv{\modB}{\phi_\alpha} \right] \\
\hat{b} &= \frac{\modB}{B_{ref}} \\
\hat{b}_\phi &= \frac{1}{B_{ref}} \pdv{\modB}{\phi_\alpha} \\
q_0 &= \frac{1}{\iota_0} \\
\hat{s} &= \rho_0\iota_0 \frac{d\iota}{d\rho} \\
\hat{\beta} &= \frac{p_0 \mu_0}{B_{ref}^2}
\hat{p} &= -\frac{a^4 \mu_0}{2 \Phi_{edge} |\Phi_{edge}|} \frac{dp}{d\rho} \\
\end{align}
$$

For comparison with VMEC quantities, note that

$$
\pdv{Q}{\rho} = \frac{ds}{d\rho} \pdv{Q}{s} = 2\rho \pdv{Q}{s}.
$$

The GIST file then contains the Fortran namelist `parameters` with:
* `s0` $= \rho_0^2$
* `my_dpdx` $= \hat{p}$
* `q0` $= q_0$
* `shat` $= \hat{s}$
* `gridpoints`
* `n_pol` $= n_{pol}$
* `Lref` $= a$
* `Bref` $= B_{ref}$
* `n0_global` $= N_{FP}$
* `beta` $= \hat{\beta}$

and the columns:
1) $\hat{g}^{11}$
2) $\hat{g}^{12}$
3) $\hat{g}^{22}$
4) $\hat{b}$
5) $\hat\Jac$
6) $L_2$
7) $-L_1$
8) $\hat{b}_\phi$

which are evaluated at `gridpoints` equidistantly spaced points in $\phi_\alpha\in\left[-n_{pol}\pi,n_{pol}\pi\right)$ for a fixed $\rho=\rho_0$ and $\alpha$.
