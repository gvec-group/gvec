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

$\Rightarrow$ `dmod_B_dr`, `dmod_B_dt`, `dmod_B_dz`, `grad_mod_B`

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

Left-handed coordinate systems $\left(s,\alpha,\theta\right)$ and $\left(x^1, x^2, x^3\right) = \left(\sqrt{s}, \frac{\sqrt{s_0}}{a}\alpha,\theta\right)$.

With respect to a left-handed $\left(\rho, \vartheta, \zeta\right)$, the fieldline coordinates are defined as

$$
\begin{align}
s &= \rho^2 &
\rho &= \sqrt{s} \\
\alpha &= \vartheta/\iota - \zeta &
\zeta &= \theta/\iota - \alpha \\
\theta &= \vartheta &
\vartheta &= \theta \\
\end{align}
$$

which can be obtained from the right-handed GVEC angles by flipping the toroidal or poloidal direction.

$$
\begin{align}
\frac{\partial s}{\partial \rho} &= 2\rho &
\frac{\partial \rho}{\partial s} &= \frac{1}{2\sqrt{s}} \\
\frac{\partial \alpha}{\partial \rho} &= -\frac{d\iota}{d\rho} \frac{\vartheta}{\iota^2} &
\frac{\partial \zeta}{\partial s} &= -\frac{d\iota}{ds} \frac{\theta}{\iota^2} \\
\frac{\partial \alpha}{\partial \vartheta} &= \frac{1}{\iota} &
\frac{\partial \zeta}{\partial \theta} &= \frac{1}{\iota} \\
\frac{\partial \alpha}{\partial \zeta} &= -1 &
\frac{\partial \zeta}{\partial \alpha} &= -1 \\
\frac{\partial \theta}{\partial \vartheta} &= 1 &
\frac{\partial \vartheta}{\partial \theta} &= 1 \\
\end{align}
$$

$$
\begin{align}
\nabla s &= 2\rho \nabla\rho \\
\nabla\alpha &= -\frac{d\iota}{d\rho} \frac{\vartheta}{\iota^2} \nabla\rho + \frac{1}{\iota} \nabla\vartheta - \nabla\zeta \\
\nabla\theta &= \nabla\vartheta
\end{align}
$$

$$
\begin{align}
\frac{\partial\modB}{\partial s} &= \frac{1}{2\sqrt{s}} \frac{\partial\modB}{\partial \rho} - \frac{d\iota}{ds}\frac{\theta}{\iota^2} \frac{\partial\modB}{\partial \zeta} \\
\frac{\partial\modB}{\partial \alpha} &= \frac{1}{\iota} \frac{\partial\modB}{\partial \zeta} \\
\frac{\partial\modB}{\partial \theta} &= \frac{\partial\modB}{\partial \vartheta} + \frac{1}{\iota} \frac{\partial\modB}{\partial \zeta}\\
\end{align}
$$
