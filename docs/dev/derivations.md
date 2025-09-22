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
