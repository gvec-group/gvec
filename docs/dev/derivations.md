# Derivations for computable Quantities

This page collects some derivations for computable quantities:

## Derivatives of `mod_B`

$$
\begin{align}
    \modB &= B_\thet B^\thet + B_\zeta B^\zeta \\
    \modB &= B^\thet B^\thet g_{\thet\thet} + 2 B^\thet B^\zeta g_{\thet\zeta} + B^\zeta B^\zeta g_{\zeta\zeta}
    \frac{\partial\modB}{\partial \alpha} &= 2 B^\thet \frac{\partial B^\thet}{\partial \alpha} g_{\thet\thet} + B^\thet B^\thet \frac{\partial g_{\thet\thet}}{\partial \alpha} \\
                                         &+ 2 \left(\frac{\partial B^\thet}{\partial \alpha} B^\zeta + B^\thet \frac{\partial B^\zeta}{\partial \alpha}\right) g_{\thet\zeta}
                                         + 2 B^\thet B^\zeta \frac{\partial g_{\thet\zeta}}{\partial \alpha} \\
                                         &+ 2 B^\zeta \frac{\partial B^\zeta}{\partial \alpha} g_{\zeta\zeta} + B^\zeta B^\zeta \frac{\partial g_{\zeta\zeta}}{\partial \alpha} \\
\end{align}
$$

for $\alpha\in\left\{\rho,\thet,\zeta\right\}$.
