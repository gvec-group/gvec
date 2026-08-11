# GENE-GIST

:::{admonition} Experimental Feature
:class: caution
This feature is experimental and still being tested.
:::

:::{versionadded} 1.3
:::

This is the interface to the gyrokinetic turbulence code *GENE* using *GIST*-files.
It can be used to generate GIST-files for specified fieldlines from a GVEC equilibrium.

The interface is installed automatically with pyGVEC and available with `pygvec to-gist`.

## Usage

`pygvec to-gist` is a pyGVEC subcommand which is installed as part of the `gvec` package.

The options for `pygvec to-gist` are:
```
$ pygvec to-gist --help
usage: pygvec to-gist [-h] [--rundir RUNDIR] (-s S | -r RHO) [-a ALPHA] [--npol NPOL]
                      [--gridpoints GRIDPOINTS] [--MNfactor MNFACTOR]
                      [-x {auto,none,pol,tor,both}] [-o OUTPUTFILE] [-p]
                      [--projectname PROJECTNAME] [-v | -q]

Produce a GENE-GIST input file from a GVEC state.

options:
  -h, --help            show this help message and exit
  --rundir RUNDIR       GVEC run directory
  -s S                  position of the target flux surface (in normalized toroidal flux, 0 < s <= 1)
  -r RHO, --rho RHO     position of the target flux surface (in square root of the normalized toroidal flux, 0 < rho <= 1)
  -a ALPHA, --alpha ALPHA
                        fieldline label as float or multiple of pi (e.g. '0.0', '2pi', 'pi/2', default 0.0)
  --npol NPOL           number of poloidal turns (default 1)
  --gridpoints GRIDPOINTS
                        number of grid points along the fieldline (default 128)
  --MNfactor MNFACTOR   multiplication factor for the maximum fourier modes for the boozer transform (default 3)
  -x {auto,none,pol,tor,both}, --flip {auto,none,pol,tor,both}
                        flip the poloidal or toroidal direction with respect to GVEC's Boozer coordinates; 'auto' determines the necessary flips to get positive toroidal and poloidal flux (default: 'auto')
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        output file name (default: '{projectname}_s{s}.gist.txt')
  -p, --plot            plot the output quantities ('{projectname}_s{s}.gist.png')
  --projectname PROJECTNAME
                        override the project name for the output files (default: use GVEC state name)
  -v, --verbose         verbosity level: -v for info, -vv for debug
  -q, --quiet           suppress output
```

A typical conversion command would then be:
```bash
pygvec to-gist -s 0.25 -p -v
```
which will look for a GVEC parameter- and statefile in the current directory.

## Derivations for GIST quantities

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

In GVEC, Boozer coordinates are right-handed and follow the directions of the logical (flux aligned) coordinates $\left(\rho,\theta,\zeta\right)$.
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
\hat{\beta} &= \frac{p_0 \mu_0}{B_{ref}^2} \\
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

### coordinate directions

In principle the directions of $\theta_B,\zeta_B$ are arbitrary.
It can however make sense to fix their directions such that certain quantities become positive.

Let us note:

$$
\begin{align}
B_{ref} &= \frac{d\Phi}{d\rho} \frac{1}{a^2\rho_0} \\
\hat\Jac &= \frac{\Jac_B}{\rho_0\iota_0 a^3} \\
\iota &= \frac{d\chi}{d\Phi} \\
\end{align}
$$

In order for $B_{ref}$ to be positive, we need to choose the direction of $\zeta_B$, such that the toroidal magnetic flux $\Phi$ is positive.

Any change in direction also changes the sign of $\iota$ and $\Jac_B$.
$\hat\Jac$ however maintains the sign.

We can therefore change the direction of $\theta_B$, such that the poloidal magnetic flux $\chi$ is also positive.
