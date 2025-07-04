# CAS3D

:::{warning}
This feature is experimental and likely does not yet produce the expected results.
:::

This is the interface to the MHD stability code *CAS3D* [^CAS3D].
It can be used to convert a GVEC equilibrium (parameterfile & statefile) into specialized netCDF file that can be read by CAS3D.

The interface is installed automatically with pyGVEC and available with `pygvec to-cas3d`.

## Usage

[Install](install.md) GVEC with python bindings as normal.
The `pygvec` executable with the `pygvec to-cas3d` subcommand will be added to the binary folder (e.g. `venv/bin/`, which is also added to `$PATH`).
Then you should simply be able to execute:
```bash
pygvec to-cas3d parameter.ini GVEC_State_0000_00000100.dat GVEC_BoozFT-CAS3D.nc --ns 3 --MN_out 10 10 --stellsym
```
which will produce the `GVEC_BoozFT-CAS3D.nc` netCDF file to be used with CAS3D.

The other options for `pygvec to-cas3d` are:
```bash
$ pygvec to-cas3d --help
usage: pygvec to-cas3d [-h] --ns NS --MN_out MN_OUT MN_OUT [--MN_booz MN_BOOZ MN_BOOZ] [--sampling SAMPLING] [--stellsym] [--pointwise POINTWISE] parameterfile statefile outputfile

Convert a GVEC equilibrium to be used in CAS3D

positional arguments:
  parameterfile         input GVEC parameter-file
  statefile             input GVEC state-file
  outputfile            output netCDF file

options:
  -h, --help            show this help message and exit
  --ns NS               number of flux surfaces (equally spaced in s=rho^2) (required)
  --MN_out MN_OUT MN_OUT
                        maximum fourier modes in the output (M, N) (required)
  --MN_booz MN_BOOZ MN_BOOZ
                        maximum fourier modes for the boozer transform (M, N)
  --sampling SAMPLING   sampling factor for the fourier transform and surface reparametrization
  --stellsym            filter the output for stellarator symmetry
  --pointwise POINTWISE
                        output pointwise data to a separate file
```

Note that per default `pygvec to-cas3d` will be parallelized with OpenMP and you can set the number of threads (for example 20) before running the converter as an environment variable, with
```bash
export OMP_NUM_THREADS=20
```

:::{note}
The Boozer Transform is limiting both the performance and memory usage of the converter. A higher resolution (e.g. `--MN_out 16 32 --sampling 4`) can require significant amounts of RAM (e..g more than 32GiB).
:::

## Interface

The netCDF export for CAS3D contains a number of flux surfaces, equidistantly spaced in the normalized toroidal flux. On each surface the magnetic field, metric tensor and second fundamental form is described using fourier coefficients. The specific quantities and their definitions are found below.

### Coordinates

CAS3D uses a flux aligned coordinate system with a radial coordinate $s\in[0,1]$, proportional to the normalized toroidal flux, and two angular coordinates $\vartheta,\zeta\in[0,1]$.
Note that the first radial position in the output is not the $s=0$ surface, but $s=10^{-8}$ to avoid the singularity at the axis.
The two angles are *Boozer-straight-fieldline-angles*, that is the components of the magnetic field $B_\vartheta(s),B_\zeta(s),\mathcal{J}B^\vartheta(s),\mathcal{J}B^\zeta(s)$ are constant on each fluxsurface.
The $\vartheta=0$ surface is on the outward side of the device and $\vartheta,\zeta$ both increase in the counter-clockwise direction when viewed from the front or above respectively.
This makes the $(s,\vartheta,\zeta)$ coordinate system left-handed ($\mathcal{J} < 0$) (and is similar to a COCOS[^COCOS] number of 3/5/13/15).

:::{warning}
The placement of the $\vartheta=0$ contour is not guaranteed by the converter right now. The expected position of that contour for a generalized axis-following frame is however not clearly defined.
:::

:::{warning}
ToDo: this derivation assumes that $\zeta=[0,1]$ on the whole device for CAS3D, that is that $\zeta=[0,\frac{1}{N_{FP}}]$ on one field period.
:::

In terms of the GVEC *Boozer-coordinates* $(\rho,\vartheta_B,\zeta_B)$, the CAS3D flux coordinates are:

$$ s = \rho^2, \qquad \vartheta = -\sigma\frac{1}{2\pi} \vartheta_B, \qquad \zeta = \sigma\frac{1}{2\pi} \zeta_B, $$

where $\sigma=1$ if the GVEC reference frame ($h$-map) is defined with counter-clockwise $\zeta$ direction and $\sigma=-1$ otherwise.
The derivatives are then:

$$ \frac{ds}{d\rho} = 2\rho, \qquad \frac{d\vartheta}{d\vartheta_B} = -\sigma\frac{1}{2\pi}, \qquad \frac{d\zeta}{d\zeta_B} = \sigma\frac{1}{2\pi}.$$

The derived quantities are then transformed as shown [here](./coordinate-conventions.md#different-conventions)

### field-periodic representation

If the number of field periods is >1, the flux surface positions $x,y,z$ in cartesian coordinates are not periodic on one field-period. A transformation to field-periodic variables $\hat{x},\hat{y},\hat{z}$ is possible, using the toroidal parameterization $\zeta\in[0,2\pi]$, see Guiliani et. al. [^xhat]. The transform forflux surface with $\vartheta$ poloidal  and $\zeta$ toroidal parameterization is then defined as:

$$
\begin{align}
\hat{x}(\vartheta,\zeta) &:=\quad x(\vartheta,\zeta)\cos(2\pi\zeta)+y(\vartheta,\zeta)\sin(2\pi\zeta), \\
\hat{y}(\vartheta,\zeta) &:=-x(\vartheta,\zeta)\cos(2\pi\zeta)+y(\vartheta,\zeta)\cos(2\pi\zeta), \\
\hat{z}(\vartheta,\zeta) &:=\quad z(\vartheta,\zeta)
\end{align}
$$

and its inverse

$$
\begin{align}
x(\vartheta,\zeta) &:= \hat{x}(\vartheta,\zeta)\cos(2\pi\zeta) - \hat{y}(\vartheta,\zeta)\sin(2\pi\zeta), \\
y(\vartheta,\zeta) &:= \hat{x}(\vartheta,\zeta)\cos(2\pi\zeta) + \hat{y}(\vartheta,\zeta)\cos(2\pi\zeta), \\
z(\vartheta,\zeta) &:= \hat{z}(\vartheta,\zeta),
\end{align}
$$

With stellarator-symmetry, $\hat{x}$ is an even periodic function on a field period and $\hat{y},\hat{z}$ are odd periodic functions on a field period.

<!--- References -->

[^CAS3D]: C. Schwab; Ideal magnetohydrodynamics: Global mode analysis of three‐dimensional plasma configurations. *Phys. Fluids B* 1 September 1993; 5 (9): 3195–3206. [DOI:10.1063/1.860656](https://doi.org/10.1063/1.860656)

[^COCOS]: O. Sauter and S. Yu. Medvedev; Tokamak Coordinate Conventions: COCOS. *Comput. Phys. Commun.* 2012; [DOI:10.1016/j.cpc.2012.09.010](https://doi.org/10.1016/j.cpc.2012.09.010)

[^xhat]: A. Guiliani, F. Wechsung, M. Landreman, G. Stadler and A. Cerfon; Direct computation of magnetic surfaces in Boozer coordinates and coil optimization for quasi-symmetry. *J. Plasma Phys.* 2022; [DOI:10.1017/S0022377822000563](https://doi.org/10.1017/S0022377822000563)
