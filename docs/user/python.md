# Python Bindings
GVEC has *Python* bindings (referred to as *pyGVEC*) to run gvec and evaluate gvec equilibria from *Python*, while relying on the compiled *Fortran* functions to perform the actual computations.
This ensures post-processing is consistent with computing the equilibrium and also improves performance.


## Installation

Please follow the instructions for installing [**gvec and its python bindings**](install.md).


## Python interface to gvec state

The low-level interface is provided with a `State` class, which can be instantiated from a given *parameter-* and *statefile*:
```python
from gvec import State

with State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    ...
```

A high-level interface for evaluations is provided in form of a [xarray](https://docs.xarray.dev/) `Dataset` that can be generated with the `Evaluations` factory:
```python
from gvec import State, Evaluations

with State("parameter.ini", "EXAMPLE_State_0001_00001000.dat") as state:
    ev = Evaluations(rho=[0.1, 0.5, 0.9], theta="int", zeta="int", state=state)
    state.compute(ev, "B")
```
Here the additional arguments configure the points in the radial, poloidal and toroidal direction respectively, with `"int"` selecting the integration points that are used internally by GVEC.
The `ev` object is an instance of the `xarray.Dataset` and the individual `xarray.DataArray`s can then be accessed using `ev.B` or `ev["B"]`.
A `xarray.Dataset` closely mirrors the structure of netCDF, grouping several variables with named dimensions and coordinates as well as metadata.
The `state.compute` method can be used to compute various quantities that are then added to the `Dataset`.
Here *pyGVEC* takes care of computing all the required intermediate quantities, which are also added to the `Dataset`.

## Troubleshooting

* `401 Error` when trying to install from the GitLab package registry:
  * `WARNING: 401 Error, Credentials not correct for https://gitlab.mpcdf.mpg.de/api/v4/projects/1395/packages/pypi/simple/gvec/`
  * you likely missing the [personal access token](https://gitlab.mpcdf.mpg.de/help/user/profile/personal_access_tokens) required for this installation method (for now)
* `Cannot open include file 'netcdf.inc'` during installation
  * `gvec/src/vmec/vmec_readin.f90(494): error #5102: Cannot open include file 'netcdf.inc'`
  * gvec likely did not recognize your hostname and cannot find netCDF. Possible fixes:
    * explicitly tell gvec which preset to use, e.g. `pip install git+ssh://git@gitlab.mpcdf.mpg.de/gvec-group/gvec.git --config-settings=cmake.define.CMAKE_HOSTNAME=tokp`
    * disable netCDF: `pip install git+ssh://git@gitlab.mpcdf.mpg.de/gvec-group/gvec.git --config-settings=cmake.define.LINK_GVEC_TO_NETCDF=Off`
* `INTEL_MKL_ERROR` when trying to use pyGVEC
  * `INTEL MKL ERROR: /usr/lib/x86_64-linux-gnu/libmkl_avx2.so: undefined symbol: mkl_sparse_optimize_bsr_trsm_i8.`
  * `Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so or libmkl_def.so.`
  * you can try `export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libmkl_def.so:/usr/lib/x86_64-linux-gnu/libmkl_avx2.so:/usr/lib/x86_64-linux-gnu/libmkl_core.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_thread.so:/usr/lib/x86_64-linux-gnu/libiomp5.so`
* `undefined reference to EVP_KDF_CTX` when trying to import pyGVEC and using `conda`
  * `/lib64/libk5crypto.so.3: undefined reference to EVP_KDF_CTX_new_id@OPENSSL_1_1_1b`
  * this can be caused by the `conda` environment conflicting with system libraries. You can try: `export LD_PRELOAD="/usr/lib64/libcrypto.so /usr/lib64/libssl.so"`

## Available Quantities for Evaluation
The following table contains the quantities that can be evaluated with the python bindings.

:::{note}
This table is not automatically generated (yet) and might be out of date.
It was last generated for `0.5.1.dev109+gbbfb3a57` on 2024-12-16.
:::

|      label       |                                    long name                                    |                           symbol                           |
| ---------------- | ------------------------------------------------------------------------------- | ---------------------------------------------------------- |
|       `B`        |                                 magnetic field                                  |                        $\mathbf{B}$                        |
|   `B_contra_t`   |                             poloidal magnetic field                             |                         $B^\theta$                         |
|  `B_contra_t_B`  |           contravariant poloidal magnetic field in Boozer coordinates           |                       $B^{\theta_B}$                       |
|   `B_contra_z`   |                             toroidal magnetic field                             |                         $B^\zeta$                          |
|  `B_contra_z_B`  |           contravariant toroidal magnetic field in Boozer coordinates           |                       $B^{\zeta_B}$                        |
|  `B_theta_avg`   |                         average poloidal magnetic field                         |                   $\overline{B_\theta}$                    |
|   `B_zeta_avg`   |                         average toroidal magnetic field                         |                    $\overline{B_\zeta}$                    |
|       `F`        |                                    MHD force                                    |                            $F$                             |
|    `F_r_avg`     |                              radial force balance                               |                    $\overline{F_\rho}$                     |
|     `I_pol`      |                                poloidal current                                 |                         $I_{pol}$                          |
|     `I_tor`      |                                toroidal current                                 |                         $I_{tor}$                          |
|       `J`        |                                 current density                                 |                        $\mathbf{J}$                        |
|   `J_contra_r`   |                      contravariant radial current density                       |                         $J^{\rho}$                         |
|   `J_contra_t`   |                     contravariant poloidal current density                      |                        $J^{\theta}$                        |
|   `J_contra_z`   |                     contravariant toroidal current density                      |                        $J^{\zeta}$                         |
|      `Jac`       |                              Jacobian determinant                               |                       $\mathcal{J}$                        |
|     `Jac_B`      |                   Jacobian determinant in Boozer coordinates                    |                      $\mathcal{J}_B$                       |
|     `Jac_P`      |                    Jacobian determinant in PEST coordinates                     |                      $\mathcal{J}_P$                       |
|     `Jac_h`      |                         reference Jacobian determinant                          |                      $\mathcal{J}_h$                       |
|     `Jac_l`      |                          logical Jacobian determinant                           |                      $\mathcal{J}_l$                       |
|       `LA`       |                          straight field line potential                          |                         $\lambda$                          |
|      `N_FP`      |                             number of field periods                             |                          $N_{FP}$                          |
|      `Phi`       |                             toroidal magnetic flux                              |                           $\Phi$                           |
|     `Phi_n`      |                        normalized toroidal magnetic flux                        |                          $\Phi_n$                          |
|       `V`        |                                  plasma volume                                  |                            $V$                             |
|     `W_MHD`      |                                total MHD energy                                 |                         $W_{MHD}$                          |
|       `X1`       |                           first reference coordinate                            |                           $X^1$                            |
|       `X2`       |                           second reference coordinate                           |                           $X^2$                            |
|      `chi`       |                             poloidal magnetic flux                              |                           $\chi$                           |
| `dB_contra_t_dr` |                radial derivative of the poloidal magnetic field                 |         $\frac{\partial B^\theta}{\partial \rho}$          |
| `dB_contra_t_dt` |               poloidal derivative of the poloidal magnetic field                |        $\frac{\partial B^\theta}{\partial \theta}$         |
| `dB_contra_t_dz` |               toroidal derivative of the poloidal magnetic field                |         $\frac{\partial B^\theta}{\partial \zeta}$         |
| `dB_contra_z_dr` |                radial derivative of the toroidal magnetic field                 |          $\frac{\partial B^\zeta}{\partial \rho}$          |
| `dB_contra_z_dt` |               poloidal derivative of the toroidal magnetic field                |         $\frac{\partial B^\zeta}{\partial \theta}$         |
| `dB_contra_z_dz` |               toroidal derivative of the toroidal magnetic field                |         $\frac{\partial B^\zeta}{\partial \zeta}$          |
|    `dJac_dr`     |                  radial derivative of the Jacobian determinant                  |        $\frac{\partial \mathcal{J}}{\partial \rho}$        |
|    `dJac_dt`     |                 poloidal derivative of the Jacobian determinant                 |       $\frac{\partial \mathcal{J}}{\partial \theta}$       |
|    `dJac_dz`     |                 toroidal derivative of the Jacobian determinant                 |       $\frac{\partial \mathcal{J}}{\partial \zeta}$        |
|   `dJac_h_dr`    |             radial derivative of the reference Jacobian determinant             |       $\frac{\partial \mathcal{J}_h}{\partial \rho}$       |
|   `dJac_h_dt`    |            poloidal derivative of the reference Jacobian determinant            |      $\frac{\partial \mathcal{J}_h}{\partial \theta}$      |
|   `dJac_h_dz`    |            toroidal derivative of the reference Jacobian determinant            |      $\frac{\partial \mathcal{J}_h}{\partial \zeta}$       |
|   `dJac_l_dr`    |              radial derivative of the logical Jacobian determinant              |       $\frac{\partial \mathcal{J}_l}{\partial \rho}$       |
|   `dJac_l_dt`    |             poloidal derivative of the logical Jacobian determinant             |      $\frac{\partial \mathcal{J}_l}{\partial \theta}$      |
|   `dJac_l_dz`    |             toroidal derivative of the logical Jacobian determinant             |      $\frac{\partial \mathcal{J}_l}{\partial \zeta}$       |
|     `dLA_dr`     |             radial derivative of the straight field line potential              |          $\frac{\partial \lambda}{\partial \rho}$          |
|    `dLA_drr`     |          second radial derivative of the straight field line potential          |        $\frac{\partial^2 \lambda}{\partial \rho^2}$        |
|    `dLA_drt`     |         radial-poloidal derivative of the straight field line potential         | $\frac{\partial^2 \lambda}{\partial \rho\partial \theta}$  |
|    `dLA_drz`     |         radial-toroidal derivative of the straight field line potential         |  $\frac{\partial^2 \lambda}{\partial \rho\partial \zeta}$  |
|     `dLA_dt`     |            poloidal derivative of the straight field line potential             |         $\frac{\partial \lambda}{\partial \theta}$         |
|    `dLA_dtt`     |         second poloidal derivative of the straight field line potential         |       $\frac{\partial^2 \lambda}{\partial \theta^2}$       |
|    `dLA_dtz`     |        poloidal-toroidal derivative of the straight field line potential        | $\frac{\partial^2 \lambda}{\partial \theta\partial \zeta}$ |
|     `dLA_dz`     |            toroidal derivative of the straight field line potential             |         $\frac{\partial \lambda}{\partial \zeta}$          |
|    `dLA_dzz`     |         second toroidal derivative of the straight field line potential         |       $\frac{\partial^2 \lambda}{\partial \zeta^2}$        |
|    `dPhi_dr`     |                         toroidal magnetic flux gradient                         |                   $\frac{d\Phi}{d\rho}$                    |
|    `dPhi_drr`    |                        toroidal magnetic flux curvature                         |                 $\frac{d^2\Phi}{d\rho^2}$                  |
|   `dPhi_n_dr`    |                   normalized toroidal magnetic flux gradient                    |                  $\frac{d\Phi_n}{d\rho}$                   |
|   `dV_dPhi_n`    |    derivative of the plasma volume w.r.t. normalized toroidal magnetic flux     |                    $\frac{dV}{d\Phi_n}$                    |
|   `dV_dPhi_n2`   | second derivative of the plasma volume w.r.t. normalized toroidal magnetic flux |                  $\frac{d^2V}{d\Phi_n^2}$                  |
|     `dX1_dr`     |               radial derivative of the first reference coordinate               |            $\frac{\partial X^1}{\partial \rho}$            |
|    `dX1_drr`     |           second radial derivative of the first reference coordinate            |          $\frac{\partial^2 X^1}{\partial \rho^2}$          |
|    `dX1_drt`     |          radial-poloidal derivative of the first reference coordinate           |   $\frac{\partial^2 X^1}{\partial \rho\partial \theta}$    |
|    `dX1_drz`     |          radial-toroidal derivative of the first reference coordinate           |    $\frac{\partial^2 X^1}{\partial \rho\partial \zeta}$    |
|     `dX1_dt`     |              poloidal derivative of the first reference coordinate              |           $\frac{\partial X^1}{\partial \theta}$           |
|    `dX1_dtt`     |          second poloidal derivative of the first reference coordinate           |         $\frac{\partial^2 X^1}{\partial \theta^2}$         |
|    `dX1_dtz`     |         poloidal-toroidal derivative of the first reference coordinate          |   $\frac{\partial^2 X^1}{\partial \theta\partial \zeta}$   |
|     `dX1_dz`     |              toroidal derivative of the first reference coordinate              |           $\frac{\partial X^1}{\partial \zeta}$            |
|    `dX1_dzz`     |          second toroidal derivative of the first reference coordinate           |         $\frac{\partial^2 X^1}{\partial \zeta^2}$          |
|     `dX2_dr`     |              radial derivative of the second reference coordinate               |            $\frac{\partial X^2}{\partial \rho}$            |
|    `dX2_drr`     |           second radial derivative of the second reference coordinate           |          $\frac{\partial^2 X^2}{\partial \rho^2}$          |
|    `dX2_drt`     |          radial-poloidal derivative of the second reference coordinate          |   $\frac{\partial^2 X^2}{\partial \rho\partial \theta}$    |
|    `dX2_drz`     |          radial-toroidal derivative of the second reference coordinate          |    $\frac{\partial^2 X^2}{\partial \rho\partial \zeta}$    |
|     `dX2_dt`     |             poloidal derivative of the second reference coordinate              |           $\frac{\partial X^2}{\partial \theta}$           |
|    `dX2_dtt`     |          second poloidal derivative of the second reference coordinate          |         $\frac{\partial^2 X^2}{\partial \theta^2}$         |
|    `dX2_dtz`     |         poloidal-toroidal derivative of the second reference coordinate         |   $\frac{\partial^2 X^2}{\partial \theta\partial \zeta}$   |
|     `dX2_dz`     |             toroidal derivative of the second reference coordinate              |           $\frac{\partial X^2}{\partial \zeta}$            |
|    `dX2_dzz`     |          second toroidal derivative of the second reference coordinate          |         $\frac{\partial^2 X^2}{\partial \zeta^2}$          |
|    `dchi_dr`     |                         poloidal magnetic flux gradient                         |                   $\frac{d\chi}{d\rho}$                    |
|    `dg_rr_dr`    |           radial derivative of the rr component of the metric tensor            |       $\frac{\partial g_{\rho\rho}}{\partial \rho}$        |
|    `dg_rr_dt`    |          poloidal derivative of the rr component of the metric tensor           |      $\frac{\partial g_{\rho\rho}}{\partial \theta}$       |
|    `dg_rr_dz`    |          toroidal derivative of the rr component of the metric tensor           |       $\frac{\partial g_{\rho\rho}}{\partial \zeta}$       |
|    `dg_rt_dr`    |           radial derivative of the rt component of the metric tensor            |      $\frac{\partial g_{\rho\theta}}{\partial \rho}$       |
|    `dg_rt_dt`    |          poloidal derivative of the rt component of the metric tensor           |     $\frac{\partial g_{\rho\theta}}{\partial \theta}$      |
|    `dg_rt_dz`    |          toroidal derivative of the rt component of the metric tensor           |      $\frac{\partial g_{\rho\theta}}{\partial \zeta}$      |
|    `dg_rz_dr`    |           radial derivative of the rz component of the metric tensor            |       $\frac{\partial g_{\rho\zeta}}{\partial \rho}$       |
|    `dg_rz_dt`    |          poloidal derivative of the rz component of the metric tensor           |      $\frac{\partial g_{\rho\zeta}}{\partial \theta}$      |
|    `dg_rz_dz`    |          toroidal derivative of the rz component of the metric tensor           |      $\frac{\partial g_{\rho\zeta}}{\partial \zeta}$       |
|    `dg_tt_dr`    |           radial derivative of the tt component of the metric tensor            |     $\frac{\partial g_{\theta\theta}}{\partial \rho}$      |
|    `dg_tt_dt`    |          poloidal derivative of the tt component of the metric tensor           |    $\frac{\partial g_{\theta\theta}}{\partial \theta}$     |
|    `dg_tt_dz`    |          toroidal derivative of the tt component of the metric tensor           |     $\frac{\partial g_{\theta\theta}}{\partial \zeta}$     |
|    `dg_tz_dr`    |           radial derivative of the tz component of the metric tensor            |      $\frac{\partial g_{\theta\zeta}}{\partial \rho}$      |
|    `dg_tz_dt`    |          poloidal derivative of the tz component of the metric tensor           |     $\frac{\partial g_{\theta\zeta}}{\partial \theta}$     |
|    `dg_tz_dz`    |          toroidal derivative of the tz component of the metric tensor           |     $\frac{\partial g_{\theta\zeta}}{\partial \zeta}$      |
|    `dg_zz_dr`    |           radial derivative of the zz component of the metric tensor            |      $\frac{\partial g_{\zeta\zeta}}{\partial \rho}$       |
|    `dg_zz_dt`    |          poloidal derivative of the zz component of the metric tensor           |     $\frac{\partial g_{\zeta\zeta}}{\partial \theta}$      |
|    `dg_zz_dz`    |          toroidal derivative of the zz component of the metric tensor           |      $\frac{\partial g_{\zeta\zeta}}{\partial \zeta}$      |
|    `diota_dr`    |                          rotational transform gradient                          |                   $\frac{d\iota}{d\rho}$                   |
|    `dnu_B_dt`    |  poloidal derivative of the Boozer potential computed from the magnetic field   |   $\left.\frac{\partial \nu_B}{\partial \theta}\right\|$   |
|    `dnu_B_dz`    |  toroidal derivative of the Boozer potential computed from the magnetic field   |   $\left.\frac{\partial \nu_B}{\partial \zeta}\right\|$    |
|     `dp_dr`      |                                pressure gradient                                |                     $\frac{dp}{d\rho}$                     |
|      `e_X1`      |                      first reference tangent basis vector                       |                     $\mathbf{e}_{X^1}$                     |
|      `e_X2`      |                      second reference tangent basis vector                      |                     $\mathbf{e}_{X^2}$                     |
|     `e_rho`      |                           radial tangent basis vector                           |                     $\mathbf{e}_\rho$                      |
|    `e_theta`     |                          poloidal tangent basis vector                          |                    $\mathbf{e}_\theta$                     |
|   `e_theta_B`    |               poloidal tangent basis vector in Boozer coordinates               |                  $\mathbf{e}_{\theta_B}$                   |
|     `e_zeta`     |                          toroidal tangent basis vector                          |                     $\mathbf{e}_\zeta$                     |
|    `e_zeta3`     |                     toroidal reference tangent basis vector                     |                   $\mathbf{e}_{\zeta^3}$                   |
|    `e_zeta_B`    |               toroidal tangent basis vector in Boozer coordinates               |                   $\mathbf{e}_{\zeta_B}$                   |
|      `g_rr`      |                        rr component of the metric tensor                        |                       $g_{\rho\rho}$                       |
|      `g_rt`      |                        rt component of the metric tensor                        |                      $g_{\rho\theta}$                      |
|      `g_rz`      |                        rz component of the metric tensor                        |                      $g_{\rho\zeta}$                       |
|      `g_tt`      |                        tt component of the metric tensor                        |                     $g_{\theta\theta}$                     |
|      `g_tz`      |                        tz component of the metric tensor                        |                     $g_{\theta\zeta}$                      |
|      `g_zz`      |                        zz component of the metric tensor                        |                      $g_{\zeta\zeta}$                      |
|    `grad_rho`    |                         radial reciprocal basis vector                          |                        $\nabla\rho$                        |
|   `grad_theta`   |                        poloidal reciprocal basis vector                         |                       $\nabla\theta$                       |
|  `grad_theta_P`  |              poloidal reciprocal basis vector in PEST coordinates               |                     $\nabla \theta_P$                      |
|   `grad_zeta`    |                        toroidal reciprocal basis vector                         |                       $\nabla\zeta$                        |
|      `iota`      |                              rotational transform                               |                          $\iota$                           |
|     `iota_0`     |               geometric contribution to the rotational transform                |                         $\iota_0$                          |
|    `iota_avg`    |                          average rotational transform                           |                       $\bar{\iota}$                        |
|   `iota_curr`    |            toroidal current contribution to the rotational transform            |                       $\iota_{tor}$                        |
|  `major_radius`  |                                  major radius                                   |                         $r_{maj}$                          |
|  `minor_radius`  |                                  minor radius                                   |                         $r_{min}$                          |
|     `mod_B`      |                          modulus of the magnetic field                          |                $\left\|\mathbf{B}\right\|$                 |
|     `mod_F`      |                            modulus of the MHD force                             |                     $\left\|F\right\|$                     |
|     `mod_J`      |                         modulus of the current density                          |                $\left\|\mathbf{J}\right\|$                 |
|   `mod_e_rho`    |                   modulus of the radial tangent basis vector                    |              $\left\|\mathbf{e}_\rho\right\|$              |
|  `mod_e_theta`   |                  modulus of the poloidal tangent basis vector                   |             $\left\|\mathbf{e}_\theta\right\|$             |
|   `mod_e_zeta`   |                  modulus of the toroidal tangent basis vector                   |             $\left\|\mathbf{e}_\zeta\right\|$              |
|  `mod_grad_rho`  |                  modulus of the radial reciprocal basis vector                  |                $\left\|\nabla\rho\right\|$                 |
| `mod_grad_theta` |                 modulus of the poloidal reciprocal basis vector                 |               $\left\|\nabla\theta\right\|$                |
| `mod_grad_zeta`  |                 modulus of the toroidal reciprocal basis vector                 |                $\left\|\nabla\zeta\right\|$                |
|      `mu0`       |                                magnetic constant                                |                          $\mu_0$                           |
|       `p`        |                                    pressure                                     |                            $p$                             |
|      `pos`       |                                 position vector                                 |                        $\mathbf{x}$                        |
|     `shear`      |                              global magnetic shear                              |                           $s_g$                            |
|    `theta_P`     |                       poloidal angle in PEST coordinates                        |                         $\theta_P$                         |
|      `xyz`       |                           cartesian vector components                           |                        $\mathbf{x}$                        |
