---
file_format: mystnb
kernelspec:
  name: python3
---
(table-of-quantities)=
# Quantities for Postprocessing
GVEC provides a number of built-in quantities for postprocessing an equilibrium. These can be used for example with the [`compute`](#gvec.core.compute.compute), [`evaluate`](#gvec.core.compute.evaluate) and [`evaluate_sfl`](#gvec.core.compute.evaluate_sfl) functions, or the [plotting utilities](../tutorials/notebooks/051_plotting.ipynb). The [`table_of_quantities`](#gvec.core.compute.table_of_quantities) function can be used to generate a table of available quantities.

See [the derivations in the developer guide](../dev/derivations.md) for more details on the definitions of the quantities listed below.

```{code-cell} ipython3
:tags: [remove-cell]

import gvec
from IPython.display import Markdown
from myst_nb import glue

def toq(label, Qs):
    txt = gvec.table_of_quantities(markdown=True, keys=Qs)
    return glue(label, txt, display=False)
```

## Solution variables
```{code-cell} ipython3
:tags: [remove-input]
Qs_solution = ["X1", "X2", "LA",
  "dX1_dr", "dX1_dt", "dX1_dz", "dX1_drr", "dX1_drt", "dX1_drz", "dX1_dtt", "dX1_dtz", "dX1_dzz",
  "dX2_dr", "dX2_dt", "dX2_dz", "dX2_drr", "dX2_drt", "dX2_drz", "dX2_dtt", "dX2_dtz", "dX2_dzz",
  "dLA_dr", "dLA_dt", "dLA_dz", "dLA_drr", "dLA_drt", "dLA_drz", "dLA_dtt", "dLA_dtz", "dLA_dzz"]
toq("toq_solution", Qs_solution)
```
```{glue:md} toq_solution
:format: myst
```

### Additional solution variables for the Boozer transform
```{note}
If a Boozer transform is performed (e.g. [`evaluate_sfl(..., sfl="boozer")`](#gvec.core.compute.evaluate_sfl)), the variables `NU_B`, `dNU_B_dr`, ... `dNU_B_dzz` are available instead of those given below.
```

```{code-cell} ipython3
:tags: [remove-input]
#Qs_boozer = ["NU_B", "dNU_B_dr", "dNU_B_dt", "dNU_B_dz",
#  "dNU_B_drr", "dNU_B_drt", "dNU_B_drz", "dNU_B_dtt", "dNU_B_dtz", "dNU_B_dzz"]
Qs_boozer = ["dNU_B_dt", "dNU_B_dz"]
toq("toq_boozer_nu", Qs_boozer)
```
```{glue:md} toq_boozer_nu
:format: myst
```

## Curvilinear coordinates
```{code-cell} ipython3
:tags: [remove-input]
Qs_metric = ["Jac_l", "Jac", "e_rho", "e_theta", "e_zeta", "grad_rho", "grad_theta", "grad_zeta",
  "g_rr", "g_rt", "g_rz", "g_tt", "g_tz", "g_zz",
  "k_rr", "k_rt", "k_rz", "k_tt", "k_tz", "k_zz",
  "II_tt", "II_tz", "II_zz",
  "mod_e_rho", "mod_e_theta", "mod_e_zeta", "mod_grad_rho", "mod_grad_theta", "mod_grad_zeta",
  "dJac_l_dr", "dJac_l_dt", "dJac_l_dz", "dJac_dr", "dJac_dt", "dJac_dz",
  "dg_rr_dr", "dg_rr_dt", "dg_rr_dz", "dg_rt_dr", "dg_rt_dt", "dg_rt_dz",
  "dg_rz_dr", "dg_rz_dt", "dg_rz_dz", "dg_tt_dr", "dg_tt_dt", "dg_tt_dz",
  "dg_tz_dr", "dg_tz_dt", "dg_tz_dz", "dg_zz_dr", "dg_zz_dt", "dg_zz_dz",
]
toq("toq_metric", Qs_metric)
```
```{glue:md} toq_metric
:format: myst
```

### Boozer coordinates
```{code-cell} ipython3
:tags: [remove-input]
Qs_metric_B = ["Jac_B", "e_rho_B", "e_theta_B", "e_zeta_B", "grad_theta_B", "grad_zeta_B",
  "g_rr_B", "g_rt_B", "g_rz_B", "g_tt_B", "g_tz_B", "g_zz_B",
  "k_tt_B", "k_tz_B", "k_zz_B", "II_tt_B", "II_tz_B", "II_zz_B",
  ]
toq("toq_metric_B", Qs_metric_B)
```
```{glue:md} toq_metric_B
:format: myst
```

### PEST coordinates
```{code-cell} ipython3
:tags: [remove-input]
Qs_metric_P = ["theta_P", "Jac_P", "e_rho_P", "e_theta_P", "e_zeta_P", "grad_theta_P",
  "g_rr_P", "g_rt_P", "g_rz_P", "g_tt_P", "g_tz_P", "g_zz_P",
  "k_tt_P", "k_tz_P", "k_zz_P", "II_tt_P", "II_tz_P", "II_zz_P",
  ]
toq("toq_metric_P", Qs_metric_P)
```
```{glue:md} toq_metric_P
:format: myst
```

## Reference coordinates
Also called the *coordinate frame* or *$h$-map*.

```{code-cell} ipython3
:tags: [remove-input]
Qs_hmap = ["Jac_h", "N_FP", "e_q1", "e_q2", "e_q3",
  "k_q1q1", "k_q1q2", "k_q1q3", "k_q2q2", "k_q2q3", "k_q3q3",
  "dJac_h_dr", "dJac_h_dt", "dJac_h_dz",
  ]
toq("toq_hmap", Qs_hmap)
```
```{glue:md} toq_hmap
:format: myst
```

## Geometric quantities
```{code-cell} ipython3
:tags: [remove-input]
Qs_geometric = ["pos", "normal", "dA", "V", "L_axis", "A_surface", "r_major", "r_minor", "aspect_ratio", "elongation",
  "dV_dPhi_n", "dV_dPhi_n2"]
toq("toq_geometric", Qs_geometric)
```
```{glue:md} toq_geometric
:format: myst
```

## Magnetic flux
```{versionchanged} 1.5
The contributions to to the rotational transform were renamed: `iota_0` → `iota_geom` and `iota_curr_0` → `iota_curr_bar`.
The old identifiers are still usable, but will be removed in a future version.
```

```{code-cell} ipython3
:tags: [remove-input]
Qs_flux = ["Phi_edge", "Phi", "dPhi_dr", "dPhi_drr", "chi", "dchi_dr", "dchi_drr", "iota", "diota_dr", "diota_drr", "iota_avg", "iota_avg2", "shear", "shear_avg", "shear_avg2", "iota_geom", "iota_curr", "iota_curr_bar", "iota_0", "iota_curr_0"]
toq("toq_flux", Qs_flux)
```
```{glue:md} toq_flux
:format: myst
```

### Current profiles
```{code-cell} ipython3
:tags: [remove-input]
Qs_current = ["I_tor", "I_pol", "dI_tor_dr"]
toq("toq_current", Qs_current)
```
```{glue:md} toq_current
:format: myst
```

## Pressure
```{code-cell} ipython3
:tags: [remove-input]
Qs_pressure = ["p", "dp_dr", "dp_drr", "beta_avg"]
toq("toq_pressure", Qs_pressure)
```
```{glue:md} toq_pressure
:format: myst
```

## Magnetic field
```{code-cell} ipython3
:tags: [remove-input]
Qs_B = ["B", "B_contra_t", "B_contra_z", "B_theta_avg", "B_zeta_avg", "mod_B", "grad_mod_B",
  "dB_contra_t_dr", "dB_contra_t_dt", "dB_contra_t_dz", "dB_contra_z_dr", "dB_contra_z_dt", "dB_contra_z_dz",
  "dmod_B_dr", "dmod_B_dt", "dmod_B_dz",
  "dB_theta_avg_dr",
  "dB_dr", "dB_dt", "dB_dz", "db_dr", "db_dt", "db_dz",
  ]
toq("toq_B", Qs_B)
```
```{glue:md} toq_B
:format: myst
```

### Boozer magnetic field
```{code-cell} ipython3
:tags: [remove-input]
Qs_B_B = ["B_contra_t_B", "B_contra_z_B", "B_rho_B", "B_theta_B", "B_zeta_B",
  "dmod_B_dr_B", "dmod_B_dt_B", "dmod_B_dz_B"]
toq("toq_B_B", Qs_B_B)
```
```{glue:md} toq_B_B
:format: myst
```

### PEST magnetic field
```{code-cell} ipython3
:tags: [remove-input]
Qs_B_P = ["B_contra_t_P", "B_rho_P", "B_theta_P", "B_zeta_P",
  "dmod_B_dr_P", "dmod_B_dt_P", "dmod_B_dz_P"]
toq("toq_B_P", Qs_B_P)
```
```{glue:md} toq_B_P
:format: myst
```

## Current density
```{code-cell} ipython3
:tags: [remove-input]
Qs_J = ["J", "J_contra_r", "J_contra_t", "J_contra_z", "mod_J"]
toq("toq_J", Qs_J)
```
```{glue:md} toq_J
:format: myst
```

### Boozer current density
```{code-cell} ipython3
:tags: [remove-input]
Qs_J_B = ["J_contra_t_B", "J_contra_z_B", "J_rho_B", "J_theta_B", "J_zeta_B"]
toq("toq_J_B", Qs_J_B)
```
```{glue:md} toq_J_B
:format: myst
```

### PEST current density
```{code-cell} ipython3
:tags: [remove-input]
Qs_J_P = ["J_contra_t_P", "J_rho_P", "J_theta_P", "J_zeta_P"]
toq("toq_J_P", Qs_J_P)
```
```{glue:md} toq_J_P
:format: myst
```

## MHD Force
```{code-cell} ipython3
:tags: [remove-input]
Qs_force = [
  "W_MHD", "F", "mod_F", "F_r_avg",
  "gradB2_volavg", "mod_F_rel", "mod_F_rel_surfavg", "mod_F_rel_volavg",
  ]
toq("toq_force", Qs_force)
```
```{glue:md} toq_force
:format: myst
```

## Derived quantities of interest
```{code-cell} ipython3
:tags: [remove-input]
Qs_derived = ["mirror_ratio", "L_gradB", "vacuum_magnetic_well_depth"]
toq("toq_derived", Qs_derived)
```
```{glue:md} toq_derived
:format: myst
```

### Mercier stability criterion
```{code-cell} ipython3
:tags: [remove-input]
Qs_Mercier = ["D_Merc", "D_Merc_Curr", "D_Merc_Geod", "D_Merc_Shear", "D_Merc_Well"]
toq("toq_Mercier", Qs_Mercier)
```
```{glue:md} toq_Mercier
:format: myst
```

### Geodesic curvature
```{code-cell} ipython3
:tags: [remove-input]
Qs_curvature = ["kappa_B", "kappa_G", "kappa_N"]
toq("toq_curvature", Qs_curvature)
```
```{glue:md} toq_curvature
:format: myst
```

## Others
```{code-cell} ipython3
:tags: [remove-input]
Qs_others = sorted(set(gvec.core.compute.QUANTITIES.keys()) - set(Qs_solution) - set(Qs_boozer) - set(Qs_metric) - set(Qs_hmap) - set(Qs_flux) - set(Qs_B) - set(Qs_pressure) - set(Qs_curvature) - set(Qs_J) - set(Qs_B_B) - set(Qs_metric_B) - set(Qs_metric_P) - set(Qs_B_P) - set(Qs_B_B) - set(Qs_J_B) - set(Qs_J_P) - set(Qs_derived) - set(Qs_Mercier) - set(Qs_current) - set(Qs_geometric) - set(Qs_force))
toq("toq_others", Qs_others)
```
```{glue:md} toq_others
:format: myst
```
