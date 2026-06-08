# SIMPLE

From the [README](https://github.com/itpplasma/SIMPLE)

> **S**ymplectic **I**ntegration **M**ethods for **P**article **L**oss **E**stimation
>
> SIMPLE computes statistical losses of guiding-center orbits for particles of given mass, charge and energy from the volume of 3D magnetic configurations. Orbits are traced via a symplectic integrator [1,2] that guarantees conservation of invariants of motion within fixed bounds over long integration periods. A classifier based on Poincarè plots [1] allows for accelerated prediction for confinement of regular (non-chaotic) orbits.
>
> The main focus of SIMPLE is to make computation of fusion alpha particle losses fast enough to be used directly in stellarator optimization codes. For that reason currently 3D configurations with nested magnetic flux surfaces are supported based on a VMEC equilibrium file in NetCDF format, and orbits are computed without taking collisions into account.

SIMPLE can use a GVEC equilibrium for the magnetic configuration, through the use of a *chartmap*-file. Such a file contains the geometry & magnetic field information in Boozer coordinates and is stored in NetCDF format. To generate a *chartmap*-file, SIMPLE provides a python script in [`tools/gvec_to_boozer_chartmap.py`](https://github.com/itpplasma/SIMPLE/blob/main/tools/gvec_to_boozer_chartmap.py) which takes a GVEC parameter & statefile as inputs. The script requires `netcdf4` and `gvec` to be installed and more information can be found using the script's `--help` option. A chartmap can then be used by specifying the `field_input` and `geom_input` parameters in the SIMPLE input file. For more details, see the SIMPLE `README`, `docs` and `examples` folder.

* Repository: <https://github.com/itpplasma/SIMPLE>
* Reference: C. G. Albert, S. V. Kasilov, and W. Kernbichler, Accelerated methods for direct computation of fusion alpha particle losses within stellarator optimization. J. Plasma Phys 86, 815860201 (2020), [DOI: 10.1017/S0022377820000203](https://doi.org/10.1017/S0022377820000203)

## Conventions

* SIMPLE can use different coordinate systems for the geometry, magnetic field and the particle tracing.
* When using a *chartmap*, the radial coordinate is proportional to the square root of the normalized toroidal flux (GVEC's $\rho$), and the angular coordinates are the Boozer angles of the GVEC equilibrium.
* SIMPLE uses a left-handed coordinate system, so there is a choice to either flip the toroidal or poloidal coordinate when generating the *chartmap* from GVEC.
* The toroidal flux in SIMPLE (`torflux`) is really the $A_\theta$ component of the vector potential at the edge and follows the signs accordingly (left-handed!).
* The starting position specified in the SIMPLE input file (`sbeg`) is always in the normalized toroidal flux coordinate (GVEC's $\rho^2$).
* SIMPLE uses CGS units internally, GVEC uses SI units.
* VMEC and GVEC use different definitions for the effective major radius, but the difference is usually small. SIMPLE mainly uses the effective major radius to compute the integration timestep, so this affects the performance, but not the results.
