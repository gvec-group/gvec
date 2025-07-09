# Other Interfaces

GVEC also has interfaces to:
* the MHD stability code *CASTOR3D* *(only $RZ\phi$ $h$-map)*
* the non-linear MHD code [JOREK](https://www.jorek.eu/) *(only $RZ\phi$ $h$-map)*
* the turbulence code [GENE](https://genecode.org/) *(only $RZ\phi$ $h$-map)*
* [HOPR](https://hopr.readthedocs.io) *(only $RZ\phi$ $h$-map)*
* the geometric plasma simulation package [Struphy](https://struphy.pages.mpcdf.de/struphy)
* the orbit-tracer [SIMPLE](https://github.com/itpplasma/SIMPLE)

## VMEC
GVEC is also compatible with the MHD equilibrium code VMEC to a certain extent.
In particular a VMEC equilibrium can be used as the initial state for a GVEC computation.
See [initialization parameters](./gvec-parameter-list.md#initialization-parameters) for details.
