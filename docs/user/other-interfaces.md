# Other Interfaces

GVEC also has interfaces to:
* the MHD stability code *CASTOR3D* *(only $RZ\phi$ $h$-map)*
* the non-linear MHD code [JOREK](https://www.jorek.eu/) *(only $RZ\phi$ $h$-map)*
* the turbulence code [GENE](https://genecode.org/) *(only $RZ\phi$ $h$-map)*
* [HOPR](https://hopr.readthedocs.io) *(only $RZ\phi$ $h$-map)*
* the geometric plasma simulation package [Struphy](https://struphy.pages.mpcdf.de/struphy)
* the ideal MHD relaxation code *MRX*

## VMEC
GVEC is also compatible with the MHD equilibrium code VMEC to a certain extent.
In particular a VMEC equilibrium can be used as the initial state for a GVEC computation.
See [initialization parameters](./gvec-parameter-list.md#initialization-parameters) for details.
Alternatively, a VMEC input file can be converted to a GVEC parameter file with the `pygvec convert-params --vmec` command, and a VMEC equilibrium can be converted to a GVEC state file with the `pygvec convert-wout` command.
Note that the equilibrium (or its inputs) cannot be converted exactly due to the different numerical methods employed by the two codes.
