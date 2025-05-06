# Other Interfaces
```{warning}
These interfaces only support the $RZ\phi$ $h$-map!
```

* the MHD stability code *CASTOR3D*
* the non-linear MHD code [JOREK](https://www.jorek.eu/)
* the turbulence code [GENE](https://genecode.org/)
* [HOPR](https://hopr.readthedocs.io)
* the geometric plasma simulation package [Struphy](https://struphy.pages.mpcdf.de/struphy) via [gvec_to_python](https://gitlab.mpcdf.mpg.de/gvec-group/gvec_to_python)

## VMEC
GVEC is also compatible with the MHD equilibrium code VMEC to a certain extend.
In particular a VMEC equilibrium can be used as the initial state for a GVEC computation.
See [initialization parameters](./gvec-parameter-list.md#initialization-parameters) for details.
