To check:

- RESTART,OUTPUT: only the root should read/write data!!  (restart.f90)
- ALL WRITE -> SWRITE + use MPIroot
- initialization on the root: calling mhd_profiles.f90, initVMEC...
- checkEvalForce!
- CMake `USE_MPI` option can not be used with
  - `COMPILE_GVEC_TO_GENE`
  - `COMPILE_GVEC_TO_HOPR`
  - `COMPILE_GVEC_TO_CASTOR3D`
  - `COMPILE_GVEC_TO_JOREK`
