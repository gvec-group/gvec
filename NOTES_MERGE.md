To check:

- RESTART,OUTPUT: only the root should read/write data!!  (restart.f90)
- ALL WRITE -> SWRITE + use MPIroot
- initialization on the root: calling mhd_profiles.f90, initVMEC...
- checkEvalForce!
