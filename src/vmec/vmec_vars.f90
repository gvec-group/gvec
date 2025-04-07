!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================


!===================================================================================================================================
!>
!!# Module **VMEC Variables**
!!
!!
!===================================================================================================================================
MODULE MODgvec_VMEC_Vars
! MODULES
USE MODgvec_Globals, ONLY: wp
USE MODgvec_rProfile_base, ONLY: c_rProfile
USE MODgvec_cubic_spline, ONLY: t_cubspl
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------


! GLOBAL VARIABLES
LOGICAL                 :: useVMEC                   !! main switch
LOGICAL                 :: switchZeta                !! switch vmec_phi = -zeta
LOGICAL                 :: switchTheta               !! switch vmec_theta = -theta
CHARACTER(LEN = 256)    :: VMECdataFile
INTEGER                 :: VMECFile_Format           !! 0: netcdf format (default), 1: nemec ascii, 2: nemec binary
INTEGER,ALLOCATABLE     :: xmAbs(:)                  !! |xm(iMode)|, 1 for m=0, 2 for even, 3 for odd
REAL(wp),ALLOCATABLE    :: Phi_prof(:)               !! TOROIDAL flux profile (called phi in VMEC)
REAL(wp),ALLOCATABLE    :: normFlux_prof(:)          !! normalized flux profile, can be either toroidal of poloidal flux)
REAL(wp),ALLOCATABLE    :: chi_prof(:)               !! POLOIDAL flux profile (called chi in VMEC)

REAL(wp),ALLOCATABLE    :: rho(:)                    !! := sqrt(phinorm) at all flux surface

TYPE(t_cubspl),ALLOCATABLE    :: Rmnc_Spl(:)           !! cubic spline fit of R cosine, array over modes
TYPE(t_cubspl),ALLOCATABLE    :: Rmns_Spl(:)           !! cubic spline fit of R sine, array over modes
TYPE(t_cubspl),ALLOCATABLE    :: lmnc_Spl(:)           !! cubic spline fit of lambda  cosine , array over modes
TYPE(t_cubspl),ALLOCATABLE    :: lmns_Spl(:)           !! cubic spline fit of lambda sine, array over modes
TYPE(t_cubspl),ALLOCATABLE    :: Zmnc_Spl(:)           !! cubic spline fit of Z cosine,array over modes
TYPE(t_cubspl),ALLOCATABLE    :: Zmns_Spl(:)           !! cubic spline fit of Z sine,array over modes
CLASS(c_rProfile), ALLOCATABLE :: vmec_Phi_profile        !! B-spline profiles in (rho^2) for Phi
CLASS(c_rProfile), ALLOCATABLE :: vmec_Chi_profile        !! B-spline profile in (rho^2) for chi
CLASS(c_rProfile), ALLOCATABLE :: vmec_iota_profile        !! B-spline profiles in (rho^2) for iota
CLASS(c_rProfile), ALLOCATABLE :: vmec_pres_profile        !! B-spline profile in (rho^2) for pressure

!===================================================================================================================================
END MODULE MODgvec_VMEC_Vars
