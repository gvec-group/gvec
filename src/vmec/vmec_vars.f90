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
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------


! GLOBAL VARIABLES 
LOGICAL                 :: useVMEC                   !! main switch
LOGICAL                 :: switchZeta                !! True: change from R,phi,Z to R,Z,phi coordinate system
LOGICAL                 :: reLambda                  !! switch for recomputing lambda
CHARACTER(LEN = 256)    :: VMECdataFile
INTEGER                 :: VMECFile_Format           !! 0: netcdf format (default), 1: nemec ascii, 2: nemec binary
INTEGER,ALLOCATABLE     :: xmAbs(:)                  !! |xm(iMode)|, 1 for m=0, 2 for even, 3 for odd
REAL(wp),ALLOCATABLE    :: Phi_prof(:)               !! TOROIDAL flux profile (called phi in VMEC)
REAL(wp),ALLOCATABLE    :: normFlux_prof(:)          !! normalized flux profile, can be either toroidal of poloidal flux) 
REAL(wp),ALLOCATABLE    :: chi_prof(:)               !! POLOIDAL flux profile (called chi in VMEC)

REAL(wp),ALLOCATABLE    :: rho(:)                    !! := sqrt(phinorm) at all flux surface 
REAL(wp),ALLOCATABLE    :: pres_Spl(:,:)             !! Spline coefficients in (rho) for Pressure, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: iota_Spl(:,:)             !! Spline coefficients in (rho) for iota, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: Phi_Spl(:,:)              !! Spline coefficients in (rho) for Phi, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: chi_Spl(:,:)              !! Spline coefficients in (rho) for chi, (1:4,nFluxVMEC)
REAL(wp),ALLOCATABLE    :: Rmnc_Spl(:,:,:)           !! modified spline coefficients R cosine, (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: Rmns_Spl(:,:,:)           !! modified spline coefficients R sine,   (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: lmnc_Spl(:,:,:)           !! modified spline coefficients of lambda cosine , (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: lmns_Spl(:,:,:)           !! modified spline coefficients of lambda sine,   (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: Zmnc_Spl(:,:,:)           !! modified spline coefficients of Z cosine, (1:4,iFlux,iMode)
REAL(wp),ALLOCATABLE    :: Zmns_Spl(:,:,:)           !! modified spline coefficients of Z sine,   (1:4,iFlux,iMode)

!===================================================================================================================================
END MODULE MODgvec_VMEC_Vars

