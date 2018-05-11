!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>

! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================


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
LOGICAL                 :: useSFL                    !! use straight-field line coordinates
LOGICAL                 :: reLambda                  !! switch for recomputing lambda
CHARACTER(LEN = 256)    :: VMECdataFile
INTEGER,ALLOCATABLE     :: xmAbs(:)                  !! |xm(iMode)|, 1 for m=0, 2 for even, 3 for odd
REAL(wp),ALLOCATABLE    :: Phi_prof(:)               !! TOROIDAL flux profile (called phi in VMEC)
REAL(wp),ALLOCATABLE    :: Phinorm_prof(:)           !! normalized TOROIDAL flux profile 
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

