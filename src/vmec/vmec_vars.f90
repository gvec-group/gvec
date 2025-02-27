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
CLASS(c_rProfile), ALLOCATABLE :: Phi_profile        !! B-spline profiles in (rho^2) for Phi
CLASS(c_rProfile), ALLOCATABLE :: Chi_profile        !! B-spline profile in (rho^2) for chi

!===================================================================================================================================
END MODULE MODgvec_VMEC_Vars

