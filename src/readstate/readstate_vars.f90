!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module ** Read State Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_ReadState_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base
USE MODgvec_c_hmap, ONLY: c_hmap
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
  INTEGER                   :: fileID_r,OutputLevel_r
  CLASS(c_hmap),ALLOCATABLE :: hmap_r                 !! container for global coordinate system
  TYPE(t_sgrid)             :: sgrid_r                !! container for the grid of X1,X2,LA
  CLASS(t_base),ALLOCATABLE :: X1_base_r              !! container for base of X1
  CLASS(t_base),ALLOCATABLE :: X2_base_r              !! container for base of X2
  CLASS(t_base),ALLOCATABLE :: LA_base_r              !! container for base of LA
  REAL(wp),ALLOCATABLE      :: X1_r(:,:)              !! spline x fourier coefs of solution X1
  REAL(wp),ALLOCATABLE      :: X2_r(:,:)              !! spline x fourier coefs of solution X2
  REAL(wp),ALLOCATABLE      :: LA_r(:,:)              !! spline x fourier coefs of solution LA 
  REAL(wp),ALLOCATABLE      :: profiles_1d(:,:)       !! spline coefficients for 1d profiles (using X1_base...needs to be improved!)
  REAL(wp)                  :: a_minor,r_major,volume !! scalars: average minor and major radius, total volume
!===================================================================================================================================
END MODULE MODgvec_ReadState_Vars

