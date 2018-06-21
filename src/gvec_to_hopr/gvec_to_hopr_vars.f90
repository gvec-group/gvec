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
!!# Module ** GVEC_TO_HOPR Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_hopr_vars
! MODULES
USE MODgvec_Globals,ONLY:wp
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base, base_new
USE MODgvec_sBase  ,ONLY: t_sbase
USE MODgvec_fbase,  ONLY: sin_cos_map 
USE MODgvec_c_hmap, ONLY: c_hmap
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
  CLASS(t_base),ALLOCATABLE         :: X1_base_r
  CLASS(t_base),ALLOCATABLE         :: X2_base_r
  CLASS(t_base),ALLOCATABLE         :: LA_base_r
  CLASS(c_hmap),ALLOCATABLE         :: hmap_r
  TYPE(t_sgrid)                     :: sgrid_r
  REAL,ALLOCATABLE                  :: X1_r(:,:),X2_r(:,:),LA_r(:,:)
  REAL,ALLOCATABLE                  :: profiles_1d(:,:)
                                         
!===================================================================================================================================
END MODULE MODgvec_gvec_to_hopr_vars

