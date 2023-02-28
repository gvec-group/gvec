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
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** HMAP new **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_hmap
! MODULES
USE MODgvec_c_hmap    , ONLY: c_hmap
IMPLICIT NONE
PUBLIC


CONTAINS

!===================================================================================================================================
!> initialize the type hmap, also readin parameters here if necessary 
!!
!===================================================================================================================================
SUBROUTINE hmap_new( sf, which_hmap,nfp)
! MODULES
USE MODgvec_Globals   , ONLY: abort,wp
USE MODgvec_hmap_RZ   , ONLY: t_hmap_RZ
USE MODgvec_hmap_RphiZ, ONLY: t_hmap_RphiZ
USE MODgvec_hmap_knot , ONLY: t_hmap_knot
USE MODgvec_hmap_cyl  , ONLY: t_hmap_cyl
USE MODgvec_hmap_frenet,ONLY: t_hmap_frenet
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: which_hmap         !! input number of field periods
  INTEGER       , INTENT(IN   ) :: nfp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_hmap), ALLOCATABLE,INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(which_hmap)
  CASE(1)
    ALLOCATE(t_hmap_RZ :: sf) 
    sf%which_hmap=which_hmap
  !CASE(2)
  !  ALLOCATE(t_hmap_RphiZ :: sf) 
  !  sf%which_hmap=which_hmap
  CASE(3)
    ALLOCATE(t_hmap_cyl :: sf) 
    sf%which_hmap=which_hmap
  CASE(10)
    ALLOCATE(t_hmap_knot :: sf) 
    sf%which_hmap=which_hmap
  CASE(20)
    ALLOCATE(sf,source=t_hmap_frenet(which_hmap=which_hmap,nfp=nfp)) 
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "this hmap choice does not exist  !")
  END SELECT 
  CALL sf%init()
  

END SUBROUTINE hmap_new


END MODULE MODgvec_hmap

