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
!!# Module ** functional **
!!
!! contains the routines to initialize and finalize the functional
!!
!===================================================================================================================================
MODULE MODgvec_functional
! MODULES
USE MODgvec_Globals    ,ONLY:wp,Unit_stdOut,abort
USE MODgvec_c_functional, ONLY: t_functional
IMPLICIT NONE

PUBLIC

!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize the type functional with number of elements
!!
!===================================================================================================================================
SUBROUTINE InitFunctional(sf, which_functional)
! MODULES
USE MODgvec_MHD3D, ONLY :t_functional_mhd3d
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: which_functional
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_functional), ALLOCATABLE,INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(which_functional)
  CASE(1)
    ALLOCATE(t_functional_mhd3d :: sf)
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "this functional choice does not exist (MHD3D=1) !")
  END SELECT 

  sf%which_functional=which_functional
  CALL sf%init()

END SUBROUTINE InitFunctional


!===================================================================================================================================
!> finalize the type functional
!!
!===================================================================================================================================
SUBROUTINE FinalizeFunctional(sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_functional), ALLOCATABLE,INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL sf%free()

END SUBROUTINE FinalizeFunctional

END MODULE MODgvec_functional

