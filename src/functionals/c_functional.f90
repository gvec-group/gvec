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
!!# Module ** c_functional **
!!
!! contains the type that points to the routines of one chosen functional
!!
!===================================================================================================================================
MODULE MOD_c_functional
! MODULES
USE MOD_Globals    ,ONLY:wp,Unit_stdOut,abort
IMPLICIT NONE

PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
TYPE, ABSTRACT :: c_functional
  CONTAINS
    PROCEDURE(i_sub_functional     ),DEFERRED :: init
    PROCEDURE(i_sub_functional     ),DEFERRED :: free

END TYPE c_functional

ABSTRACT INTERFACE
  SUBROUTINE i_sub_functional( sf)
    IMPORT c_functional
    CLASS(c_functional), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_functional

END INTERFACE
 

TYPE,ABSTRACT,EXTENDS(c_functional) :: t_functional
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: which_functional         !! points to functional (1: MHD3D) 
  !---------------------------------------------------------------------------------------------------------------------------------

END TYPE t_functional

!===================================================================================================================================


END MODULE MOD_c_functional

