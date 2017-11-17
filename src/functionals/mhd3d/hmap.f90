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
!!# Module ** hmap **
!!
!! contains only the abstract type to point to a specific map h (maps  omega_p x S^1 --> omega) 
!!
!===================================================================================================================================
MODULE MOD_hmap
! MODULES
USE MOD_Globals    ,ONLY:wp,Unit_stdOut,abort
IMPLICIT NONE

PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
TYPE, ABSTRACT :: c_hmap
  CONTAINS
    PROCEDURE(i_sub_hmap     ),DEFERRED :: init
    PROCEDURE(i_sub_hmap     ),DEFERRED :: free

END TYPE c_hmap

ABSTRACT INTERFACE
  SUBROUTINE i_sub_hmap( sf )
    IMPORT c_hmap
    CLASS(c_hmap), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_hmap

END INTERFACE
 

TYPE,ABSTRACT,EXTENDS(c_hmap) :: t_hmap
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: which_hmap         !! points to hmap (1: MHD3D) 
  !---------------------------------------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap

!===================================================================================================================================

END MODULE MOD_hmap

