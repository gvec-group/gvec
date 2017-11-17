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
!!# Module ** C_Sol_Var **
!!
!! contains only abstract type c_sol_var
!!
!===================================================================================================================================
MODULE MOD_c_sol_var
! MODULES
USE MOD_Globals,ONLY:wp
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
TYPE,ABSTRACT :: c_sol_var
  INTEGER :: nVars
  CONTAINS
  PROCEDURE(i_sub_sol_var_init  ),DEFERRED :: init
  PROCEDURE(i_sub_sol_var       ),DEFERRED :: free
  PROCEDURE(i_sub_sol_var_set_to_solvar),DEFERRED :: set_to_solvar
  PROCEDURE(i_sub_sol_var_set_to_scalar),DEFERRED :: set_to_scalar
  PROCEDURE(i_sub_sol_var_copy  ),DEFERRED :: copy
  PROCEDURE(i_fun_sol_var_norm_2),DEFERRED :: norm_2
  PROCEDURE(i_sub_sol_var_AXBY  ),DEFERRED :: AXBY
END TYPE c_sol_var

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sol_var_init( sf ,varsize)
    IMPORT c_sol_var
    INTEGER         , INTENT(IN   ) :: varsize(:)
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var_init

  SUBROUTINE i_sub_sol_var( sf ) 
    IMPORT c_sol_var
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var

  FUNCTION i_fun_sol_var_norm_2( sf ) RESULT(norm_2)
    IMPORT wp,c_sol_var
    CLASS(c_sol_var), INTENT(IN   ) :: sf
    REAL(wp)                       :: norm_2(sf%nvars)
  END FUNCTION i_fun_sol_var_norm_2

  SUBROUTINE i_sub_sol_var_copy( sf, tocopy ) 
    IMPORT c_sol_var
    CLASS(c_sol_var), INTENT(IN   ) :: tocopy
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var_copy

  SUBROUTINE i_sub_sol_var_set_to_solvar( sf, toset ,scal_in) 
    IMPORT wp,c_sol_var
    CLASS(c_sol_var), INTENT(IN   ) :: toset
    CLASS(c_sol_var), INTENT(INOUT) :: sf
    REAL(wp),INTENT(IN),OPTIONAL    :: scal_in
  END SUBROUTINE i_sub_sol_var_set_to_solvar

  SUBROUTINE i_sub_sol_var_set_to_scalar( sf, scalar ) 
    IMPORT wp,c_sol_var
    REAL(wp)        , INTENT(IN   ) :: scalar
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var_set_to_scalar

  SUBROUTINE i_sub_sol_var_AXBY( sf, aa, X, bb, Y ) 
    IMPORT wp,c_sol_var
    REAL(wp)        , INTENT(IN   ) :: aa
    CLASS(c_sol_var), INTENT(IN   ) :: X
    REAL(wp)        , INTENT(IN   ) :: bb
    CLASS(c_sol_var), INTENT(IN   ) :: Y
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var_AXBY
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MOD_c_sol_var

