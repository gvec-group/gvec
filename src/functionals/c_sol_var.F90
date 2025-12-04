!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** C_Sol_Var **
!!
!! contains only abstract type c_sol_var
!!
!===================================================================================================================================
MODULE MODgvec_c_sol_var
! MODULES
USE MODgvec_Globals,ONLY:wp
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

END MODULE MODgvec_c_sol_var
