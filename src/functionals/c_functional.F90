!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** c_functional **
!!
!! contains the type that points to the routines of one chosen functional
!!
!===================================================================================================================================
MODULE MODgvec_c_functional
! MODULES
USE MODgvec_Globals    ,ONLY:wp,Unit_stdOut,abort
IMPLICIT NONE

PRIVATE
PUBLIC t_functional
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES
TYPE, ABSTRACT :: c_functional
  CONTAINS
    PROCEDURE(i_sub_functional     ),DEFERRED :: init
    PROCEDURE(i_sub_functional     ),DEFERRED :: initSolution
    PROCEDURE(i_sub_functional_min ),DEFERRED :: minimize
    PROCEDURE(i_sub_functional     ),DEFERRED :: free

END TYPE c_functional

ABSTRACT INTERFACE
  SUBROUTINE i_sub_functional( sf)
    IMPORT c_functional
    CLASS(c_functional), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_functional

  SUBROUTINE i_sub_functional_min( sf)
    IMPORT wp,c_functional
    CLASS(c_functional), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_functional_min

END INTERFACE


TYPE,ABSTRACT,EXTENDS(c_functional) :: t_functional
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: which_functional         !! points to functional (1: MHD3D)
  !---------------------------------------------------------------------------------------------------------------------------------

END TYPE t_functional

!===================================================================================================================================


END MODULE MODgvec_c_functional
