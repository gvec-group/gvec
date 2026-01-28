!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

!===================================================================================================================================
!>
!!# Module **TEMPLATE**
!!
!! some description
!===================================================================================================================================
MODULE MODgvec_TEMPLATE

USE MODgvec_Globals, ONLY: wp

IMPLICIT NONE
PRIVATE
PUBLIC :: t_templatetype, template_subroutine

!===================================================================================================================================

! avoid global variables, unless they are parameters
! prefer parameters over macros
! don't be afraid to use strings
INTEGER, PARAMETER :: PARAM_ONE = 1
CHARACTER(LEN=3), PARAMETER :: PARAM_STR = "abc"

!===================================================================================================================================
!> template type
!!
!! some description
!! for abstract types, see the profiles directory as a template
!! avoid abstract types if they are not needed (e.g. if we only have a single implementation, even if we might want to extend it later)
!===================================================================================================================================
TYPE :: t_templatetype
  INTEGER :: var1 = 0                   !! some description
  REAL(wp) :: var2 = 0.0_wp             !! some description
  CHARACTER(LEN=10) :: var3 = "default" !! some description
  LOGICAL :: useThis = .FALSE.          !! some description

  CONTAINS

  FINAL :: template_free
  PROCEDURE :: method => TEMPLATE_METHOD
END TYPE t_templatetype

INTERFACE t_templatetype
  MODULE PROCEDURE InitTEMPLATE
END INTERFACE t_templatetype

CONTAINS

!===================================================================================================================================
!> some function
!!
!! some description
!! use functions when it is applicable (i.e. a single return value)
!===================================================================================================================================
FUNCTION xxx(input) RESULT(output)
  ! MODULES
  USE MODgvec_Globals, ONLY: wp
  IMPLICIT NONE
  ! INPUT VARIABLES -------------------------!
  REAL(wp), INTENT(IN) :: input
  ! OUTPUT VARIABLES -------------------------!
  REAL(wp) :: output
  ! LOCAL VARIABLES -------------------------!
  INTEGER :: i
  ! CODE --------------------------------------------------------------------------------------------------------------------------!

END FUNCTION xxx

!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE InitTEMPLATE()
  ! MODULES
  USE MOD_Globals,ONLY:UNIT_stdOut,fmt_sep
  USE MOD_TEMPLATE_Vars
  USE MOD_ReadInTools,ONLY:GETLOGICAL
  IMPLICIT NONE
  ! INPUT VARIABLES -------------------------!
  ! OUTPUT VARIABLES -------------------------!
  ! LOCAL VARIABLES -------------------------!
  ! CODE --------------------------------------------------------------------------------------------------------------------------!
  SWRITE(UNIT_stdOut,'(A)')'INIT TEMPLATE ...'
  useThis    = GETLOGICAL('useThis','F')

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitTEMPLATE


!===================================================================================================================================
!>
!!
!===================================================================================================================================
SUBROUTINE TEMPLATE_METHOD()
  ! MODULES
  USE MOD_Globals, ONLY:wp
  USE MOD_TEMPLATE_Vars
  IMPLICIT NONE
  ! INPUT VARIABLES -------------------------!
  ! OUTPUT VARIABLES -------------------------!
  ! LOCAL VARIABLES -------------------------!
  ! CODE --------------------------------------------------------------------------------------------------------------------------!
END SUBROUTINE TEMPLATE

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeTEMPLATE
  ! MODULES
  USE MOD_TEMPLATE_Vars
  IMPLICIT NONE
  ! INPUT VARIABLES -------------------------!
  ! OUTPUT VARIABLES -------------------------!
  ! LOCAL VARIABLES -------------------------!
  ! CODE --------------------------------------------------------------------------------------------------------------------------!

END SUBROUTINE FinalizeTEMPLATE

END MODULE MOD_TEMPLATE
