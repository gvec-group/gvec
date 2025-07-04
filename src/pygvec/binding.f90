!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

MODULE MODgvec_py_binding

IMPLICIT NONE
PUBLIC

INTERFACE  f90wrap_abort
  SUBROUTINE f90wrap_abort(ErrorMessage)
    CHARACTER(LEN=*), INTENT(IN) :: ErrorMessage
  END SUBROUTINE f90wrap_abort
END INTERFACE f90wrap_abort

CONTAINS

!================================================================================================================================!
SUBROUTINE redirect_stdout(filename)
  ! MODULES
  USE MODgvec_Globals, ONLY: Unit_stdOut, UNIT_errOut,abort
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=*), INTENT(IN) :: filename
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: ios
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CLOSE(Unit_stdOut)
  OPEN(Unit_stdOut, FILE=filename, ACTION='WRITE', FORM='FORMATTED', ACCESS='SEQUENTIAL', POSITION='APPEND', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(Unit_errOut, '(A)') 'ERROR: could not open file', filename, 'for writing'
    CALL abort(__STAMP__,"")
  END IF
END SUBROUTINE redirect_stdout

!================================================================================================================================!
SUBROUTINE flush_stdout()
  ! MODULES
  USE MODgvec_Globals, ONLY: Unit_stdOut, UNIT_errOut,abort
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: ios
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  FLUSH(Unit_stdOut, IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(Unit_errOut, '(A)') 'ERROR: could not flush standard output'
    CALL abort(__STAMP__,"")
  END IF
END SUBROUTINE flush_stdout

!================================================================================================================================!
SUBROUTINE redirect_abort()
  ! MODULES
  USE MODgvec_Globals, ONLY: RaiseExceptionPtr
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  RaiseExceptionPtr => F90WRAP_ABORT
END SUBROUTINE redirect_abort

END MODULE MODgvec_py_binding
