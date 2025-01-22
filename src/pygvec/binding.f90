!================================================================================================================================!
! Copyright (C) 2025 Robert Babin <robert.babin@ipp.mpg.de>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!================================================================================================================================!
! PyGVEC general bindings
!
! This module contains the utility functions for the Python bindings of GVEC.
!================================================================================================================================!

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

  OPEN(Unit_stdOut, FILE=filename, STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(Unit_errOut, '(A)') 'ERROR: could not open file', filename, 'for writing'
    CALL abort(__STAMP__,"")
  END IF
END SUBROUTINE redirect_stdout

!================================================================================================================================!
SUBROUTINE redirect_abort()
  ! MODULES
  USE MODgvec_Globals, ONLY: RaiseExceptionPtr
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  RaiseExceptionPtr => F90WRAP_ABORT
END SUBROUTINE redirect_abort

END MODULE MODgvec_py_binding
