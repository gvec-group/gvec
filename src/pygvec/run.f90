!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

MODULE MODgvec_py_run

#ifndef NOISOENV
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
#endif

IMPLICIT NONE
PUBLIC

LOGICAL :: initialized = .FALSE.

CONTAINS

!================================================================================================================================!
SUBROUTINE start_rungvec(parameterfile,restartfile_in,comm_in)
  ! MODULES
  USE MODgvec_Globals, ONLY: Unit_stdOut
  USE MODgvec_MPI    , ONLY: par_Init,par_finalize
  USE MODgvec_rungvec, ONLY: rungvec
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=*),INTENT(IN) :: parameterfile
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: restartfile_in
  INTEGER,INTENT(IN),OPTIONAL :: comm_in
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: comm
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  initialized = .TRUE.
  IF(PRESENT(comm_in)) THEN
    CALL par_init(comm_in)
  ELSE
    CALL par_init() !USE MPI_COMM_WORLD
  END IF

  IF(PRESENT(restartfile_in))THEN
    CALL rungvec(parameterfile,restartfile_in=restartfile_in)
  ELSE
    CALL rungvec(parameterfile)
  END IF

  CALL par_finalize()
  initialized = .FALSE.
END SUBROUTINE start_rungvec
END MODULE MODgvec_py_run
