!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!>
!!# **GVEC** Driver program
!!
!===================================================================================================================================
PROGRAM GVEC
  USE_MPI
  USE MODgvec_Globals
  USE MODgvec_MPI    ,ONLY  : par_Init,par_finalize
  USE MODgvec_rungvec, ONLY : rungvec
  USE MODgvec_cla
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  !local variables
  CHARACTER(LEN=STRLEN)   :: Parameterfile=""
  CHARACTER(LEN=STRLEN)   :: RestartFile=""
  CHARACTER(LEN=STRLEN)   :: f_str
  CHARACTER(LEN=24)       :: execname="gvec"
  LOGICAL                 :: commandFailed
  LOGICAL                 :: doRestart
  !===================================================================================================================================
  CALL enter_subregion("startup")
  CALL par_Init() !USE MPI_COMM_WORLD

  IF(MPIroot)THEN
    !USING CLAF90 module to get command line arguments!
    CALL cla_init()

    CALL cla_posarg_register('parameterfile', &
         ' path/filename of GVEC input parameter file [MANDATORY]',  cla_char,'') !
    CALL cla_posarg_register('statefile', &
         'path/filename of a GVEC State. if provided, a restart from this State is done  [OPTIONAL, DEFAULT: None provided]',  cla_char,'')
    CALL cla_validate('gvec')
    CALL cla_get('parameterfile',f_str)
    parameterfile=TRIM(f_str)
    CALL cla_get('statefile',f_str)
    RestartFile=TRIM(f_str)
    commandFailed=.FALSE.
    IF((LEN(TRIM(parameterfile)).EQ.0))THEN
      IF(.NOT.commandFailed) CALL cla_help(execname)
      commandFailed=.TRUE.
      SWRITE(UNIT_StdOut,*) " ==> input parameter filename is MANDATORY must be specified as first positional argument!!!"
    END IF
    dorestart=(LEN(TRIM(restartFile)).GT.0)

    IF(commandFailed) STOP " ...check your command line arguments!"
  END IF !MPIroot
  CALL exit_subregion("startup")
  IF(dorestart)THEN
    CALL rungvec(parameterfile,restartfile_in=restartfile)
  ELSE
    CALL rungvec(parameterfile)
  END IF

  CALL par_finalize()
END PROGRAM GVEC
