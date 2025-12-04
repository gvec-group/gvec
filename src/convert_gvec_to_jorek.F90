!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!>
!!# **GVEC_TO_JOREK**  converter program
!!
!===================================================================================================================================
PROGRAM CONVERT_GVEC_TO_JOREK
USE MODgvec_Globals
USE MODgvec_gvec_to_jorek
!$ USE omp_lib
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
REAL(wp)                :: StartTime,EndTime
!===================================================================================================================================
  __PERFINIT
  __PERFON('main')

  CALL CPU_TIME(StartTime)
!$ StartTime=OMP_GET_WTIME()

  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(5(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - CONVERT GVEC => JOREK  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
  WRITE(Unit_stdOut,'(132("="))')
!$ WRITE(UNIT_stdOut,'(A,I6)')'   Number of OpenMP threads : ',OMP_GET_MAX_THREADS()
!$ WRITE(Unit_stdOut,'(132("="))')
  CALL GET_CLA_gvec_to_jorek()

  !initialization phase
  CALL Init_gvec_to_jorek()

  CALL gvec_to_jorek_writeToFile()
  CALL Finalize_gvec_to_jorek()

  CALL CPU_TIME(EndTime)
!$ EndTime=OMP_GET_WTIME()
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' CONVERT GVEC TO JOREK FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)

  __PERFOFF('main')
  __PERFOUT('main')
END PROGRAM CONVERT_GVEC_TO_JOREK
