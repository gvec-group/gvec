!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!>
!!# **GVEC** Driver Module
!!
!===================================================================================================================================
MODULE MODgvec_rungvec

USE MODgvec_Functional, ONLY: t_functional

IMPLICIT NONE
PUBLIC

CLASS(t_functional), ALLOCATABLE :: functional
LOGICAL :: dorestart

!INTERFACE rungvec
!  MODULE PROCEDURE rungvec
!END INTERFACE
!===================================================================================================================================
CONTAINS

SUBROUTINE rungvec(parameterFile,restartfile_in)
USE_MPI
USE MODgvec_Globals
USE MODgvec_Analyze    ,ONLY : InitAnalyze,FinalizeAnalyze
USE MODgvec_Output     ,ONLY : InitOutput,FinalizeOutput
USE MODgvec_Restart    ,ONLY : InitRestart,FinalizeRestart
USE MODgvec_ReadInTools,ONLY : FillStrings,GETLOGICAL,GETINT,IgnoredStrings,FinalizeReadIn
USE MODgvec_Functional ,ONLY : t_functional, InitFunctional,FinalizeFunctional
USE MODgvec_MHD3D_Vars ,ONLY : maxIter
!$ USE omp_lib
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!INPUT VARIABLES
  CHARACTER(LEN=*),INTENT(IN)             :: Parameterfile  !! input parameters of GVEC
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL    :: RestartFile_in  !! if present, a restart will be executed
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: which_functional
INTEGER                 :: TimeArray(8)
CHARACTER(LEN=255)      :: testfile
REAL(wp)                :: StartTimeTotal,EndTimeTotal,StartTime,EndTime
!===================================================================================================================================
  __PERFINIT
  __PERFON('main')
  CALL reset_subregion()

  StartTimeTotal=GetTime()
  SWRITE(Unit_stdOut,fmt_sep)
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  SWRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
    '%%% Sys date : ',timeArray(1:3),timeArray(5:7)

  !header
  SWRITE(Unit_stdOut,fmt_sep)
  SWRITE(Unit_stdOut,'(18(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - GGGGGGGGGGGGGGG - VVVVVVVV  - - - -  VVVVVVVV - EEEEEEEEEEEEEEEEEEEEEE  - - - - CCCCCCCCCCCCCCC - -  '&
,'  - - - - - - - - - - GGGG::::::::::::G - V::::::V  - - - -  V::::::V - E::::::::::::::::::::E  - - - CCCC::::::::::::C - - - '&
,' - - - - - - - - - GGG:::::::::::::::G - V::::::V  - - - -  V::::::V - E::::::::::::::::::::E  - - CCC:::::::::::::::C - - -  '&
,'  - - - - - - -  GG:::::GGGGGGGG::::G - V::::::V  - - - -  V::::::V - EEE:::::EEEEEEEEE::::E  -  CC:::::CCCCCCCC::::C - - - - '&
,' - - - - - - - GG:::::GG  - - GGGGGG -  V:::::V  - - - -  V:::::V  - - E:::::E - - - EEEEEE  - CC:::::CC - -  CCCCCC - - - -  '&
,'  - - - - - - G:::::GG  - - - - - - - - V:::::V - - - - V:::::V - - - E:::::E - - - - - - - - C:::::CC  - - - - - - - - - - - '&
,' - - - - - - G:::::G - - - - - - - - -  V:::::V  - -  V:::::V  - - - E:::::EEEEEEEEEEE - - - C:::::C - - - - - - - - - - - -  '&
,'  - - - - - G:::::G -  GGGGGGGGGG - - - V:::::V - - V:::::V - - - - E:::::::::::::::E - - - C:::::C - - - - - - - - - - - - - '&
,' - - - - - G:::::G -  G::::::::G - - -  V:::::V   V:::::V  - - - - E:::::::::::::::E - - - C:::::C - - - - - - - - - - - - -  '&
,'  - - - - G:::::G -  GGGGG::::G - - - - V:::::V V:::::V - - - - - E:::::EEEEEEEEEEE - - - C:::::C - - - - - - - - - - - - - - '&
,' - - - - G:::::G - - -  G::::G - - - -  V:::::V:::::V  - - - - - E:::::E - - - - - - - - C:::::C - - - - - - - - - - - - - -  '&
,'  - - -  G:::::G  - -  G::::G - - - - - V:::::::::V - - - - - - E:::::E - - - EEEEEE  -  C:::::C  - -  CCCCCC - - - - - - - - '&
,' - - - - G::::::GGGGGGG::::G - - - - -  V:::::::V  - - - - - EEE:::::EEEEEEEEE::::E  - - C::::::CCCCCCC::::C - - - - - - - -  '&
,'  - - - - G:::::::::::::::G - - - - - - V:::::V - - - - - - E::::::::::::::::::::E  - - - C:::::::::::::::C - - - - - - - - - '&
,' - - - - - GG::::GGGG::::G - - - - - -  V:::V  - - - - - - E::::::::::::::::::::E  - - - - CC::::::::::::C - - - - - - - - -  '&
,'  - - - - -  GGGG  GGGGGG - - - - - - - VVV - - - - - - - EEEEEEEEEEEEEEEEEEEEEE  - - - - -  CCCCCCCCCCCC - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '
  SWRITE(Unit_stdOut,fmt_sep)
  !.only executes if compiled with OpenMP
!$ SWRITE(UNIT_stdOut,'(A,I6)')'   Number of OpenMP threads : ',OMP_GET_MAX_THREADS()
!$ SWRITE(Unit_stdOut,'(132("="))')
  !.only executes if compiled with MPI
# if MPI
  SWRITE(UNIT_stdOut,'(A,I6)')'   Number of MPI tasks : ',nRanks
  SWRITE(Unit_stdOut,fmt_sep)
# endif
#include  "configuration-cmake.f90"
  SWRITE(Unit_stdOut,fmt_sep)
  CALL FillStrings(ParameterFile) !< readin parameterfile, done on MPI root + Bcast

  testdbg =GETLOGICAL('testdbg',Proposal=.FALSE.)
  testlevel=GETINT('testlevel',Proposal=-1)
  IF(testlevel.GT.0)THEN
    testUnit=GETFREEUNIT()
    WRITE(testFile,'(A,I4.4,A)')"tests_",myRank,".out"
    OPEN(UNIT     = testUnit    ,&
         FILE     = testfile    ,&
         STATUS   = 'REPLACE'   ,&
         ACCESS   = 'SEQUENTIAL' )
  END IF

  CALL enter_subregion("initialize")
  !initialization phase
  dorestart=.FALSE.
  IF(PRESENT(RestartFile_in)) THEN
    dorestart=(LEN(TRIM(RestartFile_in)).GT.0)
  END IF
  IF(dorestart) CALL InitRestart(RestartFile_in)
  CALL InitOutput()
  CALL InitAnalyze()
  which_functional=GETINT('which_functional', Proposal=1 )
  CALL InitFunctional(functional,which_functional)

  CALL IgnoredStrings()

  CALL functional%InitSolution()
  StartTime=GetTime()
  SWRITE(Unit_stdOut,'(A,F8.2,A)') ' INITIALIZATION FINISHED! [',StartTime-StartTimeTotal,' sec ]'
  SWRITE(Unit_stdOut,fmt_sep)

  CALL exit_subregion("initialize")

  CALL functional%minimize()
  EndTime=GetTime()
  SWRITE(Unit_stdOut,'(A,2(F8.2,A))') ' FUNCTIONAL MINIMISATION FINISHED! [',EndTime-StartTime,' sec ], corresponding to [', &
       (EndTime-StartTime)/REAL(MaxIter,wp)*1.e3_wp,' msec/iteration ]'

  CALL FinalizeFunctional(functional)
  DEALLOCATE(functional)
  CALL enter_subregion("finalize")
  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  IF(dorestart) CALL FinalizeRestart()
  CALL FinalizeReadIn()
  ! do something
  IF(testlevel.GT.0)THEN
    SWRITE(UNIT_stdout,*)
    SWRITE(UNIT_stdOut,'(A)')"** TESTESTESTESTESTESTESTESTESTESTESTESTESTESTESTEST **"
    SWRITE(UNIT_stdout,*)
    n_warnings_occured=n_warnings_occured +nFailedMsg
    IF(nFailedMsg.GT.0)THEN
      SWRITE(UNIT_stdOut,'(A)')"!!!!!!!   SOME TEST(S) FAILED, see tests.out !!!!!!!!!!!!!"
    ELSE
      SWRITE(UNIT_stdOut,'(A)')"   ...   ALL IMPLEMENTED TESTS SUCCESSFULL ..."
    END IF !nFailedMsg
    SWRITE(UNIT_stdout,*)
    SWRITE(UNIT_stdOut,'(A)')"** TESTESTESTESTESTESTESTESTESTESTESTESTESTESTESTEST **"
    SWRITE(UNIT_stdout,*)
    CLOSE(testUnit)
  END IF !testlevel
  CALL exit_subregion("finalize")
  EndTimeTotal=GetTime()
  SWRITE(Unit_stdOut,fmt_sep)
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  SWRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
    '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
  SWRITE(Unit_stdOut,fmt_sep)
  IF(n_warnings_occured.EQ.0)THEN
    SWRITE(Unit_stdOut,'(A,F8.2,A)') ' GVEC SUCESSFULLY FINISHED! [',EndTimeTotal-StartTimeTotal,' sec ]'
  ELSE
    SWRITE(Unit_stdOut,'(A,F8.2,A,I8,A)') ' GVEC FINISHED! [',EndTimeTotal-StartTimeTotal,' sec ], WITH ' , n_warnings_occured , ' WARNINGS!!!!'
  END IF
  SWRITE(Unit_stdOut,fmt_sep)
  __PERFOFF('main')
  IF(MPIRoot) THEN
  __PERFOUT('main')
  END IF
END SUBROUTINE rungvec
END MODULE MODgvec_rungvec
