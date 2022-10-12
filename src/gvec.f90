!===================================================================================================================================
! Copyright (C) 2017 - 2022  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2021 - 2022  Tiago Ribeiro
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!> 
!!# **GVEC** Driver program 
!!
!===================================================================================================================================
PROGRAM GVEC
USE_MPI
USE MODgvec_MPI,        ONLY : par_Init, par_Finalize, nRanks
USE MODgvec_Globals,    ONLY : wp,fmt_sep,n_warnings_occured,testdbg,testlevel,testUnit,nfailedMsg,UNIT_stdOut,GETFREEUNIT
USE MODgvec_Globals,    ONLY : GetTime,MPIRoot
USE MODgvec_Analyze    ,ONLY: InitAnalyze,FinalizeAnalyze
USE MODgvec_Output     ,ONLY: InitOutput,FinalizeOutput
USE MODgvec_Restart    ,ONLY: InitRestart,FinalizeRestart
USE MODgvec_ReadInTools,ONLY: FillStrings,GETLOGICAL,GETINT,IgnoredStrings 
USE MODgvec_Functional ,ONLY : t_functional, InitFunctional,FinalizeFunctional
USE MODgvec_MHD3D_Vars, ONLY : maxIter
!$ USE omp_lib
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: nArgs
CHARACTER(LEN=255)      :: Parameterfile
INTEGER                 :: which_functional
REAL(wp)                :: StartTimeTotal,EndTimeTotal,StartTime,EndTime
CLASS(t_functional),ALLOCATABLE   :: functional
!===================================================================================================================================
  CALL par_Init()
  __PERFINIT
  __PERFON('main')


  StartTimeTotal=GetTime()
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.GE.1)THEN
    CALL GET_COMMAND_ARGUMENT(1,Parameterfile)
  ELSE
    STOP 'parameterfile not given, usage: "./executable parameter.ini"'
  END IF
    
  
  !header
  SWRITE(Unit_stdOut,fmt_sep)
  SWRITE(Unit_stdOut,'(18(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - GGGGGGGGGGGGGGG - VVVVVVVV  - - - -  VVVVVVVV - EEEEEEEEEEEEEEEEEEEEEE  - - - - CCCCCCCCCCCCCCC - -  '&
,'  - - - - - - - - - - GGGG::::::::::::G - V::::::V  - - - -  V::::::V - E::::::::::::::::::::E  - - - CCCC::::::::::::C - - - '&
,' - - - OMP - - - - GGG:::::::::::::::G - V::::::V  - - - -  V::::::V - E::::::::::::::::::::E  - - CCC:::::::::::::::C - - -  '&
,'  - - MPI - - -  GG:::::GGGGGGGG::::G - V::::::V  - - - -  V::::::V - EEE:::::EEEEEEEEE::::E  -  CC:::::CCCCCCCC::::C - - - - '&
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
  SWRITE(Unit_stdOut,'(132("="))')
# endif
  CALL FillStrings(ParameterFile) !<<<<

  testdbg =GETLOGICAL('testdbg',Proposal=.FALSE.)
  testlevel=GETINT('testlevel',Proposal=-1)
  IF(testlevel.GT.0)THEN
    testUnit=GETFREEUNIT()
    OPEN(UNIT     = testUnit    ,&
         FILE     = "tests.out" ,&
         STATUS   = 'REPLACE'   ,&
         ACCESS   = 'SEQUENTIAL' ) 
  END IF
  
  !initialization phase
  CALL InitRestart()
  CALL InitOutput()
  CALL InitAnalyze()
  
  which_functional=GETINT('which_functional', Proposal=1 )
  CALL InitFunctional(functional,which_functional)

  
  CALL IgnoredStrings()
  
  CALL functional%InitSolution() 
  
  StartTime=GetTime()
  CALL functional%minimize()
  EndTime=GetTime()
  SWRITE(Unit_stdOut,'(A,2(F8.2,A))') ' FUNCTIONAL MINIMISATION FINISHED! [',EndTime-StartTime,' sec ], corresponding to [', &
       (EndTime-StartTime)/REAL(MaxIter,wp)*1.e3_wp,' msec/iteration ]'

  CALL FinalizeFunctional(functional)
 
  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  CALL FinalizeRestart()


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
  EndTimeTotal=GetTime()
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
  CALL par_Finalize()

END PROGRAM GVEC


