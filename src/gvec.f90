!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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
USE MOD_Globals
USE MOD_MHDEQ      ,ONLY: InitMHDEQ,FinalizeMHDEQ
USE MOD_Analyze    ,ONLY: InitAnalyze,Analyze,FinalizeAnalyze
USE MOD_Output     ,ONLY: InitOutput,FinalizeOutput
USE MOD_Functional ,ONLY: InitFunctional,FinalizeFunctional
USE MOD_Output     ,ONLY: InitOutput,FinalizeOutput
USE MOD_ReadInTools,ONLY: GETLOGICAL,GETINT,IgnoredStrings 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: i,nArgs
CHARACTER(LEN=255)      :: Parameterfile
!===================================================================================================================================
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.EQ.1)THEN
    CALL GET_COMMAND_ARGUMENT(1,Parameterfile)
  ELSE
    STOP 'parameterfile not given, usage: "./executable parameter.ini"'
  END IF
  
  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(18(("*",A128,2X,"*",:,"\n")))')&
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
  WRITE(Unit_stdOut,'(132("="))')
  testdbg =GETLOGICAL('testdbg','F')
  testlevel=GETINT('testlevel','-1')
  
  !initialization phase
  CALL InitOutput()
  CALL InitAnalyze()
  CALL InitMHDEQ()
  
  CALL InitFunctional()
  
  CALL IgnoredStrings()
  ! do something
  IF(testlevel.GT.0)THEN
    SWRITE(UNIT_stdout,*)
    SWRITE(UNIT_stdOut,'(A)')"** TESTESTESTESTESTESTESTESTESTESTESTESTESTESTESTEST **"
    SWRITE(UNIT_stdout,*)
    IF(nFailedMsg.GT.0)THEN
      DO i=1,nFailedMsg
        SWRITE(UNIT_stdOut,'(A)')TRIM(TestFailedMsg(i))
      END DO
      SWRITE(UNIT_stdOut,'(A)')"!!!!!!!   SOME TEST(S) FAILED !!!!!!!!!!!!!"
    ELSE
      SWRITE(UNIT_stdOut,'(A)')"   ...   ALL TESTS SUCCESSFULL ..."
    END IF !nFailedMsg
    SWRITE(UNIT_stdout,*)
    SWRITE(UNIT_stdOut,'(A)')"** TESTESTESTESTESTESTESTESTESTESTESTESTESTESTESTEST **"
    SWRITE(UNIT_stdout,*)
  END IF !testlevel
  
  !finalization phase
  CALL Analyze()
  
  CALL FinalizeFunctional()
  
  CALL FinalizeMHDEQ()
  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A)') ' GVEC SUCESSFULLY FINISHED'
  WRITE(Unit_stdOut,fmt_sep)

END PROGRAM GVEC


