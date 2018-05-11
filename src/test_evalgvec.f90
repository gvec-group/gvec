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
PROGRAM TEST_EVALGVEC
USE MODgvec_Globals
USE MODgvec_Eval_GVEC    ,ONLY: InitEval_GVEC,Eval_GVEC,FinalizeEval_GVEC
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: nArgs
CHARACTER(LEN=255)      :: filename 
REAL(wp)                :: StartTime,EndTime
REAL(wp)                :: xin(3,2),xout(3,2),data_out(9,2)
REAL(wp)                :: phi_edge_axis(2) 
REAL(wp)                :: chi_edge_axis(2) 
!===================================================================================================================================
  CALL CPU_TIME(StartTime)
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.GE.1)THEN
    CALL GET_COMMAND_ARGUMENT(1,filename)
  ELSE
    STOP 'EVALGVEC: gvec filename not given, usage: "./executable gvec_file.dat"'
  END IF
    
  
  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(18(("*",A128,2X,"*",:,"\n")))')&
 '  EVAL- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' -EVAL - - - - - - - - - GGGGGGGGGGGGGGG - VVVVVVVV  - - - -  VVVVVVVV - EEEEEEEEEEEEEEEEEEEEEE  - - - - CCCCCCCCCCCCCCC - -  '&
,'  EVAL- - - - - - - - GGGG::::::::::::G - V::::::V  - - - -  V::::::V - E::::::::::::::::::::E  - - - CCCC::::::::::::C - - - '&
,' -EVAL - - - - - - GGG:::::::::::::::G - V::::::V  - - - -  V::::::V - E::::::::::::::::::::E  - - CCC:::::::::::::::C - - -  '&
,'  EVAL- - - - -  GG:::::GGGGGGGG::::G - V::::::V  - - - -  V::::::V - EEE:::::EEEEEEEEE::::E  -  CC:::::CCCCCCCC::::C - - - - '&
,' -EVAL - - - - GG:::::GG  - - GGGGGG -  V:::::V  - - - -  V:::::V  - - E:::::E - - - EEEEEE  - CC:::::CC - -  CCCCCC - - - -  '&
,'  EVAL- - - - G:::::GG  - - - - - - - - V:::::V - - - - V:::::V - - - E:::::E - - - - - - - - C:::::CC  - - - - - - - - - - - '&
,' -EVAL - - - G:::::G - - - - - - - - -  V:::::V  - -  V:::::V  - - - E:::::EEEEEEEEEEE - - - C:::::C - - - - - - - - - - - -  '&
,'  EVAL- - - G:::::G -  GGGGGGGGGG - - - V:::::V - - V:::::V - - - - E:::::::::::::::E - - - C:::::C - - - - - - - - - - - - - '&
,' -EVAL - - G:::::G -  G::::::::G - - -  V:::::V   V:::::V  - - - - E:::::::::::::::E - - - C:::::C - - - - - - - - - - - - -  '&
,'  EVAL- - G:::::G -  GGGGG::::G - - - - V:::::V V:::::V - - - - - E:::::EEEEEEEEEEE - - - C:::::C - - - - - - - - - - - - - - '&
,' -EVAL - G:::::G - - -  G::::G - - - -  V:::::V:::::V  - - - - - E:::::E - - - - - - - - C:::::C - - - - - - - - - - - - - -  '&
,'  EVAL-  G:::::G  - -  G::::G - - - - - V:::::::::V - - - - - - E:::::E - - - EEEEEE  -  C:::::C  - -  CCCCCC - - - - - - - - '&
,' -EVAL - G::::::GGGGGGG::::G - - - - -  V:::::::V  - - - - - EEE:::::EEEEEEEEE::::E  - - C::::::CCCCCCC::::C - - - - - - - -  '&
,'  EVAL- - G:::::::::::::::G - - - - - - V:::::V - - - - - - E::::::::::::::::::::E  - - - C:::::::::::::::C - - - - - - - - - '&
,' -EVAL - - GG::::GGGG::::G - - - - - -  V:::V  - - - - - - E::::::::::::::::::::E  - - - - CC::::::::::::C - - - - - - - - -  '&
,'  EVAL- - -  GGGG  GGGGGG - - - - - - - VVV - - - - - - - EEEEEEEEEEEEEEEEEEEEEE  - - - - -  CCCCCCCCCCCC - - - - - - - - - - '&
,' -EVAL - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '
  WRITE(Unit_stdOut,'(132("="))')
  
  !initialization phase
  CALL InitEval_GVEC(filename)
 
  xin(:,1)=(/0.,0.5,0.3/)
  xin(:,2)=(/1.,-0.33,-0.45/)
  CALL Eval_GVEC(2,xin,xout,data_out,phi_edge_axis,chi_edge_axis)
  WRITE(*,*)'xin',xin(:,1)
  WRITE(*,*)'xout',xout(:,1)
  WRITE(*,*)'data_out',data_out(:,1)
  WRITE(*,*)'xin',xin(:,2)
  WRITE(*,*)'xout',xout(:,2)
  WRITE(*,*)'data_out',data_out(:,2)

  CALL FinalizeEval_GVEC()

  CALL CPU_TIME(EndTime)
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' EVALGVEC FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)

END PROGRAM TEST_EVALGVEC


