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
PROGRAM TEST_GVEC_TO_HOPR
USE MODgvec_Globals
USE MODgvec_gvec_to_hopr
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: i,nArgs,SFL
CHARACTER(LEN=255)      :: filename 
REAL(wp)                :: StartTime,EndTime
REAL(wp)                :: xin(3,4),xout(3,4),data_out(9,4)
REAL(wp)                :: phi_edge_axis(2) 
REAL(wp)                :: chi_edge_axis(2) 
!===================================================================================================================================
  CALL CPU_TIME(StartTime)
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.GE.1)THEN
    CALL GET_COMMAND_ARGUMENT(1,filename)
  ELSE
    STOP ' TEST GVEC TO HOPR: gvec filename not given, usage: "./executable gvec_file.dat"'
  END IF
    
  
  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(5(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - TEST GVEC => HOPR - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
  WRITE(Unit_stdOut,'(132("="))')
  
  DO SFL=0,2
    !initialization phase
    CALL Init_gvec_to_hopr(filename,SFLcoord_in=SFL,factorSFL_in=2)
   
    WRITE(*,*)'===> SFLcoord: ',SFL
    xin(:,1)=(/0.0,0.5,0.3/)
    xin(:,2)=(/0.3,0.13,0.65/)
    xin(:,3)=(/0.6,0.43,0.15/)
    xin(:,4)=(/1.0,-0.33,-0.45/)
    CALL gvec_to_hopr(4,xin,xout,data_out,phi_edge_axis,chi_edge_axis)
    WRITE(*,*)'phi_edge_axis: ',phi_edge_axis
    WRITE(*,*)'chi_edge_axis: ',chi_edge_axis
    DO i=1,4
      WRITE(*,*)'s,thet,zeta: ',xin(:,i)
      WRITE(*,*)'x,y,z      : ',xout(:,i)
      WRITE(*,*)'pressure   : ',data_out(1,i)
      WRITE(*,*)'Bcart      : ',data_out(2:4,i)
      WRITE(*,*)'|B|        : ',SQRT(SUM(data_out(2:4,i)**2))
      WRITE(*,*)'chi,phi    : ',data_out(5:6,i)
      WRITE(*,*)'Acart      : ',data_out(7:9,i)
      WRITE(*,*)'-----------------------'
    END DO
    
    CALL Finalize_gvec_to_hopr()

  END DO
  CALL CPU_TIME(EndTime)
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' TEST GVEC TO HOPR FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)

END PROGRAM TEST_GVEC_TO_HOPR


