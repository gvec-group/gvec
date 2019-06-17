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
!!# **GVEC_TO_CASTOR3D**  converter program 
!!
!===================================================================================================================================
PROGRAM CONVERT_GVEC_TO_CASTOR3D
USE MODgvec_Globals
USE MODgvec_gvec_to_castor3d
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: err,nArgs,Ns,FourierFactor
CHARACTER(LEN=255)      :: tmp 
LOGICAL                 :: commandFailed
CHARACTER(LEN=255)      :: FileName 
CHARACTER(LEN=255)      :: FileNameOut
REAL(wp)                :: StartTime,EndTime
!===================================================================================================================================
  CALL CPU_TIME(StartTime)
  nArgs=COMMAND_ARGUMENT_COUNT()
  
  
  commandFailed=.FALSE.
  IF(nArgs.EQ.3)THEN
    CALL GET_COMMAND_ARGUMENT(1,tmp)
    READ(tmp,*,IOSTAT=err)Ns
    IF(err.NE.0) commandFailed=.TRUE.
    CALL GET_COMMAND_ARGUMENT(2,tmp)
    READ(tmp,*,IOSTAT=err)FourierFactor
    IF(err.NE.0) commandFailed=.TRUE.
    CALL GET_COMMAND_ARGUMENT(3,FileName)
  ELSE
    commandFailed=.TRUE.
  END IF
  IF(commandfailed)  STOP ' GVEC TO CASTOR3D: command not correct: "./executable Ns FourierFactor gvec_file.dat"'
  IF(Ns.LT.2) STOP ' GVEC TO CASTOR3D: choose number in radial dierction Ns>=2' 
  IF(FourierFactor.GT.16) STOP ' GVEC TO CASTOR3D: choose fourierFactor  <=16' 
    
!  ppos = INDEX(FileName,'.dat',back=.TRUE.)
!  IF ( ppos  .GT. 0 )THEN
!     FileNameOut = FileName(1:ppos-1)
!  ELSE
!    STOP ' did not find .dat extension in input file'
!  END IF
  FileNameOut='gvec2castor3d_'//TRIM(FileName)
    
  
  
  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(5(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - CONVERT GVEC => CASTOR3D  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
  WRITE(Unit_stdOut,'(132("="))')
  
  !initialization phase
  CALL Init_gvec_to_castor3d(FileName,Ns,FourierFactor)
 
  CALL gvec_to_castor3d_prepare()

  CALL gvec_to_castor3d_writeToFile(FileNameOut)
  CALL Finalize_gvec_to_castor3d()

  CALL CPU_TIME(EndTime)
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' CONVERT GVEC TO CASTOR3D FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)

END PROGRAM CONVERT_GVEC_TO_CASTOR3D


