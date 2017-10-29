!===================================================================================================================================
! Copyright (c) 2017 - 2018 Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **Output CSV**
!!
!! 
!!
!===================================================================================================================================
MODULE MOD_Output_CSV
! MODULES
USE MOD_Globals, ONLY: wp
IMPLICIT NONE
PRIVATE

INTERFACE WriteDataToCSV
  MODULE PROCEDURE WriteDataToCSV
END INTERFACE

PUBLIC::WriteDataToCSV
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write  
!!
!===================================================================================================================================
SUBROUTINE WriteDataToCSV(nVal, nPlot,VarNames,Values,FileString)
! MODULES
USE MOD_Globals,ONLY:Unit_stdOut,GETFREEUNIT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVal                    !! Number of output Values 
INTEGER,INTENT(IN)            :: nPlot                   !! number of 1D values 
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          !! Variable names, 
REAL,INTENT(IN)               :: Values(nVal,nPlot)      !! variable data
CHARACTER(LEN=*),INTENT(IN)   :: FileString              !! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVal,iPlot,ioUnit
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE DATA TO CSV FILE... "//TRIM(FileString)//".csv"
ioUnit=GETFREEUNIT()
OPEN(UNIT   = ioUnit       ,&
     FILE   = TRIM(FileString)//".csv"   ,&
     STATUS = 'NEW'    ,&
     ACCESS = 'SEQUENTIAL' ) 

WRITE(ioUnit,'(A)')   '# TITLE="Analysis,'//TRIM(FileString)//'"'
WRITE(ioUnit,'(A,I8)')'# nPlot=',nPlot

DO iVal=1,nVal-1
  WRITE(ioUnit,'(A,1X,(","),1X)',ADVANCE='NO' )  '"'//TRIM(varNames(iVal))//'"'
END DO
WRITE(ioUnit,'(A)')  '"'//TRIM(varNames(nVal))//'"'


DO iPlot=1,nPlot
  WRITE(ioUnit,'(*(e23.15,:,","))') Values(:,iPlot) 
END DO

CLOSE(ioUnit)

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToCSV
 
 

END MODULE MOD_Output_CSV
