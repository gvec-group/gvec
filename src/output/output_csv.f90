!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **Output CSV**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_Output_CSV
! MODULES
USE MODgvec_Globals, ONLY: wp
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
SUBROUTINE WriteDataToCSV(VarNames,Values,FileString,append_in,vfmt_in)
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(:)          !! Variable names,
REAL(wp),INTENT(IN)           :: Values(:,:)      !! variable data
CHARACTER(LEN=*),INTENT(IN)   :: FileString              !! Output file name
LOGICAL,INTENT(IN),OPTIONAL   :: append_in                  !! append data
CHARACTER(LEN=*),INTENT(IN),OPTIONAL   :: vfmt_in           !! value format
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nVal                    !! Number of output Values
INTEGER                       :: nPlot                   !! number of 1D values
INTEGER                        :: iVal,iPlot,ioUnit
LOGICAL                        :: append
CHARACTER(LEN=10)              :: vfmt
!===================================================================================================================================
nVal=SIZE(VarNames,1)
IF(SIZE(Values,1).NE.nVal) STOP 'number of values /= nVariables in csv output'
nPlot=SIZE(Values,2)
ioUnit=GETFREEUNIT()
IF(PRESENT(append_in))THEN
  append=append_in
ELSE
  append=.FALSE.
END IF
IF(PRESENT(vfmt_in))THEN
  vfmt=vfmt_in
ELSE
  vfmt='e23.15'
END IF
IF(append)THEN
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'     APPEND DATA TO CSV FILE "'//TRIM(FileString)//'" ...'
  OPEN(UNIT     = ioUnit       ,&
       FILE     = TRIM(FileString)   ,&
       STATUS   = 'OLD'   ,&
       POSITION = 'APPEND'   ,&
       ACCESS   = 'SEQUENTIAL' )
ELSE
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITE DATA TO CSV FILE    "'//TRIM(FileString)//'" ...'
  OPEN(UNIT     = ioUnit       ,&
       FILE     = TRIM(FileString)   ,&
       STATUS   = 'REPLACE'   ,&
       ACCESS   = 'SEQUENTIAL' )
END IF


DO iVal=1,nVal-1
  WRITE(ioUnit,'(A,1X,(","))',ADVANCE='NO' )  '"'//TRIM(varNames(iVal))//'"'
END DO
WRITE(ioUnit,'(A)')  '"'//TRIM(varNames(nVal))//'"'

DO iPlot=1,nPlot
  WRITE(ioUnit,'(*('//TRIM(vfmt)//',:,","))') Values(:,iPlot)+SIGN(1.0e-99,Values(:,iPlot))
END DO

CLOSE(ioUnit)

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToCSV



END MODULE MODgvec_Output_CSV
