!===================================================================================================================================
! Copyright (c) 2017 - 2018 Florian Hindenlang <hindenlang@gmail.com>
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/hopr)
!
! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
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
!!# Module **Globals**
!!
!! Here globally used variables /functions are defined 
!!
!===================================================================================================================================
MODULE MODgvec_Globals

#ifndef NOISOENV
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
#endif
IMPLICIT NONE

PUBLIC 

!-----------------------------------------------------------------------------------------------------------------------------------
! Select here the working precision wp
!INTEGER, PARAMETER :: wp = selected_real_kind(6,35)   !! single precision
INTEGER, PARAMETER :: wp = selected_real_kind(15,307)  !! double precision
!INTEGER, PARAMETER :: wp = selected_real_kind(33,307) !! quadruple precision
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=20)   :: fmt_sep ='(132("="))'             !! formatting of separator line: WRITE(*,fmt_sep)
REAL(wp),PARAMETER  :: PI   =ACOS(-1.0_wp)               !! pi parameter
REAL(wp),PARAMETER  :: TWOPI=2.0_wp*PI                   !! 2*pi parameter
INTEGER             :: n_warnings_occured=0              !! for final line in screen output: 0 no warnings occured
!-----------------------------------------------------------------------------------------------------------------------------------
!for testing
LOGICAL                     :: testdbg=.FALSE.           !! for debugging the tests, set true for implementing tests, false to run
INTEGER                     :: testlevel =-1             !! flag for testing routines in code: -1: off
INTEGER                     :: ntestCalled=0             !! counter for called tests
INTEGER                     :: nfailedMsg=0              !! counter for messages on failed tests 
INTEGER                     :: testUnit                  !! unit for out.test file
INTEGER                     :: ProgressBar_oldpercent    !! for progressBar
REAL(wp)                    :: ProgressBar_starttime     !! for progressBar
!-----------------------------------------------------------------------------------------------------------------------------------
#ifndef NOISOENV
INTEGER, PARAMETER          :: UNIT_stdIn  = input_unit  !! Terminal input
INTEGER, PARAMETER          :: UNIT_stdOut = output_unit !! Terminal output
INTEGER, PARAMETER          :: UNIT_errOut = error_unit  !! For error output
#else                                                     
INTEGER, PARAMETER          :: UNIT_stdIn  = 5           !! Terminal input
INTEGER, PARAMETER          :: UNIT_stdOut = 6           !! Terminal output
INTEGER, PARAMETER          :: UNIT_errOut = 0           !! For error output
#endif


INTERFACE Abort
   MODULE PROCEDURE Abort
END INTERFACE

INTERFACE ProgressBar
   MODULE PROCEDURE ProgressBar
END INTERFACE

INTERFACE GETFREEUNIT
   MODULE PROCEDURE GETFREEUNIT
END INTERFACE

INTERFACE Eval1DPoly 
  MODULE PROCEDURE Eval1DPoly
END INTERFACE

INTERFACE CROSS
   MODULE PROCEDURE CROSS
END INTERFACE

INTERFACE NORMALIZE
   MODULE PROCEDURE NORMALIZE
END INTERFACE

INTERFACE DET33
   MODULE PROCEDURE DET33
END INTERFACE

INTERFACE INV33
   MODULE PROCEDURE INV33
END INTERFACE

CONTAINS

!==================================================================================================================================
!> Terminate program correctly if an error has occurred (important in MPI mode!).
!! Uses a MPI_ABORT which terminates FLUXO if a single proc calls this routine.
!!
!==================================================================================================================================
SUBROUTINE Abort(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfo,RealInfo,ErrorCode)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      !! Source file where error has occurred
INTEGER                           :: SourceLine      !! Line in source file
CHARACTER(LEN=*)                  :: CompDate        !! Compilation date
CHARACTER(LEN=*)                  :: CompTime        !! Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    !! Error message
INTEGER,OPTIONAL                  :: IntInfo         !! Error info (integer)
REAL(wp),OPTIONAL                 :: RealInfo        !! Error info (real)
INTEGER,OPTIONAL                  :: ErrorCode       !! Error info (integer)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)                 :: IntString,RealString
#if MPI
INTEGER                           :: errOut          ! Output of MPI_ABORT
INTEGER                           :: signalout       ! Output errorcode
#endif
!==================================================================================================================================
IntString = ""
RealString = ""

IF (PRESENT(IntInfo))  WRITE(IntString,"(A,I0)")  "\nIntInfo:  ", IntInfo
IF (PRESENT(RealInfo)) WRITE(RealString,"(A,F24.19)") "\nRealInfo: ", RealInfo

WRITE(UNIT_stdOut,*) '_____________________________________________________________________________\n', &
#if MPI
                     'Program abort caused on Proc ',myRank, '\n', &
#else
                     'Program abort caused on Proc ',0, '\n', &
#endif  
                     '  in File : ',TRIM(SourceFile),' Line ',SourceLine, '\n', &
                     '  This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime), '\n', &
                     'Message: ',TRIM(ErrorMessage), &
                     TRIM(IntString), TRIM(RealString)

CALL FLUSH(UNIT_stdOut)
#if MPI
signalout=2 ! MPI_ABORT requires an output error-code /=0
IF(PRESENT(ErrorCode)) signalout=ErrorCode
CALL MPI_ABORT(MPI_COMM_WORLD,signalout,errOut)
#endif  
#if GNU
CALL BACKTRACE
#endif
ERROR STOP 2
END SUBROUTINE Abort


!==================================================================================================================================
!> Print a progress bar to screen, call either with init=T or init=F
!!
!==================================================================================================================================
SUBROUTINE ProgressBar(iter,n_iter)
! MODULES
!$ USE omp_lib
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iter,n_iter  !! iter ranges from 1...n_iter
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=8)  :: fmtstr
INTEGER           :: newpercent
REAL(wp)          :: endTime
!==================================================================================================================================
  IF(iter.LE.0)THEN !INIT
    ProgressBar_oldpercent=0
    CALL CPU_Time(Progressbar_StartTime)
!$ Progressbar_StartTime=OMP_GET_WTIME()
    SWRITE(UNIT_StdOut,'(4X,A,I8)') &
    '|       10%       20%       30%       40%       50%       60%       70%       80%       90%      100%| ... of ',n_iter
    SWRITE(UNIT_StdOut,'(4X,A1)',ADVANCE='NO')'|'
    CALL FLUSH(UNIT_stdOut)
  ELSE
    newpercent=FLOOR(REAL(iter,wp)/REAL(n_iter,wp)*(100.0_wp+1.0e-12_wp))
    WRITE(fmtstr,'(I4)')newpercent-ProgressBar_oldpercent
    IF(newpercent-ProgressBar_oldpercent.GT.0)THEN
      SWRITE(UNIT_StdOut,'('//TRIM(fmtstr)//'("."))',ADVANCE='NO')
      CALL FLUSH(UNIT_stdOut)
    END IF
    ProgressBar_oldPercent=newPercent
    IF(newpercent.EQ.100)THEN
      CALL CPU_Time(endTime)
!$  endTime=OMP_GET_WTIME()
      SWRITE(Unit_stdOut,'(A3,F8.2,A)') '| [',EndTime-ProgressBar_StartTime,' sec ]'
    END IF
  END IF 
END SUBROUTINE ProgressBar

!==================================================================================================================================
!> Get unused file unit number
!==================================================================================================================================
FUNCTION GETFREEUNIT()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER :: GetFreeUnit !! File unit number
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: connected
!==================================================================================================================================
GetFreeUnit=55
INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
IF(connected)THEN
  DO
    GetFreeUnit=GetFreeUnit+1
    INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
    IF(.NOT.connected)EXIT
  END DO
END IF
END FUNCTION GETFREEUNIT


!===================================================================================================================================
!> evalute monomial polynomial c_1+c_2*x+c_3*x^2 ...
!!
!===================================================================================================================================
FUNCTION Eval1DPoly(nCoefs,Coefs,x)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,  INTENT(IN)  :: nCoefs                   !! number of coefficients 
REAL(wp), INTENT(IN)  :: Coefs(nCoefs)            !! coefficients
REAL(wp), INTENT(IN)  :: x                        !! evaluation position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)              :: Eval1DPoly
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!===================================================================================================================================
Eval1DPoly=0.
DO i=nCoefs,1,-1
  Eval1DPoly=Eval1DPoly*x+Coefs(i)
END DO

END FUNCTION Eval1DPoly


!===================================================================================================================================
!> normalizes a nDim vector with respect to the eucledian norm
!!
!===================================================================================================================================
PURE FUNCTION NORMALIZE(v1,nVal)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nVal     !! vector size
REAL(wp),INTENT(IN) :: v1(nVal) !! vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)            :: normalize(nVal) !! result, normalized vector
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
normalize=v1/SQRT(SUM(v1*v1))
END FUNCTION NORMALIZE


!===================================================================================================================================
!> computes the cross product of to 3 dimensional vectors: cross=v1 x v2
!!
!===================================================================================================================================
PURE FUNCTION CROSS(v1,v2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: v1(3) !! first input vector 
REAL(wp),INTENT(IN) :: v2(3) !! second input vector  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)            :: cross(3)  !! result v1 x v2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
cross=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS


!===================================================================================================================================
!> compute determinant of 3x3 matrix
!!
!===================================================================================================================================
FUNCTION DET33(Mat)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)  :: Mat(3,3) !! input matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)             :: DET33 !! determinant of the input matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DET33=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
         + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
         + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION DET33


!===================================================================================================================================
!> compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
!!
!===================================================================================================================================
FUNCTION INV33(Mat, Det_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)             :: Mat(3,3) !! input matrix
REAL(wp),INTENT(IN),OPTIONAL    ::  Det_in  !! determinant of input matrix (otherwise computed here)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)             :: INV33(3,3) !! inverse matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)             :: sDet 
!===================================================================================================================================
IF(PRESENT(Det_in))THEN
  sDet=1.0_wp/Det_in
ELSE
  sDet=1.0_wp/DET33(Mat)
END IF
INV33(1,1) = ( Mat(2,2) * Mat(3,3) - Mat(2,3) * Mat(3,2) ) * sDet
INV33(1,2) = ( Mat(1,3) * Mat(3,2) - Mat(1,2) * Mat(3,3) ) * sDet
INV33(1,3) = ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * sDet
INV33(2,1) = ( Mat(2,3) * Mat(3,1) - Mat(2,1) * Mat(3,3) ) * sDet
INV33(2,2) = ( Mat(1,1) * Mat(3,3) - Mat(1,3) * Mat(3,1) ) * sDet
INV33(2,3) = ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * sDet
INV33(3,1) = ( Mat(2,1) * Mat(3,2) - Mat(2,2) * Mat(3,1) ) * sDet
INV33(3,2) = ( Mat(1,2) * Mat(3,1) - Mat(1,1) * Mat(3,2) ) * sDet
INV33(3,3) = ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * sDet

END FUNCTION INV33


END MODULE MODgvec_Globals
