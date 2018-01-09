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
!=================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **Restart**
!!
!!
!!
!===================================================================================================================================
MODULE MOD_Restart
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitRestart
  MODULE PROCEDURE InitRestart
END INTERFACE

INTERFACE WriteState
  MODULE PROCEDURE WriteStateToASCII
END INTERFACE

INTERFACE ReadState
  MODULE PROCEDURE ReadStateFromASCII
END INTERFACE

INTERFACE FinalizeRestart
  MODULE PROCEDURE FinalizeRestart
END INTERFACE

PUBLIC::InitRestart
PUBLIC::WriteState
PUBLIC::ReadState
PUBLIC::FinalizeRestart
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitRestart 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MOD_Restart_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: nArgs
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT RESTART ...'
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.EQ.2)THEN
    doRestart=.TRUE.
    CALL GET_COMMAND_ARGUMENT(2,RestartFile)
  ELSE
    doRestart=.FALSE.
  END IF

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitRestart


!===================================================================================================================================
!> write an input solution (X1,X2,LA) to an ascii .dat file 
!!
!===================================================================================================================================
SUBROUTINE WriteStateToASCII(Uin,fileID)
! MODULES
USE MOD_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MOD_Output_Vars, ONLY:ProjectName,OutputLevel
USE MOD_MHD3D_Vars, ONLY:X1_base,X2_base,LA_base
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin 
  INTEGER               , INTENT(IN   ) :: fileID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(LEN=255)  :: fileString
  INTEGER             :: ioUnit,iMode
!===================================================================================================================================
  WRITE(FileString,'(A,"_State_",I4.4,"_",I8.8,".dat")')TRIM(ProjectName),outputLevel,fileID

  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITE SOLUTION VARIABLE TO FILE    "'//TRIM(FileString)//'" ...'

  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit       ,&
     FILE     = TRIM(FileString) ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' ) 

  WRITE(ioUnit,'(A)')'## MHD3D Solution... outputLevel and fileID:'
  WRITE(ioUnit,'(I4.4,1X,I8.8)')outputLevel,fileID
  WRITE(ioUnit,'(A)')'## grid: nElems ###########################################################################################'
  WRITE(ioUnit,'(I8)')X1_base%s%grid%nElems
  WRITE(ioUnit,'(A)')'## grid: sp(0:nElems)' 
  WRITE(ioUnit,'(*(E23.15,:,","))')X1_base%s%grid%sp(:)
  WRITE(ioUnit,'(A)')'## X1_base: s%nbase,s%deg,s%continuity,f%modes,f%sin_cos,f%excl_mn_zero ###################################'
  WRITE(ioUnit,'(*(I8,:,","))')X1_base%s%nbase,X1_base%s%deg,X1_base%s%continuity,X1_base%f%modes,X1_base%f%sin_cos &
                    ,MERGE(1,0,X1_base%f%exclude_mn_zero)
  WRITE(ioUnit,'(A)')'## X2_base: s%nbase,s%deg,s%continuity,f%modes,f%sin_cos,f%excl_mn_zero ###################################'
  WRITE(ioUnit,'(*(I8,:,","))')X2_base%s%nbase,X2_base%s%deg,X2_base%s%continuity,X2_base%f%modes,X2_base%f%sin_cos &
                    ,MERGE(1,0,X2_base%f%exclude_mn_zero)
  WRITE(ioUnit,'(A)')'## LA_base: s%nbase,s%deg,s%continuity,f%modes,f%sin_cos,f%excl_mn_zero ###################################'
  WRITE(ioUnit,'(*(I8,:,","))')LA_base%s%nbase,LA_base%s%deg,LA_base%s%continuity,LA_base%f%modes,LA_base%f%sin_cos &
                    ,MERGE(1,0,LA_base%f%exclude_mn_zero)
  WRITE(ioUnit,'(A)')'## X1: m,n,X1(1:nbase,iMode) ##############################################################################'
  DO iMode=1,X1_base%f%modes
    WRITE(ioUnit,'(2(I8,","),*(E23.15,:,","))')X1_base%f%Xmn(:,iMode),Uin%X1(:,iMode)
  END DO
  WRITE(ioUnit,'(A)')'## X2: m,n,X2(1:nbase,iMode) ##############################################################################'
  DO iMode=1,X2_base%f%modes
    WRITE(ioUnit,'(2(I8,","),*(E23.15,:,","))')X2_base%f%Xmn(:,iMode),Uin%X2(:,iMode)
  END DO
  WRITE(ioUnit,'(A)')'## LA: m,n,LA(1:nbase,iMode) ##############################################################################'
  DO iMode=1,LA_base%f%modes
    WRITE(ioUnit,'(2(I8,","),*(E23.15,:,","))')LA_base%f%Xmn(:,iMode),Uin%LA(:,iMode)
  END DO
  
  CLOSE(ioUnit)
  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE WriteStateToASCII 


!===================================================================================================================================
!> read an input solution and initialize Uin (X1,X2,LA) of size X1/X2/LA_base , from an ascii .dat file 
!! if size of grid/X1/X2/LA  not equal X1/X2/X3_base
!! interpolate readin solution to the current base of Uin
!!
!===================================================================================================================================
SUBROUTINE ReadStateFromASCII()!Uin,fileString)
! MODULES
USE MOD_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MOD_Output_Vars, ONLY:ProjectName,OutputLevel
USE MOD_MHD3D_Vars, ONLY:X1_base,X2_base,LA_base,U
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin 
!  CHARACTER(LEN=255)    , INTENT(IN   ) :: fileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(LEN=255)  :: fileString
  INTEGER              :: ioUnit,iMode,nElems_r,fileID_r,OutputLevel_r
  INTEGER              :: X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  INTEGER              :: X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  INTEGER              :: LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r
  INTEGER,ALLOCATABLE  :: X1_mn_r(:,:),X2_mn_r(:,:),LA_mn_r(:,:)
  REAL(wp),ALLOCATABLE :: sp_r(:),X1_r(:,:),X2_r(:,:),LA_r(:,:)
  LOGICAL              :: sameGrid,sameX1,sameX2,sameLA
!===================================================================================================================================
  !DEBUG!
  WRITE(FileString,'(A,"_State_",I4.4,"_",I8.8,".dat")')TRIM(ProjectName),OutputLevel,99999999

  WRITE(UNIT_stdOut,'(A)')'   READ SOLUTION VARIABLE FROM FILE    "'//TRIM(FileString)//'" ...'

  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit         ,&
     FILE     = TRIM(FileString) ,&
     STATUS   = 'OLD'            ,&
     ACTION   = 'READ'           ,&
     ACCESS   = 'SEQUENTIAL' ) 

  READ(ioUnit,*) !## MHD3D Solution file
  READ(ioUnit,*) outputLevel_r,fileID_r
  READ(ioUnit,*) !## grid: nElems
  READ(ioUnit,'(I8)') nElems_r
  ALLOCATE(sp_r(0:nElems_r))

  READ(ioUnit,*) !## grid: sp(0:nElems)
  READ(ioUnit,*)sp_r(:)
  READ(ioUnit,*) !## X1_base: 
  READ(ioUnit,*) X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  READ(ioUnit,*) !## X2_base: 
  READ(ioUnit,*) X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  READ(ioUnit,*) !## LA_base: 
  READ(ioUnit,*) LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r
  ALLOCATE(X1_r(1:X1_nbase_r,1:X1_modes_r))
  ALLOCATE(X2_r(1:X2_nbase_r,1:X2_modes_r))
  ALLOCATE(LA_r(1:LA_nbase_r,1:LA_modes_r))
  ALLOCATE(X1_mn_r(2,1:X1_modes_r))
  ALLOCATE(X2_mn_r(2,1:X2_modes_r))
  ALLOCATE(LA_mn_r(2,1:LA_modes_r))
  READ(ioUnit,*) !## X1: 
  DO iMode=1,X1_modes_r
    READ(ioUnit,*)X1_mn_r(:,iMode),X1_r(:,iMode)
  END DO
  READ(ioUnit,*) !## X2: 
  DO iMode=1,X2_modes_r
    READ(ioUnit,*)X2_mn_r(:,iMode),X2_r(:,iMode)
  END DO
  READ(ioUnit,*) !## LA: 
  DO iMode=1,LA_modes_r
    READ(ioUnit,*)LA_mn_r(:,iMode),LA_r(:,iMode)
  END DO

  CLOSE(ioUnit)

  !update outputlevel
  WRITE(UNIT_stdOut,'(A,I4.4,A)')' outputLevel of restartFile: ',outputLevel_r
  outputLevel=outputLevel_r +1

  ! check if input has changed:
  ASSOCIATE(sgrid=>X1_base%s%grid)
  sameGrid=(nElems_r.EQ.sgrid%nElems)
  IF(sameGrid) sameGrid=(SUM(ABS(sgrid%sp(:)-sp_r(:))).LT.(1.0e-12_wp*nElems_r))
  END ASSOCIATE
  IF(.NOT.sameGrid)THEN
    sameX1=.FALSE.
    sameX2=.FALSE.
    sameLA=.FALSE.
  ELSE
    sameX1=( (X1_nBase_r       .EQ.          X1_base%s%nbase     ).AND. &
             (X1_deg_r         .EQ.          X1_base%s%deg       ).AND. &
             (X1_cont_r        .EQ.          X1_base%s%continuity).AND. &
             (X1_modes_r       .EQ.          X1_base%f%modes     ).AND. &
             (X1_sin_cos_r     .EQ.          X1_base%f%sin_cos   ).AND. &
             (X1_excl_mn_zero_r.EQ.MERGE(1,0,X1_base%f%exclude_mn_zero)) ) 

    sameX2=( (X2_nBase_r       .EQ.          X2_base%s%nbase     ).AND. &
             (X2_deg_r         .EQ.          X2_base%s%deg       ).AND. &
             (X2_cont_r        .EQ.          X2_base%s%continuity).AND. &
             (X2_modes_r       .EQ.          X2_base%f%modes     ).AND. &
             (X2_sin_cos_r     .EQ.          X2_base%f%sin_cos   ).AND. &
             (X2_excl_mn_zero_r.EQ.MERGE(1,0,X2_base%f%exclude_mn_zero)) ) 

    sameLA=( (LA_nBase_r       .EQ.          LA_base%s%nbase     ).AND. &
             (LA_deg_r         .EQ.          LA_base%s%deg       ).AND. &
             (LA_cont_r        .EQ.          LA_base%s%continuity).AND. &
             (LA_modes_r       .EQ.          LA_base%f%modes     ).AND. &
             (LA_sin_cos_r     .EQ.          LA_base%f%sin_cos   ).AND. &
             (LA_excl_mn_zero_r.EQ.MERGE(1,0,LA_base%f%exclude_mn_zero)) ) 
  END IF
  IF(sameX1.AND.sameX2.AND.sameLA)THEN
    U(-1)%X1(:,:)=X1_r
    U(-1)%X2(:,:)=X2_r
    U(-1)%LA(:,:)=LA_r
  ELSE
    WRITE(*,*)sameGrid,sameX1,sameX2,sameLA
    STOP 'restart from other configuration not yet implemented'
  END IF 
  

  DEALLOCATE(sp_r)
  DEALLOCATE(X1_r)
  DEALLOCATE(X2_r)
  DEALLOCATE(LA_r)
  DEALLOCATE(X1_mn_r)
  DEALLOCATE(X2_mn_r)
  DEALLOCATE(LA_mn_r)

  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE ReadStateFromASCII 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeRestart 
! MODULES
USE MOD_Restart_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
