!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **Restart**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_Restart
! MODULES
USE MODgvec_Globals, ONLY:wp,MPIroot,enter_subregion,exit_subregion
IMPLICIT NONE
PRIVATE

INTERFACE InitRestart
  MODULE PROCEDURE InitRestart
END INTERFACE

INTERFACE WriteState
  MODULE PROCEDURE WriteStateToASCII
END INTERFACE

INTERFACE RestartFromState
  MODULE PROCEDURE RestartFromState
END INTERFACE

INTERFACE FinalizeRestart
  MODULE PROCEDURE FinalizeRestart
END INTERFACE

PUBLIC::InitRestart
PUBLIC::WriteState
PUBLIC::RestartFromState
PUBLIC::FinalizeRestart
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE InitRestart(RestartFile_in)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_Restart_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL    :: Restartfile_in
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.MPIroot) RETURN
  WRITE(UNIT_stdOut,'(A)')'INIT RESTART ...'
  doRestart=.TRUE.
  IF(PRESENT(restartFile_in))THEN
    restartfile=restartfile_in
  ELSE
    restartfile=""  !nothing to save!
  END IF
  WRITE(UNIT_stdOut,'(A)')'... DONE'
  WRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitRestart


!===================================================================================================================================
!> write an input solution (X1,X2,LA) to an ascii .dat file
!!
!===================================================================================================================================
SUBROUTINE WriteStateToASCII(Uin,fileID)
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,PI,TWOPI,GETFREEUNIT
USE MODgvec_Output_Vars, ONLY:ProjectName,OutputLevel
USE MODgvec_MHD3D_Vars, ONLY:X1_base,X2_base,LA_base,sgrid,which_hmap
USE MODgvec_MHD3D_vars, ONLY: Phi_profile, chi_profile, pres_profile, iota_profile
USE MODgvec_MPI, ONLY:par_Reduce
USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
USE MODgvec_MHD3D_evalFunc, ONLY: EvalTotals
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
  INTEGER             :: ioUnit,iMode,is
  REAL(wp)            :: vol,surfAvg,a_minor,r_major
!===================================================================================================================================
  __PERFON("output_state")
  SWRITE(FileString,'(A,"_State_",I4.4,"_",I8.8,".dat")')TRIM(ProjectName),outputLevel,fileID
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITE SOLUTION VARIABLE TO FILE    "'//TRIM(FileString)//'" ...'
  !compute volume& poloidal surface average -> pi*aMinor^2=surfavg, surfavg*2*Pi*RMajor=volume
  CALL EvalTotals(Uin,vol,surfAvg)
  IF(MPIroot)THEN
    a_Minor = SQRT(surfAvg/PI)
    r_Major = vol/(TWOPI*surfAvg)

    __PERFON("write_state")
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' ...'
    ioUnit=GETFREEUNIT()
    OPEN(UNIT     = ioUnit       ,&
       FILE     = TRIM(FileString) ,&
       STATUS   = 'REPLACE'   ,&
       ACCESS   = 'SEQUENTIAL' )

    WRITE(ioUnit,'(A)')'## MHD3D Solution... outputLevel and fileID:'
    WRITE(ioUnit,'(I4.4,",",I8.8)')outputLevel,fileID
    WRITE(ioUnit,'(A)')'## grid: nElems, gridType #################################################################################'
    WRITE(ioUnit,'(*(I8,:,","))')sgrid%nElems,sgrid%grid_type
    WRITE(ioUnit,'(A)')'## grid: sp(0:nElems)'
    WRITE(ioUnit,'(*(E23.15,:,","))')X1_base%s%grid%sp(:)
    WRITE(ioUnit,'(A)')'## global: nfp,degGP,mn_nyq(2),hmap #######################################################################'
    WRITE(ioUnit,'(*(I8,:,","))')X1_base%f%nfp,X1_base%s%degGP,X1_base%f%mn_nyq,which_hmap
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
    WRITE(ioUnit,'(A)')'## at X1_base IP point positions (size nBase): spos,phi,chi,iota,pressure  ################################'
    ASSOCIATE(s_IP         => X1_base%s%s_IP, &
              nBase        => X1_base%s%nBase )
    !write Profiles at Greville interpolation points s_IP(1:nBase)
    DO is=1,nBase
      WRITE(ioUnit,'(*(E23.15,:,","))')s_IP(is),Phi_profile%eval_at_rho(s_IP(is)), &
                                                chi_profile%eval_at_rho(s_IP(is)), &
                                               iota_profile%eval_at_rho(s_IP(is)),&
                                               pres_profile%eval_at_rho(s_IP(is))
    END DO
    END ASSOCIATE
    WRITE(ioUnit,'(A)')'## a_minor,r_major,volume  ################################################################################'
    WRITE(ioUnit,'(*(E23.15,:,","))')a_Minor,r_Major,vol

    CLOSE(ioUnit)
    WRITE(UNIT_stdOut,'(A)')'...DONE.'
    __PERFOFF("write_state")
  END IF !MPIroot
  __PERFOFF("output_state")
END SUBROUTINE WriteStateToASCII


!===================================================================================================================================
!> read an input solution and initialize U(0) (X1,X2,LA) of size X1/X2/LA_base , from an ascii .dat file
!! if size of grid/X1/X2/LA  not equal X1/X2/X3_base
!! interpolate readin solution to the current base of Uin
!!
!===================================================================================================================================
SUBROUTINE RestartFromState(fileString,U_r)
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MODgvec_Output_Vars, ONLY:OutputLevel
USE MODgvec_MHD3D_Vars,  ONLY:X1_base,X2_base,LA_base,sgrid,hmap
USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base, base_new
USE MODgvec_readState_Vars, ONLY:sgrid_r,X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r,outputLevel_r
USE MODgvec_readState, ONLY: ReadState,Finalize_ReadState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=255)    , INTENT(IN   ) :: fileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: U_r
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL              :: sameGrid
  LOGICAL              :: sameX1  ,sameX2  ,sameLA, changed
!===================================================================================================================================
  IF(.NOT.MPIroot) RETURN
  WRITE(UNIT_stdOut,'(A)')'RESTARTING FROM FILE ...'
  CALL enter_subregion("restart-from-state")
  CALL ReadState(FileString,hmap_in=hmap)

  !update outputlevel
  WRITE(UNIT_stdOut,'(A,I4.4,A)')' outputLevel of restartFile: ',outputLevel_r
  outputLevel=outputLevel_r +1
  CALL sgrid%compare(sgrid_r,sameGrid)
  CALL X1_base%compare(X1_base_r,sameX1)
  CALL X2_base%compare(X2_base_r,sameX2)
  CALL LA_base%compare(LA_base_r,sameLA)

  changed=.NOT.(sameX1.AND.sameX2.AND.sameLA)

  IF(changed)THEN
    WRITE(UNIT_stdOut,'(A,4(A,L1))') '    ... restart from other configuration: \n','         sameGrid= ',sameGrid, ', sameX1= ',sameX1,', sameX2= ',sameX2,', sameLA= ',sameLA
  ELSE
    WRITE(UNIT_stdOut,'(A)') '     ... restart from same configuration ... '
  END IF
  CALL X1_base%change_base(X1_base_r,X1_r,U_r%X1)
  CALL X2_base%change_base(X2_base_r,X2_r,U_r%X2)
  CALL LA_base%change_base(LA_base_r,LA_r,U_r%LA)

  CALL Finalize_ReadState()
  CALL exit_subregion("restart-from-state")
  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE RestartFromState

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeRestart
! MODULES
USE MODgvec_Restart_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.MPIroot) RETURN
  dorestart=.FALSE.
END SUBROUTINE FinalizeRestart

END MODULE MODgvec_Restart
