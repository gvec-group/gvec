!================================================================================================================================!
! Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!================================================================================================================================!
! PyGVEC postprocessing
!
! This module contains the Python interface for postprocessing the GVEC output, starting with a parameter file and a statefile.
!================================================================================================================================!

MODULE MODpygvec_post

USE MODgvec_c_functional, ONLY: t_functional

IMPLICIT NONE

PUBLIC

CLASS(t_functional), ALLOCATABLE :: functional

CONTAINS

!================================================================================================================================!
! a simple scalar function to test the wrapper
SUBROUTINE Init(parameterfile)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Analyze,        ONLY: InitAnalyze
  USE MODgvec_Output,         ONLY: InitOutput
  ! USE MODgvec_Output_vars,    ONLY: OutputLevel
  USE MODgvec_Restart,        ONLY: InitRestart
  ! USE MODgvec_Restart,        ONLY: RestartFromState
  ! USE MODgvec_Restart_vars,   ONLY: doRestart
  ! USE MODgvec_Output_Vars,    ONLY: OutputLevel,ProjectName
  ! USE MODgvec_ReadState_Vars, ONLY: fileID_r,outputLevel_r
  ! USE MODgvec_MHD3D_Vars,     ONLY: U,F
  ! USE MODgvec_MHD3D_EvalFunc, ONLY: EvalForce
  USE MODgvec_ReadInTools,    ONLY: FillStrings,GETLOGICAL,GETINT,IgnoredStrings 
  USE MODgvec_Functional,     ONLY: InitFunctional
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=*) :: parameterfile
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: which_functional
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A)') "GVEC POST ! GVEC POST ! GVEC POST ! GVEC POST"
  WRITE(Unit_stdOut,'(132("="))')

  ! read parameter file
  CALL FillStrings(parameterfile)
  
  ! initialization phase
  CALL InitRestart()
  CALL InitOutput()
  CALL InitAnalyze()
  
  ! initialize the functional
  which_functional = GETINT('which_functional',Proposal=1)
  CALL InitFunctional(functional,which_functional)
  
  ! print the ignored parameters
  CALL IgnoredStrings()
END SUBROUTINE init

!================================================================================================================================!
SUBROUTINE ReadState(statefile)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Output_Vars,    ONLY: outputLevel,ProjectName
  USE MODgvec_ReadState_Vars, ONLY: fileID_r,outputLevel_r
  USE MODgvec_Restart,        ONLY: RestartFromState
  USE MODgvec_MHD3D_Vars,     ONLY: U
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=255) :: statefile
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  ProjectName='POST_'//TRIM(StateFile(1:INDEX(StateFile,'_State_')-1))
  CALL RestartFromState(StateFile,U(0))
  ! JacCheck=2
  ! CALL EvalForce(U(0),.TRUE.,JacCheck, F(0))
  ! CALL Analyze(FileID_r)
  outputLevel=outputLevel_r
END SUBROUTINE ReadState

!================================================================================================================================!
SUBROUTINE Finalize()
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Analyze,        ONLY: FinalizeAnalyze
  USE MODgvec_Output,         ONLY: FinalizeOutput
  ! USE MODgvec_Output_vars,    ONLY: OutputLevel
  USE MODgvec_Restart,        ONLY: FinalizeRestart
  ! USE MODgvec_Restart,        ONLY: RestartFromState
  ! USE MODgvec_Restart_vars,   ONLY: doRestart
  ! USE MODgvec_Output_Vars,    ONLY: OutputLevel,ProjectName
  ! USE MODgvec_ReadState_Vars, ONLY: fileID_r,outputLevel_r
  ! USE MODgvec_MHD3D_Vars,     ONLY: U,F
  ! USE MODgvec_MHD3D_EvalFunc, ONLY: EvalForce
  ! USE MODgvec_ReadInTools,    ONLY: FillStrings,GETLOGICAL,GETINT,IgnoredStrings 
  USE MODgvec_Functional,     ONLY: FinalizeFunctional
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL FinalizeFunctional(functional)
  DEALLOCATE(functional)
  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  CALL FinalizeRestart()
  
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A)') "GVEC POST FINISHED !"
  WRITE(Unit_stdOut,'(132("="))')
END SUBROUTINE Finalize

END MODULE MODpygvec_post