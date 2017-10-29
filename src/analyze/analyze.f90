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
!!# Module **Analyze**
!!
!! Analyze and output equilibrium data 
!!
!===================================================================================================================================
MODULE MOD_Analyze
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE Analyze
  MODULE PROCEDURE Analyze
END INTERFACE

INTERFACE FinalizeAnalyze
  MODULE PROCEDURE FinalizeAnalyze
END INTERFACE

PUBLIC::InitAnalyze
PUBLIC::Analyze
PUBLIC::FinalizeAnalyze
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitAnalyze 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut
USE MOD_Analyze_Vars
USE MOD_ReadInTools,ONLY:GETLOGICAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'INIT ANALYZE ...'
visuVMEC1D    = GETLOGICAL('visuVMEC1D','T')   
visuVMEC2D    = GETLOGICAL('visuVMEC2D','F')   

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitAnalyze


!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE Analyze()
! MODULES
USE MOD_Globals, ONLY:wp
USE MOD_Analyze_Vars
USE MOD_Output_CSV, ONLY:WriteDataToCSV
USE MOD_VMEC_Readin, ONLY:mn_mode,xm,xn,nFluxVMEC
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iMode
CHARACTER(LEN=120) :: modenames(2+mn_mode) 
REAL(wp)           :: values(2+mn_mode,nFluxVMEC) 
!===================================================================================================================================
IF(visuVMEC1D)THEN
  modenames(1)="Phi"
  modenames(2)="chi"
  DO iMode=1,mn_mode
    WRITE(modenames(2+iMode),'(A,"_",I3.3,"_",I3.3)')'Rmnc',NINT(xm(iMode)),NINT(xn(iMode))
  END DO
  values=0.0_wp
  CALL WriteDataToCSV(2+mn_mode, nFluxVMEC,modenames,Values,"Rmnc_modes")
END IF !visuVMEC1D

END SUBROUTINE Analyze 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeAnalyze 
! MODULES
USE MOD_Analyze_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE FinalizeAnalyze

END MODULE MOD_Analyze
