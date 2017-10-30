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
USE MOD_Globals,ONLY:UNIT_stdOut,fmt_sep
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
WRITE(UNIT_stdOut,fmt_sep)
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
  CALL VMEC1D_visu() 
END IF !visuVMEC1D

END SUBROUTINE Analyze 


!===================================================================================================================================
!> Visualize VMEC flux surface data for each mode, for Rmnc 
!!
!===================================================================================================================================
SUBROUTINE VMEC1D_visu()
! MODULES
USE MOD_Globals, ONLY:wp
USE MOD_Output_Vars, ONLY:ProjectName
USE MOD_Output_CSV, ONLY:WriteDataToCSV
USE MOD_VMEC_Readin
USE MOD_VMEC_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: nVal,iMode
CHARACTER(LEN=120) :: varnames(7+mn_mode) 
REAL(wp)           :: values(7+mn_mode,nFluxVMEC-1) 
!===================================================================================================================================
nVal=1
Varnames(nVal)='Phi'
values(  nVal,:)=Phi_prof(2:nFluxVMEC)

nVal=nVal+1
Varnames(nVal)='chi'
values(  nVal,:)=Chi_prof(2:nFluxVMEC)

nVal=nVal+1
Varnames(nVal)='rho'
values(  nVal,:)=rho(2:nFluxVMEC)

nVal=nVal+1
Varnames(nVal)='iota(Phi_norm)'
values(  nVal,:)=iotaf(2:nFluxVMEC)

nVal=nVal+1
Varnames(nVal)='pres(Phi_norm)'
values(  nVal,:)=presf(2:nFluxVMEC)

nVal=nVal+2
Varnames(nVal-1)='Rmnc_m_odd_n0'
Varnames(nVal)='Rmnc_m_even_n0'
values(nVal-1:nVal,:)=0.
DO iMode=1,mn_mode
  IF(NINT(xn(iMode)).EQ.0)THEN
    IF(MOD(NINT(xm(iMode)),2).NE.0)THEN
      values(nVal-1,:)= values(nVal-1,:)+Rmnc(iMode,2:nFluxVMEC)
    ELSE
      values(nVal,:)= values(nVal,:)+Rmnc(iMode,2:nFluxVMEC)
    END IF
  END IF !n=0
END DO

DO iMode=1,mn_mode
  nVal=nVal+1
  WRITE(VarNames(nVal),'(A,"_m_",I3.3,"_n_",I3.3)')'Rmnc',NINT(xm(iMode)),NINT(xn(iMode))
  values(nVal,:)=Rmnc(iMode,2:nFluxVMEC)
END DO
CALL WriteDataToCSV(nVal, nFluxVMEC-1,VarNames(1:nVal),Values(1:nVal,:),(TRIM(ProjectName)//"_Rmnc_modes"))

nval=7 !rewind
Varnames(nVal-1)='Zmns_m_odd_n0'
Varnames(nVal)='Zmns_m_even_n0'
values(nVal-1:nVal,:)=0.
DO iMode=1,mn_mode
  IF(NINT(xn(iMode)).EQ.0)THEN
    IF(MOD(NINT(xm(iMode)),2).NE.0)THEN
      values(nVal-1,:)= values(nVal-1,:)+Zmns(iMode,2:nFluxVMEC)
    ELSE
      values(nVal,:)= values(nVal,:)+Zmns(iMode,2:nFluxVMEC)
    END IF
  END IF !n=0
END DO
DO iMode=1,mn_mode
  nVal=nVal+1
  WRITE(VarNames(nVal),'(A,"_",I3.3,"_",I3.3)')'Zmns',NINT(xm(iMode)),NINT(xn(iMode))
  values(nVal,:)=Zmns(iMode,2:nFluxVMEC)
END DO
CALL WriteDataToCSV(nVal, nFluxVMEC-1,VarNames(1:nVal),Values(1:nVal,:),(TRIM(ProjectName)//"_Zmns_modes"))

END SUBROUTINE VMEC1D_visu 

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
