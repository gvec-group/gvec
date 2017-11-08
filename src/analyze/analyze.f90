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
SWRITE(UNIT_stdOut,'(A)')'INIT ANALYZE ...'
visuVMEC1D    = GETLOGICAL('visuVMEC1D','T')   
visuVMEC2D    = GETLOGICAL('visuVMEC2D','F')   

SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitAnalyze


!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE Analyze()
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
USE MOD_Output_Vars, ONLY:ProjectName
USE MOD_Output_CSV, ONLY:WriteDataToCSV
USE MOD_VMEC_Readin
USE MOD_VMEC_Vars
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iGuess,nVal,nValRewind,iMode
INTEGER,PARAMETER  :: n_Int=200
CHARACTER(LEN=120) :: varnames(  8+mn_mode) 
REAL(wp)           :: values(    8+mn_mode,nFluxVMEC) 
REAL(wp)           :: values_int(8+mn_mode,n_int) 
REAL(wp)           :: rho_int(n_int),rho_half(nFluxVMEC)
!===================================================================================================================================
!interpolation points
DO i=0,n_int-1
  rho_int(1+i)=REAL(i,wp)/REAL(n_int-1,wp)
END DO

nVal=1
Varnames(nVal)='Phi'
values(  nVal,:)=Phi_prof(:)
values_int(nVal,:)=EvalSpl(Phi_Spl)

nVal=nVal+1
Varnames(nVal)='chi'
values(  nVal,:)=Chi_prof(:)
values_int(nVal,:)=EvalSpl(chi_Spl)

nVal=nVal+1
Varnames(nVal)='rho'
values(  nVal,:)=rho(:)
values_int(nVal,:)=rho_int(:)

rho_half(1)=0.
DO i=1,nFluxVMEC-1
  rho_half(i+1)=SQRT(0.5_wp*(Phinorm_prof(i+1)+Phinorm_prof(i))) !0.5*(rho(iFlux)+rho(iFlux+1))
END DO
nVal=nVal+1
Varnames(nVal)='rho_half'
values(  nVal,:)= rho_half(:)
values_int(nVal,:)=0.

nVal=nVal+1
Varnames(nVal)='iota(Phi_norm)'
values(  nVal,:)=iotaf(:)
values_int(nVal,:)=EvalSplDeriv(chi_Spl) / ( EvalSplDeriv(Phi_Spl) )

nVal=nVal+1
Varnames(nVal)='pres(Phi_norm)'
values(  nVal,:)=presf(:)
values_int(nVal,:)=EvalSpl(pres_Spl)

nValRewind=nVal

CALL writeDataMN("Rmnc","Rmnc",Rmnc)
nval=nValRewind
CALL writeDataMN("Zmns","Zmns",Zmns)
nval=nValRewind
IF(reLambda)THEN
  CALL writeDataMN("Lmns","Lmns",Lmns)
ELSE
  CALL writeDataMN("Lmns_half","Lmns",Lmns)
END IF

!interpolated profiles
nval=nValRewind
CALL writeDataMN_int("INT_Rmnc","Rmnc",Rmnc_Spl)
nval=nValRewind
CALL writeDataMN_int("INT_Zmns","Zmns",Zmns_Spl)
nval=nValRewind
IF(reLambda)THEN
  CALL writeDataMN_int("INT_Lmns","Lmns",Lmns_Spl)
ELSE
  CALL writeDataMN_int("INT_Lmns_half","Lmns",Lmns_spl)
END IF

CONTAINS
  SUBROUTINE writeDataMN(fname,vname,xx)
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    REAL(wp),INTENT(IN)        :: xx(:,:)
    
    DO iMode=1,mn_mode
      nVal=nVal+1
      WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname),NINT(xm(iMode)),NINT(xn(iMode))/nfp
      values(nVal,:)=xx(iMode,:)
    END DO
    nVal=nVal+2
    Varnames(nVal-1)=TRIM(vname)//', m= odd, n= 000'
    Varnames(nVal)=  TRIM(vname)//', m=even, n= 000'
    values(nVal-1:nVal,:)=0.
    DO iMode=1,mn_mode
      IF(NINT(xn(iMode)).EQ.0)THEN
        IF(MOD(NINT(xm(iMode)),2).NE.0)THEN
          values(nVal-1,:)= values(nVal-1,:)+values(nVal-2-mn_mode+iMode,:)
        ELSE
          values(nVal,:)= values(nVal,:)+values(nVal-2-mn_mode+iMode,:)
        END IF
      END IF !n=0
    END DO
    CALL WriteDataToCSV(nVal, nFluxVMEC,VarNames(1:nVal),Values(1:nVal,:),(TRIM(ProjectName)//"_"//TRIM(fname)//"_modes"))

  END SUBROUTINE writeDataMN

  SUBROUTINE writeDataMN_int(fname,vname,xx_Spl)
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    REAL(wp),INTENT(IN)        :: xx_Spl(:,:,:)

    DO iMode=1,mn_mode
      nVal=nVal+1
      WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname),NINT(xm(iMode)),NINT(xn(iMode))/nfp
    END DO
    values_int(nVal-mn_mode+1:nVal,:)=EvalSplMode(xx_Spl)

    nVal=nVal+2
    Varnames(nVal-1)=TRIM(vname)//', m= odd, n= 000'
    Varnames(nVal)=  TRIM(vname)//', m=even, n= 000'
    values_int(nVal-1:nVal,:)=0.
    DO iMode=1,mn_mode
      IF(NINT(xn(iMode)).EQ.0)THEN
        IF(MOD(NINT(xm(iMode)),2).NE.0)THEN
          values_int(nVal-1,:)= values_int(nVal-1,:)+values_int(nVal-2-mn_mode+iMode,:)
        ELSE
          values_int(nVal,:)= values_int(nVal,:)+values_int(nVal-2-mn_mode+iMode,:)
        END IF
      END IF !n=0
    END DO
    CALL WriteDataToCSV(nVal, n_Int,VarNames(1:nVal),Values_int(1:nVal,:),(TRIM(ProjectName)//"_"//TRIM(fname)//"_modes"))
  END SUBROUTINE writeDataMN_int

  FUNCTION EvalSpl(xx_spl)
    REAL(wp),INTENT(IN)        :: xx_spl(:,:)
    REAL(wp)                   :: EvalSpl(n_int)
    !local
    REAL(wp)                   :: splOut(3)
    DO i=1,n_int
      CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_int(i),rho,xx_Spl(:,:),iGuess,splout) 
      EvalSpl(i)=splout(1)
    END DO
  END FUNCTION EvalSpl

  FUNCTION EvalSplDeriv(xx_spl)
    REAL(wp),INTENT(IN)        :: xx_Spl(:,:)
    REAL(wp)                   :: EvalSplDeriv(n_int)
    !local
    REAL(wp)                   :: splOut(3)
    DO i=1,n_int
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_int(i),rho,xx_Spl(:,:),iGuess,splout) 
      EvalSplDeriv(i)=splout(2)
    END DO
  END FUNCTION EvalSplDeriv

  FUNCTION EvalSplMode(xx_Spl)
    REAL(wp),INTENT(IN)        :: xx_Spl(:,:,:)
    REAL(wp)                   :: EvalSplMode(mn_mode,n_int)
    !local
    REAL(wp)                   :: rho_p,rhom,drhom,splOut(3) !for weighted spline interpolation
    DO i=1,n_int
      rho_p=rho_int(i)
      DO iMode=1,mn_mode
        SELECT CASE(xmabs(iMode))
        CASE(0)
          rhom=1.0_wp
          drhom=0.0_wp
        CASE(1)
          rhom=rho_p
          drhom=1.0_wp
        CASE(2)
          rhom=rho_p*rho_p
          drhom=2.0_wp*rho_p
        CASE DEFAULT
          rhom=rho_p**xmabs(iMode)
          drhom=REAL(xmabs(iMode),wp)*rho_p**(xmabs(iMode)-1)
        END SELECT
        CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,xx_Spl(:,:,iMode),iGuess,splout) 
        EvalSplMode(iMode,i)=rhom*splout(1)
      END DO !iMode
    END DO !rho_int(i)
  END FUNCTION EvalSplMode

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
