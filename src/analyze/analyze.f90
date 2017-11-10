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
USE MOD_ReadInTools,ONLY:GETINT
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
visuVMEC1D    = GETINT('visuVMEC1D','1')   
visuVMEC2D    = GETINT('visuVMEC2D','0')   

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
IF(visuVMEC1D.NE.0)THEN
  CALL VMEC1D_visu() 
END IF !visuVMEC1D

END SUBROUTINE Analyze 


!===================================================================================================================================
!> Visualize VMEC flux surface data for each mode, for Rmnc 
!!
!===================================================================================================================================
SUBROUTINE VMEC1D_visu()
! MODULES
USE MOD_Globals,ONLY:Pi
USE MOD_Analyze_Vars, ONLY:visuVMEC1D
USE MOD_Output_Vars, ONLY:ProjectName
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
CHARACTER(LEN=120) :: varnames(  8+2*mn_mode) 
REAL(wp)           :: values(    8+2*mn_mode,nFluxVMEC) 
REAL(wp)           :: values_int(8+2*mn_mode,n_int) 
REAL(wp)           :: rho_int(n_int),rho_half(nFluxVMEC)
LOGICAL            :: vcase(4)
CHARACTER(LEN=4)   :: vstr
!===================================================================================================================================
!visuVmec1D: all possible combinations: 1,2,3,4,12,13,14,23,24,34,123,124,234,1234
WRITE(vstr,'(I4)')visuVmec1D
vcase=.FALSE.
IF(INDEX(vstr,'1').NE.0) vcase(1)=.TRUE.
IF(INDEX(vstr,'2').NE.0) vcase(2)=.TRUE.
IF(INDEX(vstr,'3').NE.0) vcase(3)=.TRUE.
IF(INDEX(vstr,'4').NE.0) vcase(4)=.TRUE.
IF(.NOT.(ANY(vcase))) THEN
  WRITE(*,*)'visuVmec1D case not found:',visuVmec1D,' nothing visualized...'
  RETURN
END IF

!interpolation points
DO i=0,n_int-1
  rho_int(1+i)=REAL(i,wp)/REAL(n_int-1,wp)
END DO
!strech towards axis and edge
rho_int=rho_int+0.05*SIN(Pi*(2*rho_int-1))

nVal=1
Varnames(nVal)='Phi'
values(  nVal,:)=Phi_prof(:)
values_int(nVal,:)=EvalSpl(0,Phi_Spl)

nVal=nVal+1
Varnames(nVal)='chi'
values(  nVal,:)=Chi_prof(:)
values_int(nVal,:)=EvalSpl(0,chi_Spl)

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
values_int(nVal,:)=EvalSpl(0,iota_Spl) 

nVal=nVal+1
Varnames(nVal)='pres(Phi_norm)'
values(  nVal,:)=presf(:)
values_int(nVal,:)=EvalSpl(0,pres_Spl)

nValRewind=nVal


IF(vcase(1))THEN
  WRITE(*,*)'1) Visualize R,Z,lambda pointwise ...'
  nval=nValRewind
  CALL writeDataMN("Rmnc","Rmnc",0,rho,Rmnc)
  nval=nValRewind
  CALL writeDataMN("Zmns","Zmns",0,rho,Zmns)
  nval=nValRewind
  IF(reLambda)THEN
    CALL writeDataMN("Lmns","Lmns",0,rho,Lmns)
  ELSE
    CALL writeDataMN("Lmns_half","Lmns_h",0,rho_half,Lmns)
  END IF
  IF(lasym)THEN
    nval=nValRewind
    CALL writeDataMN("Rmns","Rmns",0,rho,Rmnc)
    nval=nValRewind
    CALL writeDataMN("Zmnc","Zmnc",0,rho,Zmns)
    nval=nValRewind
    IF(reLambda)THEN
      CALL writeDataMN("Lmnc","Lmnc",0,rho,Lmns)
    ELSE
      CALL writeDataMN("Lmnc_half","Lmnc_h",0,rho_half,Lmns)
    END IF
  END IF!lasym
END IF
IF(vcase(2))THEN
  WRITE(*,*)'2) Visualize R,Z,lambda interpolated...'
  nval=nValRewind
  CALL writeDataMN_int("INT_Rmnc","Rmnc",0,rho_int,Rmnc_Spl)
  nval=nValRewind
  CALL writeDataMN_int("INT_Zmns","Zmns",0,rho_int,Zmns_Spl)
  nval=nValRewind
  IF(reLambda)THEN
    CALL writeDataMN_int("INT_Lmns","Lmns",0,rho_int,Lmns_Spl)
  ELSE
    CALL writeDataMN_int("INT_Lmns_half","Lmns_h",0,rho_int,Lmns_spl)
  END IF
  IF(lasym)THEN
    nval=nValRewind
    CALL writeDataMN_int("INT_Rmns","Rmns",0,rho_int,Rmnc_Spl)
    nval=nValRewind
    CALL writeDataMN_int("INT_Zmnc","Zmnc",0,rho_int,Zmns_Spl)
    IF(reLambda)THEN
      CALL writeDataMN_int("INT_Lmnc","Lmnc",0,rho_int,Lmns_Spl)
    ELSE
      CALL writeDataMN_int("INT_Lmnc_half","Lmnc_h",0,rho_int,Lmns_spl)
    END IF
  END IF!lasym
END IF
IF(vcase(3))THEN
  WRITE(*,*)'3) Visualize dRrho,dZrho pointwise (1st order finite difference)...'
  nval=nValRewind
  CALL writeDataMN("dRmnc","dRmnc",1,rho,Rmnc)
  nval=nValRewind
  CALL writeDataMN("dZmns","dZmns",1,rho,Zmns)
  IF(lasym)THEN
    nval=nValRewind
    CALL writeDataMN("dRmns","dRmns",1,rho,Rmnc)
    nval=nValRewind
    CALL writeDataMN("dZmnc","dZmnc",1,rho,Zmns)
  END IF!lasym
END IF
IF(vcase(4))THEN
  WRITE(*,*)'4) Visualize dRrho,dZrho interpolated...'
  nval=nValRewind
  CALL writeDataMN_int("INT_dRmnc","dRmnc",1,rho_int,Rmnc_Spl)
  nval=nValRewind
  CALL writeDataMN_int("INT_dZmns","dZmns",1,rho_int,Zmns_Spl)
  IF(lasym)THEN
    !interpolated profiles
    nval=nValRewind
    CALL writeDataMN_int("INT_dRmns","dRmns",1,rho_int,Rmnc_Spl)
    nval=nValRewind
    CALL writeDataMN_int("INT_dZmnc","dZmnc",1,rho_int,Zmns_Spl)
    nval=nValRewind
  END IF!lasym
END IF

CONTAINS
  SUBROUTINE writeDataMN(fname,vname,rderiv,coord,xx)
    INTEGER,INTENT(IN)         :: rderiv !0: point values, 1: 1/2 ( (R(i+1)-R(i))/rho(i+1)-rho(i) (R(i)-R(i-1))/rho(i)-rho(i-1)) 
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    REAL(wp),INTENT(IN)        :: xx(:,:)
    REAL(wp),INTENT(IN)        :: coord(nFluxVMEC)
    !local
    REAL(wp)                   :: dxx(size(xx,1),size(xx,2)) !derivative in Rho
    IF(rderiv.EQ.1) THEN
      dxx(:,1)=(xx(:,2)-xx(:,1))/(rho(2)-rho(1))
      DO i=2,nFluxVMEC-1
        dxx(:,i)=0.5_wp*( (xx(:,i+1)-xx(:,i  ))/(rho(i+1)-rho(i  )) &
                         +(xx(:,i  )-xx(:,i-1))/(rho(i  )-rho(i-1))) !mean slope
      END DO
      dxx(:,nFluxVMEC)=(xx(:,nFluxVMEC)-xx(:,nFluxVMEC-1))/(rho(nFluxVMEC)-rho(nFluxVMEC-1))
    ELSE
      dxx=xx
    END IF
    DO iMode=1,mn_mode
      nVal=nVal+1
      WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname),NINT(xm(iMode)),NINT(xn(iMode))/nfp
      values(nVal,:)=dxx(iMode,:)
    END DO
    CALL writeNow(fname,vname,coord,values,VarNames) 

  END SUBROUTINE writeDataMN

  SUBROUTINE writeDataMN_int(fname,vname,rderiv,coord,xx_Spl)
    INTEGER,INTENT(IN)         :: rderiv !0: eval spl, 1: eval spl deriv
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    REAL(wp),INTENT(IN)        :: xx_Spl(:,:,:)
    REAL(wp),INTENT(IN)        :: coord(n_int)

    DO iMode=1,mn_mode
      nVal=nVal+1
      WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname),NINT(xm(iMode)),NINT(xn(iMode))/nfp
    END DO
    values_int(nVal-mn_mode+1:nVal,:)=EvalSplMode(rderiv,xx_Spl)

    CALL writeNow(fname,vname,coord,values_int,VarNames) 

  END SUBROUTINE writeDataMN_int

  SUBROUTINE writeNow(fname,vname,coord,values_in,VarNames_in)
    USE MOD_Output_CSV, ONLY:WriteDataToCSV
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    REAL(wp),INTENT(IN)        :: coord(:)
    REAL(wp),INTENT(INOUT)        :: values_in(:,:)
    CHARACTER(LEN=*),INTENT(INOUT):: VarNames_in(:)
    !local
    REAL(wp)                   :: minmaxval(2)
    REAL(wp) ,ALLOCATABLE      :: max_loc_val(:)
    CHARACTER(LEN=100),ALLOCATABLE :: varnames_max(:)

    minmaxval(1)=MINVAL(values_in(nVal-mn_mode:nVal,:))
    minmaxval(2)=MAXVAL(values_in(nVal-mn_mode:nVal,:))

    DO iMode=1,mn_mode
      nVal=nVal+1
      WRITE(VarNames_in(nVal),'(A)')TRIM(VarNames_in(nVal-mn_mode))//'_norm'
      values_in(nVal,:)=values_in(nVal-mn_mode,:)/(MAXVAL(ABS(values_in(nVal-mn_mode,:)))+1.0E-12)
    END DO

    nVal=nVal+2
    Varnames_in(nVal-1)=TRIM(vname)//', m= odd, n= 000'
    Varnames_in(nVal)=  TRIM(vname)//', m=even, n= 000'
    values_in(nVal-1:nVal,:)=0.
    DO iMode=1,mn_mode
      IF(NINT(xn(iMode)).EQ.0)THEN
        IF(MOD(NINT(xm(iMode)),2).NE.0)THEN
          values_in(nVal-1,:)= values_in(nVal-1,:)+values_in(nVal-2-2*mn_mode+iMode,:)
        ELSE
          values_in(nVal,:)= values_in(nVal,:)+values_in(nVal-2-2*mn_mode+iMode,:)
        END IF
      END IF !n=0
    END DO

    CALL WriteDataToCSV(VarNames_in(1:nVal),Values_in(1:nVal,:), (TRIM(ProjectName)//"_"//TRIM(fname)//"_modes"))

    ALLOCATE(max_loc_val(nVal),Varnames_max(nVal))
    DO i=1,nVal
      max_loc_val(i)=coord(MAXLOC(ABS(values_in(i,:)),1))
      Varnames_max(i)=TRIM(VarNames_in(i))//'_maxloc'
    END DO 
    CALL WriteDataToCSV(VarNames_max(:) ,RESHAPE(max_loc_val(:),(/nval,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"_modes") &
                               ,append_in=.TRUE.,vfmt_in='E10.2')
    DO i=1,nVal
      max_loc_val(i)=      MAXVAL(ABS(values_in(i,:)))+1.0E-12
      Varnames_max(i)=TRIM(VarNames_in(i))//'_maxval'
    END DO 
    CALL WriteDataToCSV(VarNames_max(:) ,RESHAPE(max_loc_val(:),(/nval,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"_modes") &
                               ,append_in=.TRUE.,vfmt_in='E10.2')
    DEALLOCATE(max_loc_val,Varnames_max)
    !write position of first flux surface
    CALL WriteDataToCSV((/'rhoFirst'/) ,RESHAPE((/rho(2)/),(/1,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"_modes") &
                               ,append_in=.TRUE.)
    !write position of first flux surface
    CALL WriteDataToCSV((/'minval_total','maxval_total'/) ,RESHAPE(minmaxval,(/2,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"_modes") &
                               ,append_in=.TRUE.)

  END SUBROUTINE writeNow

  FUNCTION EvalSpl(deriv,xx_spl)
    INTEGER,INTENT(IN)         :: deriv !0: eval spl, 1: eval spl deriv
    REAL(wp),INTENT(IN)        :: xx_spl(:,:)
    REAL(wp)                   :: EvalSpl(n_int)
    !local
    REAL(wp)                   :: splOut(3)
    DO i=1,n_int
      CALL SPLINE1_EVAL((/1,deriv,0/), nFluxVMEC,rho_int(i),rho,xx_Spl(:,:),iGuess,splout) 
      EvalSpl(i)=splout(1+deriv)
    END DO
  END FUNCTION EvalSpl

  FUNCTION EvalSplMode(rderiv,xx_Spl)
    INTEGER,INTENT(IN)         :: rderiv !0: eval spl, 1: eval spl deriv
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
        CALL SPLINE1_EVAL((/1,rderiv,0/), nFluxVMEC,rho_p,rho,xx_Spl(:,:,iMode),iGuess,splout) 
        EvalSplMode(iMode,i)=rhom*splout(1+rderiv)+REAL(rderiv,wp)*(drhom*splout(1))
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
