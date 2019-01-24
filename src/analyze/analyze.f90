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
MODULE MODgvec_Analyze
! MODULES
USE MODgvec_Globals, ONLY:wp,abort
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
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_Analyze_Vars
USE MODgvec_ReadInTools,ONLY:GETINT,GETINTARRAY,GETREALARRAY,GETLOGICAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp):: visu_minmax(3,0:1)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'INIT ANALYZE ...'
visu1D    = GETINT('visu1D',Proposal=0)   
visu2D    = GETINT('visu2D',Proposal=0)   
visu3D    = GETINT('visu3D',Proposal=0)   
SFL_theta = GETLOGICAL('SFL_theta',Proposal=.TRUE.)   


visu_minmax(1:3,0)=GETREALARRAY("visu_min",3,Proposal=(/0.0_wp,0.0_wp,0.0_wp/))
visu_minmax(1:3,1)=GETREALARRAY("visu_max",3,Proposal=(/1.0_wp,1.0_wp,1.0_wp/))

iAnalyze=-1

IF(visu1D.GT.0)THEN
  np_1d          = GETINT(     "np_1d",Proposal=5)
  IF(np_1d.LE.1) CALL abort(__STAMP__,&
     "np_1d must be >1")
END IF
IF(visu2D.GT.0)THEN
  np_visu_BC     = GETINTARRAY("np_visu_BC",2,Proposal=(/20,30/))
  IF(any(np_visu_BC.LE.1)) CALL abort(__STAMP__,&
     "all point numbers in np_visu_BC must be >1")
  visu_BC_minmax(2:3,0)=GETREALARRAY("visu_BC_min",2,Proposal=visu_minmax(2:3,0),quiet_def_in=.TRUE.)
  visu_BC_minmax(2:3,1)=GETREALARRAY("visu_BC_max",2,Proposal=visu_minmax(2:3,1),quiet_def_in=.TRUE.)
  np_visu_planes = GETINTARRAY("np_visu_planes",3, (/5,12,10/))
  IF(any(np_visu_planes.LE.1)) CALL abort(__STAMP__,&
     "all point numbers in np_visu_planes must be >1")
  visu_planes_minmax(1:3,0)=GETREALARRAY("visu_planes_min",3,Proposal=visu_minmax(1:3,0),quiet_def_in=.TRUE.)
  visu_planes_minmax(1:3,1)=GETREALARRAY("visu_planes_max",3,Proposal=visu_minmax(1:3,1),quiet_def_in=.TRUE.)
END IF
IF(visu3D.GT.0)THEN
  np_visu_3D     = GETINTARRAY("np_visu_3D",3,Proposal=(/5,12,10/))
  IF(any(np_visu_3D.LE.1)) CALL abort(__STAMP__,&
     "all point numbers in np_visu_3D must be >1")
  visu_3D_minmax(1:3,0)=GETREALARRAY("visu_3D_min",3,Proposal=visu_minmax(1:3,0),quiet_def_in=.TRUE.)
  visu_3D_minmax(1:3,1)=GETREALARRAY("visu_3D_max",3,Proposal=visu_minmax(1:3,1),quiet_def_in=.TRUE.)
END IF
SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitAnalyze


!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE Analyze(fileID_in)
! MODULES
USE MODgvec_Analyze_Vars
USE MODgvec_MHD3D_Vars, ONLY:which_init
USE MODgvec_mhd3d_visu
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN), OPTIONAL :: fileID_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: vcase(4)
CHARACTER(LEN=4)   :: vstr
INTEGER            :: FileID
!===================================================================================================================================
iAnalyze=iAnalyze+1
IF(PRESENT(FileID_in))THEN
  FileID=FileID_in
ELSE
  FileID=iAnalyze
END IF
IF(iAnalyze.EQ.0) THEN
  IF(which_init.EQ.1) THEN
    IF(visu1D.NE.0) CALL VMEC1D_visu() 
    IF(visu2D.NE.0) CALL VMEC3D_visu(np_visu_planes,visu_planes_minmax,.TRUE. ) 
    IF(visu3D.NE.0) CALL VMEC3D_visu(np_visu_3D    ,visu_3D_minmax    ,.FALSE.) 
  END IF
END IF !iAnalyze==0
IF(visu1D.NE.0)THEN
  CALL visu_1d_modes(np_1d,FileID) 
  !
END IF !visu1D
IF(visu2D.NE.0)THEN
  vcase=.FALSE.
  WRITE(vstr,'(I4)')visu2D
  IF(INDEX(vstr,'1').NE.0) vcase(1)=.TRUE.
  IF(INDEX(vstr,'2').NE.0) vcase(2)=.TRUE.
  IF(INDEX(vstr,'3').NE.0) vcase(3)=.TRUE.
  IF(INDEX(vstr,'4').NE.0) vcase(4)=.TRUE.
  IF(vcase(1))THEN
    IF(iAnalyze.EQ.0) CALL visu_BC_face(np_visu_BC(1:2),visu_BC_minmax(:,:),FileID)
  END IF
  IF(vcase(2))THEN
    CALL visu_3D(np_visu_planes,visu_planes_minmax,.TRUE.,FileID) !only planes
  END IF 
END IF !visu2d
IF(visu3D.NE.0)THEN
  vcase=.FALSE.
  WRITE(vstr,'(I4)')visu3D
  IF(INDEX(vstr,'1').NE.0) vcase(1)=.TRUE.
  IF(INDEX(vstr,'2').NE.0) vcase(2)=.TRUE.
  IF(INDEX(vstr,'3').NE.0) vcase(3)=.TRUE.
  IF(INDEX(vstr,'4').NE.0) vcase(4)=.TRUE.
  IF(vcase(1))THEN
    CALL visu_3D(np_visu_3D,visu_3D_minmax,.FALSE.,FileID) !full 3D
  END IF 
END IF !visu2d

END SUBROUTINE Analyze 


!===================================================================================================================================
!> Visualize VMEC flux surface data for each mode, for Rmnc 
!!
!===================================================================================================================================
SUBROUTINE VMEC1D_visu()
! MODULES
USE MODgvec_Globals,ONLY:Pi
USE MODgvec_Analyze_Vars, ONLY:visu1D
USE MODgvec_Output_Vars, ONLY:ProjectName
USE MODgvec_write_modes
USE MODgvec_VMEC_Readin
USE MODgvec_VMEC_Vars
USE MODgvec_VMEC, ONLY: VMEC_EvalSpl,VMEC_EvalSplMode
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,nVal,nValRewind,iMode
INTEGER,PARAMETER  :: n_Int=200
CHARACTER(LEN=120) :: varnames(  8+2*mn_mode) 
REAL(wp)           :: values(    8+2*mn_mode,nFluxVMEC) 
REAL(wp)           :: values_int(8+2*mn_mode,n_int) 
REAL(wp)           :: rho_int(n_int),rho_half(nFluxVMEC)
LOGICAL            :: vcase(4)
CHARACTER(LEN=4)   :: vstr
!===================================================================================================================================
  !visu1D: all possible combinations: 1,2,3,4,12,13,14,23,24,34,123,124,234,1234
  WRITE(vstr,'(I4)')visu1D
  vcase=.FALSE.
  IF(INDEX(vstr,'1').NE.0) vcase(1)=.TRUE.
  IF(INDEX(vstr,'2').NE.0) vcase(2)=.TRUE.
  IF(INDEX(vstr,'3').NE.0) vcase(3)=.TRUE.
  IF(INDEX(vstr,'4').NE.0) vcase(4)=.TRUE.
  IF(.NOT.(ANY(vcase))) THEN
    WRITE(*,*)'visu1D case not found:',visu1D,' nothing visualized...'
    RETURN
  END IF
  
  !interpolation points
  DO i=0,n_int-1
    rho_int(1+i)=REAL(i,wp)/REAL(n_int-1,wp)
  END DO
  !strech towards axis and edge
  rho_int=rho_int+0.05_wp*SIN(Pi*(2.0_wp*rho_int-1.0_wp))
  
  nVal=1
  Varnames(nVal)='Phi'
  values(  nVal,:)=Phi_prof(:)
  DO i=1,n_int
    values_int(nVal,i)=VMEC_EvalSpl(0,rho_int(i),Phi_Spl)
  END DO !i
  
  nVal=nVal+1
  Varnames(nVal)='chi'
  values(  nVal,:)=Chi_prof(:)
  DO i=1,n_int
    values_int(nVal,i)=VMEC_EvalSpl(0,rho_int(i),chi_Spl)
  END DO !i
  
  nVal=nVal+1
  Varnames(nVal)='rho'
  values(  nVal,:)=rho(:)
  values_int(nVal,:)=rho_int(:)
  
  rho_half(1)=0.
  DO i=1,nFluxVMEC-1
    rho_half(i+1)=SQRT(0.5_wp*(NormFlux_prof(i+1)+NormFlux_prof(i))) !0.5*(rho(iFlux)+rho(iFlux+1))
  END DO
  nVal=nVal+1
  Varnames(nVal)='rho_half'
  values(  nVal,:)= rho_half(:)
  values_int(nVal,:)=0.
  
  nVal=nVal+1
  Varnames(nVal)='iota(Phi_norm)'
  values(  nVal,:)=iotaf(:)
  DO i=1,n_int
    values_int(nVal,i)=VMEC_EvalSpl(0,rho_int(i),iota_Spl)
  END DO !i
  
  nVal=nVal+1
  Varnames(nVal)='pres(Phi_norm)'
  values(  nVal,:)=presf(:)
  DO i=1,n_int
    values_int(nVal,i)=VMEC_EvalSpl(0,rho_int(i),pres_Spl)
  END DO !i
  
  nValRewind=nVal
  
  
  IF(vcase(1))THEN
    WRITE(*,*)'1) Visualize VMEC modes R,Z,lambda interpolated...'
    nval=nValRewind
    CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Rmnc","Rmnc",0,rho_int,Rmnc_Spl)
    nval=nValRewind
    CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Zmns","Zmns",0,rho_int,Zmns_Spl)
    nval=nValRewind
    IF(reLambda.OR.(lambda_grid.EQ."full"))THEN
      CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Lmns","Lmns",0,rho_int,Lmns_Spl)
    ELSE
      CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Lmns_half","Lmns_h",0,rho_int,Lmns_spl)
    END IF
    IF(lasym)THEN
      nval=nValRewind
      CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Rmns","Rmns",0,rho_int,Rmns_Spl)
      nval=nValRewind
      CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Zmnc","Zmnc",0,rho_int,Zmnc_Spl)
      nval=nValRewind
      IF(reLambda.OR.(lambda_grid.EQ."full"))THEN
        CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Lmnc","Lmnc",0,rho_int,Lmnc_Spl)
      ELSE
        CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_Lmnc_half","Lmnc_h",0,rho_int,Lmnc_spl)
      END IF
    END IF!lasym
  END IF !vcase(1)
  IF(vcase(2))THEN
    WRITE(*,*)'2) Visualize VMEC modes dRrho,dZrho interpolated...'
    nval=nValRewind
    CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_dRmnc","dRmnc",1,rho_int,Rmnc_Spl)
    nval=nValRewind
    CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_dZmns","dZmns",1,rho_int,Zmns_Spl)
    IF(lasym)THEN
      !interpolated profiles
      nval=nValRewind
      CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_dRmns","dRmns",1,rho_int,Rmns_Spl)
      nval=nValRewind
      CALL writeDataMN_int(TRIM(ProjectName)//"_VMEC_INT_dZmnc","dZmnc",1,rho_int,Zmnc_Spl)
    END IF!lasym
  END IF !vcase(2)
  IF(vcase(3))THEN
    WRITE(*,*)'3) Visualize VMEC modes R,Z,lambda pointwise ...'
    nval=nValRewind
    CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Rmnc","Rmnc",0,rho,Rmnc)
    nval=nValRewind
    CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Zmns","Zmns",0,rho,Zmns)
    nval=nValRewind
    IF(reLambda.OR.(lambda_grid.EQ."full"))THEN
      CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Lmns","Lmns",0,rho,Lmns)
    ELSE
      CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Lmns_half","Lmns_h",0,rho_half,Lmns)
    END IF
    IF(lasym)THEN
      nval=nValRewind
      CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Rmns","Rmns",0,rho,Rmns)
      nval=nValRewind
      CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Zmnc","Zmnc",0,rho,Zmnc)
      nval=nValRewind
      IF(reLambda.OR.(lambda_grid.EQ."full"))THEN
        CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Lmnc","Lmnc",0,rho,Lmnc)
      ELSE
        CALL writeDataMN(TRIM(ProjectName)//"_VMEC_Lmnc_half","Lmnc_h",0,rho_half,Lmnc)
      END IF
    END IF!lasym
  END IF !vcase(3)
  IF(vcase(4))THEN
    WRITE(*,*)'4) Visualize VMEC modes dRrho,dZrho pointwise (1st order finite difference)...'
    nval=nValRewind
    CALL writeDataMN(TRIM(ProjectName)//"_VMEC_dRmnc","dRmnc",1,rho,Rmnc)
    nval=nValRewind
    CALL writeDataMN(TRIM(ProjectName)//"_VMEC_dZmns","dZmns",1,rho,Zmns)
    IF(lasym)THEN
      nval=nValRewind
      CALL writeDataMN(TRIM(ProjectName)//"_VMEC_dRmns","dRmns",1,rho,Rmns)
      nval=nValRewind
      CALL writeDataMN(TRIM(ProjectName)//"_VMEC_dZmnc","dZmnc",1,rho,Zmnc)
    END IF!lasym
  END IF !vcase(4)
  
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
    CALL write_modes(TRIM(fname)//'.csv',vname,nVal,mn_mode,INT(xm),INT(xn),coord,rho(2),values,VarNames) 

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
      DO i=1,n_int
        values_int(    nVal,i)  =VMEC_EvalSplMode((/iMode/),rderiv,coord(i),xx_Spl)
      END DO
    END DO

    CALL write_modes(TRIM(fname)//'.csv',vname,nVal,mn_mode,INT(xm),INT(xn),coord,rho(2),values_int,VarNames) 

  END SUBROUTINE writeDataMN_int


END SUBROUTINE VMEC1D_visu 

!===================================================================================================================================
!> Visualize VMEC flux surface data  in planes or 3D, number of radial posisiotns fixed to nFluxVMEC+1, only R,Z,Lambda
!!
!===================================================================================================================================
SUBROUTINE VMEC3D_visu(np_in,minmax,only_planes)
! MODULES
USE MODgvec_Analyze_Vars,ONLY:SFL_theta
USE MODgvec_Globals,ONLY:TWOPI,PI,UNIT_stdOut
USE MODgvec_VMEC_Readin
USE MODgvec_VMEC_Vars
USE MODgvec_Output_Vars,ONLY:Projectname
USE MODgvec_Output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_Newton,     ONLY: NewtonRoot1D_FdF
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: np_in(3)     !! (1) #points in s & theta per element,(2:3) #elements in  theta,zeta
  REAL(wp)      , INTENT(IN   ) :: minmax(3,0:1)  !! minimum /maximum range in s,theta,zeta [0,1] 
  LOGICAL       , INTENT(IN   ) :: only_planes  !! true: visualize only planes, false:  full 3D
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: i_s,j_s,i_m,i_n,iMode,nplot(3),mn_IP(2)
  INTEGER,PARAMETER  :: nVal=6
  REAL(wp)           :: coord_visu(     3,nFluxVMEC,np_in(1),np_in(3),np_in(2))
  REAL(wp)           :: var_visu(    nVal,nFluxVMEC,np_in(1),np_in(3),np_in(2))
  REAL(wp)           :: thet(np_in(1),np_in(2)),zeta(np_in(3)),R,Z,LA,sinmn(mn_mode),cosmn(mn_mode)
  REAL(wp)           :: xIP(2),theta_star
  REAL(wp)           :: cosmn_nyq,sqrtg,absB
  CHARACTER(LEN=40)  :: VarNames(nVal)          !! Names of all variables that will be written out
  CHARACTER(LEN=255) :: filename
!===================================================================================================================================
  IF(only_planes)THEN
    SWRITE(UNIT_stdOut,'(A)') 'Start VMEC visu planes...'
  ELSE
    SWRITE(UNIT_stdOut,'(A)') 'Start VMEC visu 3D...'
  END IF
  IF((minmax(1,1)-minmax(1,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
     'WARNING visu3D, nothing to visualize since s-range is <=0, s_min= ',minmax(1,0),', s_max= ',minmax(1,1)
    RETURN
  ELSEIF((minmax(2,1)-minmax(2,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visu3D, nothing to visualize since theta-range is <=0, theta_min= ',minmax(2,0),', theta_max= ',minmax(2,1)
    RETURN
  ELSEIF((minmax(3,1)-minmax(3,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visu3D, nothing to visualize since zeta-range is <=0, zeta_min= ',minmax(3,0),', zeta_max= ',minmax(3,1)
    RETURN
  END IF
  VarNames(1)="Phi"
  VarNames(2)="iota"
  VarNames(3)="pressure"
  VarNames(4)="lambda"
  VarNames(5)="Jacobian"
  VarNames(6)="|B|"
  mn_IP(1:2)=np_in(2:3)
  ASSOCIATE(n_s=>nFluxVMEC )
  DO i_m=1,mn_IP(1)
    DO j_s=1,np_in(1)
      thet(j_s,i_m)=TWOPI*(minmax(2,0)+(minmax(2,1)-minmax(2,0)) &
                            *REAL((j_s-1)+(i_m-1)*(np_in(1)-1),wp)/REAL((np_in(1)-1)*mn_IP(1),wp))
    END DO !j_s
  END DO
  IF(ABS((minMax(2,1)-minmax(2,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
    thet(np_in(1),mn_IP(1))=thet(1,1)
  END IF
  DO i_n=1,mn_IP(2)
    zeta(i_n)=TWOPI*(minmax(3,0)+(minmax(3,1)-minmax(3,0))*REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp))
  END DO
  IF(ABS((minMax(3,1)-minmax(3,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
    zeta(mn_IP(2))=zeta(1)
  END IF
  
  DO i_s=1,n_s
    var_visu(1,i_s,:,:,:)=Phi_Prof(i_s)
    var_visu(2,i_s,:,:,:)=iotaf(i_s)
    var_visu(3,i_s,:,:,:)=presf(i_s)
  END DO !i_s
  DO i_n=1,mn_IP(2)
    xIP(2)=zeta(i_n)
    DO i_m=1,mn_IP(1)
      DO j_s=1,np_in(1)
        !xIP(1)=thet(j_s,i_m) 
        DO i_s=1,n_s
          !SFL
          IF(SFL_theta)THEN
            theta_star=thet(j_s,i_m)
            XIP(1)=NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1_wp*PI,theta_star   , theta_star,FRdFR)
          ELSE
            XIP(1)=thet(j_s,i_m)
          END IF
          
          DO iMode=1,mn_mode
            sinmn(iMode)=SIN(xm(iMode)*xIP(1)-xn(iMode)*xIP(2))
            cosmn(iMode)=COS(xm(iMode)*xIP(1)-xn(iMode)*xIP(2))
          END DO !iMode
          R=0.0_wp
          Z=0.0_wp
          LA=0.0_wp
          DO iMode=1,mn_mode
            R =R +Rmnc(iMode,i_s)*cosmn(iMode)
            Z =Z +Zmns(iMode,i_s)*sinmn(iMode)
            LA=LA+Lmns(iMode,i_s)*sinmn(iMode)
          END DO !iMode
          IF(lasym)THEN
            DO iMode=1,mn_mode
              R =R +Rmns(iMode,i_s)*sinmn(iMode)
              Z =Z +Zmnc(iMode,i_s)*cosmn(iMode)
              LA=LA+Lmnc(iMode,i_s)*cosmn(iMode)
            END DO !iMode
          END IF !lasym
          coord_visu(1,i_s,j_s,i_n,i_m) = R*COS(xIP(2))
          coord_visu(2,i_s,j_s,i_n,i_m) =-R*SIN(xIP(2)) !vmec data(R,phi,Z) was flipped in input to (R,Z,phi)!!
          coord_visu(3,i_s,j_s,i_n,i_m) = Z
          var_visu(  4,i_s,j_s,i_n,i_m) = LA
          sqrtg=0.
          absB =0.
          DO iMode=1,mn_mode_nyq
            cosmn_nyq=COS(xm_nyq(iMode)*xIP(1)-xn_nyq(iMode)*xIP(2))
            sqrtg=sqrtg+gmnc(iMode,i_s)*cosmn_nyq
            absB =absB +bmnc(iMode,i_s)*cosmn_nyq
          END DO !iMode
          var_visu(  5,i_s,j_s,i_n,i_m) = sqrtg*2.0*sqrt(phi_Prof(i_s)/phi_prof(n_s))  !VMEC: s=Phi_norm, but should match GVEC s~sqrt(phi_norm)
          var_visu(  6,i_s,j_s,i_n,i_m) = absB
        END DO !i_s=1,n_s
      END DO !j_s=1,np_in(1)
    END DO !i_n
  END DO !i_m
  !overwrite data on axis with theta average of first flux surface (index 2)
  DO i_n=1,mn_IP(2)
    var_visu(4,1,:,i_n,:) =SUM(var_visu(4,2,:,i_n,:))/REAL(mn_IP(1)*np_in(1))
    var_visu(5,1,:,i_n,:) =SUM(var_visu(5,2,:,i_n,:))/REAL(mn_IP(1)*np_in(1))
    var_visu(6,1,:,i_n,:) =SUM(var_visu(6,2,:,i_n,:))/REAL(mn_IP(1)*np_in(1))
  END DO !i_m

  IF(only_planes)THEN
    nplot(1:2)=(/n_s,np_in(1) /)-1
    WRITE(filename,'(A,A)')TRIM(Projectname),"_visu_planes_VMEC.vtu"
    CALL WriteDataToVTK(2,3,nVal,nplot(1:2),mn_IP(1)*mn_IP(2),VarNames, &
                        coord_visu(:,:,:,:,:), &
                          var_visu(:,:,:,:,:),TRIM(filename))
  ELSE
    !3D
    nplot(1:3)=(/n_s,np_in(1),mn_IP(2)/)-1
    WRITE(filename,'(A,A)')TRIM(Projectname),"_visu_3D_VMEC.vtu"
    CALL WriteDataToVTK(3,3,nVal,nplot,mn_IP(1),VarNames, &
                        coord_visu(:,:,:,:,:), &
                          var_visu(:,:,:,:,:),TRIM(filename))
  END IF

  END ASSOCIATE

  SWRITE(UNIT_stdOut,'(A)') '... DONE.'

!for iteration on theta^*
CONTAINS 

  FUNCTION FRdFR(theta_iter)
    !uses current zeta=XIP(2),i_s where newton is called
    IMPLICIT NONE
    REAL(wp) :: theta_iter
    REAL(wp) :: FRdFR(2) !output
    !--------------------------------------------------- 
    DO iMode=1,mn_mode
      sinmn(iMode)=SIN(xm(iMode)*theta_iter-xn(iMode)*xIP(2))
      cosmn(iMode)=COS(xm(iMode)*theta_iter-xn(iMode)*xIP(2))
    END DO !iMode
    FRdFR(1)=0. !LA
    FRdFR(2)=0. !dLA/dtheta
    DO iMode=1,mn_mode
      FRdFR(1)=FRdFR(1)          +Lmns(iMode,i_s)*sinmn(iMode)
      FRdFR(2)=FRdFR(2)+xm(iMode)*Lmns(iMode,i_s)*cosmn(iMode)
    END DO !iMode
    IF(lasym)THEN
      DO iMode=1,mn_mode
        FRdFR(1)=FRdFR(1)          +Lmnc(iMode,i_s)*cosmn(iMode)
        FRdFR(2)=FRdFR(2)-xm(iMode)*Lmnc(iMode,i_s)*sinmn(iMode)
      END DO !iMode
    END IF !lasym
    FRdFR(1)=FRdFR(1)+theta_iter
    FRdFR(2)=FRdFR(2)+1.0_wp
  END FUNCTION FRdFR

END SUBROUTINE VMEC3D_visu 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeAnalyze 
! MODULES
USE MODgvec_Analyze_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE FinalizeAnalyze

END MODULE MODgvec_Analyze
