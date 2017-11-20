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
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_MHD3D_visu
! MODULES
USE MOD_Globals,ONLY: wp,Unit_stdOut,abort
IMPLICIT NONE
PUBLIC



!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE visu_BC_face(mn_IP )
! MODULES
USE MOD_Globals, ONLY:TWOPI
USE MOD_MHD3D_vars
USE MOD_output_vtk, ONLY: WriteDataToVTK
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: mn_IP(2) !! muber of points in theta,zeta direction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i_m,i_n,nplot(2)
  REAL(wp) :: x(2),q(3)
  REAL(wp) :: X1_visu,X2_visu
  REAL(wp) :: x_visu(3,mn_IP(1),mn_IP(2),1)
  INTEGER,PARAMETER  :: nVal=1
  REAL(wp) :: var_visu(nVal,mn_IP(1),mn_IP(2),1)
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
!===================================================================================================================================
DO i_n=1,mn_IP(2)
  DO i_m=1,mn_IP(1)
    x=TWOPI*(/REAL(i_m-1,wp)/REAL(mn_IP(1)-1,wp),REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp)/)
    X1_visu=X1_base%f%evalDOF_x(x,0,X1_b)
    X2_visu=X2_base%f%evalDOF_x(x,0,X2_b)
    q=(/X1_visu,X2_visu,x(2)/)
    x_visu(  :,i_m,i_n,1)=hmap%eval(q)
    var_visu(1,i_m,i_n,1)=LA_base%f%evalDOF_x(x,0,LA_b)
  END DO !i_m
END DO !i_n
VarNames(1)="lambda"
nplot(:)=mn_IP-1
CALL WriteDataToVTK(2,3,nVal,nplot,1,VarNames,x_visu,var_visu,"visu_BC.vtu")

END SUBROUTINE visu_BC_face
!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE visu_planes(n_s,mn_IP )
! MODULES
USE MOD_Globals, ONLY:TWOPI
USE MOD_MHD3D_vars
USE MOD_output_vtk, ONLY: WriteDataToVTK
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: n_s,mn_IP(2) !! muber of points in theta,zeta direction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode,i_s,i_m,i_n,iElem,nElems,nplot(2)
  REAL(wp) :: s,x(2),q(3)
  REAL(wp) :: X1_s(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu
  INTEGER,PARAMETER  :: nVal=1
  REAL(wp) :: x_visu(     3,n_s,mn_IP(1),mn_IP(2)*sgrid%nElems)
  REAL(wp) :: var_visu(nVal,n_s,mn_IP(1),mn_IP(2)*sgrid%nElems)
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
!===================================================================================================================================
nElems=sgrid%nElems
DO iElem=1,nElems
  DO i_s=1,n_s
    s=sgrid%sp(iElem-1)+REAL(i_s-1,wp)/REAL(n_s-1,wp)*sgrid%ds(iElem)
    DO iMode=1,X1_base%f%modes
      X1_s(iMode)=X1_base%s%evalDOF_s(s,0,U(0)%X1(:,iMode))
    END DO
    DO iMode=1,X2_base%f%modes
      X2_s(iMode)=X2_base%s%evalDOF_s(s,0,U(0)%X2(:,iMode))
    END DO
    DO iMode=1,LA_base%f%modes
      LA_s(iMode)=LA_base%s%evalDOF_s(s,0,U(0)%LA(:,iMode))
    END DO
    DO i_n=1,mn_IP(2)
      DO i_m=1,mn_IP(1)
        x=TWOPI*(/REAL(i_m-1,wp)/REAL(mn_IP(1)-1,wp),REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp)/)
        X1_visu=X1_base%f%evalDOF_x(x,0,X1_s)
        X2_visu=X2_base%f%evalDOF_x(x,0,X2_s)
        var_visu(1,i_s,i_m,1+(i_n-1)+(iElem-1)*mn_IP(2) )=LA_base%f%evalDOF_x(x,0,LA_s)
        q=(/X1_visu,X2_visu,x(2)/)
        x_visu(:,i_s,i_m,1+(i_n-1)+(iElem-1)*mn_IP(2) )=hmap%eval(q)
      END DO !i_m
    END DO !i_n
  END DO !i_s
END DO !iElem
VarNames(1)="lambda"

nplot(:)=(/n_s,mn_IP(1)/)-1

CALL WriteDataToVTK(2,3,nVal,nplot,(mn_IP(2)*nElems),VarNames,x_visu,var_visu,"visu_planes.vtu")

WRITE(*,*)'    ... min x,y,z',MINVAL(x_visu(1,:,:,:)),MINVAL(x_visu(2,:,:,:)),MINVAL(x_visu(3,:,:,:))
WRITE(*,*)'    ... max x,y,z',MAXVAL(x_visu(1,:,:,:)),MAXVAL(x_visu(2,:,:,:)),MAXVAL(x_visu(3,:,:,:))

END SUBROUTINE visu_planes

!===================================================================================================================================
!> Visualize 
!!
!===================================================================================================================================
SUBROUTINE visu_1d_modes(n_s)
! MODULES
USE MOD_Globals,      ONLY:Eval1DPoly
USE MOD_Analyze_Vars, ONLY:visu1D
USE MOD_Output_Vars,  ONLY:ProjectName
USE MOD_base,         ONLY: t_base
USE MOD_fbase,        ONLY: sin_cos_map
USE MOD_MHD3D_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: n_s  !! number of visualization points per element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i_s,iElem,i,j,iMode,nvisu,max_modes
INTEGER            :: nVal,nValRewind,addval
CHARACTER(LEN=120),ALLOCATABLE :: varnames(:) 
REAL(wp)          ,ALLOCATABLE :: values_visu(:,:)
REAL(wp)          ,ALLOCATABLE :: s_visu(:)
LOGICAL            :: vcase(4)
CHARACTER(LEN=4)   :: vstr
CHARACTER(LEN=40)  :: fname
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

max_modes=MAXVAL((/X1_base%f%modes,X2_base%f%modes,LA_base%f%modes/))
nvisu   =sgrid%nElems*n_s+1

addval = 5
ALLOCATE(varnames(   addval+2*max_modes+2))
ALLOCATE(values_visu(addval+2*max_modes+2,nvisu))
ALLOCATE(s_visu(nvisu))

s_visu(1)=0.0_wp
j=1
DO iElem=1,sgrid%nElems
  DO i_s=1,n_s
    j=j+1
    s_visu(j)=sgrid%sp(iElem-1)+REAL(i_s,wp)/REAL(n_s,wp)*sgrid%ds(iElem)
  END DO
END DO


nVal=1
Varnames(   nVal)='Phi'
values_visu(nVal,:)=0.0_wp !TODO

nVal=nVal+1
Varnames(   nVal)='chi'
values_visu(nVal,:)=0.0_wp !TODO

nVal=nVal+1
Varnames(   nVal)='rho'
values_visu(nVal,:)=s_visu(:)

nVal=nVal+1
Varnames(nVal)='iota(Phi_norm)'
DO j=1,nvisu
  values_visu(  nVal,j)=Eval1DPoly(n_iota_coefs,iota_coefs,s_visu(j))
END DO

nVal=nVal+1
Varnames(nVal)='mass(Phi_norm)'
DO j=1,nvisu
  values_visu(  nVal,j)=Eval1DPoly(n_mass_coefs,mass_coefs,s_visu(j))
END DO

nValRewind=nVal

IF(vcase(2))THEN
  WRITE(*,*)'2) Visualize modes in 1D: R,Z,lambda interpolated...'
  nval=nValRewind
  fname="INT_X1"//TRIM(sin_cos_map(X1_base%f%sin_cos))
  CALL writeDataMN_visu(fname,"X1mn",0,s_visu,X1_base,U(0)%X1)
  nval=nValRewind
  fname="INT_X2"//TRIM(sin_cos_map(X2_base%f%sin_cos))
  CALL writeDataMN_visu(fname,"X2mn",0,s_visu,X2_base,U(0)%X2)
  nval=nValRewind
  fname="INT_LA"//TRIM(sin_cos_map(LA_base%f%sin_cos))
  CALL writeDataMN_visu(fname,"LAmn",0,s_visu,LA_base,U(0)%LA)
END IF
IF(vcase(4))THEN
  WRITE(*,*)'4) Visualize modes in 1D: dRrho,dZrho interpolated...'
  nval=nValRewind
  fname="INT_dX1"//TRIM(sin_cos_map(X1_base%f%sin_cos))
  CALL writeDataMN_visu(fname,"dX1mn",DERIV_S,s_visu,X1_base,U(0)%X1)
  nval=nValRewind
  fname="INT_dX2"//TRIM(sin_cos_map(X2_base%f%sin_cos))
  CALL writeDataMN_visu(fname,"dX2mn",DERIV_S,s_visu,X2_base,U(0)%X2)
END IF

!
DEALLOCATE(varnames)
DEALLOCATE(values_visu)
DEALLOCATE(s_visu)

CONTAINS

  SUBROUTINE writeDataMN_visu(fname,vname,rderiv,coord,base_in,xx_in)
    INTEGER,INTENT(IN)         :: rderiv !0: eval spl, 1: eval spl deriv
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    REAL(wp),INTENT(INOUT)     :: xx_in(:,:)
    TYPE(t_base),INTENT(IN)    :: base_in
    REAL(wp),INTENT(IN)        :: coord(nvisu)

    DO iMode=1,base_in%f%modes
      nVal=nVal+1
      WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname),base_in%f%Xmn(1,iMode),base_in%f%Xmn(2,iMode)/nfp
      DO j=1,nvisu
        values_visu(nVal,j)=base_in%s%evalDOF_s(s_visu(j),rderiv,xx_in(:,iMode))
      END DO !j
    END DO

    CALL writeNow(fname,vname,base_in%f%modes,base_in%f%Xmn(1,:),base_in%f%Xmn(2,:),coord,values_visu,VarNames) 

  END SUBROUTINE writeDataMN_visu

  SUBROUTINE writeNow(fname,vname,modes,xm,xn,coord,values_in,VarNames_in)
    USE MOD_Output_CSV, ONLY:WriteDataToCSV
    CHARACTER(LEN=*),INTENT(IN):: fname
    CHARACTER(LEN=*),INTENT(IN):: vname
    INTEGER         ,INTENT(IN):: modes
    INTEGER         ,INTENT(IN):: xm(1:modes)
    INTEGER         ,INTENT(IN):: xn(1:modes)
    REAL(wp),INTENT(IN)        :: coord(:)
    REAL(wp),INTENT(INOUT)        :: values_in(:,:)
    CHARACTER(LEN=*),INTENT(INOUT):: VarNames_in(:)
    !local
    REAL(wp)                   :: minmaxval(2)
    REAL(wp) ,ALLOCATABLE      :: max_loc_val(:)
    CHARACTER(LEN=100),ALLOCATABLE :: varnames_max(:)

    minmaxval(1)=MINVAL(values_in(nVal-modes:nVal,:))
    minmaxval(2)=MAXVAL(values_in(nVal-modes:nVal,:))

    DO iMode=1,modes
      nVal=nVal+1
      WRITE(VarNames_in(nVal),'(A)')TRIM(VarNames_in(nVal-modes))//'_norm'
      values_in(nVal,:)=values_in(nVal-modes,:)/(MAXVAL(ABS(values_in(nVal-modes,:)))+1.0E-12)
    END DO

    nVal=nVal+2
    Varnames_in(nVal-1)=TRIM(vname)//', m= odd, n= 000'
    Varnames_in(nVal)=  TRIM(vname)//', m=even, n= 000'
    values_in(nVal-1:nVal,:)=0.
    DO iMode=1,modes
      IF(xn(iMode).EQ.0)THEN
        IF(MOD(xm(iMode),2).NE.0)THEN
          values_in(nVal-1,:)= values_in(nVal-1,:)+values_in(nVal-2-2*modes+iMode,:)
        ELSE
          values_in(nVal,:)= values_in(nVal,:)+values_in(nVal-2-2*modes+iMode,:)
        END IF
      END IF !n=0
    END DO

    CALL WriteDataToCSV(VarNames_in(1:nVal),Values_in(1:nVal,:), (TRIM(ProjectName)//"_"//TRIM(fname)//"modes"))

    ALLOCATE(max_loc_val(nVal),Varnames_max(nVal))
    DO i=1,nVal
      max_loc_val(i)=coord(MAXLOC(ABS(values_in(i,:)),1))
      Varnames_max(i)=TRIM(VarNames_in(i))//'_maxloc'
    END DO 
    CALL WriteDataToCSV(VarNames_max(:) ,RESHAPE(max_loc_val(:),(/nval,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"modes") &
                               ,append_in=.TRUE.,vfmt_in='E10.2')
    DO i=1,nVal
      max_loc_val(i)=      MAXVAL(ABS(values_in(i,:)))+1.0E-12
      Varnames_max(i)=TRIM(VarNames_in(i))//'_maxval'
    END DO 
    CALL WriteDataToCSV(VarNames_max(:) ,RESHAPE(max_loc_val(:),(/nval,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"modes") &
                               ,append_in=.TRUE.,vfmt_in='E10.2')
    DEALLOCATE(max_loc_val,Varnames_max)
    !write position of first flux surface
    CALL WriteDataToCSV((/'rhoFirst'/) ,RESHAPE((/sgrid%sp(1)/),(/1,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"modes") &
                               ,append_in=.TRUE.)
    !write position of first flux surface
    CALL WriteDataToCSV((/'minval_total','maxval_total'/) ,RESHAPE(minmaxval,(/2,1/)) &
                               ,(TRIM(ProjectName)//"_"//TRIM(fname)//"modes") &
                               ,append_in=.TRUE.)

  END SUBROUTINE writeNow

END SUBROUTINE visu_1d_modes

END MODULE MOD_MHD3D_visu

