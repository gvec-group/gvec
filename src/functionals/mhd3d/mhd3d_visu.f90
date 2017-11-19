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


END MODULE MOD_MHD3D_visu

