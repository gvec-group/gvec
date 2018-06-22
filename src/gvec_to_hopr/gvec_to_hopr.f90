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
!!# Module **gvec_to_hopr**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_hopr
! MODULES
USE MODgvec_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE Init_gvec_to_hopr
  MODULE PROCEDURE Init_gvec_to_hopr
END INTERFACE
!
INTERFACE gvec_to_hopr
  MODULE PROCEDURE gvec_to_hopr
END INTERFACE

INTERFACE Finalize_gvec_to_hopr
  MODULE PROCEDURE Finalize_gvec_to_hopr
END INTERFACE

PUBLIC::Init_gvec_to_hopr
PUBLIC::gvec_to_hopr
PUBLIC::Finalize_gvec_to_hopr
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE Init_gvec_to_hopr(fileName) 
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_gvec_to_hopr_Vars
USE MODgvec_readState,ONLY: ReadState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileName !< name of GVEC file
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT GVEC TO HOPR ...'

  CALL ReadState(fileName)

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE Init_gvec_to_hopr


!===================================================================================================================================
!> Evaluate gvec state at a list of s,theta,zeta positions
!!
!===================================================================================================================================
SUBROUTINE gvec_to_hopr(nNodes,xIn,xOut,data_out,phi_axis_edge,chi_axis_edge)
! MODULES
USE MODgvec_Globals, ONLY: CROSS
USE MODgvec_gvec_to_hopr_vars
USE MODgvec_readState_vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER          :: nNodes
REAL,INTENT( IN) :: xIn(3,nNodes)  !!s=sqrt(psi_norm),theta,zeta positions for evaluation, psi_norm is normalized toroidal flux
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: xOut(3,nNodes)  !! x,y,z cartesian coordinates
REAL,INTENT(OUT) :: data_out(9,nNodes)  !! pressure,Bcart(3),chi,phi,Acart(3)
REAL,INTENT(OUT) :: phi_axis_edge(2)
REAL,INTENT(OUT) :: chi_axis_edge(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,iMode
REAL    :: spos,thet,zeta,X1_int,X2_int,LA_int
REAL(wp),DIMENSION(1:X1_base_r%f%modes) :: X1_s,dX1ds_s
REAL(wp),DIMENSION(1:X2_base_r%f%modes) :: X2_s,dX2ds_s
REAL(wp),DIMENSION(1:LA_base_r%f%modes) :: LA_s
REAL    :: dX1ds   ,dX2ds
REAL    :: dX1dthet,dX2dthet
REAL    :: dX1dzeta,dX2dzeta
REAL    :: dLAdthet,dLAdzeta
REAL    :: phi_int,chi_int,pres_int,iota_int
REAL    :: phiPrime_int,ChiPrime_int         !prime refers to d/ds , where s=sqrt(phi_norm)
REAL    :: sqrtG
REAL    :: Bcart(3),Acart(3),qvec(3)
REAL    :: e_s(3),e_thet(3),e_zeta(3)
REAL    :: grad_s(3),grad_thet(3),grad_zeta(3)
!===================================================================================================================================
phi_axis_edge(1)= X1_base_r%s%evalDOF_s(1.0e-08, 0,profiles_1d(:,1))
chi_axis_edge(1)= X1_base_r%s%evalDOF_s(1.0e-08, 0,profiles_1d(:,2))
phi_axis_edge(2)= X1_base_r%s%evalDOF_s(1.0, 0,profiles_1d(:,1))
chi_axis_edge(2)= X1_base_r%s%evalDOF_s(1.0, 0,profiles_1d(:,2))
DO iNode=1,nNodes
  spos=xIn(1,iNode)
  thet=xIn(2,iNode)
  zeta=xIn(3,iNode)

  phi_int      = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,1))
  chi_int      = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,2))
  iota_int     = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))
  pres_int     = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,4))
  PhiPrime_int = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  ChiPrime_int = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2))

  DO iMode=1,X1_base_r%f%modes
    X1_s(   iMode) =X1_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode))
    dX1ds_s(iMode) =X1_base_r%s%evalDOF_s(spos,DERIV_S,X1_r(:,iMode))
  END DO
  DO iMode=1,X2_base_r%f%modes
    X2_s(   iMode) =X2_base_r%s%evalDOF_s(spos,      0,X2_r(:,iMode))
    dX2ds_s(iMode) =X2_base_r%s%evalDOF_s(spos,DERIV_S,X2_r(:,iMode))
  END DO
  DO iMode=1,LA_base_r%f%modes
    LA_s(   iMode) =LA_base_r%s%evalDOF_s(spos,      0,LA_r(:,iMode))
  END DO
  X1_int     = X1_base_r%f%evalDOF_x((/thet,zeta/),         0,X1_s)
  dX1ds      = X1_base_r%f%evalDOF_x((/thet,zeta/),         0,dX1ds_s)
  dX1dthet   = X1_base_r%f%evalDOF_x((/thet,zeta/),DERIV_THET,X1_s)
  dX1dzeta   = X1_base_r%f%evalDOF_x((/thet,zeta/),DERIV_ZETA,X1_s)
  X2_int     = X2_base_r%f%evalDOF_x((/thet,zeta/),         0,X2_s)
  dX2ds      = X2_base_r%f%evalDOF_x((/thet,zeta/),         0,dX2ds_s)
  dX2dthet   = X2_base_r%f%evalDOF_x((/thet,zeta/),DERIV_THET,X2_s)
  dX2dzeta   = X2_base_r%f%evalDOF_x((/thet,zeta/),DERIV_ZETA,X2_s)
  LA_int     = LA_base_r%f%evalDOF_x((/thet,zeta/),         0,LA_s)
  dLAdthet   = LA_base_r%f%evalDOF_x((/thet,zeta/),DERIV_THET,LA_s)
  dLAdzeta   = LA_base_r%f%evalDOF_x((/thet,zeta/),DERIV_ZETA,LA_s)

  qvec=(/X1_int,X2_int,zeta/)
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,0.0_wp/))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 

  Bcart(:)=  (  e_thet(:)*(iota_int-dLAdzeta )  &
              + e_zeta(:)*(1.0_wp+dLAdthet) )*(PhiPrime_int/sqrtG)
  grad_s    = CROSS(e_thet,e_zeta) !/sqrtG
  grad_thet = CROSS(e_zeta,e_s   ) !/sqrtG
  grad_zeta = CROSS(e_s   ,e_thet) !/sqrtG

  Acart(:)=  ( phi_int*grad_thet(:)-(LA_int*PhiPrime_int)*grad_s(:)  -chi_int*grad_zeta)/sqrtG

  xOut(:,iNode)=hmap_r%eval(qvec)

  data_out(  1,iNode)=pres_int
  data_out(2:4,iNode)=Bcart(:)
  data_out(  5,iNode)=chi_int
  data_out(  6,iNode)=phi_int
  data_out(7:9,iNode)=Acart(:)
END DO

END SUBROUTINE gvec_to_hopr


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE Finalize_gvec_to_hopr 
! MODULES
USE MODgvec_ReadState,ONLY:Finalize_ReadState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL Finalize_ReadState()
END SUBROUTINE Finalize_gvec_to_hopr

END MODULE MODgvec_gvec_to_hopr
