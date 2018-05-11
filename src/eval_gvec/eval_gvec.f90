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
!!# Module **Eval_GVEC**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_Eval_GVEC
! MODULES
USE MODgvec_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitEval_GVEC
  MODULE PROCEDURE InitEval_GVEC
END INTERFACE
!
INTERFACE Eval_GVEC
  MODULE PROCEDURE Eval_GVEC
END INTERFACE

INTERFACE FinalizeEval_GVEC
  MODULE PROCEDURE FinalizeEval_GVEC
END INTERFACE

PUBLIC::InitEval_GVEC
PUBLIC::Eval_GVEC
PUBLIC::FinalizeEval_GVEC
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitEval_GVEC(fileName) 
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_Eval_GVEC_Vars
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
  SWRITE(UNIT_stdOut,'(A)')'INIT EVAL GVEC ...'

  CALL ReadState(fileName)

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitEval_GVEC


!===================================================================================================================================
!> read an input solution and initialize U(0) (X1,X2,LA) of size X1/X2/LA_base , from an ascii .dat file 
!! if size of grid/X1/X2/LA  not equal X1/X2/X3_base
!! interpolate readin solution to the current base of Uin
!!
!===================================================================================================================================
SUBROUTINE ReadState(fileString)
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MODgvec_Eval_GVEC_Vars
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base, base_new
USE MODgvec_fbase,  ONLY: sin_cos_map 
USE MODgvec_hmap,  ONLY: hmap_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=*)    , INTENT(IN   ) :: fileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: fileID_r,OutputLevel_r
  INTEGER              :: ioUnit,iMode,is,nElems_r,grid_type_r,nfp_r,degGP_r,mn_nyq_r(2),which_hmap_r 
  INTEGER              :: X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  INTEGER              :: X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  INTEGER              :: LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r
  INTEGER,ALLOCATABLE  :: X1_mn_r(:,:),X2_mn_r(:,:),LA_mn_r(:,:)
  REAL(wp),ALLOCATABLE :: sp_r(:),profiles_IP(:,:)
  INTEGER              :: X1_mn_max_r(2),X2_mn_max_r(2),LA_mn_max_r(2)
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)')'   READ SOLUTION VARIABLE FROM FILE    "'//TRIM(FileString)//'" ...'

  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit         ,&
     FILE     = TRIM(FileString) ,&
     STATUS   = 'OLD'            ,&
     ACTION   = 'READ'           ,&
     ACCESS   = 'SEQUENTIAL' ) 

  READ(ioUnit,*) !## MHD3D Solution file
  READ(ioUnit,*) outputLevel_r,fileID_r
  READ(ioUnit,*) !## grid: nElems, grid_type
  READ(ioUnit,*) nElems_r,grid_type_r
  ALLOCATE(sp_r(0:nElems_r))

  READ(ioUnit,*) !## grid: sp(0:nElems)
  READ(ioUnit,*)sp_r(:)
  READ(ioUnit,*) !## global: nfp, degGP, mn_nyq
  READ(ioUnit,*) nfp_r, degGP_r,mn_nyq_r,which_hmap_r
  READ(ioUnit,*) !## X1_base: 
  READ(ioUnit,*) X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  READ(ioUnit,*) !## X2_base:                 
  READ(ioUnit,*) X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  READ(ioUnit,*) !## LA_base:                 
  READ(ioUnit,*) LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r

  ALLOCATE(X1_r(1:X1_nbase_r,1:X1_modes_r))
  ALLOCATE(X2_r(1:X2_nbase_r,1:X2_modes_r))
  ALLOCATE(LA_r(1:LA_nbase_r,1:LA_modes_r))

  ALLOCATE(profiles_IP(1:X1_nbase_r,5))
  ALLOCATE(X1_mn_r(2,1:X1_modes_r))
  ALLOCATE(X2_mn_r(2,1:X2_modes_r))
  ALLOCATE(LA_mn_r(2,1:LA_modes_r))
  READ(ioUnit,*) !## X1: 
  DO iMode=1,X1_modes_r
    READ(ioUnit,*)X1_mn_r(:,iMode),X1_r(:,iMode)
  END DO
  READ(ioUnit,*) !## X2: 
  DO iMode=1,X2_modes_r
    READ(ioUnit,*)X2_mn_r(:,iMode),X2_r(:,iMode)
  END DO
  READ(ioUnit,*) !## LA: 
  DO iMode=1,LA_modes_r
    READ(ioUnit,*)LA_mn_r(:,iMode),LA_r(:,iMode)
  END DO
  READ(ioUnit,*) !## profiles at X1_base IP points : spos,phi,chi,iota,pressure 
  DO is=1,X1_nbase_r
    READ(ioUnit,*)profiles_IP(is,:)
  END DO

  CLOSE(ioUnit)

  CALL hmap_new(hmap_r,which_hmap_r)
  ! check if input has changed:

  CALL sgrid_r%init(nElems_r,grid_type_r)

  !needed to build base of restart file
  X1_mn_max_r = (/MAXVAL(X1_mn_r(1,:)),MAXVAL(X1_mn_r(2,:))/nfp_r/)
  X2_mn_max_r = (/MAXVAL(X2_mn_r(1,:)),MAXVAL(X2_mn_r(2,:))/nfp_r/)
  LA_mn_max_r = (/MAXVAL(LA_mn_r(1,:)),MAXVAL(LA_mn_r(2,:))/nfp_r/)

  CALL base_new(X1_base_r,X1_deg_r,X1_cont_r,sgrid_r,degGP_r,X1_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(X1_sin_cos_r),(X1_excl_mn_zero_r.EQ.1))
  CALL base_new(X2_base_r,X2_deg_r,X2_cont_r,sgrid_r,degGP_r,X2_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(X2_sin_cos_r),(X2_excl_mn_zero_r.EQ.1))
  CALL base_new(LA_base_r,LA_deg_r,LA_cont_r,sgrid_r,degGP_r,LA_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(LA_sin_cos_r),(LA_excl_mn_zero_r.EQ.1))

  ALLOCATE(profiles_1d(1:X1_nbase_r,4))
  !convert to spline DOF
  profiles_1d(:,1) =X1_base_r%s%initDOF( profiles_IP(:,1+1) ) !phi
  profiles_1d(:,2) =X1_base_r%s%initDOF( profiles_IP(:,1+2) ) !chi
  profiles_1d(:,3) =X1_base_r%s%initDOF( profiles_IP(:,1+3) ) !iota
  profiles_1d(:,4) =X1_base_r%s%initDOF( profiles_IP(:,1+4) ) !pressure

  DEALLOCATE(sp_r,profiles_IP)
  DEALLOCATE(X1_mn_r)
  DEALLOCATE(X2_mn_r)
  DEALLOCATE(LA_mn_r)


  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE ReadState


!===================================================================================================================================
!> Evaluate gvec state at a list of s,theta,zeta positions
!!
!===================================================================================================================================
SUBROUTINE Eval_GVEC(nNodes,xIn,xOut,data_out,phi_axis_edge,chi_axis_edge)
! MODULES
USE MODgvec_Eval_GVEC_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER          :: nNodes
REAL,INTENT( IN) :: xIn(3,nNodes)  !!s,theta,zeta positions for evaluation
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: xOut(3,nNodes)  !! x,y,z coordinates
REAL,INTENT(OUT) :: data_out(9,nNodes)  !! pressure,Bcart(3),chi,phi,Acart(3)
REAL,INTENT(OUT) :: phi_axis_edge(2)
REAL,INTENT(OUT) :: chi_axis_edge(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode
REAL    :: spos,thet,zeta,R,Z,lam
REAL    :: dRds,dZds
REAL    :: dRdthet,dRdzeta
REAL    :: dZdthet,dZdzeta
REAL    :: dlamdthet,dlamdzeta
REAL    :: phi_s,chi_s,phiPrime_s,ChiPrime_s,pres_s,iota_s
REAL    :: sqrtG
REAL    :: xp(3),Bcart(3),Acart(3),q(3),e_s(3),e_thet(3),e_zeta(3)
!===================================================================================================================================
phi_axis_edge(1)= X1_base_r%s%evalDOF_s(1.0e-08, 0,profiles_1d(:,1))
chi_axis_edge(1)= X1_base_r%s%evalDOF_s(1.0e-08, 0,profiles_1d(:,2))
phi_axis_edge(2)= X1_base_r%s%evalDOF_s(1.0, 0,profiles_1d(:,1))
chi_axis_edge(2)= X1_base_r%s%evalDOF_s(1.0, 0,profiles_1d(:,2))
DO iNode=1,nNodes
  xp(1)  =MAX(1.0e-08,xIn(1,iNode))
  xp(2:3) = xIn(2:3,iNode)
  spos=xp(1)
  thet=xp(2)
  zeta=xp(3)
  phi_s  = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,1))
  chi_s  = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,2))
  iota_s = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))
  pres_s = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,4))
  PhiPrime_s = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  ChiPrime_s = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2))

  R          =X1_base_r%evalDOF_x(xp,(/      0,         0/),X1_r)
  dRds       =X1_base_r%evalDOF_x(xp,(/DERIV_S,         0/),X1_r)
  dRdthet    =X1_base_r%evalDOF_x(xp,(/      0,DERIV_THET/),X1_r)
  dRdzeta    =X1_base_r%evalDOF_x(xp,(/      0,DERIV_ZETA/),X1_r)
  Z          =X2_base_r%evalDOF_x(xp,(/      0,         0/),X2_r)
  dZds       =X2_base_r%evalDOF_x(xp,(/DERIV_S,         0/),X2_r)
  dZdthet    =X2_base_r%evalDOF_x(xp,(/      0,DERIV_THET/),X2_r)
  dZdzeta    =X2_base_r%evalDOF_x(xp,(/      0,DERIV_ZETA/),X2_r)
  lam        =LA_base_r%evalDOF_x(xp,(/      0,         0/),LA_r)
  dlamdthet  =LA_base_r%evalDOF_x(xp,(/      0,DERIV_THET/),LA_r)
  dlamdzeta  =LA_base_r%evalDOF_x(xp,(/      0,DERIV_ZETA/),LA_r)

  q=(/R,Z,zeta/)
  e_s    = hmap_r%eval_dxdq(q,(/dRds   ,dZds    ,0.0_wp/))
  e_thet = hmap_r%eval_dxdq(q,(/dRdthet,dZdthet ,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(q,(/dRdzeta ,dZdzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(q)*(dRds*dZdthet -dZds*dRdthet) 

  Bcart(:)=  (  e_thet(:)*(iota_s-dlamdzeta )  &
              + e_zeta(:)*(1.0_wp+dlamdthet) )*(PhiPrime_s/sqrtG)
  Acart(:)=  ( phi_s*e_thet(:)-e_s(:)*(lam*PhiPrime_s)  -chi_s*e_zeta)/sqrtG

  xOut(:,iNode)=hmap_r%eval(q)

  data_out(  1,iNode)=pres_s
  data_out(2:4,iNode)=Bcart(:)
  data_out(  5,iNode)=chi_s
  data_out(  6,iNode)=phi_s
  data_out(7:9,iNode)=Acart(:)
END DO

END SUBROUTINE Eval_GVEC


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeEval_GVEC 
! MODULES
USE MODgvec_Eval_GVEC_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  DEALLOCATE(X1_r)
  DEALLOCATE(X2_r)
  DEALLOCATE(LA_r)
  DEALLOCATE(profiles_1d)
  CALL sgrid_r%free()
!  CALL hmap_r%free()
  CALL X1_base_r%free()
  CALL X2_base_r%free()
  CALL LA_base_r%free()
!  DEALLOCATE(hmap_r)
  DEALLOCATE(X1_base_r)
  DEALLOCATE(X2_base_r)
  DEALLOCATE(LA_base_r)

END SUBROUTINE FinalizeEval_GVEC

END MODULE MODgvec_Eval_GVEC
