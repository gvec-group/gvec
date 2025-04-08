!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
USE MODgvec_Globals, ONLY:wp,UNIT_stdOut,fmt_sep
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
SUBROUTINE Init_gvec_to_hopr(fileName,SFLcoord_in,factorSFL_in)
! MODULES
USE MODgvec_gvec_to_hopr_Vars
USE MODgvec_readState         ,ONLY: ReadState,eval_phiPrime_r,eval_iota_r
USE MODgvec_ReadState_vars    ,ONLY: hmap_r,X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r
USE MODgvec_transform_sfl     ,ONLY: transform_sfl_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileName !< name of GVEC file
INTEGER,OPTIONAL,INTENT(IN) :: SFLcoord_in  !< which SFL coordinate system: =0: none(GVEC), =1: PEST =2: BOOZER
INTEGER,OPTIONAL,INTENT(IN) :: factorSFL_in !< increase number of fourier modes to represent transformed coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: mn_max(2),factorSFL
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT GVEC TO HOPR ...'

  CALL ReadState(fileName)


  IF(PRESENT(SFLcoord_in))THEN
    SFLcoord=SFLcoord_in
  ELSE
    SFLcoord=0
  END IF
  IF((SFLcoord.LT.0).OR.(SFLcoord.GT.2))THEN
    STOP "init_gvec_to_hopr: SFLcoord must be 0,1,2"
  END IF

  IF(SFLcoord.NE.0)THEN
    IF(PRESENT(factorSFL_in))THEN
      factorSFL=factorSFL_in
    ELSE
      factorSFL=4 !default
    END IF
    mn_max(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
    mn_max(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))
    mn_max=mn_max*factorSFL !*SFLfactor on modes
    CALL transform_sfl_new(trafoSFL,mn_max,SFLcoord,X1_base_r%s%deg,X1_base_r%s%continuity, &
                           X1_base_r%s%degGP,X1_base_r%s%grid ,hmap_r,X1_base_r,X2_base_r,LA_base_r,eval_phiPrime_r,eval_iota_r)
    CALL trafoSFL%buildTransform(X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r)
  END IF


  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE Init_gvec_to_hopr

SUBROUTINE gvec_to_hopr(nNodes,xIn,xOut,data_out,phi_axis_edge,chi_axis_edge)
! MODULES
USE MODgvec_Globals, ONLY: CROSS
USE MODgvec_gvec_to_hopr_vars
USE MODgvec_ReadState_Vars    ,ONLY: profiles_1d,sbase_prof !for profiles
USE MODgvec_ReadState_vars    ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_ReadState_vars    ,ONLY: LA_r,X1_r,X2_r
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER ,INTENT(IN) :: nNodes
REAL(wp),INTENT(IN) :: xIn(3,nNodes)  !!s=sqrt(psi_norm),theta,zeta positions for evaluation, psi_norm is normalized toroidal flux
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: xOut(3,nNodes)  !! x,y,z cartesian coordinates
REAL(wp),INTENT(OUT) :: data_out(9,nNodes)  !! pressure,Bcart(3),chi,phi,Acart(3)
REAL(wp),INTENT(OUT) :: phi_axis_edge(2)
REAL(wp),INTENT(OUT) :: chi_axis_edge(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! local Variables
!===================================================================================================================================
phi_axis_edge(1)= sbase_prof%evalDOF_s(1.0e-08_wp, 0,profiles_1d(:,1))
chi_axis_edge(1)= sbase_prof%evalDOF_s(1.0e-08_wp, 0,profiles_1d(:,2))
phi_axis_edge(2)= sbase_prof%evalDOF_s(1.0_wp, 0,profiles_1d(:,1))
chi_axis_edge(2)= sbase_prof%evalDOF_s(1.0_wp, 0,profiles_1d(:,2))

SELECT CASE(SFLcoord)
CASE(0) !DEFAULT, NO SFL coordinate
  CALL gvec_to_hopr_SFL(nNodes,xIn,X1_base_r,X1_r,X2_base_r,X2_r,LA_base_r,LA_r,xOut,data_out)
CASE(1) !PEST coordinates
  CALL gvec_to_hopr_SFL(nNodes,xIn,trafoSFL%X1sfl_base,trafoSFL%X1sfl,trafoSFL%X2sfl_base,trafoSFL%X2sfl,LA_base_r,LA_r,&
                        xOut,data_out) !LA only used as placeholder
CASE(2) !BOOZER coordinates
  CALL gvec_to_hopr_SFL(nNodes,xIn,trafoSFL%X1sfl_base,trafoSFL%X1sfl,trafoSFL%X2sfl_base,trafoSFL%X2sfl, &
                        trafoSFL%GZsfl_base,trafoSFL%GZsfl,xOut,data_out)
END SELECT

END SUBROUTINE gvec_to_hopr

!===================================================================================================================================
!> Evaluate gvec state at a list of s,theta,zeta positions
!!
!===================================================================================================================================
SUBROUTINE gvec_to_hopr_SFL(nNodes,xIn,X1_base_in,X1_in,X2_base_in,X2_in,LG_base_in,LG_in,xOut,data_out)
! MODULES
USE MODgvec_Globals, ONLY: CROSS,ProgressBar
USE MODgvec_gvec_to_hopr_vars
USE MODgvec_ReadState_Vars, ONLY: profiles_1d,hmap_r,sbase_prof !for profiles
USE MODgvec_Base,   ONLY: t_base
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER          :: nNodes
REAL(wp),INTENT( IN) :: xIn(3,nNodes)  !!s=sqrt(psi_norm),theta,zeta positions for evaluation, psi_norm is normalized toroidal flux
CLASS(t_base) ,INTENT(IN) :: X1_base_in,X2_base_in,LG_base_in
REAL(wp)      ,INTENT(IN) :: X1_in(1:X1_base_in%s%nBase,1:X1_base_in%f%modes)
REAL(wp)      ,INTENT(IN) :: X2_in(1:X2_base_in%s%nBase,1:X2_base_in%f%modes)
REAL(wp)      ,INTENT(IN) :: LG_in(1:LG_base_in%s%nBase,1:LG_base_in%f%modes) ! is either LA if SFLcoord=0/1 or G if SFLcoord=2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: xOut(3,nNodes)  !! x,y,z cartesian coordinates
REAL(wp),INTENT(OUT) :: data_out(9,nNodes)  !! pressure,Bcart(3),chi,phi,Acart(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode
REAL    :: spos,xp(2),X1_int,X2_int,LA_int
REAL(wp),DIMENSION(1:X1_base_in%f%modes) :: X1_s,dX1ds_s
REAL(wp),DIMENSION(1:X2_base_in%f%modes) :: X2_s,dX2ds_s
REAL(wp),DIMENSION(1:LG_base_in%f%modes) :: LG_s,dGds_s
REAL(wp):: dX1ds   ,dX2ds
REAL(wp):: dX1dthet,dX2dthet
REAL(wp):: dX1dzeta,dX2dzeta
REAL(wp):: dLAdthet,dLAdzeta
REAL(wp):: G_int,dGds,dGdthet,dGdzeta
REAL(wp):: phi_int,chi_int,pres_int,iota_int
REAL(wp):: dPhids_int,dChids_int         !prime refers to d/ds , where s=sqrt(phi_norm)
REAL(wp):: sqrtG
REAL(wp):: Bcart(3),Acart(3),qvec(3)
REAL(wp):: e_s(3),e_thet(3),e_zeta(3)
REAL(wp):: grad_s(3),grad_thet(3),grad_zeta(3)
!===================================================================================================================================
LA_int  = 0.0_wp !only changed for SFLcoord=0
dLAdthet= 0.0_wp !only changed for SFLcoord=0
dLAdzeta= 0.0_wp !only changed for SFLcoord=0
G_int   = 0.0_wp !only changed for SFLcoords=2
dGds    = 0.0_wp !only changed for SFLcoords=2
dGdthet = 0.0_wp !only changed for SFLcoords=2
dGdzeta = 0.0_wp !only changed for SFLcoords=2

!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)                                   &
!$OMP   FIRSTPRIVATE(LA_int,dLAdthet,dLAdzeta,G_int,dGds,dGdthet,dGdzeta)         &
!$OMP   PRIVATE(iNode,spos,X1_int,dX1ds,dX1dthet,dX1dzeta,X2_int,dX2ds,dX2dthet,dX2dzeta,      &
!$OMP           phi_int,chi_int,iota_int,pres_int,dPhids_int,dChids_int,X1_s,dX1ds_s,X2_s,dX2ds_s,LG_s,dGds_s,          &
!$OMP           xp,qvec,e_s,e_thet,e_zeta,sqrtG,Bcart,Acart,grad_s,grad_thet,grad_zeta)           &
!$OMP   SHARED(nNodes,SFLcoord,xIn,sbase_prof,profiles_1d,hmap_r,X1_base_in,X2_base_in,LG_base_in,X1_in,X2_in,LG_in,xOut,data_out)
DO iNode=1,nNodes
  spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,xIn(1,iNode)))
!  thet=xIn(2,iNode)
!  zeta=xIn(3,iNode)
  xp(1:2)=xIn(2:3,iNode)

  phi_int      = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,1))
  chi_int      = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,2))
  iota_int     = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,3))
  pres_int     = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,4))
  dPhids_int = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  !dChids_int = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2))
  dChids_int  = dPhids_int*iota_int

  X1_s(   :) = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,      0,X1_in(:,:))
  dX1ds_s(:) = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,DERIV_S,X1_in(:,:))
  X2_s(   :) = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,      0,X2_in(:,:))
  dX2ds_s(:) = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,DERIV_S,X2_in(:,:))

  IF(SFLcoord.EQ.0)THEN !GVEC coordinates
    LG_s(  :) = LG_base_in%s%evalDOF2D_s(spos,LG_base_in%f%modes,       0,LG_in(:,:))
  ELSEIF(SFLcoord.EQ.2)THEN !BOOZER
    LG_s(  :) = LG_base_in%s%evalDOF2D_s(spos,LG_base_in%f%modes,       0,LG_in(:,:))
    dGds_s(:) = LG_base_in%s%evalDOF2D_s(spos,LG_base_in%f%modes, DERIV_S,LG_in(:,:))
  END IF

  X1_int     = X1_base_in%f%evalDOF_x(xp,         0,X1_s)
  dX1ds      = X1_base_in%f%evalDOF_x(xp,         0,dX1ds_s)
  dX1dthet   = X1_base_in%f%evalDOF_x(xp,DERIV_THET,X1_s)
  dX1dzeta   = X1_base_in%f%evalDOF_x(xp,DERIV_ZETA,X1_s)
  X2_int     = X2_base_in%f%evalDOF_x(xp,         0,X2_s)
  dX2ds      = X2_base_in%f%evalDOF_x(xp,         0,dX2ds_s)
  dX2dthet   = X2_base_in%f%evalDOF_x(xp,DERIV_THET,X2_s)
  dX2dzeta   = X2_base_in%f%evalDOF_x(xp,DERIV_ZETA,X2_s)
  IF(SFLcoord.EQ.0)THEN !GVEC coordinates (else=0)
    LA_int     = LG_base_in%f%evalDOF_x(xp,         0,LG_s)
    dLAdthet   = LG_base_in%f%evalDOF_x(xp,DERIV_THET,LG_s)
    dLAdzeta   = LG_base_in%f%evalDOF_x(xp,DERIV_ZETA,LG_s)
  END IF

  IF(SFLcoord.EQ.2)THEN !BOOZER coordinates (else=0)
    G_int   = LG_base_in%f%evalDOF_x(xp,         0, LG_s)
    dGds    = LG_base_in%f%evalDOF_x(xp,         0, dGds_s)
    dGdthet = LG_base_in%f%evalDOF_x(xp,DERIV_THET, LG_s)
    dGdzeta = LG_base_in%f%evalDOF_x(xp,DERIV_ZETA, LG_s)
  END IF
  qvec=(/X1_int,X2_int,xp(2)-G_int/) !(X1,X2,zeta)
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,0.0_wp -dGds   /))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp -dGdthet/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp -dGdzeta/))
  sqrtG  = SUM(e_s*(CROSS(e_thet,e_zeta)))

  Bcart(:)=  (  e_thet(:)*(iota_int-dLAdzeta )  &
              + e_zeta(:)*(1.0_wp+dLAdthet) )*(dPhids_int/sqrtG)
  grad_s    = CROSS(e_thet,e_zeta) !/sqrtG
  grad_thet = CROSS(e_zeta,e_s   ) !/sqrtG
  grad_zeta = CROSS(e_s   ,e_thet) !/sqrtG

  Acart(:)=  ( phi_int*grad_thet(:)-(LA_int*dPhids_int)*grad_s(:)  -chi_int*grad_zeta)/sqrtG

  xOut(:,iNode)=hmap_r%eval(qvec)

  data_out(  1,iNode)=pres_int
  data_out(2:4,iNode)=Bcart(:)
  data_out(  5,iNode)=chi_int
  data_out(  6,iNode)=phi_int
  data_out(7:9,iNode)=Acart(:)
END DO
!$OMP END PARALLEL DO

END SUBROUTINE gvec_to_hopr_SFL


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE Finalize_gvec_to_hopr
! MODULES
USE MODgvec_gvec_to_hopr_Vars
USE MODgvec_ReadState,ONLY:Finalize_ReadState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'FINALIZE GVEC_TO_HOPR ...'
  CALL Finalize_ReadState()
  IF(SFLcoord.NE.0) THEN
    CALL trafoSFL%free()
    DEALLOCATE(trafoSFL)
  END IF
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE Finalize_gvec_to_hopr

END MODULE MODgvec_gvec_to_hopr
