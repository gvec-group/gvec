!===================================================================================================================================
! Copyright (C) 2017 - 2019  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **Transform SFL**
!!
!! Transform to Straight-field line angles, PEST / BOOZER 
!!
!===================================================================================================================================
MODULE MODgvec_Transform_SFL
! MODULES
USE MODgvec_Globals, ONLY:wp,abort
IMPLICIT NONE
PRIVATE

INTERFACE BuildTransform_SFL
  MODULE PROCEDURE BuildTransform_SFL
END INTERFACE

INTERFACE FinalizeTransform_SFL
  MODULE PROCEDURE FinalizeTransform_SFL
END INTERFACE


PUBLIC :: BuildTransform_SFL
PUBLIC :: FinalizeTransform_SFL
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Builds X1 and X2 in SFL coordinates
!!
!===================================================================================================================================
SUBROUTINE BuildTransform_SFL(Ns,mn_max,whichSFL)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut
USE MODgvec_base,ONLY:t_base
USE MODgvec_transform_sfl_vars
USE MODgvec_readstate_vars,ONLY: X1_base_r,X1_r,X2_base_r,X2_r,LA_base_r,LA_r
USE MODgvec_get_boozer
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER     ,INTENT(IN) :: Ns                                               !< number of new radial points 
  INTEGER     ,INTENT(IN) :: mn_max(2)                                        !< maximum number for new variables in SFL coordinates
  INTEGER     ,INTENT(IN) :: whichSFL                                         !< either =1: PEST, =2:Boozer 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                    :: fac_nyq              !< for number of integr. points  (=3...4 at least)
  REAL(wp),ALLOCATABLE       :: Gthet(:,:) !< for Boozer
!===================================================================================================================================

!possibility to increase the resolution of the grid for SFL coordinates
CALL sgrid_sfl%init(MAX(Ns-1,X1_base_r%s%grid%nElems),X1_base_r%s%grid%grid_type)
!CALL sgrid_sfl%copy(X1_base_r%s%grid)

whichSFLcoord=whichSFL
SELECT CASE(whichSFLcoord)
CASE(1) !PEST
  fac_nyq=4
  CALL Transform_Angles(LA_base_r,LA_r,X1_base_r,X1_r,mn_max,fac_nyq,sgrid_sfl,"X1",X1sfl_base,X1sfl)
  CALL Transform_Angles(LA_base_r,LA_r,X2_base_r,X2_r,mn_max,fac_nyq,sgrid_sfl,"X2",X2sfl_base,X2sfl)
CASE(2) !BOOZER
  !uses restart data for X1,X2,LA and profiles!
  fac_nyq=4
  CALL Get_Boozer(mn_max,fac_nyq,sgrid_sfl,GZ_base,Gthet,GZ)
  
  fac_nyq=4
  CALL Transform_Angles(GZ_base,Gthet,X1_base_r,X1_r,mn_max,fac_nyq,sgrid_sfl,"X1",X1sfl_base,X1sfl,B_in=GZ)
  CALL Transform_Angles(GZ_base,Gthet,X2_base_r,X2_r,mn_max,fac_nyq,sgrid_sfl,"X2",X2sfl_base,X2sfl,B_in=GZ)
  CALL Transform_Angles(GZ_base,Gthet,GZ_base  ,GZ  ,mn_max,fac_nyq,sgrid_sfl,"GZ",GZsfl_base,GZsfl,B_in=GZ)

  DEALLOCATE(Gthet)

CASE DEFAULT
  SWRITE(UNIT_stdOut,*)'This input for SFL coordinate transform is not valid: whichSFL=',whichSFL
  STOP
END SELECT


END SUBROUTINE BuildTransform_SFL


!===================================================================================================================================
!> Transform a function from VMEC angles q(s,theta,zeta) to new angles q*(s,theta*,zeta*) 
!> by projection onto the modes of the new angles: sigma_mn(theta*,zeta*)
!> using the same representation in s (for now!)
!> Here, new angles are theta*=theta+A(theta,zeta), zeta*=zeta+B(theta,zeta), 
!> with A,B periodic functions and zero average and same base 
!> Note that in this routine, the integral is transformed back to (theta,zeta)
!> q*_mn = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) dtheta* dzeta*
!>       = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) [(1+dA/dtheta)*(1+dB/dzeta)-(dA/dzeta*dB/dzeta)] dtheta dzeta
!!
!===================================================================================================================================
SUBROUTINE Transform_Angles(AB_base_in,A_in,q_base_in,q_in,mn_max,fac_nyq,sgrid_in,q_name,q_base_out,q_out,B_in)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,Progressbar
USE MODgvec_base   ,ONLY: t_base,base_new
USE MODgvec_sGrid  ,ONLY: t_sgrid
USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_Base),INTENT(IN) :: AB_base_in                                    !< basis of A and B
  REAL(wp)     ,INTENT(IN) :: A_in(1:AB_base_in%s%nBase,1:AB_base_in%f%modes) !< coefficients of thet*=thet+A(s,theta,zeta)  
  CLASS(t_Base),INTENT(IN) :: q_base_in                                     !< basis of function f 
  REAL(wp)     ,INTENT(IN) :: q_in(1:q_base_in%s%nBase,1:q_base_in%f%modes) !< coefficients of f 
  INTEGER      ,INTENT(IN) :: mn_max(2)                                     !< maximum number for new variables in SFL coordinates
  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
                                                                            !< n_IP=fac_nyq*max(mn_max_in)
  CLASS(t_sgrid),INTENT(IN),TARGET :: sgrid_in                              !< change grid for q_base_out
  REAL(wp)    ,INTENT(IN),OPTIONAL :: B_in(1:AB_base_in%s%nBase,1:AB_base_in%f%modes) !< coefficients of zeta*=zeta+B(s,theta,zeta)
  CHARACTER(LEN=*),INTENT(IN):: q_name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base),ALLOCATABLE,INTENT(INOUT) :: q_base_out          !< new fourier basis of function q in new angles
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT) :: q_out(:,:)          !< coefficients of q in new angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,i_mn,mn_IP
  INTEGER               :: mn_nyq(2)
  REAL(wp)              :: spos,dthet_dzeta 
  REAL(wp)              :: check(1:3)
  LOGICAL               :: docheck=.TRUE. 
  LOGICAL               :: Bpresent
  REAL(wp)              :: A_s(1:AB_base_in%f%modes) 
  REAL(wp)              :: B_s(1:AB_base_in%f%modes) 
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes) 
  CLASS(t_fBase),ALLOCATABLE :: AB_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE :: q_fbase_nyq
  REAL(wp), ALLOCATABLE :: q_IP(:)       ! q evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: f_IP(:)       ! =q*(1+dlambda/dtheta) evaluated at integration points
  REAL(wp), ALLOCATABLE :: A_IP(:),dAdthet_IP(:),dAdzeta_IP(:) 
  REAL(wp), ALLOCATABLE :: B_IP(:),dBdthet_IP(:),dBdzeta_IP(:) 
  REAL(wp), ALLOCATABLE :: modes_IP(:,:) ! mn modes of q evaluated at theta*,zeta* for all integration points
!===================================================================================================================================
  Bpresent=PRESENT(B_in)
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  IF(mn_max(2).EQ.0) mn_nyq(2)=1 !exception: 2D configuration

  SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6),A,L)')'TRANSFORM '//TRIM(q_name)//' TO NEW ANGLE COORDINATES, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq, ', B_in= ',Bpresent
                     
  !initialize
  
  ! extended base for q in the new angles, and on the new grid
  CALL base_new(q_base_out,  q_base_in%s%deg,        &
                             q_base_in%s%continuity, &
                             sgrid_in, & !q_base_in%s%grid,       &
                             q_base_in%s%degGP,      &
                             mn_max,mn_nyq,          & 
                             q_base_in%f%nfp,        &
                 sin_cos_map(q_base_in%f%sin_cos),   &
                             .FALSE.)!m=n=0 should be always there, because of coordinate transform

  !total number of integration points
  mn_IP = q_base_out%f%mn_IP

  SWRITE(UNIT_StdOut,*)'        ...Init 1 Done'
  !same base for q, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new( q_fbase_nyq,  q_base_in%f%mn_max,  mn_nyq, &
                                q_base_in%f%nfp, &
                    sin_cos_map(q_base_in%f%sin_cos), &
                                q_base_in%f%exclude_mn_zero)

  SWRITE(UNIT_StdOut,*)'        ...Init 2 Done'
  !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new(AB_fbase_nyq,  AB_base_in%f%mn_max,  mn_nyq, &
                                AB_base_in%f%nfp, &
                    sin_cos_map(AB_base_in%f%sin_cos), &
                                AB_base_in%f%exclude_mn_zero)
  
  SWRITE(UNIT_StdOut,*)'        ...Init 3 Done'
  nBase=q_base_out%s%nBase

  ALLOCATE(modes_IP(1:q_base_out%f%modes,1:mn_IP))
  ALLOCATE(f_IP(1:mn_IP),q_IP(1:mn_IP))
  ALLOCATE(A_IP(1:mn_IP),dAdthet_IP(1:mn_IP))
  IF(Bpresent)THEN
    ALLOCATE(dAdzeta_IP(1:mn_IP),B_IP(1:mn_IP),dBdthet_IP(1:mn_IP),dBdzeta_IP(1:mn_IP))
  END IF

  ALLOCATE(q_out(nBase,1:q_base_out%f%modes))

  check(1)=HUGE(1.)
  check(2:3)=0.

  dthet_dzeta  =q_base_out%f%d_thet*q_base_out%f%d_zeta !integration weights

  CALL ProgressBar(0,nBase)!init
  DO is=1,nBase
    spos=q_base_out%s%s_IP(is) !interpolation points for q_in
    !evaluate lambda at spos
    A_s(:)= AB_base_in%s%evalDOF2D_s(spos,AB_base_in%f%modes,   0,A_in(:,:))
    !evaluate lambda at spos
    ! TEST EXACT CASE: A_s=0.
    !evaluate lambda at integration points
    A_IP       = AB_fbase_nyq%evalDOF_IP(         0,A_s(:))
    dAdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET,A_s(:))
    IF(Bpresent)THEN
      B_s(:)     = AB_base_in%s%evalDOF2D_s(spos,AB_base_in%f%modes,   0,B_in(:,:))
      B_IP       = AB_fbase_nyq%evalDOF_IP(         0,B_s(:))
      dBdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET,B_s(:))
      dBdzeta_IP = AB_fbase_nyq%evalDOF_IP(DERIV_ZETA,B_s(:))

      dAdzeta_IP = AB_fbase_nyq%evalDOF_IP(DERIV_ZETA,A_s(:))
    END IF !Bpresent
    !evaluate q_in at spos
    q_in_s(:)= q_base_in%s%evalDOF2D_s(spos,q_base_in%f%modes,   0,q_in(:,:))
    q_IP     = q_fbase_nyq%evalDOF_IP(0,q_in_s(:))
    ! f_IP=q*(1+dlambda/dtheta) 
    IF(.NOT.Bpresent)THEN
      f_IP  = q_IP(:)*(1.0_wp + dAdthet_IP(:))
    ELSE
      f_IP  = q_IP(:)*((1.0_wp + dAdthet_IP(:))*(1.0_wp + dBdzeta_IP(:))-dAdzeta_IP(:)*dBdthet_IP(:))
    END IF
    
    !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
    IF(.NOT.Bpresent)THEN
      DO i_mn=1,mn_IP
        modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_fbase_nyq%x_IP(1,i_mn)+A_IP(i_mn),q_fbase_nyq%x_IP(2,i_mn)/))
      END DO
    ELSE
      DO i_mn=1,mn_IP
        modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_fbase_nyq%x_IP(1,i_mn)+A_IP(i_mn),q_fbase_nyq%x_IP(2,i_mn)+B_IP(i_mn)/))
      END DO
    END IF
    !projection: integrate (sum over mn_IP), includes normalization of base!
    q_out(is,:)=(dthet_dzeta*q_base_out%f%snorm_base(:))*(MATMUL(modes_IP(:,1:mn_IP),f_IP(1:mn_IP)))  

    !CHECK at all IP points:
    IF(doCheck)THEN
      f_IP=MATMUL(q_out(is,:),modes_IP(:,:))
  
      f_IP = ABS(f_IP - q_IP)/SUM(ABS(q_IP))*REAL(mn_IP,wp)
      check(1)=MIN(check(1),MINVAL(f_IP))
      check(2)=MAX(check(2),MAXVAL(f_IP))
      check(3)=check(3)+MAXVAL(f_IP)/REAL(nBase,wp)
!      IF(spos.LT.0.9) CYCLE
!      WRITE(*,*)'     ------  spos= ',spos
!      WRITE(*,*)'check |q_in-q_out|/(surfavg|q_in|) (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)
!      WRITE(*,*)'max,min q_in',MAXVAL(q_IP),MINVAL(q_IP) 
!      WRITE(*,*)'max,min,sum q_out |modes|',MAXVAL(ABS(q_out(is,:))),MINVAL(ABS(q_out(is,:))),SUM(ABS(q_out(is,:)))
    END IF !doCheck

    CALL ProgressBar(is,nBase)
  END DO !is

  !finalize
  CALL AB_fbase_nyq%free()
  CALL q_fbase_nyq%free()
  DEALLOCATE(AB_fbase_nyq,q_fbase_nyq,modes_IP,f_IP)
  DEALLOCATE(A_IP,dAdthet_IP)
  IF(Bpresent) DEALLOCATE(dAdzeta_IP,B_IP,dBdthet_IP,dBdzeta_IP)

  IF(doCheck) WRITE(UNIT_StdOut,'(2A,3E11.3)')'   CHECK ERROR of '//TRIM(q_name)//':', &
                                             ' |q_in-q_out|/(surfavg|q_in|) (min/max/avg)',check(1:3)
  !transform back to corresponding representation of DOF in s
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_out%s%initDOF( q_out(:,iMode) )
  END DO
  SWRITE(UNIT_StdOut,'(A)') '...DONE.'

END SUBROUTINE Transform_Angles

!===================================================================================================================================
!> Transform a function from VMEC angles q(s,theta,zeta) to q*(s,theta*,zeta*) 
!> by projection onto the modes of the new angles: sigma_mn(theta*,zeta*)
!> using the same representation in s (for now!)
!> Here, PEST angles are theta*=theta+lambda(theta,zeta), zeta*=zeta
!> Note that in this routine, an equidistant grid in theta*,zeta* is used and found via Newton-iteration, to eqvaluate q and 
!> integrate directly
!> q*_mn = iint_0^2pi q(theta*,zeta*) sigma_mn(theta*,zeta*) dtheta* dzeta*
!!
!===================================================================================================================================
SUBROUTINE Transform_to_PEST_2(LA_base_in,LA,q_base_in,q_in,mn_max,fac_nyq,sgrid_in,q_base_out,q_out)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,PI
USE MODgvec_base   ,ONLY: t_base,base_new
USE MODgvec_sGrid  ,ONLY: t_sgrid
USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map
USE MODgvec_Newton ,ONLY: NewtonRoot1D_FdF
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_Base),INTENT(IN) :: LA_base_in                                    !< basis of Lambda function
  REAL(wp)     ,INTENT(IN) :: LA(1:LA_base_in%s%nBase,1:LA_base_in%f%modes) !< coefficients of lambda(s,theta,zeta) variable 
  CLASS(t_Base),INTENT(IN) :: q_base_in                                     !< basis of function f 
  REAL(wp)     ,INTENT(IN) :: q_in(1:q_base_in%s%nBase,1:q_base_in%f%modes) !< coefficients of f 
  INTEGER      ,INTENT(IN) :: mn_max(2)                                     !< maximum number for new variables in SFL coordinates
  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
                                                                            !< n_IP=fac_nyq*max(mn_max_in)
  CLASS(t_sgrid),INTENT(IN),TARGET :: sgrid_in                              !< grid for q_base_out
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base),ALLOCATABLE,INTENT(INOUT) :: q_base_out          !< new fourier basis of function q in PEST angles (same sbase)
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT) :: q_out(:,:)          !< fourier coefficients of q in PEST angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,i_mn,mn_IP
  INTEGER               :: mn_nyq(2)
  REAL(wp)              :: spos,dthet_dzeta 
  REAL(wp)              :: xIP(2),theta_star 
  REAL(wp)              :: check(1:3)
  LOGICAL               :: doCheck=.TRUE.
  REAL(wp)              :: LA_s(1:LA_base_in%f%modes) 
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes) 
  CLASS(t_fBase),ALLOCATABLE :: q_fbase_nyq
  REAL(wp), ALLOCATABLE :: q_IP(:)       ! q evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: f_IP(:)       ! =q*(1+dlambda/dtheta) evaluated at integration points
!===================================================================================================================================
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  IF(mn_max(2).EQ.0) mn_nyq(2)=1 !exception: 2D configuration

  WRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'Transform variable to PEST coordinates, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
                     
  !initialize
  
  ! extended base for q in the new angles, and on the new grid
  CALL base_new(q_base_out,  q_base_in%s%deg,        &
                             q_base_in%s%continuity, &
                             sgrid_in, & !q_base_in%s%grid,       &
                             q_base_in%s%degGP,      &
                             mn_max,mn_nyq,          & 
                             q_base_in%f%nfp,        &
                 sin_cos_map(q_base_in%f%sin_cos),   &
                             q_base_in%f%exclude_mn_zero)
!  CALL fbase_new(q_fbase_out,  mn_max,mn_nyq,          & 
!                               q_base_in%f%nfp,        &
!                   sin_cos_map(q_base_in%f%sin_cos),   &
!                               q_base_in%f%exclude_mn_zero)

  !total number of integration points
  mn_IP = q_base_out%f%mn_IP

  !WRITE(UNIT_StdOut,*)'        ...Init 1 Done'
  !same base for q, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new( q_fbase_nyq, q_base_in%f%mn_max,  mn_nyq, &
                                q_base_in%f%nfp, &
                    sin_cos_map(q_base_in%f%sin_cos), &
                                q_base_in%f%exclude_mn_zero)

  !WRITE(UNIT_StdOut,*)'        ...Init 2 Done'

  nBase=q_base_out%s%nBase

  ALLOCATE(f_IP(1:mn_IP),q_IP(1:mn_IP))
  ALLOCATE(q_out(nBase,1:q_base_out%f%modes))

  check(1)=HUGE(1.)
  check(2:3)=0.

  dthet_dzeta  =q_base_out%f%d_thet*q_base_out%f%d_zeta !integration weights

  DO is=1,nBase
    IF(MOD(is,MAX(1,nBase/100)).EQ.0) THEN
      SWRITE(UNIT_stdOut,'(8X,I4,A4,I4,A13,A1)',ADVANCE='NO')is, ' of ',nBase,' evaluated...',ACHAR(13)
    END IF
    spos=q_base_out%s%s_IP(is) !interpolation points for q_in
    !evaluate lambda at spos
    LA_s(:)= LA_base_in%s%evalDOF2D_s(spos,LA_base_in%f%modes,     0,LA(:,:))
    ! TEST EXACT CASE: LA_s=0.
    !evaluate q_in at spos
    q_in_s(:)= q_base_in%s%evalDOF2D_s(spos,q_base_in%f%modes,     0,q_in(:,:))
    !find theta(theta*) which fulfills theta_star=theta+lambda(theta,zeta)
    DO i_mn=1,mn_IP
      xIP(2)     = q_base_out%f%x_IP(2,i_mn)
      theta_star = q_base_out%f%x_IP(1,i_mn)
      xIP(1)     = NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1_wp*PI,theta_star   , theta_star,FRdFR)
      q_IP(i_mn) =q_base_in%f%evalDOF_x(xIP, 0,q_in_s )
    END DO
    !projection: integrate (sum over mn_IP), includes normalization of base!
    q_out(is,:)=(dthet_dzeta*q_base_out%f%snorm_base(:))*(MATMUL(q_IP,q_base_out%f%base_IP(1:mn_IP,:)))  


    !CHECK at all IP (theta,zeta) points:

    IF(doCheck)THEN
      f_IP=MATMUL(q_base_out%f%base_IP,q_out(is,:))
      
      f_IP = ABS(f_IP - q_IP)
      check(1)=MIN(check(1),MINVAL(f_IP))
      check(2)=MAX(check(2),MAXVAL(f_IP))
      check(3)=check(3)+MAXVAL(f_IP)/REAL(nBase,wp)
      !WRITE(*,*)'     ------  spos= ',spos
      !WRITE(*,*)'check |q_in-q_out| (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)
    END IF !doCheck

  END DO !is

  !finalize
  CALL q_fbase_nyq%free()
  DEALLOCATE(q_fbase_nyq,f_IP)

  IF(doCheck) WRITE(UNIT_StdOut,'(A,3E11.3)')'   check |q_in-q_out| (min/max/avg)',check(1:3)
  !transform back to corresponding representation of DOF in s
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_out%s%initDOF( q_out(:,iMode) )
  END DO

!for iteration on theta^*
CONTAINS 

  FUNCTION FRdFR(theta_iter)
    !uses current zeta=XIP(2) where newton is called, and LA_s from subroutine above
    IMPLICIT NONE
    REAL(wp) :: theta_iter
    REAL(wp) :: FRdFR(2) !output
    !--------------------------------------------------- 
    FRdFR(1)=theta_iter+LA_base_in%f%evalDOF_x((/theta_iter,xIP(2)/),0,LA_s)  !theta_iter+lambda
    FRdFR(2)=1.0_wp+LA_base_in%f%evalDOF_x((/theta_iter,xIP(2)/),DERIV_THET,LA_s) !1+dlambda/dtheta
  END FUNCTION FRdFR

END SUBROUTINE Transform_to_PEST_2

!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE FinalizeTransform_SFL
! MODULES
USE MODgvec_Transform_SFL_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL sgrid_sfl%free()
  CALL X1sfl_base%free()
  CALL X2sfl_base%free()
  CALL GZsfl_base%free()
  CALL GZ_base%free()
  SDEALLOCATE(X1sfl_base)
  SDEALLOCATE(X2sfl_base)
  SDEALLOCATE(GZsfl_base)
  SDEALLOCATE(GZ_base)
  SDEALLOCATE(X1sfl)
  SDEALLOCATE(X2sfl)
  SDEALLOCATE(GZsfl)
  SDEALLOCATE(GZ)


END SUBROUTINE FinalizeTransform_SFL

END MODULE MODgvec_Transform_SFL
