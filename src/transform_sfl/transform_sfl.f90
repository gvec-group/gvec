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
SUBROUTINE BuildTransform_SFL(X1_base_in,X1_in,X2_base_in,X2_in,LA_base_in,LA_in,mn_max,whichSFL)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut
USE MODgvec_base,ONLY:t_base
USE MODgvec_transform_sfl_vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  TYPE(t_Base),INTENT(IN) :: X1_base_in                                    !< basis of X1 function, in GVEC coordinates
  TYPE(t_Base),INTENT(IN) :: X2_base_in                                    !< basis of X2 function, in GVEC coordinates
  TYPE(t_Base),INTENT(IN) :: LA_base_in                                    !< basis of Lambda function, in GVEC coordinates
  REAL(wp)    ,INTENT(IN) :: X1_in(1:X1_base_in%s%nBase,1:X1_base_in%f%modes) !< coefficients 
  REAL(wp)    ,INTENT(IN) :: X2_in(1:X2_base_in%s%nBase,1:X2_base_in%f%modes) !< coefficients 
  REAL(wp)    ,INTENT(IN) :: LA_in(1:LA_base_in%s%nBase,1:LA_base_in%f%modes) !< coefficients 
  INTEGER     ,INTENT(IN) :: mn_max(2)                                        !< maximum number for new variables in SFL coordinates
  INTEGER     ,INTENT(IN) :: whichSFL                                         !< either =1: PEST, =2:Boozer 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                    :: fac_nyq              !< for number of integr. points  (=3...4 at least)
!===================================================================================================================================
whichSFLcoord=whichSFL
SELECT CASE(whichSFLcoord)
CASE(1) !PEST
  fac_nyq=4
  CALL Transform_to_PEST(LA_base_in,LA_in,X1_base_in,X1_in,mn_max,fac_nyq,X1sfl_base,X1sfl)
  CALL Transform_to_PEST(LA_base_in,LA_in,X2_base_in,X2_in,mn_max,fac_nyq,X2sfl_base,X2sfl)
CASE(2) !BOOZER
  SWRITE(UNIT_stdOut,*)'BOOZER SFL coordinate input is not yet implemented'
  STOP
CASE DEFAULT
  SWRITE(UNIT_stdOut,*)'This input for SFL coordinate transform is not valid: whichSFL=',whichSFL
  STOP
END SELECT


END SUBROUTINE BuildTransform_SFL


!===================================================================================================================================
!> Transform a function from VMEC angles q(s,theta,zeta) to q*(s,theta*,zeta*) 
!> by projection onto the modes of the new angles: sigma_mn(theta*,zeta*)
!> using the same representation in s (for now!)
!> Here, PEST angles are theta*=theta+lambda(theta,zeta), zeta*=zeta
!> Note that in this routine, the integral is transformed back to (theta,zeta)
!> q*_mn = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) dtheta* dzeta*
!>       = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) (1+dlambda/dtheta) dtheta dzeta
!!
!===================================================================================================================================
SUBROUTINE Transform_to_PEST(LA_base_in,LA,q_base_in,q_in,mn_max,fac_nyq,q_base_out,q_out)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut
USE MODgvec_base,ONLY:t_base,base_new
USE MODgvec_fbase,ONLY:t_fbase,fbase_new,sin_cos_map
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base),ALLOCATABLE,INTENT(INOUT) :: q_base_out          !< new fourier basis of function q in PEST angles (same sbase)
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT) :: q_out(:,:)          !< fourier coefficients of q in PEST angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,i_mn,mn_IP
  INTEGER               :: mn_nyq(2)
  REAL(wp)              :: spos,dthet_dzeta 
  REAL(wp)              :: check(1:3)
  LOGICAL               :: docheck=.TRUE. 
  REAL(wp)              :: LA_s(1:LA_base_in%f%modes) 
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes) 
  CLASS(t_fBase),ALLOCATABLE :: LA_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE :: q_fbase_nyq
  REAL(wp), ALLOCATABLE :: q_IP(:)       ! q evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: f_IP(:)       ! =q*(1+dlambda/dtheta) evaluated at integration points
  REAL(wp), ALLOCATABLE :: LA_IP(:)      ! lambda evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: modes_IP(:,:) ! mn modes of q evaluated at theta*,zeta* for all integration points
  REAL(wp), ALLOCATABLE :: q_gIP(:,:)    ! interpolated values of q_out
!===================================================================================================================================
!keep the same base in radial direction
  nBase=q_base_in%s%nBase
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  IF(mn_max(2).EQ.0) mn_nyq(2)=1 !exception: 2D configuration

  WRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'Transform variable to PEST coordinates, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
                     
  !initialize
  
  ! extended base for q in the new angles

  CALL base_new(q_base_out,  q_base_in%s%deg,        &
                             q_base_in%s%continuity, &
                             q_base_in%s%grid,       &
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
  CALL fbase_new( q_fbase_nyq,  q_base_in%f%mn_max,  mn_nyq, &
                                q_base_in%f%nfp, &
                    sin_cos_map(q_base_in%f%sin_cos), &
                                q_base_in%f%exclude_mn_zero)

  !WRITE(UNIT_StdOut,*)'        ...Init 2 Done'
  !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new(LA_fbase_nyq, LA_base_in%f%mn_max,  mn_nyq, &
                                LA_base_in%f%nfp, &
                    sin_cos_map(LA_base_in%f%sin_cos), &
                                LA_base_in%f%exclude_mn_zero)
  
  !WRITE(UNIT_StdOut,*)'        ...Init 3 Done'

  ALLOCATE(q_gIP(1:q_base_out%f%modes,1:nBase))
  ALLOCATE(modes_IP(1:q_base_out%f%modes,1:mn_IP))
  ALLOCATE(f_IP(1:mn_IP),LA_IP(1:mn_IP),q_IP(1:mn_IP))

  check(1)=HUGE(1.)
  check(2:3)=0.

  dthet_dzeta  =q_base_out%f%d_thet*q_base_out%f%d_zeta !integration weights

  DO is=nBase,1,-1
    IF(MOD(is,MAX(1,nBase/100)).EQ.0) THEN
      SWRITE(UNIT_stdOut,'(8X,I4,A4,I4,A13,A1)',ADVANCE='NO')is, ' of ',nBase,' evaluated...',ACHAR(13)
    END IF
    spos=q_base_in%s%s_IP(is) !interpolation points for q_in
    !evaluate lambda at spos
    DO iMode=1,LA_base_in%f%modes
      LA_s(  iMode)= LA_base_in%s%evalDOF_s(spos,      0,LA(:,iMode))
    END DO
    !IF(spos.LT.1.0e-8_wp) LA_s=0.
    ! TEST EXACT CASE: LA_s=0.
    !evaluate q_in at spos
    DO iMode=1,q_base_in%f%modes
      q_in_s(iMode)= q_base_in%s%evalDOF_s(spos,      0,q_in(:,iMode))
    END DO
    !evaluate lambda at integration points
    LA_IP = LA_fbase_nyq%evalDOF_IP(0,LA_s(:))
    q_IP  = q_fbase_nyq%evalDOF_IP(0,q_in_s(:))
    ! f_IP=q*(1+dlambda/dtheta) 
    f_IP  = q_IP(:)*(1.0_wp + LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:)))
    
    !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
    DO i_mn=1,mn_IP
      modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_fbase_nyq%x_IP(1,i_mn)+LA_IP(i_mn),q_fbase_nyq%x_IP(2,i_mn)/))
    END DO
    !projection: integrate (sum over mn_IP), includes normalization of base!
    q_gIP(:,is)=(dthet_dzeta*q_base_out%f%snorm_base(:))*(MATMUL(modes_IP(:,1:mn_IP),f_IP(1:mn_IP)))  

    !CHECK at all IP points:
    IF(doCheck)THEN
      f_IP=MATMUL(q_gIP(:,is),modes_IP(:,:))
  
      f_IP = ABS(f_IP - q_IP)
      check(1)=MIN(check(1),MINVAL(f_IP))
      check(2)=MAX(check(2),MAXVAL(f_IP))
      check(3)=check(3)+MAXVAL(f_IP)/REAL(nBase,wp)
      !WRITE(*,*)'     ------  spos= ',spos
      !WRITE(*,*)'check |q_in-q_out| (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)
    END IF !doCheck

  END DO !is

  !finalize
  CALL LA_fbase_nyq%free()
  CALL q_fbase_nyq%free()
  DEALLOCATE(LA_fbase_nyq,q_fbase_nyq,modes_IP,f_IP,LA_IP)

  IF(doCheck) WRITE(UNIT_StdOut,'(A,3E11.3)')'   check |q_in-q_out| (min/max/avg)',check(1:3)
  !transform back to corresponding representation of DOF in s
  ALLOCATE(q_out(nBase,1:q_base_out%f%modes))
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_in%s%initDOF( q_gIP(iMode,:) )
  END DO
  DEALLOCATE(q_gIP)

END SUBROUTINE Transform_to_PEST

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
SUBROUTINE Transform_to_PEST_2(LA_base_in,LA,q_base_in,q_in,mn_max,fac_nyq,q_base_out,q_out)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,PI
USE MODgvec_base,ONLY:t_base,base_new
USE MODgvec_fbase,ONLY:t_fbase,fbase_new,sin_cos_map
USE MODgvec_Newton,         ONLY: NewtonRoot1D_FdF
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
  REAL(wp), ALLOCATABLE :: q_gIP(:,:)    ! interpolated values of q_out
!===================================================================================================================================
  !keep the same base in radial direction
  nBase=q_base_in%s%nBase
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  IF(mn_max(2).EQ.0) mn_nyq(2)=1 !exception: 2D configuration

  WRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'Transform variable to PEST coordinates, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
                     
  !initialize
  
  ! extended base for q in the new angles
  CALL base_new(q_base_out,  q_base_in%s%deg,        &
                             q_base_in%s%continuity, &
                             q_base_in%s%grid,       &
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

  ALLOCATE(q_gIP(1:q_base_out%f%modes,1:nBase))
  ALLOCATE(f_IP(1:mn_IP),q_IP(1:mn_IP))

  check(1)=HUGE(1.)
  check(2:3)=0.

  dthet_dzeta  =q_base_out%f%d_thet*q_base_out%f%d_zeta !integration weights

  DO is=1,nBase
    IF(MOD(is,MAX(1,nBase/100)).EQ.0) THEN
      SWRITE(UNIT_stdOut,'(8X,I4,A4,I4,A13,A1)',ADVANCE='NO')is, ' of ',nBase,' evaluated...',ACHAR(13)
    END IF
    spos=q_base_in%s%s_IP(is) !interpolation points for q_in
    !evaluate lambda at spos
    DO iMode=1,LA_base_in%f%modes
      LA_s(  iMode)= LA_base_in%s%evalDOF_s(spos,      0,LA(:,iMode))
    END DO
    ! TEST EXACT CASE: LA_s=0.
    !evaluate q_in at spos
    DO iMode=1,q_base_in%f%modes
      q_in_s(iMode)= q_base_in%s%evalDOF_s(spos,      0,q_in(:,iMode))
    END DO
    !find theta(theta*) which fulfills theta_star=theta+lambda(theta,zeta)
    DO i_mn=1,mn_IP
      xIP(2)     = q_base_out%f%x_IP(2,i_mn)
      theta_star = q_base_out%f%x_IP(1,i_mn)
      xIP(1)     = NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1_wp*PI,theta_star   , theta_star,FRdFR)
      q_IP(i_mn) =q_base_in%f%evalDOF_x(xIP, 0,q_in_s )
    END DO
    !projection: integrate (sum over mn_IP), includes normalization of base!
    q_gIP(:,is)=(dthet_dzeta*q_base_out%f%snorm_base(:))*(MATMUL(q_IP,q_base_out%f%base_IP(1:mn_IP,:)))  


    !CHECK at all IP (theta,zeta) points:

    IF(doCheck)THEN
      f_IP=MATMUL(q_base_out%f%base_IP,q_gIP(:,is))
      
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
  ALLOCATE(q_out(nBase,1:q_base_out%f%modes))
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_out%s%initDOF( q_gIP(iMode,:) )
  END DO
  DEALLOCATE(q_gIP)

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
  CALL X1sfl_base%free()
  CALL X2sfl_base%free()
  CALL GZsfl_base%free()
  SDEALLOCATE(X1sfl)
  SDEALLOCATE(X2sfl)
  SDEALLOCATE(GZsfl)


END SUBROUTINE FinalizeTransform_SFL

END MODULE MODgvec_Transform_SFL
