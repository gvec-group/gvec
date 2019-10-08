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

INTERFACE test_sfl
  MODULE PROCEDURE test_sfl
END INTERFACE

INTERFACE Transform_To_PEST
  MODULE PROCEDURE Transform_To_PEST
END INTERFACE


PUBLIC :: test_sfl
PUBLIC :: Transform_To_PEST
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!>
!!
!===================================================================================================================================
SUBROUTINE test_sfl(LA_base_in,LA,q_base_in,q_in)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut
USE MODgvec_base,ONLY:t_base
USE MODgvec_fbase,ONLY:t_fbase,fbase_new,sin_cos_map
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  TYPE(t_Base),INTENT(IN) :: LA_base_in                                    !< basis of Lambda function
  REAL(wp)    ,INTENT(IN) :: LA(1:LA_base_in%s%nBase,1:LA_base_in%f%modes) !< coefficients of lambda(s,theta,zeta) variable 
  TYPE(t_Base),INTENT(IN) :: q_base_in                                     !< basis of function f 
  REAL(wp)    ,INTENT(IN) :: q_in(1:q_base_in%s%nBase,1:q_base_in%f%modes) !< coefficients of f 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CLASS(t_fBase),ALLOCATABLE :: q_fbase_out          !< only fourier basis of function q in PEST angles
  REAL(wp)     ,ALLOCATABLE  :: q_out(:,:)           !< fourier coefficients of q in PEST angles 
!===================================================================================================================================
WRITE(*,*)'TESTING transform_to_PEST:'

CALL Transform_to_PEST(LA_base_in,LA,q_base_in,q_in,3,3,q_fbase_out,q_out)


CALL q_fbase_out%free()
DEALLOCATE(q_fbase_out)
DEALLOCATE(q_out)

END SUBROUTINE test_sfl


!===================================================================================================================================
!> Transform a function from VMEC angles q(s,theta,zeta) to q*(s,theta*,zeta*) 
!> by projection onto the modes of the new angles: sigma_mn(theta*,zeta*)
!> using the same representation in s (for now!)
!> Here, PEST angles are theta*=theta+lambda(theta,zeta), zeta*=zeta
!> q*_mn = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) dtheta* dzeta*
!>       = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) (1+dlambda/dtheta) dtheta dzeta
!!
!===================================================================================================================================
SUBROUTINE Transform_to_PEST(LA_base_in,LA,q_base_in,q_in,fac_modes,fac_nyq,q_fbase_out,q_out)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,TWOPI
USE MODgvec_base,ONLY:t_base
USE MODgvec_fbase,ONLY:t_fbase,fbase_new,sin_cos_map
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_Base),INTENT(IN) :: LA_base_in                                    !< basis of Lambda function
  REAL(wp)     ,INTENT(IN) :: LA(1:LA_base_in%s%nBase,1:LA_base_in%f%modes) !< coefficients of lambda(s,theta,zeta) variable 
  CLASS(t_Base),INTENT(IN) :: q_base_in                                     !< basis of function f 
  REAL(wp)     ,INTENT(IN) :: q_in(1:q_base_in%s%nBase,1:q_base_in%f%modes) !< coefficients of f 
  INTEGER      ,INTENT(IN) :: fac_modes                                     !< factor on number of modes for output (=2...3 )
  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
                                                                            !< n_IP=fac_nyq*fac_modes*max(mn_max_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase),ALLOCATABLE,INTENT(INOUT) :: q_fbase_out          !< new fourier basis of function q in PEST angles (same sbase)
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT)  :: q_out(:,:)           !< fourier coefficients of q in PEST angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,i_mn,mn_IP
  INTEGER               :: mn_max(2), mn_nyq(2)
  REAL(wp)              :: spos,dthet_dzeta 
  REAL(wp)              :: check(1:3)
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
  ! double the modes, keep same fourier basis properties (symmetry, zero mode)
  mn_max(1:2)=fac_modes*q_base_in%f%mn_max
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  IF(mn_max(2).EQ.0) mn_nyq(2)=1 !exception: 2D configuration

  WRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'Transform variable to PEST coordinates, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
                     
  !initialize
  
  ! extended base for q in the new angles
  CALL fbase_new(q_fbase_out,  mn_max,mn_nyq,          & 
                               q_base_in%f%nfp,        &
                   sin_cos_map(q_base_in%f%sin_cos),   &
                               q_base_in%f%exclude_mn_zero)

  !total number of integration points
  mn_IP = q_fbase_out%mn_IP

  WRITE(UNIT_StdOut,*)'        ...Init 1 Done'
  !same base for q, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new( q_fbase_nyq, q_base_in%f%mn_max,  mn_nyq, &
                                q_base_in%f%nfp, &
                    sin_cos_map(q_base_in%f%sin_cos), &
                                q_base_in%f%exclude_mn_zero)

  WRITE(UNIT_StdOut,*)'        ...Init 2 Done'
  !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new(LA_fbase_nyq, LA_base_in%f%mn_max,  mn_nyq, &
                                LA_base_in%f%nfp, &
                    sin_cos_map(LA_base_in%f%sin_cos), &
                                LA_base_in%f%exclude_mn_zero)
  
  WRITE(UNIT_StdOut,*)'        ...Init 3 Done'

  ALLOCATE(q_gIP(1:q_fbase_out%modes,1:nBase))
  ALLOCATE(modes_IP(1:q_fbase_out%modes,1:mn_IP))
  ALLOCATE(f_IP(1:mn_IP),LA_IP(1:mn_IP))

  check(1)=HUGE(1.)
  check(2:3)=0.

  dthet_dzeta  =q_fbase_out%d_thet*q_fbase_out%d_zeta !integration weights

  DO is=nBase,1,-1
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
    !evaluate lambda at integration points
    LA_IP = LA_fbase_nyq%evalDOF_IP(0,LA_s(:))
    q_IP  = q_fbase_nyq%evalDOF_IP(0,q_in_s(:))
    ! f_IP=q*(1+dlambda/dtheta) 
    f_IP  = q_IP(:)*(1.0_wp + LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:)))
    
    !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
    DO i_mn=1,mn_IP
      modes_IP(:,i_mn)= q_fbase_out%eval(0,(/q_fbase_nyq%x_IP(1,i_mn)+LA_IP(i_mn),q_fbase_nyq%x_IP(2,i_mn)/))
    END DO

    q_gIP(:,is)=(dthet_dzeta*q_fbase_out%snorm_base(:))*(MATMUL(modes_IP(:,1:mn_IP),f_IP(1:mn_IP)))  !integrate (sum over mn_IP), includes normalization of base!


    !CHECK at all IP points:
    f_IP=MATMUL(TRANSPOSE(modes_IP(:,:)),q_gIP(:,is))
!    IF(q_base_in%f%mn_zero_mode.NE.-1) WRITE(*,*)'check zero modes  q_in,qout' , &
!                                       q_in_s(q_base_in%f%mn_zero_mode),q_gIP(q_fbase_out%mn_zero_mode,is)
!    WRITE(*,*)'check bounds LA_s  (min/max/avg)',MINVAL(LA_IP),MAXVAL(LA_IP),SUM(LA_IP)/REAL(mn_IP,wp)
!    WRITE(*,*)'check bounds q_in  (min/max/avg)',MINVAL(q_IP),MAXVAL(q_IP),SUM(q_IP)/REAL(mn_IP,wp)
!    WRITE(*,*)'check bounds q_out (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)
    f_IP = ABS(f_IP - q_IP)
    check(1)=MIN(check(1),MINVAL(f_IP))
    check(2)=MAX(check(2),MAXVAL(f_IP))
    check(3)=check(3)+MAXVAL(f_IP)/REAL(nBase,wp)
    WRITE(*,*)'--------  spos= ',spos
    WRITE(*,*)'check |q_in-q_out| (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)

  END DO !is

  !finalize
  CALL LA_fbase_nyq%free()
  CALL q_fbase_nyq%free()
  DEALLOCATE(LA_fbase_nyq,q_fbase_nyq,modes_IP,f_IP,LA_IP)

  WRITE(UNIT_StdOut,'(A,3E11.3)')'   check |q_in-q_out| (min/max/avg)',check(1:3)
  !transform back to corresponding representation of DOF in s
  ALLOCATE(q_out(nBase,1:q_fbase_out%modes))
  DO iMode=1,q_fbase_out%modes
    q_out(:,iMode)=q_base_in%s%initDOF( q_gIP(iMode,:) )
  END DO
  DEALLOCATE(q_gIP)

END SUBROUTINE Transform_to_PEST


END MODULE MODgvec_Transform_SFL
