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
USE MODgvec_base   ,ONLY: t_base,t_sgrid
USE MODgvec_hmap,  ONLY: c_hmap
IMPLICIT NONE
PRIVATE

TYPE :: t_transform_sfl
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER                     :: whichSFLcoord !! 
  INTEGER                     :: fac_nyq,mn_max(2),mn_nyq(2),deg,continuity,degGP,nfp
  INTEGER                     :: X1sfl_sin_cos,X2sfl_sin_cos,GZ_sin_cos
  LOGICAL                     :: relambda      !! =T: recompute lambda from mapping on full fourier series, =F: use LA from eq.
  TYPE(t_sgrid)               :: sgrid_sfl     !! grid for SFL coordinates
  CLASS(c_hmap),  POINTER     :: hmap          !! pointer to hmap class

  CLASS(t_base),  ALLOCATABLE :: X1sfl_base    !! container for base of variable X1 in SFL coordinates
  CLASS(t_base),  ALLOCATABLE :: X2sfl_base    !! container for base of variable X2 in SFL coordinates
  CLASS(t_base),  ALLOCATABLE :: GZ_base       !! container for base of variable and Gzeta (transforms to BOOZER!)
  CLASS(t_base),  ALLOCATABLE :: GZsfl_base    !! container for base of variable Gzeta in SFL coordinates
  REAL(wp),       ALLOCATABLE :: X1sfl(:,:)    !! data (1:nBase,1:modes) of X1 in SFL coords.
  REAL(wp),       ALLOCATABLE :: X2sfl(:,:)    !! data (1:nBase,1:modes) of X2 in SFL coords.
  REAL(wp),       ALLOCATABLE :: Gthet(:,:)    !! data (1:nBase,1:modes) of Gthet in GVEC coords. (for BOOZER)
  REAL(wp),       ALLOCATABLE :: GZ(:,:)       !! data (1:nBase,1:modes) of GZ in GVEC coords. (for BOOZER)
  REAL(wp),       ALLOCATABLE :: Gtsfl(:,:)    !! data (1:nBase,1:modes) of Gt in SFL coords.  (for BOOZER)
  REAL(wp),       ALLOCATABLE :: GZsfl(:,:)    !! data (1:nBase,1:modes) of GZ in SFL coords.  (for BOOZER)
  PROCEDURE(i_func_evalprof), POINTER, NOPASS  :: eval_phiPrime
  PROCEDURE(i_func_evalprof), POINTER, NOPASS  :: eval_iota
  CONTAINS
  !PROCEDURE(i_func_evalprof),DEFERRED :: evalphiPrime
  PROCEDURE :: init       => transform_sfl_init
  PROCEDURE :: BuildTransform => BuildTransform_SFL
  PROCEDURE :: get_boozer  => get_boozer_sinterp
  PROCEDURE :: free        => transform_sfl_free
END TYPE t_transform_sfl

ABSTRACT INTERFACE
  FUNCTION i_func_evalprof(spos)
    IMPORT wp
    REAL(wp),INTENT(IN):: spos
    REAL(wp)           :: i_func_evalprof
  END FUNCTION i_func_evalprof
END INTERFACE

INTERFACE transform_sfl_new
  MODULE PROCEDURE transform_sfl_new
END INTERFACE

PUBLIC :: t_transform_sfl,transform_sfl_new
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Allocate class and call init
!!
!===================================================================================================================================
SUBROUTINE transform_sfl_new(sf,mn_max_in, whichSFL,relambda,deg_in,continuity_in,degGP_in,grid_in,  &
                             hmap_in,X1_base_in,X2_base_in,LA_base_in,eval_phiPrime_in,eval_iota_in)
  ! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER     ,INTENT(IN) :: mn_max_in(2)                                        !< maximum number for new variables in SFL coordinates
  INTEGER     ,INTENT(IN) :: whichSFL                                         !< either =1: PEST, =2:Boozer 
  LOGICAL     ,INTENT(IN) :: relambda                                         !< recompute lambda...
  INTEGER     ,INTENT(IN) :: deg_in,continuity_in,degGP_in                    !< for output base (X1,X2,G)
  CLASS(t_sgrid),INTENT(IN),TARGET :: grid_in                            !! grid information
  CLASS(t_base),INTENT(IN),TARGET :: X1_base_in,X2_base_in,LA_base_in           !< base classes belong to solution U_in
  CLASS(c_hmap),INTENT(IN),TARGET :: hmap_in
  PROCEDURE(i_func_evalprof)     :: eval_phiPrime_in,eval_iota_in  !!procedure pointers to profile evaluation functions.
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  CLASS(t_transform_sfl), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
  !===================================================================================================================================
  ALLOCATE(t_transform_sfl :: sf)
  sf%hmap => hmap_in
  sf%eval_phiPrime=>eval_PhiPrime_in
  sf%eval_iota=>eval_iota_in
  !TEST
  !WRITE(*,*)'DEBUG,phiprime= ? ',sf%eval_phiPrime(0.0_wp),sf%eval_phiPrime(1.0_wp)
  !WRITE(*,*)'DEBUG,iota= ? ',sf%eval_iota(0.0_wp),sf%eval_iota(1.0_wp)
  
  !pass any grid here
  CALL sf%sgrid_sfl%copy(grid_in)
  sf%mn_max=mn_max_in; sf%deg=deg_in; sf%continuity=continuity_in ; sf%degGP = degGP_in
  sf%nfp = X1_base_in%f%nfp
  sf%whichSFLcoord=whichSFL
  sf%relambda=.TRUE. !relambda !if =True, J^s=0 will be recomputed, for exact integrability condition of boozer transform  (but slower!)
  sf%fac_nyq=4  !hard coded for now
  ! use maximum number of integration points from maximum mode number in both directions
  sf%mn_nyq(1:2)=sf%fac_nyq*MAXVAL(sf%mn_max)+1 
  IF(sf%mn_max(2).EQ.0) sf%mn_nyq(2)=1 !exception: 2D configuration
  sf%X1sfl_sin_cos=X1_base_in%f%sin_cos
  sf%X2sfl_sin_cos=X2_base_in%f%sin_cos
  sf%GZ_sin_cos   =LA_base_in%f%sin_cos
  CALL sf%init()
END SUBROUTINE transform_sfl_new

!===================================================================================================================================
!> get_new 
!!
!===================================================================================================================================
SUBROUTINE transform_SFL_init(sf)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut
USE MODgvec_base   ,ONLY: t_base,base_new
USE MODgvec_fbase  ,ONLY: sin_cos_map
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_transform_sfl), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! extended base for q in the new angles, and on the new grid
CALL base_new(sf%X1sfl_base,  sf%deg, sf%continuity, sf%sgrid_sfl, sf%degGP,      &
               sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%X1sfl_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform
CALL base_new(sf%X2sfl_base,   sf%deg, sf%continuity, sf%sgrid_sfl,sf%degGP,      &
              sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%X2sfl_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform
ALLOCATE(sf%X1sfl(sf%X1sfl_base%s%nBase,sf%X1sfl_base%f%modes)); sf%X1sfl=0.0_wp
ALLOCATE(sf%X2sfl(sf%X2sfl_base%s%nBase,sf%X2sfl_base%f%modes)); sf%X2sfl=0.0_wp

SELECT CASE(sf%whichSFLcoord)
CASE(1) !PEST
 ! nothing to initialize additionally
CASE(2) !BOOZER
  CALL base_new(sf%GZ_base, sf%deg, sf%continuity, sf%sgrid_sfl,sf%degGP,      &
                 sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%GZ_sin_cos),.TRUE.) !exclude m=n=0
                 
  ALLOCATE(sf%Gthet(sf%GZ_base%s%nBase,sf%GZ_base%f%modes)); sf%Gthet=0.0_wp
  ALLOCATE(sf%GZ(   sf%GZ_base%s%nBase,sf%GZ_base%f%modes)); sf%GZ=0.0_wp

  CALL base_new(sf%GZsfl_base, sf%deg, sf%continuity, sf%sgrid_sfl,sf%degGP,      &
  sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%GZ_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform

  ALLOCATE(sf%Gtsfl(sf%GZsfl_base%s%nBase,sf%GZsfl_base%f%modes));sf%Gtsfl=0.0_wp
  ALLOCATE(sf%GZsfl(sf%GZsfl_base%s%nBase,sf%GZsfl_base%f%modes));sf%GZsfl=0.0_wp

CASE DEFAULT
SWRITE(UNIT_stdOut,*)'This input for SFL coordinate transform is not valid: whichSFL=',sf%whichSFLcoord
CALL abort(__STAMP__, &
           "wrong input for SFL coordinate transform")
END SELECT

sf%initialized=.TRUE.
END SUBROUTINE transform_sfl_init




!===================================================================================================================================
!> Builds X1 and X2 in SFL coordinates
!!
!===================================================================================================================================
SUBROUTINE BuildTransform_SFL(sf,X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in)
! MODULES
USE MODgvec_base   ,ONLY: t_base,base_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------

! INPUT VARIABLES
  CLASS(t_base),INTENT(IN) :: X1_base_in,X2_base_in,LA_base_in           !< base classes belong to solution U_in
  REAL(wp),INTENT(IN)      :: X1_in(1:X1_base_in%s%nbase,1:X1_base_in%f%modes)
  REAL(wp),INTENT(IN)      :: X2_in(1:X2_base_in%s%nbase,1:X2_base_in%f%modes)
  REAL(wp),INTENT(IN)      :: LA_in(1:LA_base_in%s%nbase,1:LA_base_in%f%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_transform_sfl), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(sf%whichSFLcoord)
CASE(1) !PEST
  CALL Transform_Angles_sinterp(LA_base_in,LA_in,X1_base_in,X1_in,"X1",sf%X1sfl_base,sf%X1sfl)
  CALL Transform_Angles_sinterp(LA_base_in,LA_in,X2_base_in,X2_in,"X2",sf%X2sfl_base,sf%X2sfl)
CASE(2) !BOOZER
  CALL sf%Get_Boozer(X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in) ! fill sf%GZ,sf%Gthet

  CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,sf%GZ_base,sf%GZ,"GZ",sf%GZsfl_base,sf%GZsfl,B_in=sf%GZ)
  CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,sf%GZ_base,sf%Gthet,"Gt",sf%Gzsfl_base,sf%Gtsfl,B_in=sf%GZ)
  CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,X1_base_in,X1_in,"X1",sf%X1sfl_base,sf%X1sfl,B_in=sf%GZ)
  CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,X2_base_in,X2_in,"X2",sf%X2sfl_base,sf%X2sfl,B_in=sf%GZ)
END SELECT

END SUBROUTINE BuildTransform_SFL


!===================================================================================================================================
!> Transform a function from VMEC angles q(s,theta,zeta) to new angles q*(s,theta*,zeta*) 
!> by projection onto the modes of the new angles: sigma_mn(theta*,zeta*)
!> using a given in s 
!> Here, new angles are theta*=theta+A(theta,zeta), zeta*=zeta+B(theta,zeta), 
!> with A,B periodic functions and zero average and same base 
!> Note that in this routine, the integral is transformed back to (theta,zeta)
!> q*_mn = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) dtheta* dzeta*
!>       = iint_0^2pi q(theta,zeta) sigma_mn(theta*,zeta*) [(1+dA/dtheta)*(1+dB/dzeta)-(dA/dzeta*dB/dzeta)] dtheta dzeta
!!
!===================================================================================================================================
SUBROUTINE Transform_Angles_sinterp(AB_base_in,A_in,q_base_in,q_in,q_name,q_base_out,q_out,B_in)
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
  CHARACTER(LEN=*),INTENT(IN):: q_name
  CLASS(t_base),INTENT(IN) :: q_base_out                                    !< new fourier basis of function q in new angles, defined mn_max,mn_nyq 
  REAL(wp)    ,INTENT(IN),OPTIONAL :: B_in(1:AB_base_in%s%nBase,1:AB_base_in%f%modes) !< coefficients of zeta*=zeta+B(s,theta,zeta)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES      
  REAL(wp) ,INTENT(INOUT) :: q_out(q_base_out%s%nBase,1:q_base_out%f%modes)          !< coefficients of q in new angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,i,is,iMode,i_mn,mn_IP,mn_max(2),mn_nyq(2)
  INTEGER               :: BCtype_axis(0:4),BC_type_axis_imode
  REAL(wp)              :: spos,dthet_dzeta 
  REAL(wp)              :: check(1:7,q_base_out%s%nBase)
  LOGICAL               :: docheck=.TRUE. 
  LOGICAL               :: Bpresent
  REAL(wp)              :: A_s(1:AB_base_in%f%modes) 
  REAL(wp)              :: B_s(1:AB_base_in%f%modes) 
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes) 
  REAL(wp), ALLOCATABLE :: q_IP(:),q_m(:,:)   ! q evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: f_IP(:)       ! =q*(1+dlambda/dtheta) evaluated at integration points
  REAL(wp), ALLOCATABLE :: modes_IP(:,:) ! mn modes of q evaluated at theta*,zeta* for all integration points
  CLASS(t_fBase),ALLOCATABLE        :: q_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE        :: AB_fbase_nyq
  REAL(wp),DIMENSION(:),ALLOCATABLE :: A_IP,dAdthet_IP,B_IP,dBdthet_IP,dBdzeta_IP,dAdzeta_IP
!===================================================================================================================================
  Bpresent=PRESENT(B_in)
  mn_max(1:2) =q_base_out%f%mn_max
  mn_nyq(1:2) =q_base_out%f%mn_nyq

  SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6),A,L)')'TRANSFORM '//TRIM(q_name)//' TO NEW ANGLE COORDINATES, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq, ', B_in= ',Bpresent
                     
  __PERFON('transform_angles')
  __PERFON('init')
  !initialize


  !total number of integration points
  mn_IP = q_base_out%f%mn_IP
  nBase = q_base_out%s%nBase

  SWRITE(UNIT_StdOut,*)'        ...Init q_out Base Done'


  !same base for X1, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new( q_fbase_nyq,  q_base_in%f%mn_max,  mn_nyq, &
                                q_base_in%f%nfp, &
                    sin_cos_map(q_base_in%f%sin_cos), &
                                q_base_in%f%exclude_mn_zero)
  SWRITE(UNIT_StdOut,*)'        ...Init q_nyq Base Done'
  !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new(AB_fbase_nyq,  AB_base_in%f%mn_max,  mn_nyq, &
                                AB_base_in%f%nfp, &
                    sin_cos_map(AB_base_in%f%sin_cos), &
                                AB_base_in%f%exclude_mn_zero)
  
  SWRITE(UNIT_StdOut,*)'        ...Init AB_nyq Base Done'

  IF(.NOT.Bpresent) THEN
    ALLOCATE(A_IP(1:mn_IP),dAdthet_IP(1:mn_IP))
  ELSE ! Bpresent
    ALLOCATE(A_IP(1:mn_IP),dAdthet_IP(1:mn_IP),dAdzeta_IP(1:mn_IP),B_IP(1:mn_IP),dBdthet_IP(1:mn_IP),dBdzeta_IP(1:mn_IP))
  END IF !.NOT.Bpresent

  ALLOCATE(f_IP(1:mn_IP),q_IP(1:mn_IP),modes_IP(1:q_base_out%f%modes,1:mn_IP))

  ALLOCATE(q_m(1:q_base_out%f%modes,nBase))
  __PERFOFF('init')

  dthet_dzeta  =q_base_out%f%d_thet*q_base_out%f%d_zeta !integration weights

  CALL ProgressBar(0,nBase)!init
  DO is=1,nBase
    __PERFON('eval_data_s')
    spos=MIN(MAX(1.0e-12_wp,q_base_out%s%s_IP(is)),1.0_wp-1.0e-12_wp) !interpolation points for q_in
    !evaluate q_in at spos
    q_in_s(:)= q_base_in%s%evalDOF2D_s(spos,q_base_in%f%modes,   0,q_in(:,:))
    !evaluate A at spos
    A_s(:)= AB_base_in%s%evalDOF2D_s(spos,AB_base_in%f%modes,   0,A_in(:,:))
    IF(Bpresent)THEN
      B_s(:)     = AB_base_in%s%evalDOF2D_s(spos,AB_base_in%f%modes,   0,B_in(:,:))
    END IF
    __PERFOFF('eval_data_s')
    __PERFON('eval_data_f')
    !evaluate lambda at spos
    ! TEST EXACT CASE: A_s=0.
    
    IF(.NOT.Bpresent)THEN
      q_IP       =  q_fbase_nyq%evalDOF_IP(         0, q_in_s(  :))
      A_IP       = AB_fbase_nyq%evalDOF_IP(         0, A_s(  :))
      dAdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET, A_s(  :))
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)        &
!$OMP   SHARED(mn_IP,q_base_out,q_IP,A_IP,dAdthet_IP,f_IP,modes_IP)
      DO i_mn=1,mn_IP
        f_IP(i_mn) = q_IP(i_mn)*(1.0_wp + dAdthet_IP(i_mn))
        !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
        modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_base_out%f%x_IP(1,i_mn)+A_IP(i_mn),q_base_out%f%x_IP(2,i_mn)/))
      END DO !i_mn=1,mn_IP
!$OMP END PARALLEL DO 

    ELSE !Bpresent
      q_IP       =  q_fbase_nyq%evalDOF_IP(         0, q_in_s(  :))
      A_IP       = AB_fbase_nyq%evalDOF_IP(         0, A_s(  :))
      dAdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET, A_s(  :))
      dAdzeta_IP = AB_fbase_nyq%evalDOF_IP(DERIV_ZETA, A_s(  :))
      B_IP       = AB_fbase_nyq%evalDOF_IP(         0, B_s(  :))
      dBdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET, B_s(  :))
      dBdzeta_IP = AB_fbase_nyq%evalDOF_IP(DERIV_ZETA, B_s(  :))
      
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)        &
!$OMP   SHARED(mn_IP,q_base_out,q_IP,A_IP,dAdthet_IP,B_IP,dBdthet_IP,dBdzeta_IP,dAdzeta_IP,f_IP,modes_IP)
      DO i_mn=1,mn_IP
        f_IP(i_mn) = q_IP(i_mn)*((1.0_wp + dAdthet_IP(i_mn))*(1.0_wp + dBdzeta_IP(i_mn))-dAdzeta_IP(i_mn)*dBdthet_IP(i_mn))
        !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
        modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_base_out%f%x_IP(1,i_mn)+A_IP(i_mn),q_base_out%f%x_IP(2,i_mn)+B_IP(i_mn)/))
      END DO !i_mn=1,mn_IP
!$OMP END PARALLEL DO 
    END IF !Bpresent
    __PERFON('project')
    __MATVEC_N(q_m(:,is),modes_IP(:,:),f_IP(:)) 

    __PERFOFF('project')
    __PERFOFF('eval_data_f')

    !projection: integrate (sum over mn_IP), includes normalization of base!
    !q_m(:,is)=(dthet_dzeta*q_base_out%f%snorm_base(:))*(MATMUL(modes_IP(:,1:mn_IP),f_IP(1:mn_IP)))  

    q_m(:,is)=dthet_dzeta*q_base_out%f%snorm_base(:)*q_m(:,is)

    __PERFON('check')
    !CHECK at all IP points
    IF(doCheck) THEN

!      __MATVEC_N(f_IP,q_base_out%f%base_IP,q_m(:,is)) !other points
!      check(6)=MIN(check(6),MINVAL(f_IP))
!      check(7)=MAX(check(7),MAXVAL(f_IP))
      check(6,is)=MINVAL(q_IP)
      check(7,is)=MAXVAL(q_IP)
      !f_IP=MATMUL(q_m(:,is),modes_IP(:,:))
      __MATVEC_T(f_IP,modes_IP,q_m(:,is)) !same points
      check(4,is)=MINVAL(f_IP)
      check(5,is)=MAXVAL(f_IP)

  
      !f_IP = ABS(f_IP - q_IP)/SUM(ABS(q_IP))*REAL(mn_IP,wp)
      f_IP=ABS(f_ip- q_IP)
      check(1,is)=MINVAL(f_IP)
      check(2,is)=MAXVAL(f_IP)
      check(3,is)=SUM(f_IP)/REAL(mn_IP,wp)
      !WRITE(*,*)'     ------  spos= ',spos
      !WRITE(*,*)'check |q_in-q_out|/(surfavg|q_in|) (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)
      !WRITE(*,*)'max,min,sum q_out |modes|',MAXVAL(ABS(q_out(is,:))),MINVAL(ABS(q_out(is,:))),SUM(ABS(q_out(is,:)))
    END IF !doCheck

    __PERFOFF('check')
    CALL ProgressBar(is,nBase)
  END DO !is

  !BCtype_axis(MN_ZERO    )= BC_TYPE_OPEN   !do nothing
  !BCtype_axis(M_ZERO     )= BC_TYPE_OPEN   !do nothing
  !BCtype_axis(M_ODD_FIRST)= BC_TYPE_DIRICHLET !=0
  !BCtype_axis(M_ODD      )= BC_TYPE_DIRICHLET !=0
  !BCtype_axis(M_EVEN     )= BC_TYPE_DIRICHLET !=0
  !Apply automatic, smooth axis BC:
  BCtype_axis=0
  !transform back to corresponding representation of DOF in s
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode)
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_out%s%initDOF( q_m(iMode,:) )
    BC_type_axis_imode=BCtype_axis(q_base_out%f%zero_odd_even(iMode))
    IF(BC_type_axis_imode.EQ.0) & !AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
        BC_type_axis_imode=-1*MIN(q_base_out%s%deg+1,q_base_out%f%Xmn(1,iMode))
    CALL q_base_out%s%applyBCtoDOF(q_out(:,iMode),(/BC_type_axis_imode,BC_TYPE_OPEN/),(/0.,0./))
  END DO
!$OMP END PARALLEL DO 

  !finalize
  CALL q_fbase_nyq%free()
  CALL AB_fbase_nyq%free()
  DEALLOCATE( q_fbase_nyq, AB_fbase_nyq, A_IP, dAdthet_IP)
  IF(Bpresent) DEALLOCATE(dAdzeta_IP, B_IP,dBdthet_IP, dBdzeta_IP )

  DEALLOCATE(modes_IP,q_IP,f_IP,q_m)

  IF(doCheck) THEN

    DO i=4,1,-1
      is=MAX(1,nBase/i)
      WRITE(UNIT_StdOut,'(A,E11.4)')'at rho=',q_base_out%s%s_IP(is)
      WRITE(UNIT_StdOut,'(A,2E21.11)') '   MIN/MAX of input  '//TRIM(q_name)//':', check(6:7,is)
      WRITE(UNIT_StdOut,'(A,2E21.11)') '   MIN/MAX of output '//TRIM(q_name)//':', check(4:5,is)
      WRITE(UNIT_StdOut,'(2A,3E11.4)')'    ERROR of '//TRIM(q_name)//':', &
                                      ' |q_in-q_out| (min/max/avg)',check(1:2,is),check(3,is)
    END DO                                  
  END IF

  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
  __PERFOFF('transform_angles')
END SUBROUTINE Transform_Angles_sinterp

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
SUBROUTINE Transform_Angles_sproject(AB_base_in,A_in,q_base_in,q_in,q_name,q_base_out,q_out,B_in)
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
  CHARACTER(LEN=*),INTENT(IN):: q_name
  CLASS(t_Base),INTENT(IN) :: q_base_out          !< new fourier basis of function q in new angles, defines mn_max and mn_nyq
  REAL(wp)    ,INTENT(IN),OPTIONAL :: B_in(1:AB_base_in%s%nBase,1:AB_base_in%f%modes) !< coefficients of zeta*=zeta+B(s,theta,zeta)

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES

  REAL(wp) ,INTENT(INOUT) :: q_out(q_base_out%s%nBase,1:q_base_out%f%modes)          !< coefficients of q in new angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iBase,nBase,iGP,nGP,iElem,nElems,degGP,deg,iMode,i_mn,mn_IP,modes,mn_max(2),mn_nyq(2)
  INTEGER               :: BCtype_axis(0:4),BC_type_axis_imode
  REAL(wp)              :: spos,dthet_dzeta 
  REAL(wp)              :: check(1:7)
  LOGICAL               :: docheck=.TRUE. 
  LOGICAL               :: Bpresent
  REAL(wp)              :: A_s(1:AB_base_in%f%modes) 
  REAL(wp)              :: B_s(1:AB_base_in%f%modes) 
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes) 
  REAL(wp), ALLOCATABLE :: q_IP(:),q_m(:,:)   ! q evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: f_IP(:)       ! =q*(1+dlambda/dtheta) evaluated at integration points
  REAL(wp), ALLOCATABLE :: modes_IP(:,:) ! mn modes of q evaluated at theta*,zeta* for all integration points
  CLASS(t_fBase),ALLOCATABLE        :: q_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE        :: AB_fbase_nyq
  REAL(wp),DIMENSION(:),ALLOCATABLE :: A_IP,dAdthet_IP,B_IP,dBdthet_IP,dBdzeta_IP,dAdzeta_IP
!===================================================================================================================================
  Bpresent=PRESENT(B_in)
  mn_max(1:2) =q_base_out%f%mn_max
  mn_nyq(1:2) =q_base_out%f%mn_nyq


  SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6),A,L)')'TRANSFORM '//TRIM(q_name)//' TO NEW ANGLE COORDINATES, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq, ', B_in= ',Bpresent
                     
  __PERFON('transform_angles')
  __PERFON('init')


  !total number of integration points
  mn_IP = q_base_out%f%mn_IP
  modes = q_base_out%f%modes
  nBase = q_base_out%s%nBase
  nGP   = q_base_out%s%nGP
  nElems= q_base_out%s%grid%nElems
  degGP = q_base_out%s%degGP
  deg   = q_base_out%s%deg

  SWRITE(UNIT_StdOut,*)'        ...Init q_out Base Done'


  !same base for X1, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new( q_fbase_nyq,  q_base_in%f%mn_max,  mn_nyq, &
                                q_base_in%f%nfp, &
                    sin_cos_map(q_base_in%f%sin_cos), &
                                q_base_in%f%exclude_mn_zero)
  SWRITE(UNIT_StdOut,*)'        ...Init q_nyq Base Done'
  !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new(AB_fbase_nyq,  AB_base_in%f%mn_max,  mn_nyq, &
                                AB_base_in%f%nfp, &
                    sin_cos_map(AB_base_in%f%sin_cos), &
                                AB_base_in%f%exclude_mn_zero)
  
  SWRITE(UNIT_StdOut,*)'        ...Init AB_nyq Base Done'

  IF(.NOT.Bpresent) THEN
    ALLOCATE(A_IP(1:mn_IP),dAdthet_IP(1:mn_IP))
  ELSE ! Bpresent
    ALLOCATE(A_IP(1:mn_IP),dAdthet_IP(1:mn_IP),dAdzeta_IP(1:mn_IP),B_IP(1:mn_IP),dBdthet_IP(1:mn_IP),dBdzeta_IP(1:mn_IP))
  END IF !.NOT.Bpresent

  ALLOCATE(f_IP(1:mn_IP),q_IP(1:mn_IP),modes_IP(1:modes,1:mn_IP))

  ALLOCATE(q_m(1:modes,1:nGP))
  __PERFOFF('init')

  check(1:7)=(/HUGE(1.),0.,0.,HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.)/)

  dthet_dzeta  =q_base_out%f%d_thet*q_base_out%f%d_zeta !integration weights

  CALL ProgressBar(0,nGP)!init
  DO iGP=1,nGP
    __PERFON('eval_data_s')
    spos=q_base_out%s%s_GP(iGP)
    !evaluate q_in at spos
    q_in_s(:)= q_base_in%s%evalDOF2D_s(spos,q_base_in%f%modes,   0,q_in(:,:))
    !evaluate A at spos
    A_s(:)= AB_base_in%s%evalDOF2D_s(spos,AB_base_in%f%modes,   0,A_in(:,:))
    IF(Bpresent)THEN
      B_s(:)     = AB_base_in%s%evalDOF2D_s(spos,AB_base_in%f%modes,   0,B_in(:,:))
    END IF
    __PERFOFF('eval_data_s')
    __PERFON('eval_data_f')
    !evaluate lambda at spos
    ! TEST EXACT CASE: A_s=0.
    
    IF(.NOT.Bpresent)THEN
      q_IP       =  q_fbase_nyq%evalDOF_IP(         0, q_in_s(  :))
      A_IP       = AB_fbase_nyq%evalDOF_IP(         0, A_s(  :))
      dAdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET, A_s(  :))
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)        &
!$OMP   SHARED(mn_IP,q_base_out,q_IP,A_IP,dAdthet_IP,f_IP,modes_IP)
      DO i_mn=1,mn_IP
        f_IP(i_mn) = q_IP(i_mn)*(1.0_wp + dAdthet_IP(i_mn))
        !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
        modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_base_out%f%x_IP(1,i_mn)+A_IP(i_mn),q_base_out%f%x_IP(2,i_mn)/))
      END DO !i_mn=1,mn_IP
!$OMP END PARALLEL DO 

    ELSE !Bpresent
      q_IP       =  q_fbase_nyq%evalDOF_IP(         0, q_in_s(  :))
      A_IP       = AB_fbase_nyq%evalDOF_IP(         0, A_s(  :))
      dAdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET, A_s(  :))
      dAdzeta_IP = AB_fbase_nyq%evalDOF_IP(DERIV_ZETA, A_s(  :))
      B_IP       = AB_fbase_nyq%evalDOF_IP(         0, B_s(  :))
      dBdthet_IP = AB_fbase_nyq%evalDOF_IP(DERIV_THET, B_s(  :))
      dBdzeta_IP = AB_fbase_nyq%evalDOF_IP(DERIV_ZETA, B_s(  :))
      
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)        &
!$OMP   SHARED(mn_IP,q_base_out,q_IP,A_IP,dAdthet_IP,B_IP,dBdthet_IP,dBdzeta_IP,dAdzeta_IP,f_IP,modes_IP)
      DO i_mn=1,mn_IP
        f_IP(i_mn) = q_IP(i_mn)*((1.0_wp + dAdthet_IP(i_mn))*(1.0_wp + dBdzeta_IP(i_mn))-dAdzeta_IP(i_mn)*dBdthet_IP(i_mn))
        !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
        modes_IP(:,i_mn)= q_base_out%f%eval(0,(/q_base_out%f%x_IP(1,i_mn)+A_IP(i_mn),q_base_out%f%x_IP(2,i_mn)+B_IP(i_mn)/))
      END DO !i_mn=1,mn_IP
!$OMP END PARALLEL DO 
    END IF !Bpresent
    __PERFON('project')
    __MATVEC_N(q_m(:,iGP),modes_IP(:,:),f_IP(:)) 

    __PERFOFF('project')
    __PERFOFF('eval_data_f')

    !projection: integrate (sum over mn_IP), includes normalization of base!
    !q_m(:,is)=(dthet_dzeta*q_base_out%f%snorm_base(:))*(MATMUL(modes_IP(:,1:mn_IP),f_IP(1:mn_IP)))  

    q_m(:,iGP)=dthet_dzeta*q_base_out%f%snorm_base(:)*q_m(:,iGP)

    __PERFON('check')
    !CHECK at all IP points, only every 4th radial point:
    IF(doCheck) THEN

!      __MATVEC_N(f_IP,q_base_out%f%base_IP,q_m(:,is)) !other points
!      check(6)=MIN(check(6),MINVAL(f_IP))
!      check(7)=MAX(check(7),MAXVAL(f_IP))

      !f_IP=MATMUL(q_m(:,is),modes_IP(:,:))
      __MATVEC_T(f_IP,modes_IP,q_m(:,iGP)) !same points
      check(6)=MIN(check(6),MINVAL(f_IP))
      check(7)=MAX(check(7),MAXVAL(f_IP))
  
      f_IP = ABS(f_IP - q_IP)/SUM(ABS(q_IP))*REAL(mn_IP,wp)
      check(1)=MIN(check(1),MINVAL(f_IP))
      check(2)=MAX(check(2),MAXVAL(f_IP))
      check(3)=check(3)+MAXVAL(f_IP)/REAL(nBase,wp)
      check(4)=MIN(check(4),MINVAL(q_IP))
      check(5)=MAX(check(5),MAXVAL(q_IP))
!      WRITE(*,*)'     ------  spos= ',spos
!      WRITE(*,*)'check |q_in-q_out|/(surfavg|q_in|) (min/max/avg)',MINVAL(f_IP),MAXVAL(f_IP),SUM(f_IP)/REAL(mn_IP,wp)
!      WRITE(*,*)'max,min,sum q_out |modes|',MAXVAL(ABS(q_out(is,:))),MINVAL(ABS(q_out(is,:))),SUM(ABS(q_out(is,:)))
    END IF !doCheck
    q_m(:,iGP)=q_m(:,iGP)*q_base_out%s%w_GP(iGP)
    __PERFOFF('check')
    CALL ProgressBar(iGP,nGP)
  END DO !iGP

  !BCtype_axis(MN_ZERO    )= BC_TYPE_OPEN   !do nothing
  !BCtype_axis(M_ZERO     )= BC_TYPE_OPEN   !do nothing
  !BCtype_axis(M_ODD_FIRST)= BC_TYPE_DIRICHLET !=0
  !BCtype_axis(M_ODD      )= BC_TYPE_DIRICHLET !=0
  !BCtype_axis(M_EVEN     )= BC_TYPE_DIRICHLET !=0
  !Apply automatic, smooth axis BC:
  BCtype_axis=0
  !transform back to corresponding representation of DOF in s
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,iElem,iGP,iBase)
  DO iMode=1,modes
    q_out(:,iMode)=0.0_wp
    DO iElem=1,nElems
      iGP=(iElem-1)*(degGP+1)+1
      ibase=q_base_out%s%base_offset(iElem)
      q_out(iBase:iBase+deg,iMode) = q_out(iBase:iBase+deg,iMode) & 
                                    + MATMUL(q_m(iMode,iGP:iGP+degGP),q_base_out%s%base_GP(0:degGP,0:deg,iElem))
    END DO !iElem
    CALL q_base_out%s%mass%solve_inplace(1,q_out(:,iMode))
    BC_type_axis_imode=BCtype_axis(q_base_out%f%zero_odd_even(iMode))
    IF(BC_type_axis_imode.EQ.0) & !AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
        BC_type_axis_imode=-1*MIN(q_base_out%s%deg+1,q_base_out%f%Xmn(1,iMode))
    CALL q_base_out%s%applyBCtoDOF(q_out(:,iMode),(/BC_type_axis_imode,BC_TYPE_OPEN/),(/0.,0./))
  END DO !iMode
!$OMP END PARALLEL DO 

  !finalize
  CALL q_fbase_nyq%free()
  CALL AB_fbase_nyq%free()
  DEALLOCATE( q_fbase_nyq, AB_fbase_nyq, A_IP, dAdthet_IP)
  IF(Bpresent) DEALLOCATE(dAdzeta_IP, B_IP,dBdthet_IP, dBdzeta_IP )

  DEALLOCATE(modes_IP,q_IP,f_IP,q_m)

  IF(doCheck) THEN
    WRITE(UNIT_StdOut,'(A,2E21.11)') '   MIN/MAX of input  '//TRIM(q_name)//':', check(4:5)
    WRITE(UNIT_StdOut,'(A,2E21.11)') '   MIN/MAX of output '//TRIM(q_name)//':', check(6:7)
    WRITE(UNIT_StdOut,'(2A,3E11.4)')'       CHECK ERROR of '//TRIM(q_name)//':', &
                                             ' |q_in-q_out|/(surfavg|q_in|) (min/max/avg)',check(1:2), &
                                                                            check(3)/(0.25*REAL(nBase,wp))
  END IF

  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
  __PERFOFF('transform_angles')
END SUBROUTINE Transform_Angles_sproject

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
SUBROUTINE Transform_to_PEST_2(LA_base_in,LA,q_base_in,q_in,q_base_out,q_out)
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
    CLASS(t_Base),INTENT(IN) :: q_base_out          !< new fourier basis of function q in new angles, defines mn_max and mn_nyq
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) ,INTENT(INOUT) :: q_out(q_base_out%s%nBase,1:q_base_out%f%modes)          !< coefficients of q in new angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,i_mn,mn_IP,mn_max(2),mn_nyq(2)
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
  mn_max(1:2) =q_base_out%f%mn_max
  mn_nyq(1:2) =q_base_out%f%mn_nyq
  WRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'Transform variable to PEST coordinates, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
                     
 

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
!> Builds the boozer transform coordinate, from restart data
!! theta^B = theta + lambda + iota(s)*G(s,theta,zeta) 
!! zeta^B  = zeta G(s,theta,zeta) 
!!
!! since in Boozer, the covariant magnetic field components are the current profiles,
!! B = Itor(s) grad theta^B + Ipol(s) grad zeta^B + X grad s
!!   = Itor(s) grad (theta+lambda+iota*G) + Ipol(s) grad (zeta + G) + X grad s
!!   = (Itor*(1+dlambda/dtheta) + (Itor*iota+Ipol)*dG/dtheta) grad theta + (Itor*(dlambda/dzeta)+(Itor*iota+Ipol)*dG/dzeta)
!!=> dG/dtheta = (B_theta - Itor - Itor*dlambda/dtheta ) / (Itor*iota+Ipol)
!!=> dG/dzeta  = (B_zeta  - Ipol - Itor*dlambda/dzeta  ) / (Itor*iota+Ipol)
!! There is a integrability condition for G: 
!!     d/dzeta(dG/dtheta)-d/dthet(dG/dzeta)=d/dzeta(dB_theta/dtheta)-d/dthet(dB_theta/dzeta)=0
!! which is equivalent to impose J^s=0. 
!! now if lambda is recomputed via a projection of J^s=0 onto the same fourier series as G, the compatibility condition is 
!! EXACTLY(!) fullfilled. 
!===================================================================================================================================
SUBROUTINE Get_Boozer_sinterp(sf,X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in)
  ! MODULES
  USE MODgvec_Globals,ONLY: UNIT_stdOut,ProgressBar
  USE MODgvec_LinAlg
  USE MODgvec_base   ,ONLY: t_base,base_new
  USE MODgvec_sGrid  ,ONLY: t_sgrid
  USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map 
  USE MODgvec_lambda_solve, ONLY:  Lambda_setup_and_solve
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_base),INTENT(IN) :: X1_base_in,X2_base_in,LA_base_in   !< base classes for input U_in
    REAL(wp),INTENT(IN):: X1_in(1:X1_base_in%s%nbase,1:X1_base_in%f%modes)
    REAL(wp),INTENT(IN):: X2_in(1:X2_base_in%s%nbase,1:X2_base_in%f%modes)
    REAL(wp),INTENT(IN):: LA_in(1:LA_base_in%s%nbase,1:LA_base_in%f%modes)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES  
    CLASS(t_transform_sfl), INTENT(INOUT) :: sf !! self !hmap, profiles, G_base_out,Gthet,GZ !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER               :: mn_max(2),mn_nyq(2),nBase,is,iMode,modes,i_mn,mn_IP
    INTEGER               :: nfp,BCtype_axis(0:4),BC_type_axis_imode
    REAL(wp)              :: spos,dthet_dzeta,dPhids_int,iota_int,dChids_int
    REAL(wp)              :: b_thet,b_zeta,qloc(3),q_thet(3),q_zeta(3) 
    REAL(wp)              :: detJ,Itor,Ipol,stmp
  
    REAL(wp)                          ::  X1_s(  1:X1_base_in%f%modes) 
    REAL(wp)                          :: dX1ds_s(1:X1_base_in%f%modes) 
    REAL(wp)                          ::  X2_s(  1:X2_base_in%f%modes) 
    REAL(wp)                          :: dX2ds_s(1:X2_base_in%f%modes) 
    REAL(wp),ALLOCATABLE              :: LA_s(:)
    REAL(wp),DIMENSION(sf%GZ_base%f%modes) :: GZ_m,GZ_n
    REAL(wp),DIMENSION(sf%GZ_base%f%mn_IP) :: Bcov_thet_IP,Bcov_zeta_IP
    REAL(wp),DIMENSION(sf%GZ_base%f%mn_IP) :: dLAdthet_IP,dLAdzeta_IP
    REAL(wp),DIMENSION(sf%GZ_base%f%mn_IP) :: LA_IP,fm_IP,fn_IP,ft_IP,gam_tt,gam_tz,gam_zz 
    REAL(wp),DIMENSION(sf%GZ_base%f%mn_IP) :: X1_IP,dX1ds_IP,dX1dthet_IP,dX1dzeta_IP
    REAL(wp),DIMENSION(sf%GZ_base%f%mn_IP) :: X2_IP,dX2ds_IP,dX2dthet_IP,dX2dzeta_IP
    CLASS(t_fBase),ALLOCATABLE             :: X1_fbase_nyq
    CLASS(t_fBase),ALLOCATABLE             :: X2_fbase_nyq
    CLASS(t_fBase),ALLOCATABLE             :: LA_fbase_nyq
  !===================================================================================================================================
    nfp = X1_base_in%f%nfp
    mn_max(1:2)=sf%GZ_base%f%mn_max
    mn_nyq(1:2)=sf%GZ_base%f%mn_nyq
    SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'GET BOOZER ANGLE TRANSFORM, nfp=',nfp, &
                                ', mn_max_in=',LA_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
    __PERFON('get_boozer')
    __PERFON('init')
    
    mn_IP        = sf%GZ_base%f%mn_IP  !total number of integration points
    modes        = sf%GZ_base%f%modes  !number of modes in output
    nBase        = sf%GZ_base%s%nBase  !number of radial points in output
    dthet_dzeta  = sf%GZ_base%f%d_thet*sf%GZ_base%f%d_zeta !integration weights
  
     !transpose basis functions and include norm for projection
    SWRITE(UNIT_StdOut,*)'        ...Init G_out Base Done'
    !same base for X1, but with new mn_nyq (for pre-evaluation of basis functions)
    CALL fbase_new( X1_fbase_nyq, X1_base_in%f%mn_max,  mn_nyq, &
                                  X1_base_in%f%nfp, &
                      sin_cos_map(X1_base_in%f%sin_cos), &
                                  X1_base_in%f%exclude_mn_zero)
    SWRITE(UNIT_StdOut,*)'        ...Init X1_nyq Base Done'
  
    CALL fbase_new( X2_fbase_nyq, X2_base_in%f%mn_max,  mn_nyq, &
                                  X2_base_in%f%nfp, &
                      sin_cos_map(X2_base_in%f%sin_cos), &
                                  X2_base_in%f%exclude_mn_zero)
    SWRITE(UNIT_StdOut,*)'        ...Init X2_nyq Base Done'
     IF(.NOT.sf%relambda)THEN
      !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
      CALL fbase_new(LA_fbase_nyq,  LA_base_in%f%mn_max,  mn_nyq, &
                                  LA_base_in%f%nfp, &
                      sin_cos_map(LA_base_in%f%sin_cos), &
                                  LA_base_in%f%exclude_mn_zero)
      SWRITE(UNIT_StdOut,*)'        ...Init LA_nyq Base Done'
    END IF
    ALLOCATE(LA_s(1:MERGE(sf%GZ_base%f%modes,LA_base_in%f%modes,sf%relambda)))
    !!!ALLOCATE(LA_s(1:LA_fbase_nyq%modes))
    GZ_m=0.0_wp; GZ_n=0.0_wp
  
    __PERFOFF('init')
  
  
    CALL ProgressBar(0,nBase) !INIT
    DO is=nBase,1,-1
      __PERFON('eval_data')
      spos=MIN(MAX(1.0e-08_wp,sf%GZ_base%s%s_IP(is)),1.0_wp-1.0e-12_wp) !interpolation points for q_in
  
      dPhids_int  = sf%eval_phiPrime(spos)
      iota_int    = sf%eval_iota(spos)
      dChids_int  = dPhids_int*iota_int 
  
      !interpolate radially
      X1_s(:)    = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,      0,X1_in(:,:))
      dX1ds_s(:) = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,DERIV_S,X1_in(:,:))
  
      X2_s(:)    = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,      0,X2_in(:,:))
      dX2ds_s(:) = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,DERIV_S,X2_in(:,:))

      IF(.NOT.sf%relambda) THEN
        LA_s(:)    = LA_base_in%s%evalDOF2D_s(spos,LA_base_in%f%modes,      0,LA_in(:,:))
      END IF
  
      !evaluate at integration points
      X1_IP       = X1_fbase_nyq%evalDOF_IP(         0, X1_s(  :))
      dX1ds_IP    = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s(:))
      dX1dthet_IP = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s(  :))
      dX1dzeta_IP = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s(  :))
  
      X2_IP       = X2_fbase_nyq%evalDOF_IP(         0, X2_s(  :))
      dX2ds_IP    = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s(:))
      dX2dthet_IP = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s(  :))
      dX2dzeta_IP = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s(  :))
  

  
      __PERFOFF('eval_data')
      __PERFON('eval_bsub')
      __PERFON('eval_metrics')
      
  !$OMP PARALLEL DO &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE)  &
  !$OMP   PRIVATE(i_mn,qloc,q_thet,q_zeta,detJ)  &
  !$OMP   SHARED(sf,mn_IP,X1_IP,X2_IP, &
  !$OMP          dX1dthet_IP,dX2dthet_IP,dX1dzeta_IP,dX2dzeta_IP,             &
  !$OMP          dX1ds_IP,dX2ds_IP,gam_tt,gam_tz,gam_zz)
      !evaluate metrics on (theta,zeta)
      DO i_mn=1,mn_IP
        qloc(  1:3) = (/ X1_IP(     i_mn), X2_IP(     i_mn),sf%GZ_base%f%x_IP(2,i_mn)/)
        q_thet(1:3) = (/dX1dthet_IP(i_mn),dX2dthet_IP(i_mn),0.0_wp/) !dq(1:2)/dtheta
        q_zeta(1:3) = (/dX1dzeta_IP(i_mn),dX2dzeta_IP(i_mn),1.0_wp/) !dq(1:2)/dzeta
        detJ        = sf%hmap%eval_Jh(qloc)*( dX1ds_IP(i_mn)*dX2dthet_IP(i_mn) &
                      -dX2ds_IP(i_mn)*dX1dthet_IP(i_mn) )
        gam_tt(i_mn)  = sf%hmap%eval_gij(q_thet,qloc,q_thet)/detJ   !g_theta,theta
        gam_tz(i_mn)  = sf%hmap%eval_gij(q_thet,qloc,q_zeta)/detJ   !g_theta,zeta =g_zeta,theta
        gam_zz(i_mn)  = sf%hmap%eval_gij(q_zeta,qloc,q_zeta)/detJ   !g_zeta,zeta
      END DO !i_mn
  !$OMP END PARALLEL DO 
      __PERFOFF('eval_metrics')

      IF(sf%relambda)THEN
      __PERFON('new_lambda')
        CALL Lambda_setup_and_solve(sf%GZ_base%f,dPhids_int,dchids_int,gam_tt,gam_tz,gam_zz,LA_s)
        !CALL Lambda_setup_and_solve(LA_fbase_nyq,dPhids_int,dchids_int,gam_tt,gam_tz,gam_zz,LA_s)
        LA_IP(:)    = sf%GZ_base%f%evalDOF_IP(         0,LA_s(:))
        dLAdthet_IP = sf%GZ_base%f%evalDOF_IP(DERIV_THET,LA_s(:))
        dLAdzeta_IP = sf%GZ_base%f%evalDOF_IP(DERIV_ZETA,LA_s(:))
        __PERFOFF('new_lambda')
      ELSE       
        LA_IP(:)    = LA_fbase_nyq%evalDOF_IP(         0,LA_s(:))
        dLAdthet_IP = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:))
        dLAdzeta_IP = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s(:))
      END IF
      

      Itor=0.0_wp;Ipol=0.0_wp
  !$OMP PARALLEL DO &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE)  &
  !$OMP   PRIVATE(i_mn,b_thet,b_zeta)  &
  !$OMP   REDUCTION(+:Itor,Ipol) &
  !$OMP   SHARED(mn_IP,dchids_int,dPhids_int,dLAdzeta_IP,dLAdthet_IP,gam_tt,gam_tz,gam_zz,Bcov_thet_IP,Bcov_zeta_IP)
      !evaluate B_theta,B_zeta (and integrate for currents)
      DO i_mn=1,mn_IP
        b_thet = dchids_int- dPhids_int*dLAdzeta_IP(i_mn)    !b_theta
        b_zeta = dPhids_int*(1.0_wp   + dLAdthet_IP(i_mn))    !b_zeta
  
        Bcov_thet_IP(i_mn) = (gam_tt(i_mn)*b_thet + gam_tz(i_mn)*b_zeta)
        Bcov_zeta_IP(i_mn) = (gam_tz(i_mn)*b_thet + gam_zz(i_mn)*b_zeta) 
        Itor=Itor+Bcov_thet_IP(i_mn) 
        Ipol=Ipol+Bcov_zeta_IP(i_mn) 
      END DO !i_mn
  !$OMP END PARALLEL DO 
      Itor=(Itor/REAL(mn_IP,wp)) !Itor=zero mode of Bcov_thet
      Ipol=(Ipol/REAL(mn_IP,wp)) !Ipol=zero mode of Bcov_thet
  
  !    Itor=(1.0_wp/REAL(mn_IP,wp))*SUM(Bcov_thet_IP(:)) 
  !    Ipol=(1.0_wp/REAL(mn_IP,wp))*SUM(Bcov_zeta_IP(:))
  
      __PERFOFF('eval_bsub')
      __PERFON('project')
  
      __PERFON('project_G')
  
      stmp=1.0_wp/(Itor*iota_int+Ipol)
  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i_mn)        &
  !$OMP   SHARED(mn_IP,Itor,Ipol,stmp,dLAdthet_IP,Bcov_thet_IP,fm_IP)
      DO i_mn=1,mn_IP
        fm_IP(i_mn)  = (Bcov_thet_IP(i_mn)-Itor-Itor*dLAdthet_IP(i_mn))*stmp
      END DO
  !$OMP END PARALLEL DO 
  
      !projection: only onto base_dthet
      CALL sf%GZ_base%f%projectIPtoDOF(.FALSE.,1.0_wp,DERIV_THET,fm_IP(:),GZ_m(:))
  
      IF(sf%GZ_base%f%mn_max(2).GT.0) THEN !3D case
  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i_mn)        &
  !$OMP   SHARED(mn_IP,Itor,Ipol,stmp,dLAdzeta_IP,Bcov_zeta_IP,fn_IP)
        DO i_mn=1,mn_IP
          fn_IP(i_mn)= (Bcov_zeta_IP(i_mn)-Ipol-Itor*dLAdzeta_IP(i_mn))*stmp
        END DO
  !$OMP END PARALLEL DO 
  
        !projection onto base_dzeta
        CALL sf%GZ_base%f%projectIPtoDOF(.FALSE.,1.0_wp,DERIV_ZETA,fn_IP(:),GZ_n(:))
      END IF !3D case (n_max >0)
  
      ! only if n=0, use formula from base_dthet projected G, else use base_dzeta projected G
      DO iMode=1,modes
        ASSOCIATE(m=>sf%GZ_base%f%Xmn(1,iMode),n=>sf%GZ_base%f%Xmn(2,iMode))
        IF(m.NE.0) GZ_m(iMode)=GZ_m(iMode)*(dthet_dzeta*sf%GZ_base%f%snorm_base(iMode))/REAL(m*m,wp)
        IF(n.NE.0) GZ_n(iMode)=GZ_n(iMode)*(dthet_dzeta*sf%GZ_base%f%snorm_base(iMode))/REAL(n*n,wp)
        IF(n.EQ.0)THEN 
          sf%GZ(is,iMode)=GZ_m(iMode)
        ELSE 
          sf%GZ(is,iMode)=GZ_n(iMode)
        END IF
        IF((m.NE.0).AND.(n.NE.0))THEN
          !compare G_mn results:,
          !WRITE(*,"(A,I3,X,A,I3,X,2(A,X,E11.4,X),A,E11.4)")'DEBUG m=',m,'n=',n,'GZ_m',GZ_m(iMode),'GZ_n',GZ_n(iMode),'GZ_m - GZ_n=',GZ_m(iMode)-GZ_n(iMode)
        END IF
        END ASSOCIATE !m,n
      END DO
      !write(*,*)'DEBUG ===',is
  
      __PERFOFF('project_G')
  
      !interpolate G and compute Gthet=iota*G+lambda
      ft_IP(:) =LA_IP(:)+iota_int*sf%GZ_base%f%evalDOF_IP(0,sf%GZ(is,:))
      !project onto Gthet base
      CALL sf%GZ_base%f%projectIPtoDOF(.FALSE.,1.0_wp, 0,ft_IP(:),sf%Gthet(is,:))
      DO iMode=1,modes
        sf%Gthet(is,iMode)=sf%Gthet(is,iMode)*dthet_dzeta*sf%GZ_base%f%snorm_base(iMode)
      END DO
  
      __PERFOFF('project')
      CALL ProgressBar(is,nBase)
    END DO !is
    !transform back to corresponding representation of DOF in s
    
    !BCtype_axis(MN_ZERO    )= BC_TYPE_DIRICHLET !=0 (should not be here!)
    !BCtype_axis(M_ZERO     )= BC_TYPE_NEUMANN  ! derivative zero
    !BCtype_axis(M_ODD_FIRST)= BC_TYPE_DIRICHLET !=0
    !BCtype_axis(M_ODD      )= BC_TYPE_DIRICHLET !=0
    !BCtype_axis(M_EVEN     )= BC_TYPE_DIRICHLET !=0
    DEALLOCATE(LA_s)
    !Apply automatic, smooth axis BC:
    BCtype_axis=0
  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode)
    DO iMode=1,modes
      sf%GZ(   :,iMode)= sf%GZ_base%s%initDOF( sf%GZ(   :,iMode) )
      sf%Gthet(:,iMode)= sf%GZ_base%s%initDOF( sf%Gthet(:,iMode) )
      BC_type_axis_imode=BCtype_axis(sf%GZ_base%f%zero_odd_even(iMode))
      IF(BC_type_axis_imode.EQ.0) & !AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
        BC_type_axis_imode=-1*MIN(sf%GZ_base%s%deg+1,sf%GZ_base%f%Xmn(1,iMode))
      CALL sf%GZ_base%s%applyBCtoDOF(sf%GZ(   :,iMode),(/BC_Type_axis_imode,BC_TYPE_OPEN/),(/0.,0./))
      CALL sf%GZ_base%s%applyBCtoDOF(sf%Gthet(:,iMode),(/BC_Type_axis_imode,BC_TYPE_OPEN/),(/0.,0./))
    END DO
  !$OMP END PARALLEL DO 
  
  
    SWRITE(UNIT_StdOut,'(A)') '...DONE.'
    __PERFOFF('get_boozer')
  END SUBROUTINE Get_Boozer_sinterp
  


!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE transform_SFL_free(sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_transform_sfl), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL sf%sgrid_sfl%free()
  IF(ALLOCATED(sf%X1sfl_base))THEN
    CALL sf%X1sfl_base%free()
    DEALLOCATE(sf%X1sfl_base)
  END IF
  IF(ALLOCATED(sf%X2sfl_base))THEN
    CALL sf%X2sfl_base%free()
    DEALLOCATE(sf%X2sfl_base)
  END IF
  IF(ALLOCATED(sf%GZsfl_base))THEN
    CALL sf%GZsfl_base%free()
    DEALLOCATE(sf%GZsfl_base)
  END IF
  IF(ALLOCATED(sf%GZ_base))THEN
    CALL sf%GZ_base%free()
    DEALLOCATE(sf%GZ_base)
  END IF
  SDEALLOCATE(sf%X1sfl)
  SDEALLOCATE(sf%X2sfl)
  SDEALLOCATE(sf%GZsfl)
  SDEALLOCATE(sf%Gthet)
  SDEALLOCATE(sf%GZ)

  sf%initialized=.FALSE.

END SUBROUTINE transform_SFL_free

END MODULE MODgvec_Transform_SFL
