!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
USE MODgvec_Globals, ONLY:wp,abort,MPIroot
USE MODgvec_base   ,ONLY: t_base
USE MODgvec_fbase   ,ONLY: t_fbase
USE MODgvec_sGrid   ,ONLY: t_sgrid
USE MODgvec_SFL_boozer, ONLY: t_sfl_boozer
USE MODgvec_hmap,  ONLY: PP_T_HMAP
USE MODgvec_newton,  ONLY: c_newton_Root1D_FdF
IMPLICIT NONE
PRIVATE

TYPE :: t_transform_sfl
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER                     :: whichSFLcoord !!
  INTEGER                     :: fac_nyq,mn_max(2),deg,continuity,degGP,nfp
  INTEGER                     :: mn_nyq(2)     !! number of integration points in fbase for X1sfl,X2sfl,Gtsfl,GZsfl
  INTEGER                     :: mn_nyq_booz(2)!! number of integration points in fbase for boozer transform (fbase of GZ_base)
  INTEGER                     :: X1sfl_sin_cos,X2sfl_sin_cos,GZ_sin_cos
  LOGICAL                     :: booz_relambda !! =T: recompute lambda from mapping on full fourier series, =F: use LA from eq.
  INTEGER                     :: to_angle_method !! =1: "interpolate": root search of interpolation points (mn_nyq=2*(m_max,n_max)+1) in sfl angles,
                                                 !! =2: "integrate": project to sfl angles with trapezoidal rule (integral with transform, mn_nyq=4*(m_max,n_max)+1)
  TYPE(t_sgrid)               :: sgrid_sfl     !! grid for SFL coordinates
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP),  POINTER     :: hmap          !! pointer to hmap class
#else
  CLASS(PP_T_HMAP),  POINTER    :: hmap          !! pointer to hmap class
#endif

  CLASS(t_base),  ALLOCATABLE :: X1sfl_base    !! container for base of variable X1 in SFL coordinates
  CLASS(t_base),  ALLOCATABLE :: X2sfl_base    !! container for base of variable X2 in SFL coordinates
  CLASS(t_base),  ALLOCATABLE :: GZ_base       !! container for base of variable  Gthet and Gzeta (transforms to BOOZER!)
  CLASS(t_base),  ALLOCATABLE :: GZsfl_base    !! container for base of variable Gtheta and Gzeta in SFL coordinates
  REAL(wp),       ALLOCATABLE :: X1sfl(:,:)    !! data (1:nBase,1:modes) of X1 in SFL coords.
  REAL(wp),       ALLOCATABLE :: X2sfl(:,:)    !! data (1:nBase,1:modes) of X2 in SFL coords.
  REAL(wp),       ALLOCATABLE :: Gthet(:,:)    !! data (1:nBase,1:modes) of Gthet in GVEC coords. (for BOOZER)
  REAL(wp),       ALLOCATABLE :: GZ(:,:)       !! data (1:nBase,1:modes) of GZ in GVEC coords. (for BOOZER)
  REAL(wp),       ALLOCATABLE :: Gtsfl(:,:)    !! data (1:nBase,1:modes) of Gt in SFL coords.  (for BOOZER)
  REAL(wp),       ALLOCATABLE :: GZsfl(:,:)    !! data (1:nBase,1:modes) of GZ in SFL coords.  (for BOOZER)
  TYPE(t_sfl_boozer),ALLOCATABLE :: booz      !! subclass needed for boozer transform
  PROCEDURE(i_func_evalprof), POINTER, NOPASS  :: eval_phiPrime
  PROCEDURE(i_func_evalprof), POINTER, NOPASS  :: eval_iota
  CONTAINS
  !PROCEDURE(i_func_evalprof),DEFERRED :: evalphiPrime
  PROCEDURE :: init       => transform_sfl_init
  PROCEDURE :: BuildTransform => BuildTransform_SFL
  PROCEDURE :: free        => transform_sfl_free
END TYPE t_transform_sfl

TYPE, EXTENDS(c_newton_Root1D_FdF) :: t_newton_Root1D_FdF_pest
  TYPE(t_fbase), POINTER :: LA_fbase_in
  REAL(wp), POINTER :: LA_in(:)
  REAL(wp) :: zeta
  CONTAINS
  PROCEDURE :: FRdFR => pest_newton_FRdFR
END TYPE t_newton_Root1D_FdF_pest

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

!INTERFACE sfl_boozer_new
!  MODULE PROCEDURE sfl_boozer_new
!END INTERFACE


PUBLIC :: t_transform_sfl,transform_sfl_new, find_pest_angles, get_pest_newton
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Allocate class and call init
!!
!===================================================================================================================================
SUBROUTINE transform_sfl_new(sf,mn_max_in, whichSFL,deg_in,continuity_in,degGP_in,grid_in,  &
                             hmap_in,X1_base_in,X2_base_in,LA_base_in,eval_phiPrime_in,eval_iota_in,booz_relambda)
  ! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER     ,INTENT(IN) :: mn_max_in(2)                                        !< maximum number for new variables in SFL coordinates
  INTEGER     ,INTENT(IN) :: whichSFL                                         !< either =1: PEST, =2:Boozer
  INTEGER     ,INTENT(IN) :: deg_in,continuity_in,degGP_in                    !< for output base (X1,X2,G)
  CLASS(t_sgrid),INTENT(IN),TARGET :: grid_in                            !! grid information
  CLASS(t_base),INTENT(IN),TARGET :: X1_base_in,X2_base_in,LA_base_in           !< base classes belong to solution U_in
#ifdef PP_WHICH_HMAP
  TYPE( PP_T_HMAP),INTENT(IN),TARGET :: hmap_in
#else
  CLASS(PP_T_HMAP),INTENT(IN),TARGET :: hmap_in
#endif
  PROCEDURE(i_func_evalprof)     :: eval_phiPrime_in,eval_iota_in  !!procedure pointers to profile evaluation functions.
  LOGICAL,INTENT(IN),OPTIONAL :: booz_relambda !! for boozer transform, recompute lambda (recommended)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  TYPE(t_transform_sfl), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
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
  IF(PRESENT(booz_relambda))THEN
    sf%booz_relambda=booz_relambda
  ELSE
    sf%booz_relambda=.TRUE. !relambda !if =True, J^s=0 will be recomputed, for exact integrability condition of boozer transform  (but slower!)
  END IF
  sf%fac_nyq=4  !hard coded for now
  sf%mn_nyq_booz(1:2)=sf%fac_nyq*MAXVAL(sf%mn_max)+1  ! for boozer transform
  ! use maximum number of integration points from maximum mode number in both directions
  sf%to_angle_method=1 !HARD CODED TO NEWTON ROOT SEARCH for interpolation points, faster than projection
  SELECT CASE(sf%to_angle_method)
  CASE(1) !INTERPOLATION
    sf%mn_nyq(1:2)=2*sf%mn_max+1  !only interpolation (=Fourier transform)
  CASE(2) !PROJECTION (previous way)
    sf%mn_nyq(1:2)=sf%mn_nyq_booz !projection with trapezoidal rule
  END SELECT
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
USE MODgvec_SFL_Boozer,ONLY: sfl_boozer_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_transform_sfl), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: irho
REAL(wp),ALLOCATABLE :: rho_pos(:),iota(:),phiPrime(:)
!===================================================================================================================================
__PERFON('transform_SFL_init')
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
  sf%mn_max,sf%mn_nyq_booz,sf%nfp,sin_cos_map(sf%GZ_sin_cos),.TRUE.) !exclude m=n=0

  ALLOCATE(sf%Gthet(sf%GZ_base%s%nBase,sf%GZ_base%f%modes)); sf%Gthet=0.0_wp
  ALLOCATE(sf%GZ(   sf%GZ_base%s%nBase,sf%GZ_base%f%modes)); sf%GZ=0.0_wp

  CALL base_new(sf%GZsfl_base, sf%deg, sf%continuity, sf%sgrid_sfl,sf%degGP,      &
  sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%GZ_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform

  ALLOCATE(rho_pos(1:sf%GZsfl_base%s%nBase),iota(1:sf%GZsfl_base%s%nBase),phiPrime(1:sf%GZsfl_base%s%nBase))
  DO irho=1,sf%GZsfl_base%s%nBase
    rho_pos(irho)=MIN(MAX(1.0e-4_wp,sf%GZsfl_base%s%s_IP(irho)),1.0_wp-1.0e-12_wp)
    iota(irho)=sf%eval_iota(rho_pos(irho))
    phiPrime(irho)=sf%eval_phiPrime(rho_pos(irho))
  END DO
  CALL sfl_boozer_new(sf%booz,sf%mn_max,sf%mn_nyq_booz,sf%nfp,sin_cos_map(sf%GZ_sin_cos),sf%hmap,sf%GZsfl_base%s%nBase, &
                      rho_pos,iota,phiPrime,relambda_in=sf%booz_relambda)
  DEALLOCATE(rho_pos,iota,phiPrime)
  ALLOCATE(sf%Gtsfl(sf%GZsfl_base%s%nBase,sf%GZsfl_base%f%modes));sf%Gtsfl=0.0_wp
  ALLOCATE(sf%GZsfl(sf%GZsfl_base%s%nBase,sf%GZsfl_base%f%modes));sf%GZsfl=0.0_wp

CASE DEFAULT
  CALL abort(__STAMP__, &
           "whichSFLcoord for coordinate transform not found, expecting 1(PEST) or 2(Boozer)", &
           TypeInfo="InvalidParameterError")
END SELECT
sf%initialized=.TRUE.
__PERFOFF('transform_SFL_init')
END SUBROUTINE transform_sfl_init


!===================================================================================================================================
!> Builds X1 and X2 in SFL coordinates
!!
!===================================================================================================================================
SUBROUTINE BuildTransform_SFL(sf,X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in)
! MODULES
USE MODgvec_base   ,ONLY: t_base,base_new
USE MODgvec_SFL_Boozer,ONLY: find_boozer_Angles
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
  INTEGER :: irho,nrho,mnIP
  REAL(wp),ALLOCATABLE :: Gt_rho(:,:),thetzeta_trafo(:,:,:)
  REAL(wp):: spos
!===================================================================================================================================
__PERFON('BuildTransform_SFL')
nrho = sf%X1sfl_base%s%nBase
mnIP = sf%X1sfl_base%f%mn_IP
ALLOCATE(thetzeta_trafo(2,mnIP,nrho))
SELECT CASE(sf%whichSFLcoord)
CASE(1) !PEST
  !interpolate lambda at rho positions (for find_pest_angles)

  SELECT CASE(sf%to_angle_method)
  CASE(1)
    ALLOCATE(Gt_rho(LA_base_in%f%modes,nrho))
    DO irho=1,nrho
      spos=MIN(MAX(1.0e-4_wp,sf%X1sfl_base%s%s_IP(irho)),1.0_wp-1.0e-12_wp)
      Gt_rho(:,irho)=LA_base_in%s%evalDOF2D_s(spos,LA_base_in%f%modes,0,LA_in(:,:))
    END DO
    CALL find_pest_angles(nrho,LA_base_in%f,Gt_rho,sf%X1sfl_base%f%mn_IP,sf%X1sfl_base%f%x_IP,thetzeta_trafo)
    DEALLOCATE(Gt_rho)
    CALL transform_Angles_3d(X1_base_in,X1_in,"X1sfl",sf%X1sfl_base,sf%X1sfl,thetzeta_trafo)
    CALL transform_Angles_3d(X2_base_in,X2_in,"X2sfl",sf%X2sfl_base,sf%X2sfl,thetzeta_trafo)
  CASE(2)
    CALL Transform_Angles_sinterp(LA_base_in,LA_in,X1_base_in,X1_in,"X1",sf%X1sfl_base,sf%X1sfl)
    CALL Transform_Angles_sinterp(LA_base_in,LA_in,X2_base_in,X2_in,"X2",sf%X2sfl_base,sf%X2sfl)
  END SELECT
CASE(2) !BOOZER
  CALL sf%booz%Get_Boozer(X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in) ! fill sf%booz%lambda,sf%booz%nu for find_angles
  ALLOCATE(Gt_rho(sf%GZ_base%f%modes,nrho))
  DO irho=1,nrho
    Gt_rho(:,irho)=sf%booz%lambda(:,irho)+sf%booz%iota(irho)*sf%booz%nu(:,irho)
  END DO
  CALL to_spline_with_BC(sf%GZ_base,Gt_rho,sf%Gthet)
  CALL to_spline_with_BC(sf%GZ_base,sf%booz%nu,sf%GZ)
  DEALLOCATE(Gt_rho)
  SELECT CASE(sf%to_angle_method)
  CASE(1)
    CALL sf%booz%find_angles(sf%X1sfl_base%f%mn_IP,sf%X1sfl_base%f%x_IP, thetzeta_trafo)
    CALL transform_Angles_3d(X1_base_in,X1_in   ,"X1sfl",sf%X1sfl_base,sf%X1sfl,thetzeta_trafo)
    CALL transform_Angles_3d(X2_base_in,X2_in   ,"X2sfl",sf%X2sfl_base,sf%X2sfl,thetzeta_trafo)
    CALL Transform_Angles_3d(sf%GZ_base,sf%Gthet,"Gtsfl",sf%GZsfl_base,sf%Gtsfl,thetzeta_trafo)
    CALL Transform_Angles_3d(sf%GZ_base,sf%GZ   ,"GZsfl",sf%GZsfl_base,sf%GZsfl,thetzeta_trafo)
  CASE(2)
    CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,sf%GZ_base,sf%GZ,"GZ",sf%GZsfl_base,sf%GZsfl,B_in=sf%GZ)
    CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,sf%GZ_base,sf%Gthet,"Gt",sf%Gzsfl_base,sf%Gtsfl,B_in=sf%GZ)
    CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,X1_base_in,X1_in,"X1",sf%X1sfl_base,sf%X1sfl,B_in=sf%GZ)
    CALL Transform_Angles_sinterp(sf%GZ_base,sf%Gthet,X2_base_in,X2_in,"X2",sf%X2sfl_base,sf%X2sfl,B_in=sf%GZ)
  END SELECT

END SELECT

DEALLOCATE(thetzeta_trafo)
__PERFOFF('BuildTransform_SFL')
END SUBROUTINE BuildTransform_SFL


!===================================================================================================================================
!> Transform a function from the GVEC angles q(s,theta,zeta) to new angles q*(s,theta*,zeta*)
!! by using interpolation in angular direction (fourier transform)
!! and spline interpolation in radial direction (at s_IP points of output base)
!! the interpolation points are given by thetazeta_IP,
!! which are the angle positions of an equidistant interpolation grid in PEST/Boozer angles
!!
!===================================================================================================================================
SUBROUTINE Transform_Angles_3d(q_base_in,q_in,q_name,q_base_out,q_out,thetazeta_IP)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,Progressbar
USE MODgvec_base   ,ONLY: t_base
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_Base),INTENT(IN) :: q_base_in                                     !< basis of function f
  REAL(wp)     ,INTENT(IN) :: q_in(1:q_base_in%s%nBase,1:q_base_in%f%modes) !< coefficients of f
  CHARACTER(LEN=*),INTENT(IN):: q_name
  CLASS(t_base),INTENT(IN) :: q_base_out                                    !< new fourier basis of function q in new angles, defined mn_max,mn_nyq
  REAL(wp)     ,INTENT(IN) :: thetazeta_IP(2,q_base_out%f%mn_IP,q_base_out%s%nBase) !< theta zeta evaluation points corresponding to equispaced integration points in boozer/pest angles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) ,INTENT(INOUT) :: q_out(q_base_out%s%nBase,1:q_base_out%f%modes)          !< spline/fourier coefficients of q in new angles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,irho,mn_IP,mn_max(2),mn_nyq(2)
  REAL(wp)              :: spos
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes),q_IP(q_base_out%f%mn_IP)
  REAL(wp)              :: q_m(q_base_out%f%modes,q_base_out%s%nBase) !output fourier on radial interpolations
!===================================================================================================================================
  mn_max(1:2) =q_base_out%f%mn_max
  mn_nyq(1:2) =q_base_out%f%mn_nyq
  SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'TRANSFORM '//TRIM(q_name)//' TO NEW ANGLE COORDINATES, nfp=',q_base_in%f%nfp, &
                              ', mn_max_in=',q_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq

  __PERFON('transform_angles')
  !total number of integration points
  mn_IP = q_base_out%f%mn_IP
  nBase = q_base_out%s%nBase

  CALL ProgressBar(0,nBase)!init

  DO irho=1,nBase
    spos=q_base_out%s%s_IP(irho) !interpolation points for q_in
    !evaluate q_in at spos
    q_in_s(:)  = q_base_in%s%evalDOF2D_s(spos,q_base_in%f%modes,   0,q_in(:,:))
    q_IP       = q_base_in%f%evalDOF_xn(mn_IP,thetazeta_IP(:,:,irho),0, q_in_s(:))
    q_m(:,irho)= q_base_out%f%initDOF(q_IP(:))
    CALL ProgressBar(irho,nBase)
  END DO !is

  CALL to_spline_with_BC(q_base_out,q_m,q_out)

  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
  __PERFOFF('transform_angles')
END SUBROUTINE Transform_Angles_3d

!===================================================================================================================================
!> Helper routine to go from spline interpolation points to spline coefficients and apply smooth axis boundary condition.
!!
!===================================================================================================================================
SUBROUTINE to_spline_with_BC(q_base_out,q_m,q_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_Base),INTENT(IN) :: q_base_out                                  !< basis of function f
  REAL(wp)     ,INTENT(IN) :: q_m(1:q_base_out%f%modes,1:q_base_out%s%nBase) !< coefficients of f, sampled at s_IP points
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) ,INTENT(OUT)    :: q_out(q_base_out%s%nBase,1:q_base_out%f%modes)          !< spline/fourier coefficients of q in new angles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iMode
  INTEGER               :: BCtype_axis
  REAL(wp)              :: BCval_axis,BCval_edge
!===================================================================================================================================
  !transform back to corresponding representation of DOF in s
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,BCtype_axis,BCval_axis,BCval_edge)
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_out%s%initDOF( q_m(iMode,:) )
    BCval_edge=q_m(iMode,q_base_out%s%nBase)
    !NEW m-dependent smooth BC at axis, for m>deg, switch off all DOF up to deg
    BCtype_axis=-1*MIN(q_base_out%s%deg+1,q_base_out%f%Xmn(1,iMode)) !
    SELECT CASE(q_base_out%f%zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      BCval_axis=q_m(iMode,1)
    CASE DEFAULT
      BCval_axis=0.0_wp
    END SELECT
    CALL q_base_out%s%applyBCtoDOF(q_out(:,iMode),(/BCtype_axis,BC_TYPE_DIRICHLET/),(/BCval_axis,BCval_edge/))
  END DO
!$OMP END PARALLEL DO
END SUBROUTINE to_spline_with_BC

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
USE MODgvec_Globals,ONLY: UNIT_stdOut,Progressbar,testlevel
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
  INTEGER               :: nBase,i,is,i_mn,mn_IP,mn_max(2),mn_nyq(2)
  REAL(wp)              :: spos,dthet_dzeta
  REAL(wp)              :: check(1:7,q_base_out%s%nBase)
  LOGICAL               :: docheck
  LOGICAL               :: Bpresent
  REAL(wp)              :: A_s(1:AB_base_in%f%modes)
  REAL(wp)              :: B_s(1:AB_base_in%f%modes)
  REAL(wp)              :: q_in_s(1:q_base_in%f%modes)
  REAL(wp), ALLOCATABLE :: q_IP(:),q_m(:,:)   ! q evaluated at spos and all integration points
  REAL(wp), ALLOCATABLE :: f_IP(:)       ! =q*(1+dlambda/dtheta) evaluated at integration points
  REAL(wp), ALLOCATABLE :: modes_IP(:,:) ! mn modes of q evaluated at theta*,zeta* for all integration points
  TYPE(t_fBase),ALLOCATABLE        :: q_fbase_nyq
  TYPE(t_fBase),ALLOCATABLE        :: AB_fbase_nyq
  REAL(wp),DIMENSION(:),ALLOCATABLE :: A_IP,dAdthet_IP,B_IP,dBdthet_IP,dBdzeta_IP,dAdzeta_IP
!===================================================================================================================================
  docheck=(testlevel.GT.0)
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
    spos=MIN(MAX(1.0e-4_wp,q_base_out%s%s_IP(is)),1.0_wp-1.0e-12_wp) !interpolation points for q_in
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

    !CHECK at all IP points
    IF(doCheck) THEN
      __PERFON('check')
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
      __PERFOFF('check')
    END IF !doCheck

    CALL ProgressBar(is,nBase)
  END DO !is

  CALL to_spline_with_BC(q_base_out,q_m,q_out)

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
!> on one flux surface, find for a list of in thet*_j,zeta*_j, the corresponding (thet_j,zeta_j) positions, given
!> Here, new PEST angles are
!> theta*=theta+lambda(theta,zeta)
!>  zeta*=zeta,
!> so a 1D root search in theta is is enough
!!
!===================================================================================================================================
SUBROUTINE find_pest_angles(nrho,fbase_in,LA_in,tz_dim,tz_pest,thetzeta_out)
  ! MODULES
  USE MODgvec_Globals,ONLY: UNIT_stdOut,ProgressBar,testlevel
  USE MODgvec_fbase  ,ONLY: t_fbase
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER      ,INTENT(IN) :: nrho   !! number of surfaces, (second dimension  of LA_in and nu_in modes)
    TYPE(t_fbase),INTENT(IN) :: fbase_in     !< same basis of lambda and nu
    REAL(wp)     ,INTENT(IN) :: LA_in(1:fbase_in%modes,nrho) !< fourier coefficients of thet*=thet+LA(theta,zeta)+iota*nu(theta,zeta)
    INTEGER      ,INTENT(IN) :: tz_dim                 !< size of the list in thetstar,zetastar
    REAL(wp)     ,INTENT(IN) :: tz_pest(2,tz_dim) !< theta,zeta positions in pest angle (same for all rho)

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)    ,INTENT(OUT) :: thetzeta_out(2,tz_dim,nrho)  !! theta,zeta position in original angles, for given pest angles
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER    :: irho,j
    REAL(wp)   :: theta_star,zeta
    REAL(wp)   :: check(tz_dim),maxerr(nrho)
    LOGICAL    :: docheck
  !===================================================================================================================================
   __PERFON('find_pest_angles')
   docheck=(testlevel.GT.0)
   SWRITE(UNIT_StdOut,'(A,2(I8,A))')'Find pest angles via Newton on  nrho=',nrho,' times ntheta_zeta= ',tz_dim, " points"
  CALL ProgressBar(0,nrho)!init
  DO irho=1,nrho
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP   PRIVATE(j,theta_star,zeta) FIRSTPRIVATE(irho) &
!$OMP   SHARED(tz_dim,tz_pest,thetzeta_out,fbase_in,LA_in)
    DO j=1,tz_dim
      theta_star=tz_pest(1,j); zeta=tz_pest(2,j)
      thetzeta_out(1,j,irho)=get_pest_newton(theta_star,zeta,fbase_in,LA_in(:,irho))
      thetzeta_out(2,j,irho)=zeta
    END DO! j
!$OMP END PARALLEL DO
    CALL ProgressBar(irho,nrho)
  END DO !irho

  IF(docheck)THEN
    __PERFON("check")
    DO irho=1,nrho
      check=fbase_in%evalDOF_xn(tz_dim,thetzeta_out(:,:,irho),0,LA_in(:,irho))
      maxerr(irho)=maxval(abs(check+(thetzeta_out(1,:,irho)-tz_pest(1,:))))
    END DO

    IF(ANY(maxerr.GT.1.0e-12))THEN
      WRITE(UNIT_stdout,*)'CHECK PEST THETA*',maxerr
      CALL abort(__STAMP__, &
          "Find_pest_Angles: Error in theta*")
    END IF
    __PERFOFF("check")
  END IF! docheck
  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
  __PERFOFF('find_pest_angles')

END SUBROUTINE find_pest_angles


!===================================================================================================================================
!> This function returns the result of the 1D newton root search for the pest theta angle
!!
!===================================================================================================================================
FUNCTION get_pest_newton(theta_star,zeta,LA_fbase_in,LA_in) RESULT(thet_out)
  USE MODgvec_fbase  ,ONLY: t_fbase
  USE MODgvec_Newton ,ONLY: NewtonRoot1D_FdF
  USE MODgvec_Globals,ONLY: PI
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp)     ,INTENT(IN) :: theta_star !< initial guess = thet*
    REAL(wp)     ,INTENT(IN) :: zeta
    TYPE(t_fbase),INTENT(IN), TARGET ::LA_fbase_in     !<  basis of lambda
    REAL(wp)     ,INTENT(IN), TARGET :: LA_in(1:LA_fbase_in%modes) !< fourier coefficients of thet*=thet+LA(theta,zeta)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)              :: thet_out !< theta position in original coordinates
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    TYPE(t_newton_Root1D_FdF_pest) :: fobj
  !===================================================================================================================================
    fobj%zeta = zeta
    fobj%LA_fbase_in => LA_fbase_in
    fobj%LA_in => LA_in

    thet_out = NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1_wp*PI, &
                                theta_star, theta_star,fobj) !start, rhs,func
END FUNCTION get_pest_newton


!===================================================================================================================================
!> Function for 1D newton root search for PEST
!!
!===================================================================================================================================
FUNCTION pest_newton_FRdFR(sf, x) RESULT(A_FRdFR)
  !uses current zeta where newton is called, and A from subroutine above
  CLASS(t_newton_Root1D_FdF_pest), INTENT(IN) :: sf
  REAL(wp), INTENT(IN) :: x ! theta_iter
  REAL(wp) :: A_FRdFR(2) !output function and derivative
  !---------------------------------------------------
  A_FRdFR(1)=x      + sf%LA_fbase_in%evalDOF_x((/x,sf%zeta/),         0,sf%LA_in) !theta_iter+lambda = thet* (rhs)
  A_FRdFR(2)=1.0_wp + sf%LA_fbase_in%evalDOF_x((/x,sf%zeta/),DERIV_THET,sf%LA_in) !1+dlambda/dtheta
END FUNCTION pest_newton_FRdFR


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
