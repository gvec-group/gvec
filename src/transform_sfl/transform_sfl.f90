!===================================================================================================================================
! Copyright (C) 2017 - 2024  Florian Hindenlang <hindenlang@gmail.com>
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
USE MODgvec_Globals, ONLY:wp,abort,MPIroot
USE MODgvec_base   ,ONLY: t_base
USE MODgvec_fbase   ,ONLY: t_fbase
USE MODgvec_sGrid   ,ONLY: t_sgrid
USE MODgvec_hmap,  ONLY: c_hmap
USE MODgvec_SFL_boozer, ONLY: t_sfl_boozer
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
  CLASS(t_base),  ALLOCATABLE :: GZ_base       !! container for base of variable  Gthet and Gzeta (transforms to BOOZER!)
  CLASS(t_base),  ALLOCATABLE :: GZsfl_base    !! container for base of variable Gtheta and Gzeta in SFL coordinates
  REAL(wp),       ALLOCATABLE :: X1sfl(:,:)    !! data (1:nBase,1:modes) of X1 in SFL coords.
  REAL(wp),       ALLOCATABLE :: X2sfl(:,:)    !! data (1:nBase,1:modes) of X2 in SFL coords.
  REAL(wp),       ALLOCATABLE :: Gthet(:,:)    !! data (1:nBase,1:modes) of Gthet in GVEC coords. (for BOOZER)
  REAL(wp),       ALLOCATABLE :: GZ(:,:)       !! data (1:nBase,1:modes) of GZ in GVEC coords. (for BOOZER)
  REAL(wp),       ALLOCATABLE :: Gtsfl(:,:)    !! data (1:nBase,1:modes) of Gt in SFL coords.  (for BOOZER)
  REAL(wp),       ALLOCATABLE :: GZsfl(:,:)    !! data (1:nBase,1:modes) of GZ in SFL coords.  (for BOOZER)
  CLASS(t_sfl_boozer),ALLOCATABLE :: booz      !! subclass needed for boozer transform
  PROCEDURE(i_func_evalprof), POINTER, NOPASS  :: eval_phiPrime
  PROCEDURE(i_func_evalprof), POINTER, NOPASS  :: eval_iota
  CONTAINS
  !PROCEDURE(i_func_evalprof),DEFERRED :: evalphiPrime
  PROCEDURE :: init       => transform_sfl_init
  PROCEDURE :: BuildTransform => BuildTransform_SFL
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




!INTERFACE sfl_boozer_new
!  MODULE PROCEDURE sfl_boozer_new
!END INTERFACE


PUBLIC :: t_transform_sfl,transform_sfl_new, find_pest_angles
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Allocate class and call init
!!
!===================================================================================================================================
SUBROUTINE transform_sfl_new(sf,mn_max_in, whichSFL,deg_in,continuity_in,degGP_in,grid_in,  &
                             hmap_in,X1_base_in,X2_base_in,LA_base_in,eval_phiPrime_in,eval_iota_in)
  ! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER     ,INTENT(IN) :: mn_max_in(2)                                        !< maximum number for new variables in SFL coordinates
  INTEGER     ,INTENT(IN) :: whichSFL                                         !< either =1: PEST, =2:Boozer 
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

! extended base for q in the new angles, and on the new grid
CALL base_new(sf%X1sfl_base,  sf%deg, sf%continuity, sf%sgrid_sfl, sf%degGP,      &
               sf%mn_max,2*sf%mn_max+1,sf%nfp,sin_cos_map(sf%X1sfl_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform
CALL base_new(sf%X2sfl_base,   sf%deg, sf%continuity, sf%sgrid_sfl,sf%degGP,      &
              sf%mn_max,2*sf%mn_max+1,sf%nfp,sin_cos_map(sf%X2sfl_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform
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
  sf%mn_max,2*sf%mn_max+1,sf%nfp,sin_cos_map(sf%GZ_sin_cos), .FALSE.)!m=n=0 should be always there, because of coordinate transform

  ALLOCATE(rho_pos(1:sf%GZsfl_base%s%nBase),iota(1:sf%GZsfl_base%s%nBase),phiPrime(1:sf%GZsfl_base%s%nBase))
  DO irho=1,sf%GZsfl_base%s%nBase
    rho_pos(irho)=MIN(MAX(1.0e-4_wp,sf%GZsfl_base%s%s_IP(irho)),1.0_wp-1.0e-12_wp)
    iota(irho)=sf%eval_iota(rho_pos(irho))  
    phiPrime(irho)=sf%eval_phiPrime(rho_pos(irho))
  END DO
  CALL sfl_boozer_new(sf%booz,sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%GZ_sin_cos),sf%hmap,sf%GZsfl_base%s%nBase,rho_pos,iota,phiPrime)
  DEALLOCATE(rho_pos,iota,phiPrime)
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
nrho = sf%X1sfl_base%s%nBase
mnIP = sf%X1sfl_base%f%mn_IP
ALLOCATE(thetzeta_trafo(2,mnIP,nrho))
SELECT CASE(sf%whichSFLcoord)
CASE(1) !PEST
  ALLOCATE(Gt_rho(LA_base_in%f%modes,nrho))
  DO irho=1,nrho
    spos=MIN(MAX(1.0e-4_wp,sf%X1sfl_base%s%s_IP(irho)),1.0_wp-1.0e-12_wp)
    Gt_rho(:,irho)=LA_base_in%s%evalDOF2D_s(spos,LA_base_in%f%modes,0,LA_in(:,:))
  END DO
  CALL find_pest_angles(nrho,LA_base_in%f,Gt_rho,sf%X1sfl_base%f%mn_IP,sf%X1sfl_base%f%x_IP,thetzeta_trafo)
  DEALLOCATE(Gt_rho)
  CALL transform_Angles_3d(X1_base_in,X1_in,"X1sfl",sf%X1sfl_base,sf%X1sfl,thetzeta_trafo)
  CALL transform_Angles_3d(X2_base_in,X2_in,"X2sfl",sf%X2sfl_base,sf%X2sfl,thetzeta_trafo)
  !TODO INTERPOLATE SPLINE IN RHO AND APPLY BC!!! -> PUT ALL THIS INTO NEW transform_angles routine using full interpolation
WRITE(*,*)'TEST_rootsearch, PEST:'
CASE(2) !BOOZER
  CALL sf%booz%Get_Boozer(X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in) ! fill sf%booz%lambda,sf%booz%nu
  !!!!
  !INIT Gthet,GZ as splines... needed?
  CALL find_boozer_Angles(nrho,sf%booz%iota,sf%booz%nu_fbase,sf%booz%lambda(:,:),sf%booz%nu(:,:), &
       sf%X1sfl_base%f%mn_IP,sf%X1sfl_base%f%x_IP, thetzeta_trafo)
  ALLOCATE(Gt_rho(sf%GZ_base%f%modes,nrho))
  DO irho=1,nrho
    Gt_rho(:,irho)=sf%booz%lambda(:,irho)+sf%booz%iota(irho)*sf%booz%nu(:,irho)
  END DO
  CALL to_spline_with_BC(sf%GZ_base,Gt_rho,sf%Gthet)
  CALL to_spline_with_BC(sf%GZ_base,sf%booz%nu,sf%GZ)
  DEALLOCATE(Gt_rho)
  CALL transform_Angles_3d(X1_base_in,X1_in   ,"X1sfl",sf%X1sfl_base,sf%X1sfl,thetzeta_trafo)
  CALL transform_Angles_3d(X2_base_in,X2_in   ,"X2sfl",sf%X2sfl_base,sf%X2sfl,thetzeta_trafo)
  CALL Transform_Angles_3d(sf%GZ_base,sf%Gthet,"Gtsfl",sf%GZsfl_base,sf%Gtsfl,thetzeta_trafo)
  CALL Transform_Angles_3d(sf%GZ_base,sf%GZ   ,"GZsfl",sf%GZsfl_base,sf%GZsfl,thetzeta_trafo)
END SELECT

DEALLOCATE(thetzeta_trafo)

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
    q_m(:,irho)=q_base_out%f%initDOF(q_IP(:))
    CALL ProgressBar(irho,nBase)
  END DO !is

  CALL to_spline_with_BC(q_base_out,q_m,q_out)

  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
  __PERFOFF('transform_angles')
END SUBROUTINE Transform_Angles_3d

!===================================================================================================================================
!> Helper routine to go from spline interpolation points to spline coefficients and apply axis boundary condition
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
  REAL(wp) ,INTENT(INOUT) :: q_out(q_base_out%s%nBase,1:q_base_out%f%modes)          !< spline/fourier coefficients of q in new angles 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iMode
  INTEGER               :: BCtype_axis(0:4)
!===================================================================================================================================
  BCtype_axis(MN_ZERO    )= BC_TYPE_OPEN   !do nothing
  BCtype_axis(M_ZERO     )= BC_TYPE_OPEN   !do nothing
  BCtype_axis(M_ODD_FIRST)= BC_TYPE_DIRICHLET !=0
  BCtype_axis(M_ODD      )= BC_TYPE_DIRICHLET !=0
  BCtype_axis(M_EVEN     )= BC_TYPE_DIRICHLET !=0
  !transform back to corresponding representation of DOF in s
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode)
  DO iMode=1,q_base_out%f%modes
    q_out(:,iMode)=q_base_out%s%initDOF( q_m(iMode,:) )
    CALL q_base_out%s%applyBCtoDOF(q_out(:,iMode),(/BCtype_axis(q_base_out%f%zero_odd_even(iMode)),BC_TYPE_OPEN/)  &
                                              ,(/0.,0./))
  END DO
!$OMP END PARALLEL DO 
END SUBROUTINE to_spline_with_BC

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
  USE MODgvec_Globals,ONLY: UNIT_stdOut
  USE MODgvec_fbase  ,ONLY: t_fbase
  USE MODgvec_Newton ,ONLY: NewtonRoot2D
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER      ,INTENT(IN) :: nrho   !! number of surfaces, (second dimension  of LA_in and nu_in modes)
    CLASS(t_fBase),INTENT(IN) ::fbase_in     !< same basis of lambda and nu
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
    LOGICAL    :: docheck=.TRUE.
  !===================================================================================================================================
   __PERFON('find_pest_angles')
  SWRITE(UNIT_StdOut,'(A)') 'FIND PEST ANGLES VIA NEWTON'
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(irho,j,theta_star,zeta)
  DO irho=1,nrho
    DO j=1,tz_dim
      theta_star=tz_pest(1,j); zeta=tz_pest(2,j)          
      thetzeta_out(1,j,irho)=get_pest_newton(theta_star,zeta,fbase_in,LA_in(:,irho))
      thetzeta_out(2,j,irho)=zeta
    END DO! j
  END DO !irho
!$OMP END PARALLEL DO
  IF(docheck)THEN
    DO irho=1,nrho
      check=fbase_in%evalDOF_xn(tz_dim,thetzeta_out(:,:,irho),0,LA_in(:,irho))
      maxerr(irho)=maxval(abs(check+(thetzeta_out(1,:,irho)-tz_pest(1,:))))
    END DO
    WRITE(*,*)'CHECK PEST THETA*',MAXVAL(maxerr)
    IF(ANY(maxerr.GT.1.0e-12))THEN
      CALL abort(__STAMP__, & 
          "Find_pest_Angles: Error in theta*")
    END IF
  END IF
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
    CLASS(t_fBase),INTENT(IN) ::LA_fbase_in     !<  basis of lambda
    REAL(wp)     ,INTENT(IN) :: LA_in(1:LA_fbase_in%modes) !< fourier coefficients of thet*=thet+LA(theta,zeta)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES      
    REAL(wp)              :: thet_out !< theta position in original coordinates
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  !===================================================================================================================================
    thet_out=NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1_wp*PI, &
                                         theta_star, theta_star,A_FRdFR) !start, rhs,func
  CONTAINS 
!for newton root search 
  FUNCTION A_FRdFR(theta_iter)
    !uses current zeta where newton is called, and A from subroutine above
    IMPLICIT NONE
    REAL(wp) :: theta_iter
    REAL(wp) :: A_FRdFR(2) !output function and derivative
    !--------------------------------------------------- 
    A_FRdFR(1)=theta_iter+LA_fbase_in%evalDOF_x((/theta_iter,zeta/),         0,LA_in(:))  !theta_iter+lambda = thet* (rhs)
    A_FRdFR(2)=1.0_wp    +LA_fbase_in%evalDOF_x((/theta_iter,zeta/),DERIV_THET,LA_in(:)) !1+dlambda/dtheta  
  END FUNCTION A_FRdFR
  
END FUNCTION get_pest_newton


  


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
