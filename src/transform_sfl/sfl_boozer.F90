!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"
!===================================================================================================================================
!>
!!# Module **SFL boozer**
!!
!! Transform to Straight-field line  BOOZER coordinates
!!
!===================================================================================================================================
MODULE MODgvec_SFL_Boozer
! MODULES
USE MODgvec_Globals, ONLY:wp,abort,MPIroot
USE MODgvec_fbase   ,ONLY: t_fbase
USE MODgvec_hmap,  ONLY: PP_T_HMAP,PP_T_HMAP_AUXVAR
USE MODgvec_Newton, ONLY: c_newton_Root2D
IMPLICIT NONE
PRIVATE

!===================================================================================================================================
!> Class for the computation of the boozer transform, on a given radial grid positions rho, with iota(rho) phiPrime(rho) values.
!!
!! theta^Boozer=theta+lambda+iota(rho)*nu(rho,theta,zeta)
!! zeta^Boozer =zeta                  +nu(rho,theta,zeta)
!===================================================================================================================================
TYPE :: t_sfl_boozer
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER  :: nrho       !! number of rho positions
  LOGICAL  :: relambda   !! if =True, J^s=0 will be recomputed, for exact integrability condition of boozer transform  (but slower!)
  TYPE(t_fbase), ALLOCATABLE :: nu_fbase

  REAL(wp),ALLOCATABLE::rho_pos(:),iota(:),phiPrime(:) !! rho positions, iota and phiPrime at these rho positions
  ! computed in the boozer transform
  REAL(wp),ALLOCATABLE::lambda(:,:),nu(:,:)   !! Fourier modes for all rho positions of lambda (recomputed on the fourier space of nu) and nu for boozer transform , (iMode,irho)
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP),  POINTER     :: hmap          !! pointer to hmap class
  TYPE(PP_T_HMAP_AUXVAR),ALLOCATABLE   :: hmap_xv(:) !! auxiliary variables for hmap
#else
  CLASS(c_hmap),  POINTER     :: hmap          !! pointer to hmap class
  CLASS(c_hmap_auxvar),ALLOCATABLE   :: hmap_xv(:) !! auxiliary variables for hmap
#endif
  CONTAINS
  PROCEDURE :: get_boozer  => get_boozer_sinterp
  PROCEDURE :: free        => sfl_boozer_free
  PROCEDURE :: find_angles => self_find_boozer_angles
  PROCEDURE :: find_angles_irho => self_find_boozer_angles_irho
END TYPE t_sfl_boozer

TYPE, EXTENDS(c_newton_Root2D) :: t_newton_Root2D_boozer
  TYPE(t_fbase), POINTER :: AB_fbase_in
  REAL(wp), POINTER :: A_in(:), B_in(:)  ! len: modes
  REAL(wp) :: x0(2)
  CONTAINS
  PROCEDURE :: FR  => get_booz_newton_FR
  PROCEDURE :: dFR => get_booz_newton_dFR
END TYPE t_newton_Root2D_boozer

INTERFACE sfl_boozer_new
  MODULE PROCEDURE sfl_boozer_new
END INTERFACE

PUBLIC :: t_sfl_boozer,sfl_boozer_new,find_boozer_angles
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize sfl boozer class
!!
!===================================================================================================================================
SUBROUTINE sfl_boozer_new(sf,mn_max,mn_nyq,nfp,sin_cos,hmap_in,nrho,rho_pos,iota,phiPrime,relambda_in)
  ! MODULES
  USE MODgvec_fbase   ,ONLY: fbase_new
  USE MODgvec_hmap,  ONLY: hmap_new_auxvar
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  INTEGER,INTENT(IN) :: mn_max(2)  !! maximum Fourier modes in theta and zeta
  INTEGER,INTENT(IN) :: mn_nyq(2)  !! number of equidistant integration points (trapezoidal rule) in m and n
  INTEGER,INTENT(IN) :: nfp        !! number of field periods
  CHARACTER(LEN=8)   :: sin_cos      !! can be either only sine: " _sin_" only cosine: " _cos_" or full: "_sin_cos_"
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP),INTENT(IN),TARGET :: hmap_in
#else
  CLASS(c_hmap)  ,INTENT(IN),TARGET :: hmap_in
#endif
  INTEGER,INTENT(IN) :: nrho       !! number of rho positions
  REAL(wp),INTENT(IN) :: rho_pos(nrho),iota(nrho),phiPrime(nrho)  !! rho positions, iota and phiPrime at these rho positions
  LOGICAL, INTENT(IN),OPTIONAL :: relambda_in  !! DEFAULT=TRUE: lambda is recomputed on the given fourier resolution, RECOMMENDED
                                   !!   for exact integrability condition of boozer transform, but slower.
                                   !! FALSE: lambda from equilibrium solution is taken.
  ! OUTPUT VARIABLES
  TYPE(t_sfl_boozer), ALLOCATABLE,INTENT(INOUT) :: sf !! self
  !=================================================================================================================================
  ALLOCATE(sf)
  sf%nrho = nrho
  ALLOCATE(sf%rho_pos(nrho),sf%iota(nrho),sf%phiPrime(nrho))
  sf%rho_pos = rho_pos
  IF(ANY((sf%rho_pos.LT.1e-4_wp).OR.(sf%rho_pos.GT. 1.0_wp))) CALL abort(__STAMP__, &
         "sfl_boozer_new: rho_pos must be >=1e-4 and <=1.0", &
         TypeInfo="InvalidParameterError")
  sf%iota = iota
  sf%phiPrime = phiPrime
  IF(PRESENT(relambda_in)) THEN
    sf%relambda=relambda_in
  ELSE
    sf%relambda = .TRUE. ! default
  END IF
  CALL fbase_new(sf%nu_fbase,mn_max,mn_nyq,nfp,sin_cos,.TRUE.)
  sf%hmap => hmap_in
  CALL hmap_new_auxvar(sf%hmap,sf%nu_fbase%x_IP(2,:),sf%hmap_xv,.TRUE.)
  ALLOCATE(sf%lambda(sf%nu_fbase%modes,nrho),sf%nu(sf%nu_fbase%modes,nrho))
  sf%initialized=.TRUE.
END SUBROUTINE sfl_boozer_new

!===================================================================================================================================
!> finalize sfl boozer class
!!
!===================================================================================================================================
SUBROUTINE sfl_boozer_free(sf)
  ! MODULES
  IMPLICIT NONE
  CLASS(t_sfl_boozer), INTENT(INOUT) :: sf !! self
  !=================================================================================================================================
  DEALLOCATE(sf%rho_pos,sf%lambda,sf%nu,sf%iota,sf%phiPrime)
  CALL sf%nu_fbase%free()
  DEALLOCATE(sf%nu_fbase)
  DEALLOCATE(sf%hmap_xv)
  NULLIFY(sf%hmap)
  sf%initialized=.FALSE.
END SUBROUTINE sfl_boozer_free



!===================================================================================================================================
!> Builds the boozer transform coordinate
!! theta^B = theta + lambda + iota(s)*nu(s,theta,zeta)
!! zeta^B  = zeta +nu(s,theta,zeta)
!!
!! since in Boozer, the covariant magnetic field components are the current profiles,
!! B = Itor(s) grad theta^B + Ipol(s) grad zeta^B + X grad s
!!   = Itor(s) grad (theta+lambda+iota*nu) + Ipol(s) grad (zeta + nu) + X grad s
!!   = (Itor*(1+dlambda/dtheta) + (Itor*iota+Ipol)*dnu/dtheta) grad theta + (Itor*(dlambda/dzeta)+(Itor*iota+Ipol)*dnu/dzeta)
!!=> dnu/dtheta = (B_theta - Itor - Itor*dlambda/dtheta ) / (Itor*iota+Ipol)
!!=> dnu/dzeta  = (B_zeta  - Ipol - Itor*dlambda/dzeta  ) / (Itor*iota+Ipol)
!! There is a integrability condition for nu:
!!     d/dzeta(dnu/dtheta)-d/dthet(dnu/dzeta)=d/dzeta(dB_theta/dtheta)-d/dthet(dB_theta/dzeta)=0
!! which is equivalent to impose J^s=0.
!! now if lambda is recomputed via a projection of J^s=0 onto the same fourier series as nu, the compatibility condition is
!! EXACTLY(!) fullfilled.
!===================================================================================================================================
SUBROUTINE Get_Boozer_sinterp(sf,X1_base_in,X2_base_in,LA_base_in,X1_in,X2_in,LA_in)
  ! MODULES
  USE MODgvec_Globals,ONLY: UNIT_stdOut,ProgressBar
  USE MODgvec_base,ONLY: t_base
  USE MODgvec_fbase,ONLY: t_fbase,fbase_new,sin_cos_map
  USE MODgvec_LinAlg
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
    CLASS(t_sfl_boozer), INTENT(INOUT) :: sf !!!-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER               :: mn_max(2),mn_nyq(2),irho,iMode,modes,i_mn,mn_IP
    INTEGER               :: nfp
    REAL(wp)              :: spos,dthet_dzeta,dPhids_int,iota_int,dChids_int
    REAL(wp)              :: b_thet,b_zeta
    REAL(wp)              :: detJ,Itor,Ipol,stmp

    REAL(wp)                          ::  X1_s(  1:X1_base_in%f%modes)
    REAL(wp)                          :: dX1ds_s(1:X1_base_in%f%modes)
    REAL(wp)                          ::  X2_s(  1:X2_base_in%f%modes)
    REAL(wp)                          :: dX2ds_s(1:X2_base_in%f%modes)
    REAL(wp),ALLOCATABLE              :: LA_s(:,:)
    REAL(wp),DIMENSION(sf%nu_fbase%modes) :: nu_m,nu_n
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: Bcov_thet_IP,Bcov_zeta_IP
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: dLAdthet_IP,dLAdzeta_IP
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: LA_IP,fm_IP,fn_IP,gam_tt,gam_tz,gam_zz
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: X1_IP,dX1ds_IP,dX1dthet,dX1dzeta
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: X2_IP,dX2ds_IP,dX2dthet,dX2dzeta
    TYPE(t_fbase),ALLOCATABLE             :: X1_fbase_nyq
    TYPE(t_fbase),ALLOCATABLE             :: X2_fbase_nyq
    TYPE(t_fbase),ALLOCATABLE             :: LA_fbase_nyq
  !===================================================================================================================================
    nfp = X1_base_in%f%nfp
    IF(nfp.NE.sf%nu_fbase%nfp) CALL abort(__STAMP__, &
                   'GET BOOZER ANGLE TRANSFORM,  sf%nu_fbase%nfp not the same as in X1')
    mn_max(1:2)=sf%nu_fbase%mn_max
    mn_nyq(1:2)=sf%nu_fbase%mn_nyq
    SWRITE(UNIT_StdOut,'(A,A,A,I4,3(A,2I6))')'GET BOOZER ANGLE TRANSFORM (', &
                                             TRIM(MERGE('RECOMPUTE      ','USE EQUILIBRIUM',sf%relambda)),' LAMBDA!), nfp=',nfp, &
                                            ', mn_max_in=',LA_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
    __PERFON('get_boozer')
    __PERFON('init')

    mn_IP        = sf%nu_fbase%mn_IP  !total number of integration points
    modes        = sf%nu_fbase%modes  !number of modes in output
    dthet_dzeta  = sf%nu_fbase%d_thet*sf%nu_fbase%d_zeta !integration weights

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
      ALLOCATE(LA_s(1:LA_base_in%f%modes,sf%nrho))
      SWRITE(UNIT_StdOut,*)'        ...Init LA_nyq Base Done'
    END IF

    !!!ALLOCATE(LA_s(1:LA_fbase_nyq%modes))
    nu_m=0.0_wp; nu_n=0.0_wp

    __PERFOFF('init')


    CALL ProgressBar(0,sf%nrho) !INIT
    DO irho=1,sf%nrho
      __PERFON('eval_data')
      spos=sf%rho_pos(irho)

      dPhids_int  = sf%phiPrime(irho)
      iota_int    = sf%iota(irho)
      dChids_int  = dPhids_int*iota_int

      !interpolate radially
      X1_s(:)    = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,      0,X1_in(:,:))
      dX1ds_s(:) = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,DERIV_S,X1_in(:,:))

      X2_s(:)    = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,      0,X2_in(:,:))
      dX2ds_s(:) = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,DERIV_S,X2_in(:,:))

      IF(.NOT.sf%relambda) THEN
        LA_s(:,irho) = LA_base_in%s%evalDOF2D_s(spos,LA_base_in%f%modes,      0,LA_in(:,:))
      END IF

      !evaluate at integration points
      X1_IP    = X1_fbase_nyq%evalDOF_IP(         0, X1_s(  :))
      dX1ds_IP = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s(:))
      dX1dthet = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s(  :))
      dX1dzeta = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s(  :))

      X2_IP    = X2_fbase_nyq%evalDOF_IP(         0, X2_s(  :))
      dX2ds_IP = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s(:))
      dX2dthet = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s(  :))
      dX2dzeta = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s(  :))


      __PERFOFF('eval_data')
      __PERFON('eval_bsub')
      __PERFON('eval_metrics')

  !$OMP PARALLEL DO &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE)  &
  !$OMP   PRIVATE(i_mn,detJ)  &
  !$OMP   SHARED(sf,mn_IP,dX1ds_IP,dX2ds_IP,dX1dthet,dX2dthet,dX1dzeta,dX2dzeta,X1_IP,X2_IP,gam_tt,gam_tz,gam_zz)
      !evaluate metrics on (theta,zeta)
      DO i_mn=1,mn_IP

        detJ        =  ( dX1ds_IP(i_mn)*dX2dthet(i_mn) -dX2ds_IP(i_mn)*dX1dthet(i_mn) ) &
                     * sf%hmap%eval_Jh_aux(X1_IP(i_mn),X2_IP(i_mn),sf%hmap_xv(i_mn)) !Jp*Jh
        gam_tt(i_mn)  = sf%hmap%eval_gij_aux(dX1dthet(i_mn),dX2dthet(i_mn),0.0_wp, &
                                              X1_IP  (i_mn), X2_IP(  i_mn),        &
                                             dX1dthet(i_mn),dX2dthet(i_mn),0.0_wp, &
                                             sf%hmap_xv(i_mn) )/detJ   !g_theta,theta
        gam_tz(i_mn)  = sf%hmap%eval_gij_aux(dX1dthet(i_mn),dX2dthet(i_mn),0.0_wp, &
                                              X1_IP  (i_mn), X2_IP(  i_mn),        &
                                             dX1dzeta(i_mn),dX2dzeta(i_mn),1.0_wp, &
                                             sf%hmap_xv(i_mn) )/detJ   !g_zeta,theta
        gam_zz(i_mn)  = sf%hmap%eval_gij_aux(dX1dzeta(i_mn),dX2dzeta(i_mn),1.0_wp, &
                                              X1_IP  (i_mn), X2_IP(  i_mn),        &
                                             dX1dzeta(i_mn),dX2dzeta(i_mn),1.0_wp, &
                                             sf%hmap_xv(i_mn) )/detJ   !g_zeta,zeta
      END DO !i_mn
  !$OMP END PARALLEL DO

      __PERFOFF('eval_metrics')

      IF(sf%relambda)THEN
      __PERFON('new_lambda')
        CALL Lambda_setup_and_solve(sf%nu_fbase,dPhids_int,dchids_int,gam_tt,gam_tz,gam_zz,sf%lambda(:,irho))
        !CALL Lambda_setup_and_solve(LA_fbase_nyq,dPhids_int,dchids_int,gam_tt,gam_tz,gam_zz,LA_s)
        LA_IP(:)    = sf%nu_fbase%evalDOF_IP(         0,sf%lambda(:,irho))
        dLAdthet_IP = sf%nu_fbase%evalDOF_IP(DERIV_THET,sf%lambda(:,irho))
        dLAdzeta_IP = sf%nu_fbase%evalDOF_IP(DERIV_ZETA,sf%lambda(:,irho))
        __PERFOFF('new_lambda')
      ELSE
        LA_IP(:)    = LA_fbase_nyq%evalDOF_IP(         0,LA_s(:,irho))
        dLAdthet_IP = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:,irho))
        dLAdzeta_IP = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s(:,irho))
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

      stmp=1.0_wp/(Itor*iota_int+Ipol)
  !$OMP PARALLEL DO        &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i_mn)        &
  !$OMP   SHARED(mn_IP,Itor,Ipol,stmp,dLAdthet_IP,Bcov_thet_IP,fm_IP)
      DO i_mn=1,mn_IP
        fm_IP(i_mn)  = (Bcov_thet_IP(i_mn)-Itor-Itor*dLAdthet_IP(i_mn))*stmp
      END DO
  !$OMP END PARALLEL DO

      !projection: only onto base_dthet
      CALL sf%nu_fbase%projectIPtoDOF(.FALSE.,1.0_wp,DERIV_THET,fm_IP(:),nu_m(:))

      IF(sf%nu_fbase%mn_max(2).GT.0) THEN !3D case
  !$OMP PARALLEL DO        &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i_mn)        &
  !$OMP   SHARED(mn_IP,Itor,Ipol,stmp,dLAdzeta_IP,Bcov_zeta_IP,fn_IP)
        DO i_mn=1,mn_IP
          fn_IP(i_mn)= (Bcov_zeta_IP(i_mn)-Ipol-Itor*dLAdzeta_IP(i_mn))*stmp
        END DO
  !$OMP END PARALLEL DO

        !projection onto base_dzeta
        CALL sf%nu_fbase%projectIPtoDOF(.FALSE.,1.0_wp,DERIV_ZETA,fn_IP(:),nu_n(:))
      END IF !3D case (n_max >0)

      ! only if n=0, use formula from base_dthet projected G, else use base_dzeta projected G
      DO iMode=1,modes
        ASSOCIATE(m=>sf%nu_fbase%Xmn(1,iMode),n=>sf%nu_fbase%Xmn(2,iMode))
        IF(m.NE.0) nu_m(iMode)=nu_m(iMode)*(dthet_dzeta*sf%nu_fbase%snorm_base(iMode))/REAL(m*m,wp)
        IF(n.NE.0) nu_n(iMode)=nu_n(iMode)*(dthet_dzeta*sf%nu_fbase%snorm_base(iMode))/REAL(n*n,wp)
        IF(n.EQ.0)THEN
          sf%nu(iMode,irho)=nu_m(iMode)
        ELSE
          sf%nu(iMode,irho)=nu_n(iMode)
        END IF
        IF((m.NE.0).AND.(n.NE.0))THEN
          !compare G_mn results:,
          !WRITE(*,"(A,I3,X,A,I3,X,2(A,X,E11.4,X),A,E11.4)")'DEBUG m=',m,'n=',n,'nu_m',nu_m(iMode),'nu_n',nu_n(iMode),'nu_m - nu_n=',nu_m(iMode)-nu_n(iMode)
        END IF
        END ASSOCIATE !m,n
      END DO
      !write(*,*)'DEBUG ===',is



      __PERFOFF('project')
      CALL ProgressBar(irho,sf%nrho)
    END DO !is
    CALL X1_fbase_nyq%free() ; DEALLOCATE(X1_fbase_nyq)
    CALL X2_fbase_nyq%free() ; DEALLOCATE(X2_fbase_nyq)
    IF(.NOT. sf%relambda) THEN
      CALL sf%nu_fbase%change_base(LA_fbase_nyq,2,LA_s,sf%lambda) !save lambda to sfl_boozer structure!
      CALL LA_fbase_nyq%free() ; DEALLOCATE(LA_fbase_nyq) ; DEALLOCATE(LA_s)
    END IF
    SWRITE(UNIT_StdOut,'(A)') '...DONE.'
    __PERFOFF('get_boozer')
  END SUBROUTINE Get_Boozer_sinterp


!===================================================================================================================================
!> interface to find_boozer_angles from the class t_sfl_boozer
!!
!===================================================================================================================================
SUBROUTINE self_find_boozer_angles(sf,tz_dim,tz_boozer,thetzeta_out)
  ! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sfl_boozer), INTENT(IN) :: sf
  INTEGER             ,INTENT(IN) :: tz_dim                !< size of the list of in thetstar,zetastar
  REAL(wp)            ,INTENT(IN) :: tz_boozer(2,tz_dim) !< theta,zeta positions in boozer angle (same for all rho)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)  ,INTENT(OUT)   :: thetzeta_out(2,tz_dim,sf%nrho)  !! theta,zeta position in original angles, for given boozer angles
!===================================================================================================================================
  CALL find_boozer_angles(sf%nrho,sf%iota,sf%nu_fbase,sf%lambda,sf%nu,tz_dim,tz_boozer,thetzeta_out)
END SUBROUTINE self_find_boozer_angles


!===================================================================================================================================
!> interface to find_boozer_angles from the class t_sfl_boozer
!!
!===================================================================================================================================
SUBROUTINE self_find_boozer_angles_irho(sf,irho,tz_dim,tz_boozer,thetzeta_out)
  ! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sfl_boozer), INTENT(IN) :: sf
  INTEGER             ,INTENT(IN) :: irho                !< index of the flux surface to find boozer angles on, within 1:nrho
  INTEGER             ,INTENT(IN) :: tz_dim              !< size of the list of in thetstar,zetastar
  REAL(wp)            ,INTENT(IN) :: tz_boozer(2,tz_dim) !< theta,zeta positions in boozer angle (same for all rho)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)  ,INTENT(OUT)   :: thetzeta_out(2,tz_dim)     !< theta,zeta position in original angles, for given boozer angles
!===================================================================================================================================
  CALL find_boozer_angles(1,sf%iota(irho:irho),sf%nu_fbase,sf%lambda(:,irho:irho),sf%nu(:,irho:irho),tz_dim,tz_boozer,thetzeta_out)
END SUBROUTINE self_find_boozer_angles_irho


!===================================================================================================================================
!> on one flux surface, find for an given list of  (thet*_j,zeta*_j), the corresponding (thet_j,zeta_j) positions, given
!> Here, new boozer angles are
!> theta*=theta+Gt(theta,zeta)
!>  zeta*=zeta+nu(theta,zeta),
!> with Gt=lambda+iota*nu and nu periodic functions and zero average and same base
!> Note that in this routine, we will use a 2d root search with a newton method, setting
!> [f1,f2]^T = [thet+A(thet,zeta)-thet*=0,  zeta+B(thet,zeta)-zeta*=0]^T
!> that includes the derivatives (Jacobian), so that the newton step needs to the solved:
!> -[f1]    [ 1+dA/dthet    dA/dzeta] [dthet]
!>  |  | =  |                       | |     |
!> -[f2]    [  dB/dthet   1+dB/dzeta] [dzeta]
!!
!===================================================================================================================================
SUBROUTINE find_boozer_angles(nrho,iota,fbase_in,LA_in,nu_in,tz_dim,tz_boozer,thetzeta_out)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,PI,ProgressBar,testlevel
USE MODgvec_fbase  ,ONLY: t_fbase
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER      ,INTENT(IN) :: nrho   !! number of surfaces, (second dimension  of LA_in and nu_in modes)
  REAL(wp)     ,INTENT(IN) :: iota(nrho)  !! iota at the rho positions.
  TYPE(t_fbase),INTENT(IN) ::fbase_in     !< same basis of lambda and nu
  REAL(wp)     ,INTENT(IN) :: LA_in(1:fbase_in%modes,nrho) !< fourier coefficients of thet*=thet+LA(theta,zeta)+iota*nu(theta,zeta)
  REAL(wp)     ,INTENT(IN) :: nu_in(1:fbase_in%modes,nrho) !< coefficients of zeta*=zeta+nu(theta,zeta)
  INTEGER      ,INTENT(IN) :: tz_dim                !< size of the list of in thetstar,zetastar
  REAL(wp)     ,INTENT(IN) :: tz_boozer(2,tz_dim) !< theta,zeta positions in boozer angle (same for all rho)

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) ,INTENT(OUT)   :: thetzeta_out(2,tz_dim,nrho)  !! theta,zeta position in original angles, for given boozer angles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: irho,j
  REAL(wp)              :: bounds(2),x0(2)
  REAL(wp)              :: check(tz_dim),maxerr(2,nrho)
  REAL(wp)              :: Gt(1:fbase_in%modes)  !! transform in theta: lambda + iota*nu
  LOGICAL               :: docheck
!===================================================================================================================================
  __PERFON('find_boozer_angles')
  SWRITE(UNIT_StdOut,'(A,2(I8,A))')'Find boozer angles via 2D Newton on  nrho=',nrho,' times ntheta_zeta= ',tz_dim, " points"
  docheck=(testlevel.GT.0)
  bounds=(/PI, PI/fbase_in%nfp/)
  CALL ProgressBar(0,nrho)!init
  DO irho=1,nrho
    Gt(:)=LA_in(:,irho)+iota(irho)*nu_in(:,irho)
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP   PRIVATE(j,x0) FIRSTPRIVATE(irho,bounds) &
!$OMP   SHARED(tz_dim,tz_boozer,thetzeta_out,fbase_in,Gt,nu_in)
    DO j=1,tz_dim
        x0=tz_boozer(:,j)
        thetzeta_out(:,j,irho)=get_booz_newton(x0,bounds,fbase_in,Gt,nu_in(:,irho))
    END DO !j
!$OMP END PARALLEL DO
    CALL ProgressBar(irho,nrho)
  END DO !irho

  IF(docheck)THEN
    DO irho=1,nrho
      check=fbase_in%evalDOF_xn(tz_dim,thetzeta_out(:,:,irho),0,(LA_in(:,irho)+iota(irho)*nu_in(:,irho)))
      maxerr(1,irho)=maxval(abs(check+(thetzeta_out(1,:,irho)-tz_boozer(1,:))))
      check=fbase_in%evalDOF_xn(tz_dim,thetzeta_out(:,:,irho),0,nu_in(:,irho))
      maxerr(2,irho)=maxval(abs(check+(thetzeta_out(2,:,irho)-tz_boozer(2,:))))
    END DO

    IF(ANY(maxerr(:,:).GT.1.0e-12)) THEN
        WRITE(UNIT_stdout,*)'CHECK BOOZER THETA*',MAXVAL(maxerr(1,:))
        WRITE(UNIT_stdout,*)'CHECK BOOZER ZETA*', maxerr
        CALL abort(__STAMP__, &
        "find_boozer_angles: Error in transform")
    END IF
  END IF  !docheck
  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
__PERFOFF('find_boozer_angles')

END SUBROUTINE find_boozer_angles

!===================================================================================================================================
!> This function returns the result of the 2D newton root search for the boozer angle
!!
!===================================================================================================================================
FUNCTION get_booz_newton(x0,bounds,AB_fbase_in,A_in,B_in) RESULT(x_out)
  USE MODgvec_fbase  ,ONLY: t_fbase
  USE MODgvec_Newton ,ONLY: NewtonRoot2D
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!INPUT VARIABLES
  REAL(wp),INTENT(IN) :: x0(2),bounds(2)
  TYPE(t_fbase),INTENT(IN), TARGET :: AB_fbase_in
  REAL(wp),INTENT(IN), TARGET :: A_in(1:AB_fbase_in%modes),B_in(1:AB_fbase_in%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
  REAL(wp) :: x_out(2)
!-----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
  TYPE(t_newton_Root2D_boozer) :: fobj
!===================================================================================================================================

  !                                     a     b       maxstep  , xinit    ,funcs, funcs_jac
  fobj%AB_fbase_in => AB_fbase_in
  fobj%A_in => A_in
  fobj%B_in => B_in
  fobj%x0 = x0
  x_out = NewtonRoot2D(1.0e-12_wp,x0-bounds,x0+bounds,0.1_wp*bounds,x0,fobj)
END FUNCTION get_booz_newton

!===================================================================================================================================
!> Target function for finding the logical angle for given boozer angles
!!
!===================================================================================================================================
FUNCTION get_booz_newton_FR(sf, x) RESULT(FF)
  IMPLICIT NONE
  CLASS(t_newton_Root2D_boozer), INTENT(IN) :: sf
  REAL(wp), INTENT(IN) :: x(2) ! xiter
  REAL(wp) :: FF(2) !two functions of x1,x2 to find root of
  REAL(wp),DIMENSION(sf%AB_fbase_in%modes) :: base_x

  base_x = sf%AB_fbase_in%eval(0, x) !base evaluation

  FF(1) = x(1) - sf%x0(1) + DOT_PRODUCT(base_x, sf%A_in)
  FF(2) = x(2) - sf%x0(2) + DOT_PRODUCT(base_x, sf%B_in)
END FUNCTION get_booz_newton_FR

!===================================================================================================================================
!> Derivative of the target function for finding the logical angle for given boozer angles
!!
!===================================================================================================================================
FUNCTION get_booz_newton_dFR(sf, x) RESULT(dFF)
  IMPLICIT NONE
  CLASS(t_newton_Root2D_boozer), INTENT(IN) :: sf
  REAL(wp), INTENT(IN) :: x(2)
  REAL(wp) :: dFF(2,2) !jacobian
  REAL(wp),DIMENSION(sf%AB_fbase_in%modes) :: base_dthet, base_dzeta

  base_dthet = sf%AB_fbase_in%eval(DERIV_THET, x) !dbase/dtheta
  base_dzeta = sf%AB_fbase_in%eval(DERIV_ZETA, x) !dbase/dtheta

  dFF(1,:) = (/1.0_wp + DOT_PRODUCT(base_dthet, sf%A_in),          DOT_PRODUCT(base_dzeta, sf%A_in)/)
  dFF(2,:) = (/         DOT_PRODUCT(base_dthet, sf%B_in), 1.0_wp + DOT_PRODUCT(base_dzeta, sf%B_in)/)
END FUNCTION get_booz_newton_dFR

END MODULE MODgvec_SFL_Boozer
