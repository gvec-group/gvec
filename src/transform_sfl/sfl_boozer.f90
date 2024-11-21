!===================================================================================================================================
! Copyright (C) 2024  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **SFL boozer**
!!
!! Transform to Straight-field line  BOOZER coordinates 
!!
!===================================================================================================================================
MODULE MODgvec_SFL_Boozer
! MODULES
USE MODgvec_Globals, ONLY:wp,abort,MPIroot
USE MODgvec_fbase   ,ONLY: t_fbase
USE MODgvec_hmap,  ONLY: c_hmap
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
  CLASS(t_fbase), ALLOCATABLE :: nu_fbase
  CLASS(c_hmap),  POINTER     :: hmap          !! pointer to hmap class
  REAL(wp),ALLOCATABLE::rho_pos(:),iota(:),phiPrime(:) !! rho positions, iota and phiPrime at these rho positions
  ! computed in the boozer transform
  REAL(wp),ALLOCATABLE::lambda(:,:),nu(:,:)   !! Fourier modes for all rho positions of lambda (recomputed on the fourier space of nu) and nu for boozer transform , (iMode,irho)
  CONTAINS
  PROCEDURE :: get_boozer => get_boozer_sinterp
  PROCEDURE :: free       => sfl_boozer_free
  !PROCEDURE :: find_angles => find_boozer_angles
END TYPE t_sfl_boozer


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
SUBROUTINE sfl_boozer_new(sf,mn_max,mn_nyq,nfp,sin_cos,hmap_in,nrho,rho_pos,iota,phiPrime)
  ! MODULES
  USE MODgvec_fbase   ,ONLY: fbase_new
  
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  INTEGER,INTENT(IN) :: mn_max(2)  !! maximum Fourier modes in theta and zeta
  INTEGER,INTENT(IN) :: mn_nyq(2)  !! number of equidistant integration points (trapezoidal rule) in m and n
  INTEGER,INTENT(IN) :: nfp        !! number of field periods
  CHARACTER(LEN=8)   :: sin_cos      !! can be either only sine: " _sin_" only cosine: " _cos_" or full: "_sin_cos_"
  CLASS(c_hmap),INTENT(IN),TARGET :: hmap_in
  INTEGER,INTENT(IN) :: nrho       !! number of rho positions
  REAL(wp),INTENT(IN) :: rho_pos(nrho),iota(nrho),phiPrime(nrho)  !! rho positions, iota and phiPrime at these rho positions
  ! OUTPUT VARIABLES
  CLASS(t_sfl_boozer), ALLOCATABLE,INTENT(INOUT) :: sf !! self
  !=================================================================================================================================
  ALLOCATE(t_sfl_boozer::sf)
  sf%nrho = nrho
  ALLOCATE(sf%rho_pos(nrho),sf%iota(nrho),sf%phiPrime(nrho))
  sf%rho_pos = rho_pos
  IF(ANY((sf%rho_pos.LT.1e-4_wp).OR.(sf%rho_pos.GT. 1.0_wp))) CALL abort(__STAMP__, &
         "sfl_boozer_new: rho_pos must be >1e-4 and <=1.0")
  sf%iota = iota
  sf%phiPrime = phiPrime
  sf%relambda = .TRUE. ! default
  CALL fbase_new(sf%nu_fbase,mn_max,mn_nyq,nfp,sin_cos,.TRUE.)
  sf%hmap => hmap_in
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
    REAL(wp)              :: b_thet,b_zeta,qloc(3),q_thet(3),q_zeta(3) 
    REAL(wp)              :: detJ,Itor,Ipol,stmp
  
    REAL(wp)                          ::  X1_s(  1:X1_base_in%f%modes) 
    REAL(wp)                          :: dX1ds_s(1:X1_base_in%f%modes) 
    REAL(wp)                          ::  X2_s(  1:X2_base_in%f%modes) 
    REAL(wp)                          :: dX2ds_s(1:X2_base_in%f%modes) 
    REAL(wp),ALLOCATABLE              :: LA_s(:)
    REAL(wp),DIMENSION(sf%nu_fbase%modes) :: nu_m,nu_n
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: Bcov_thet_IP,Bcov_zeta_IP
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: dLAdthet_IP,dLAdzeta_IP
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: LA_IP,fm_IP,fn_IP,gam_tt,gam_tz,gam_zz 
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: X1_IP,dX1ds_IP,dX1dthet_IP,dX1dzeta_IP
    REAL(wp),DIMENSION(sf%nu_fbase%mn_IP) :: X2_IP,dX2ds_IP,dX2dthet_IP,dX2dzeta_IP
    CLASS(t_fBase),ALLOCATABLE             :: X1_fbase_nyq
    CLASS(t_fBase),ALLOCATABLE             :: X2_fbase_nyq
    CLASS(t_fBase),ALLOCATABLE             :: LA_fbase_nyq
  !===================================================================================================================================
    nfp = X1_base_in%f%nfp
    IF(nfp.NE.sf%nu_fbase%nfp) CALL abort(__STAMP__, &
                   'GET BOOZER ANGLE TRANSFORM,  sf%nu_fbase%nfp not the same as in X1')
    mn_max(1:2)=sf%nu_fbase%mn_max
    mn_nyq(1:2)=sf%nu_fbase%mn_nyq
    SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'GET BOOZER ANGLE TRANSFORM, nfp=',nfp, &
                                ', mn_max_in=',LA_base_in%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
    __PERFON('get_boozer')
    __PERFON('init')
    
    mn_IP        = sf%nu_fbase%mn_IP  !total number of integration points
    modes        = sf%nu_fbase%modes  !number of modes in output
    dthet_dzeta  = sf%nu_fbase%d_thet*sf%nu_fbase%d_zeta !integration weights
  
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
      ALLOCATE(LA_s(1:LA_base_in%f%modes))
      SWRITE(UNIT_StdOut,*)'        ...Init LA_nyq Base Done'
    END IF
    
    !!!ALLOCATE(LA_s(1:LA_fbase_nyq%modes))
    nu_m=0.0_wp; nu_n=0.0_wp
  
    __PERFOFF('init')
  
  
    CALL ProgressBar(0,sf%nrho) !INIT
    DO irho=sf%nrho,1,-1
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
        qloc(  1:3) = (/ X1_IP(     i_mn), X2_IP(     i_mn),sf%nu_fbase%x_IP(2,i_mn)/)
        q_thet(1:3) = (/dX1dthet_IP(i_mn),dX2dthet_IP(i_mn),0.0_wp/) !dq(1:2)/dtheta
        q_zeta(1:3) = (/dX1dzeta_IP(i_mn),dX2dzeta_IP(i_mn),1.0_wp/) !dq(1:2)/dzeta
        detJ        =  sf%hmap%eval_Jh(qloc)*( dX1ds_IP(i_mn)*dX2dthet_IP(i_mn) &
                      -dX2ds_IP(i_mn)*dX1dthet_IP(i_mn) )
        gam_tt(i_mn)  = sf%hmap%eval_gij(q_thet,qloc,q_thet)/detJ   !g_theta,theta
        gam_tz(i_mn)  = sf%hmap%eval_gij(q_thet,qloc,q_zeta)/detJ   !g_theta,zeta =g_zeta,theta
        gam_zz(i_mn)  = sf%hmap%eval_gij(q_zeta,qloc,q_zeta)/detJ   !g_zeta,zeta
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
  
      __PERFON('project_nu')
  
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
  
      __PERFOFF('project_nu')
  
      
  
      __PERFOFF('project')
      CALL ProgressBar(irho,sf%nrho)
    END DO !is
    IF(.NOT. sf%relambda) DEALLOCATE(LA_s)
    SWRITE(UNIT_StdOut,'(A)') '...DONE.'
    __PERFOFF('get_boozer')
  END SUBROUTINE Get_Boozer_sinterp

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
!> 
!!
!===================================================================================================================================
SUBROUTINE find_boozer_angles(nrho,iota,fbase_in,LA_in,nu_in,tz_dim,tz_boozer,thetzeta_out)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,PI
USE MODgvec_fbase  ,ONLY: t_fbase
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER      ,INTENT(IN) :: nrho   !! number of surfaces, (second dimension  of LA_in and nu_in modes)
  REAL(wp)     ,INTENT(IN) :: iota(nrho)  !! iota at the rho positions.
  CLASS(t_fBase),INTENT(IN) ::fbase_in     !< same basis of lambda and nu
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
  REAL(wp)              :: bounds(2)
  REAL(wp)              :: check(tz_dim),maxerr(2,nrho)
  REAL(wp)              :: Gt(1:fbase_in%modes,nrho)  !! transform in theta: lambda + iota*nu
  LOGICAL               :: docheck=.TRUE.
!===================================================================================================================================

  SWRITE(UNIT_StdOut,'(A,I4)')'Find boozer angles via 2D Newton on  ',tz_dim, " points"
                     
  __PERFON('find_boozer_angles')
  bounds=(/PI, PI/fbase_in%nfp/)
  DO irho=1,nrho
    Gt(:,irho)=LA_in(:,irho)+iota(irho)*nu_in(:,irho)
  END DO
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(irho,j)
  DO irho=1,nrho
    DO j=1,tz_dim
        !x0=tz_boozer(:,i,j) 
        thetzeta_out(:,j,irho)=get_booz_newton(tz_boozer(:,j),bounds,fbase_in,Gt(:,irho),nu_in(:,irho))
    END DO !j 
  END DO !irho
!$OMP END PARALLEL DO

  IF(docheck)THEN
    DO irho=1,nrho
      check=fbase_in%evalDOF_xn(tz_dim,thetzeta_out(:,:,irho),0,Gt(:,irho))
      maxerr(1,irho)=maxval(abs(check+(thetzeta_out(1,:,irho)-tz_boozer(1,:))))
      check=fbase_in%evalDOF_xn(tz_dim,thetzeta_out(:,:,irho),0,nu_in(:,irho))
      maxerr(2,irho)=maxval(abs(check+(thetzeta_out(2,:,irho)-tz_boozer(2,:))))
    END DO
  
    WRITE(*,*)'CHECK BOOZER THETA*',MAXVAL(maxerr(1,:))
    IF(ANY(maxerr(1,:).GT.1.0e-12)) CALL abort(__STAMP__, &
        "find_boozer_angles: Error in  theta*")
    WRITE(*,*)'CHECK BOOZER ZETA*', maxerr
    IF(ANY(maxerr(2,:).GT.1.0e-12)) CALL abort(__STAMP__, &
        "find_boozer_angles: Error in zeta*")
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
  CLASS(t_fBase),INTENT(IN) :: AB_fbase_in
  REAL(wp),INTENT(IN) :: A_in(1:AB_fbase_in%modes),B_in(1:AB_fbase_in%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
  REAL(wp) :: x_out(2)
!-----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
  REAL(wp),DIMENSION(AB_fbase_in%modes)::base_x,base_dthet,base_dzeta  !used inside sub-functions
!===================================================================================================================================

  !                                     a     b       maxstep  , xinit    ,funcs, funcs_jac
  x_out = NewtonRoot2D(1.0e-12_wp,x0-bounds,x0+bounds,0.1_wp*bounds,x0,ABtrafo,ABtrafo_jac)

  CONTAINS 
  !for newton root search 
  FUNCTION ABtrafo(xiter) RESULT(FF) 
    !uses current x0=(zetastar,thetastar) , and A,B and derivatives from subroutine above
    IMPLICIT NONE
    REAL(wp) :: xiter(2)
    REAL(wp) :: FF(2) !two functions of x1,x2 to find root of
    
    base_x =AB_fbase_in%eval(0,xiter) !base evaluation

    FF(1)=xiter(1)-x0(1)+ DOT_PRODUCT(base_x,A_in)
    FF(2)=xiter(2)-x0(2)+ DOT_PRODUCT(base_x,B_in)
  END FUNCTION ABtrafo

  FUNCTION ABtrafo_jac(xiter) RESULT(dFF) 
    !uses current x0=(zetastar,thetastar) , and A,B and derivatives from subroutine above
    IMPLICIT NONE
    REAL(wp) :: xiter(2)
    REAL(wp) :: dFF(2,2) !jacobian 
    base_dthet =AB_fbase_in%eval(DERIV_THET,xiter) !dbase/dtheta
    base_dzeta =AB_fbase_in%eval(DERIV_ZETA,xiter) !dbase/dtheta

    dFF(1,:)=(/1.0_wp+ DOT_PRODUCT(base_dthet,A_in),        DOT_PRODUCT(base_dzeta,A_in)/)
    dFF(2,:)=(/        DOT_PRODUCT(base_dthet,B_in),1.0_wp+ DOT_PRODUCT(base_dzeta,B_in)/)
  END FUNCTION ABtrafo_jac

END FUNCTION get_booz_newton


END MODULE MODgvec_SFL_Boozer
