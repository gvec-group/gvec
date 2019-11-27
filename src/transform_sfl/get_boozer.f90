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
!!# Module **GET BOOZER**
!!
!! compute the angular transform to boozer coordinates, G(s,theta,zeta):
!! theta^B = theta + lambda + iota(s)*G(s,theta,zeta) 
!! zeta^B  = zeta  + G(s,theta,zeta) 
!!
!===================================================================================================================================
MODULE MODgvec_get_boozer
! MODULES
USE MODgvec_Globals, ONLY:wp,abort
IMPLICIT NONE
PRIVATE

INTERFACE Get_Boozer
  MODULE PROCEDURE Get_Boozer
END INTERFACE


PUBLIC :: Get_Boozer
!===================================================================================================================================

CONTAINS

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
!===================================================================================================================================
SUBROUTINE Get_Boozer(mn_max,fac_nyq,sgrid_in,G_base_out,Gthet,GZ)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,TWOPI,PI,ProgressBar
USE MODgvec_LinAlg
USE MODgvec_base   ,ONLY: t_base,base_new
USE MODgvec_sGrid  ,ONLY: t_sgrid
USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map
USE MODgvec_ReadState_vars  ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_ReadState_vars  ,ONLY: LA_r,X1_r,X2_r 
USE MODgvec_ReadState_Vars  ,ONLY: profiles_1d,hmap_r,sbase_prof !for profiles

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER      ,INTENT(IN) :: mn_max(2)                                     !< maximum number for new variables in SFL coordinates
  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
                                                                            !< n_IP=fac_nyq*max(mn_max_in)
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: sgrid_in                          !< change grid for G_base_out
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base),ALLOCATABLE,INTENT(INOUT) :: G_base_out      !< new fourier basis of function Gthet,Gzeta 
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT) :: Gthet(:,:)  !< coefficients of Gthet=LA+iota*G
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT) :: GZ(:,:)      !< coefficients of G (zeta^B=zeta+G
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,modes,i_mn,mn_IP
  INTEGER               :: mn_nyq(2),nfp
  REAL(wp)              :: spos,dthet_dzeta,dPhids_int,iota_int,dChids_int
  REAL(wp)              :: b_thet,b_zeta,qloc(3),q_thet(3),q_zeta(3) 
  REAL(wp)              :: J_h,g_tt,g_tz,g_zz,sdetJ,Itor,Ipol,stmp 

  REAL(wp)                          ::  X1_s(  1:X1_base_r%f%modes) 
  REAL(wp)                          :: dX1ds_s(1:X1_base_r%f%modes) 
  REAL(wp)                          ::  X2_s(  1:X2_base_r%f%modes) 
  REAL(wp)                          :: dX2ds_s(1:X2_base_r%f%modes) 
  REAL(wp)                          :: LA_s(   1:LA_base_r%f%modes) 
  REAL(wp),DIMENSION(:)  ,ALLOCATABLE :: Bcov_thet_IP,Bcov_zeta_IP
  REAL(wp),DIMENSION(:)  ,ALLOCATABLE :: dLAdthet_IP,dLAdzeta_IP
  REAL(wp),DIMENSION(:)  ,ALLOCATABLE :: LA_IP,fm_IP,fn_IP,ft_IP 
  CLASS(t_fBase),ALLOCATABLE        :: X1_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE        :: X2_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE        :: LA_fbase_nyq
  REAL(wp),DIMENSION(:),ALLOCATABLE :: X1_IP,dX1ds_IP,dX1dthet_IP,dX1dzeta_IP
  REAL(wp),DIMENSION(:),ALLOCATABLE :: X2_IP,dX2ds_IP,dX2dthet_IP,dX2dzeta_IP
!===================================================================================================================================
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  nfp = X1_base_r%f%nfp
  IF(mn_max(2).EQ.0) mn_nyq(2)=1 !exception: 2D configuration

  SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'GET BOOZER ANGLE TRANSFORM, nfp=',nfp, &
                              ', mn_max_in=',LA_base_r%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
  __PERFON('get_boozer')
  __PERFON('init')
  !initialize
  
  ! extended base for q in the new angles, and on the new grid
  CALL base_new(G_base_out,  X1_base_r%s%deg,        &
                             X1_base_r%s%continuity, &
                             sgrid_in, & !X1_base_r%s%grid,       &
                             X1_base_r%s%degGP,      &
                             mn_max,mn_nyq,nfp,      & 
                 sin_cos_map(LA_base_r%f%sin_cos),   &
!'_sin_   ', &
                             LA_base_r%f%exclude_mn_zero) !exclude m=n=0

  
  mn_IP        = G_base_out%f%mn_IP  !total number of integration points
  modes        = G_base_out%f%modes  !number of modes in output
  nBase        = G_base_out%s%nBase  !number of radial points in output
  dthet_dzeta  = G_base_out%f%d_thet*G_base_out%f%d_zeta !integration weights


  !transpose basis functions and include norm for projection

  SWRITE(UNIT_StdOut,*)'        ...Init G_out Base Done'
  !same base for X1, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new( X1_fbase_nyq, X1_base_r%f%mn_max,  mn_nyq, &
                                X1_base_r%f%nfp, &
                    sin_cos_map(X1_base_r%f%sin_cos), &
                                X1_base_r%f%exclude_mn_zero)
  SWRITE(UNIT_StdOut,*)'        ...Init X1_nyq Base Done'

  CALL fbase_new( X2_fbase_nyq, X2_base_r%f%mn_max,  mn_nyq, &
                                X2_base_r%f%nfp, &
                    sin_cos_map(X2_base_r%f%sin_cos), &
                                X2_base_r%f%exclude_mn_zero)
  SWRITE(UNIT_StdOut,*)'        ...Init X2_nyq Base Done'

  !same base for lambda, but with new mn_nyq (for pre-evaluation of basis functions)
  CALL fbase_new(LA_fbase_nyq,  LA_base_r%f%mn_max,  mn_nyq, &
                                LA_base_r%f%nfp, &
                    sin_cos_map(LA_base_r%f%sin_cos), &
                                LA_base_r%f%exclude_mn_zero)
  SWRITE(UNIT_StdOut,*)'        ...Init LA_nyq Base Done'

  ALLOCATE( X1_IP(1:mn_IP),dX1ds_IP(1:mn_IP), dX1dthet_IP(1:mn_IP),dX1dzeta_IP(1:mn_IP),&
            X2_IP(1:mn_IP),dX2ds_IP(1:mn_IP), dX2dthet_IP(1:mn_IP),dX2dzeta_IP(1:mn_IP) )

  ALLOCATE(dLAdthet_IP(1:mn_IP), dLAdzeta_IP(1:mn_IP),Bcov_thet_IP(1:mn_IP),Bcov_zeta_IP(1:mn_IP))
     
  ALLOCATE(LA_IP(1:mn_IP),fm_IP(mn_IP),fn_IP(mn_IP),ft_IP(mn_IP))
  
  ALLOCATE(Gthet(nBase,1:modes))
  ALLOCATE(GZ(   nBase,1:modes))

  __PERFOFF('init')


  CALL ProgressBar(0,nBase) !INIT
  DO is=1,nBase
    __PERFON('eval_data')
    spos=MIN(MAX(1.0e-08_wp,G_base_out%s%s_IP(is)),1.0_wp-1.0e-12_wp) !interpolation points for q_in

    dPhids_int  = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
    iota_int    = sbase_prof%evalDOF_s(spos,        0,profiles_1d(:,3))
    dChids_int  = dPhids_int*iota_int 

    !interpolate radially
    X1_s(:)    = X1_base_r%s%evalDOF2D_s(spos,X1_base_r%f%modes,      0,X1_r(:,:))
    dX1ds_s(:) = X1_base_r%s%evalDOF2D_s(spos,X1_base_r%f%modes,DERIV_S,X1_r(:,:))

    X2_s(:)    = X2_base_r%s%evalDOF2D_s(spos,X2_base_r%f%modes,      0,X2_r(:,:))
    dX2ds_s(:) = X2_base_r%s%evalDOF2D_s(spos,X2_base_r%f%modes,DERIV_S,X2_r(:,:))

    LA_s(:)    = LA_base_r%s%evalDOF2D_s(spos,LA_base_r%f%modes,      0,LA_r(:,:))

    !evaluate at integration points
    X1_IP       = X1_fbase_nyq%evalDOF_IP(         0, X1_s(  :))
    dX1ds_IP    = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s(:))
    dX1dthet_IP = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s(  :))
    dX1dzeta_IP = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s(  :))

    X2_IP       = X2_fbase_nyq%evalDOF_IP(         0, X2_s(  :))
    dX2ds_IP    = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s(:))
    dX2dthet_IP = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s(  :))
    dX2dzeta_IP = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s(  :))

    LA_IP(:)    = LA_fbase_nyq%evalDOF_IP(         0,LA_s(:))
    dLAdthet_IP = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:))
    dLAdzeta_IP = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s(:))

    __PERFOFF('eval_data')
    __PERFON('eval_bsub')
    
    Itor=0.0_wp;Ipol=0.0_wp
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)  &
!$OMP   PRIVATE(i_mn,b_thet,b_zeta,qloc,q_thet,q_zeta,J_h,g_tt,g_tz,g_zz,sdetJ)  &
!$OMP   REDUCTION(+:Itor,Ipol) &
!$OMP   SHARED(mn_IP,hmap_r,dchids_int,dPhids_int,dLAdzeta_IP,dLAdthet_IP,X1_IP,X2_IP, &
!$OMP          G_base_out,dX1dthet_IP,dX2dthet_IP,dX1dzeta_IP,dX2dzeta_IP,             &
!$OMP          dX1ds_IP,dX2ds_IP,Bcov_thet_IP,Bcov_zeta_IP)
    !evaluate (theta*,zeta*) modes of q_in at (theta,zeta)
    DO i_mn=1,mn_IP
      b_thet = dchids_int- dPhids_int*dLAdzeta_IP(i_mn)    !b_theta
      b_zeta = dPhids_int*(1.0_wp   + dLAdthet_IP(i_mn))    !b_zeta

      qloc(  1:3) = (/ X1_IP(     i_mn), X2_IP(     i_mn),G_base_out%f%x_IP(2,i_mn)/)
      q_thet(1:3) = (/dX1dthet_IP(i_mn),dX2dthet_IP(i_mn),0.0_wp/) !dq(1:2)/dtheta
      q_zeta(1:3) = (/dX1dzeta_IP(i_mn),dX2dzeta_IP(i_mn),1.0_wp/) !dq(1:2)/dzeta

      J_h         = hmap_r%eval_Jh(qloc)
      g_tt        = hmap_r%eval_gij(q_thet,qloc,q_thet)   !g_theta,theta
      g_tz        = hmap_r%eval_gij(q_thet,qloc,q_zeta)   !g_theta,zeta =g_zeta,theta
      g_zz        = hmap_r%eval_gij(q_zeta,qloc,q_zeta)   !g_zeta,zeta

      sdetJ       = 1.0_wp/(J_h*( dX1ds_IP(i_mn)*dX2dthet_IP(i_mn) &
                                 -dX2ds_IP(i_mn)*dX1dthet_IP(i_mn) ))

      Bcov_thet_IP(i_mn) = (g_tt*b_thet + g_tz*b_zeta)*sdetJ  
      Bcov_zeta_IP(i_mn) = (g_tz*b_thet + g_zz*b_zeta)*sdetJ  
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
    CALL G_base_out%f%projectIPtoDOF(.FALSE.,1.0_wp,DERIV_THET,fm_IP(:),GZ(is,:))

!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i_mn)        &
!$OMP   SHARED(mn_IP,Itor,Ipol,stmp,dLAdzeta_IP,Bcov_zeta_IP,fn_IP)
    DO i_mn=1,mn_IP
      fn_IP(i_mn)= (Bcov_zeta_IP(i_mn)-Ipol-Itor*dLAdzeta_IP(i_mn))*stmp
    END DO
!$OMP END PARALLEL DO 

    !linearly superimpose: G_mn = alpha* 1/(4pi^2 m^2) Gthet_mn + (1-alpha) 1/(4pi^2 n^2) Gzeta_mn
    !              alpha = m^2/(m^2+n^2) => =0 for m=0 =1 for n=0
    !       =>  (1-alpha)= n^2/(m^2+n^2) => =1 for m=0 =0 for n=0
    !       =>    G_mn = 1/(4pi^2(m^2+n^2)) (Gthet_mn+Gzeta_mn)
    !GZ(is,:)=GZ(is,:)+MATMUL(TRANSPOSE(G_base_out%f%base_dzeta_IP),fn_IP)

    !add projection onto base_dzeta
    CALL G_base_out%f%projectIPtoDOF(.TRUE.,1.0_wp,DERIV_ZETA,fn_IP(:),GZ(is,:))

    ! include norm and weighting 1/(m^2+n^2)
    DO iMode=1,modes
      GZ(is,iMode)=GZ(is,iMode)*(dthet_dzeta*G_base_out%f%snorm_base(iMode))/REAL(SUM(G_base_out%f%Xmn(1:2,iMode)**2),wp)
    END DO

    __PERFOFF('project_G')

    !interpolate G and compute Gthet=iota*G+lambda
    ft_IP(:) =LA_IP(:)+iota_int*G_base_out%f%evalDOF_IP(0,GZ(is,:))
    !project onto Gthet base
    CALL G_base_out%f%projectIPtoDOF(.FALSE.,1.0_wp, 0,ft_IP(:),Gthet(is,:))
    DO iMode=1,modes
      Gthet(is,iMode)=Gthet(is,iMode)*dthet_dzeta*G_base_out%f%snorm_base(iMode)
    END DO

    __PERFOFF('project')
    CALL ProgressBar(is,nBase)
  END DO !is
  !transform back to corresponding representation of DOF in s
  DO iMode=1,modes
    GZ(   :,iMode)= G_base_out%s%initDOF( GZ(   :,iMode) )
    Gthet(:,iMode)= G_base_out%s%initDOF( Gthet(:,iMode) )
  END DO

  DEALLOCATE(LA_IP,dLAdthet_IP,dLAdzeta_IP,&
             fm_IP,fn_IP,ft_IP,   &
             Bcov_thet_IP,Bcov_zeta_IP)

  SWRITE(UNIT_StdOut,'(A)') '...DONE.'
  __PERFOFF('get_boozer')
END SUBROUTINE Get_Boozer


END MODULE MODgvec_get_boozer
