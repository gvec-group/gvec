!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **MHD3D Evalfunc**
!!
!! Evaluate the MHD3D functional and its derivative 
!!
!===================================================================================================================================
MODULE MOD_MHD3D_evalFunc
! MODULES
USE MOD_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep
IMPLICIT NONE
PUBLIC

!evaluations at radial gauss points, size(1:base%s%nGP)                                       
REAL(wp),ALLOCATABLE :: pres_GP(:)      !! mass profile 
REAL(wp),ALLOCATABLE :: iota_GP(:)      !! iota profile 
REAL(wp),ALLOCATABLE :: PhiPrime_GP(:)  !! s derivative of toroidal flux : Phi'(s)

!evaluations at all integration points, size(1:base%f%mn_IP,1:base%s%nGP)                                       
REAL(wp),ALLOCATABLE :: X1_IP_GP(:,:)   !! evaluation of X1
REAL(wp),ALLOCATABLE :: X2_IP_GP(:,:)   !! evaluation of X2
REAL(wp),ALLOCATABLE :: J_h(:,:)        !! Jacobian of the mapping h (X1,X2,zeta) -->(x,y,z) 
REAL(wp),ALLOCATABLE :: J_p(:,:)        !! Jacobian of poloidal mapping: dX1_ds*dX2_dtheta - dX2_ds*dX1_theta
REAL(wp),ALLOCATABLE :: detJ(:,:)       !! global Jacobian: detJ=sqrt(det g)=J_h*J_p ) 
REAL(wp),ALLOCATABLE :: dX1_ds(:,:)     !! radial derivative of X1
REAL(wp),ALLOCATABLE :: dX2_ds(:,:)     !! radial derivative of X2
REAL(wp),ALLOCATABLE :: dX1_dthet(:,:)  !! theta  derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dthet(:,:)  !! theta  derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dthet(:,:)  !! theta  derivative of lambda
REAL(wp),ALLOCATABLE :: dX1_dzeta(:,:)  !! zeta   derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dzeta(:,:)  !! zeta   derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dzeta(:,:)  !! zeta   derivative of lambda
REAL(wp),ALLOCATABLE :: b_thet(   :,:)  !! b_thet=(iota-dlamba_dzeta,1+dlambda_dtheta), normalized contravariant magnetic field 
REAL(wp),ALLOCATABLE :: b_zeta(   :,:)  !! b_zeta=1+dlambda_dtheta, normalized contravariant magnetic field 
REAL(wp),ALLOCATABLE :: sJ_bcov_thet(:,:)  !! covariant normalized magnetic field, scaled with 1/J:  
REAL(wp),ALLOCATABLE :: sJ_bcov_zeta(:,:)  !! sJ_bcov_alpha=1/detJ (g_{alpha,theta} b_theta + g_{alpha,zeta) b_zeta)
!REAL(wp),ALLOCATABLE :: g_tt(     :,:)  !! metric tensor g_(theta,theta)
!REAL(wp),ALLOCATABLE :: g_tz(     :,:)  !! metric tensor g_(theta,zeta )=g_(zeta,theta)
!REAL(wp),ALLOCATABLE :: g_zz(     :,:)  !! metric tensor g_(zeta ,zeta )

!private module variables, set in init
INTEGER                ,PRIVATE :: nElems
INTEGER                ,PRIVATE :: nGP
INTEGER                ,PRIVATE :: degGP
INTEGER                ,PRIVATE :: mn_IP
REAL(wp)               ,PRIVATE :: dthet_dzeta
REAL(wp),ALLOCATABLE   ,PRIVATE :: s_GP(:),w_GP(:),zeta_IP(:)
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitializeMHD3D_evalFunc() 
! MODULES
USE MOD_MHD3D_Vars,ONLY:X1_base
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D_EVALFUNC...'
  !same for all basis
  nElems  = X1_base%s%grid%nElems 
  nGP     = X1_base%s%nGP  
  degGP   = X1_base%s%degGP
  mn_IP   = X1_base%f%mn_IP
  dthet_dzeta  =X1_base%f%d_thet*X1_base%f%d_zeta
  ALLOCATE(s_GP(1:nGP),w_GP(1:nGP),zeta_IP(1:mn_IP))
  s_GP    = X1_base%s%s_GP(:)
  w_GP    = X1_base%s%w_GP(:)
  zeta_IP = X1_base%f%x_IP(2,:)

  ALLOCATE(pres_GP(           nGP) )
  ALLOCATE(iota_GP(           nGP) )
  ALLOCATE(PhiPrime_GP(       nGP) )
  ALLOCATE(J_h(         mn_IP,nGP) )
  ALLOCATE(J_p(         mn_IP,nGP) )
  ALLOCATE(detJ(        mn_IP,nGP) )
  ALLOCATE(X1_IP_GP(    mn_IP,nGP) )
  ALLOCATE(X2_IP_GP(    mn_IP,nGP) )
  ALLOCATE(dX1_ds(      mn_IP,nGP) )
  ALLOCATE(dX2_ds(      mn_IP,nGP) )
  ALLOCATE(dX1_dthet(   mn_IP,nGP) )
  ALLOCATE(dX2_dthet(   mn_IP,nGP) )
  ALLOCATE(dLA_dthet(   mn_IP,nGP) )
  ALLOCATE(dX1_dzeta(   mn_IP,nGP) )
  ALLOCATE(dX2_dzeta(   mn_IP,nGP) )
  ALLOCATE(dLA_dzeta(   mn_IP,nGP) )
  ALLOCATE(b_thet(      mn_IP,nGP) )
  ALLOCATE(b_zeta(      mn_IP,nGP) )
  ALLOCATE(sJ_bcov_thet(mn_IP,nGP) )
  ALLOCATE(sJ_bcov_zeta(mn_IP,nGP) )
!  ALLOCATE(g_tt(        mn_IP,nGP) )
!  ALLOCATE(g_tz(        mn_IP,nGP) )
!  ALLOCATE(g_zz(        mn_IP,nGP) )
 

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
  

END SUBROUTINE InitializeMHD3D_evalFunc

!===================================================================================================================================
!> Evaluate auxiliary variables at input state, writes onto module variables!!!
!!
!===================================================================================================================================
SUBROUTINE EvalAux(Uin)
! MODULES
USE MOD_Globals         , ONLY: n_warnings_occured
USE MOD_MHD3D_Profiles  , ONLY: Eval_iota,Eval_pres,Eval_phiPrime
USE MOD_MHD3D_vars      , ONLY: X1_base,X2_base,LA_base,hmap
USE MOD_sol_var_MHD3D   , ONLY: t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin   !! input solution 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iGP,i_mn,IP_GP(2)
  REAL(wp)  :: qloc(3),q_thet(3),q_zeta(3),g_tt_loc,g_tz_loc,g_zz_loc
!===================================================================================================================================

  !1D data at gauss-points
  DO iGP=1,nGP
    iota_GP(    iGP) = Eval_iota(    s_GP(iGP))
    pres_GP(    iGP) = Eval_pres(    s_GP(iGP))
    PhiPrime_GP(iGP) = Eval_PhiPrime(s_GP(iGP))
  END DO !iGP

  !2D data: interpolation points x gauss-points
  X1_IP_GP   = X1_base%evalDOF(     0   ,Uin%X1)
  X2_IP_GP   = X2_base%evalDOF(     0   ,Uin%X2)
  dX1_ds    = X1_base%evalDOF(DERIV_S   ,Uin%X1)
  dX2_ds    = X2_base%evalDOF(DERIV_S   ,Uin%X2)
  dX1_dthet = X1_base%evalDOF(DERIV_THET,Uin%X1)
  dX2_dthet = X2_base%evalDOF(DERIV_THET,Uin%X2)
  dLA_dthet = LA_base%evalDOF(DERIV_THET,Uin%LA)
  dX1_dzeta = X1_base%evalDOF(DERIV_ZETA,Uin%X1)
  dX2_dzeta = X2_base%evalDOF(DERIV_ZETA,Uin%X2)
  dLA_dzeta = LA_base%evalDOF(DERIV_ZETA,Uin%LA)

  q_thet(3)=0.0_wp !dq(3)/dtheta
  q_zeta(3)=1.0_wp !dq(3)/zeta
  DO iGP=1,nGP
    DO i_mn=1,mn_IP
      J_p(  i_mn,iGP) = ( dX1_ds(i_mn,iGP)*dX2_dthet(i_mn,iGP) &
                         -dX2_ds(i_mn,iGP)*dX1_dthet(i_mn,iGP) )

      qloc(  1:3) = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)

      J_h( i_mn,iGP) = hmap%eval_Jh(qloc)
      detJ(i_mn,iGP) = J_p(i_mn,iGP)*J_h(i_mn,iGP)
    END DO !i_mn
  END DO !iGP
  IF(MINVAL(detJ) .LT.1.0e-12) THEN
    n_warnings_occured=n_warnings_occured+1
    IP_GP= MINLOC(detJ(:,:))
    WRITE(UNIT_stdOut,'(4X,A8,I8,4(A,E11.3))')'WARNING ',n_warnings_occured, &
                                                 ' : min(J)= ',MINVAL(detJ),' at s= ',s_GP(IP_GP(2)), &
                                                                       ' theta= ',X1_base%f%x_IP(1,IP_GP(1)), &
                                                                        ' zeta= ',X1_base%f%x_IP(2,IP_GP(1)) 
    IP_GP= MAXLOC(detJ(:,:))
    WRITE(UNIT_stdOut,'(4X,16X,4(A,E11.3))')'     ...max(J)= ',MAXVAL(detJ),' at s= ',s_GP(IP_GP(2)), &
                                                                       ' theta= ',X1_base%f%x_IP(1,IP_GP(1)), &
                                                                        ' zeta= ',X1_base%f%x_IP(2,IP_GP(1)) 
    CALL abort(__STAMP__, &
        'EvalAux: Jacobian smaller that  1.0e-12!!!' )
  END IF
  DO iGP=1,nGP
    DO i_mn=1,mn_IP
      b_thet(i_mn,iGP) = (iota_GP(iGP)- dLA_dzeta(i_mn,iGP))    !b_theta
      b_zeta(i_mn,iGP) = (1.0_wp      + dLA_dthet(i_mn,iGP))    !b_zeta

      qloc(  1:3) = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:2) = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP)/) !dq(1:2)/dtheta
      q_zeta(1:2) = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP)/) !dq(1:2)/dzeta

      g_tt_loc               = hmap%eval_gij(q_thet,qloc,q_thet)   !g_theta,theta
      g_tz_loc               = hmap%eval_gij(q_thet,qloc,q_zeta)   !g_theta,zeta =g_zeta,theta
      g_zz_loc               = hmap%eval_gij(q_zeta,qloc,q_zeta)   !g_zeta,zeta
      sJ_bcov_thet(i_mn,iGP) = (g_tt_loc*b_thet(i_mn,iGP) + g_tz_loc*b_zeta(i_mn,iGP))/detJ(i_mn,iGP)
      sJ_bcov_zeta(i_mn,iGP) = (g_tz_loc*b_thet(i_mn,iGP) + g_zz_loc*b_zeta(i_mn,iGP))/detJ(i_mn,iGP)
    END DO !i_mn
  END DO !iGP

END SUBROUTINE EvalAux

!===================================================================================================================================
!> Evaluate 3D MHD energy
!! NOTE: set callEvalaux TRUE if not called before for the same Uin !!
!!
!===================================================================================================================================
FUNCTION EvalEnergy(Uin,callEvalAux) RESULT(W_MHD3D)
! MODULES
USE MOD_MHD3D_Vars, ONLY: s2mu_0,sgammM1
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin   !! input solution 
  LOGICAL               , INTENT(IN   ) :: callEvalAux !! call auxiliary variable computation on current Uin
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                              :: W_MHD3D   !! total integral of MHD3D energy 
                                                     !! W_MHD3D= int ( p/(gamma-1) + 1/(2mu_0) |B|^2) detJ ds dtheta dzeta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iGP
  REAL(wp) :: Wmag_GP(1:nGP)    !! magnetic energy at gauss points 
                                !! = 1/(dtheta*dzeta) * ( int [1/detJ * b_alpha*g_{alpha,beta}*b_beta]_iGP dtheta dzeta )
  REAL(wp) :: Vprime_GP(1:nGP)  !! =  1/(dtheta*dzeta) *( int detJ|_iGP ,dtheta dzeta)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'COMPUTE ENERGY...'
  IF(callEvalAux) CALL EvalAux(Uin)
  DO iGP=1,nGP
    Wmag_GP(iGP)   = SUM( b_thet(:,iGP)*sJ_bcov_thet(:,iGP)  &
                         +b_zeta(:,iGP)*sJ_bcov_zeta(:,iGP)  )
    Vprime_GP(iGP) = SUM(detJ(:,iGP))
  END DO !iGP

  W_MHD3D= dthet_dzeta* (  s2mu_0 *SUM(PhiPrime_GP(:)**2*Wmag_GP(:)*w_GP(:)) &
                         + sgammM1*SUM(    pres_GP(:) *Vprime_GP(:)*w_GP(:)) )
  
  SWRITE(UNIT_stdOut,'(A,E21.11)')'... DONE: ',W_MHD3D
  SWRITE(UNIT_stdOut,fmt_sep)
END FUNCTION EvalEnergy

!===================================================================================================================================
!> Evaluate the variation of the Energy = Force... F_j=-(D_U W(U))_j = -DW(u_h)*testfunc_j
!! NOTE: set callEvalaux TRUE if not called before for the same Uin !!
!!
!===================================================================================================================================
SUBROUTINE EvalForce(Uin,callEvalAux,F_MHD3D)
! MODULES
USE MOD_MHD3D_Vars, ONLY: X1_base,X2_base,LA_base,hmap,s2mu_0
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin   !! input solution 
  LOGICAL               , INTENT(IN   ) :: callEvalAux !! call auxiliary variable computation on current Uin
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: F_MHD3D   !! variation of the energy projected onto the basis functions of Uin 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: ibase,nBase,iMode,modes,iGP,i_mn,deg,iElem
  REAL(wp)  :: qloc(3),q_thet(3),q_zeta(3),y1,y1_thet(3),y1_zeta(3),y2,y2_thet(3),y2_zeta(3)
  REAL(wp)  :: F_GP   (1:nGP)    
  REAL(wp)  :: F_ds_GP(1:nGP)    
  REAL(wp)  :: dW(1:mn_IP, 1:nGP)        != p(s)+|Phi'(s)|^2 (b^alpha *g_{alpha,beta} *b^beta)/(2mu_0 *detJ^2)
  REAL(wp)  :: btt_sJ(1:mn_IP, 1:nGP)    != b^theta*b^theta/detJ 
  REAL(wp)  :: btz_sJ(1:mn_IP, 1:nGP)    != b^theta*b^zeta /detJ 
  REAL(wp)  :: bzz_sJ(1:mn_IP, 1:nGP)    != b^zeta *b^zeta /detJ 
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'COMPUTE FORCE...'
  IF(callEvalAux) CALL EvalAux(Uin)

  !additional auxiliary variables for X1 and X2 force
  DO iGP=1,nGP
    DO i_mn=1,mn_IP
      btt_sJ(i_mn,iGP)= b_thet( i_mn,iGP)/detJ(i_mn,iGP) ! *b_thet(i_mn,iGP)
      bzz_sJ(i_mn,iGP)= b_zeta( i_mn,iGP)/detJ(i_mn,iGP)  !*b_zeta(i_mn,iGP)
      dW(    i_mn,iGP)=  btt_sJ(i_mn,iGP)*sJ_bcov_thet(i_mn,iGP) &
                        +bzz_sJ(i_mn,iGP)*sJ_bcov_zeta(i_mn,iGP)
      btz_sJ(i_mn,iGP)= bzz_sJ( i_mn,iGP)*b_thet(i_mn,iGP)
      btt_sJ(i_mn,iGP)= btt_sJ( i_mn,iGP)*b_thet(i_mn,iGP)
      bzz_sJ(i_mn,iGP)= bzz_sJ( i_mn,iGP)*b_zeta(i_mn,iGP)
    END DO
    dW(:,iGP)=s2mu_0*(PhiPrime_GP(iGP)**2*dW(:,iGP))+pres_GP(iGP)
  END DO
  q_thet(3)=0.0_wp
  q_zeta(3)=1.0_wp

  ASSOCIATE(F_X1=>F_MHD3D%X1)
  nBase = X1_Base%s%nBase 
  modes = X1_Base%f%modes
  deg   = X1_base%s%deg
  F_X1(:,:)=0.0_wp
  Y1_thet(1:3)=(/0.0_wp,0.0_wp,0.0_wp/) ! test (Y1_dthet,0,0)
  Y1_zeta(1:3)=(/0.0_wp,0.0_wp,1.0_wp/) ! test (Y1_dzeta,0,1)
  DO iMode=1,modes
    DO iGP=1,nGP
      F_GP(iGP)=0.0_wp
      DO i_mn=1,mn_IP
        qloc(1:3)   =(/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
        q_thet(1:2) =(/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP)/)
        q_zeta(1:2) =(/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP)/)
        Y1           =X1_base%f%base_IP(i_mn,iMode)
        Y1_thet(1)   =X1_base%f%base_dthet_IP(i_mn,iMode)
        Y1_zeta(1)   =X1_base%f%base_dzeta_IP(i_mn,iMode)
        F_GP(iGP) = F_GP(iGP) &
          +dW(    i_mn,iGP)*(  J_p(i_mn,iGP)*hmap%eval_Jh_dq1(qloc(:))*Y1          &
                             - J_h(i_mn,iGP)*dX2_ds(i_mn,iGP)         *Y1_thet(1) )&
          -btt_sJ(i_mn,iGP)*( 2.0_wp*hmap%eval_gij(    q_thet,qloc,Y1_thet)        &
                                    +hmap%eval_gij_dq1(q_thet,qloc, q_thet)*Y1 )   &
          -bzz_sJ(i_mn,iGP)*( 2.0_wp*hmap%eval_gij(    q_zeta,qloc,Y1_zeta)        & 
                                    +hmap%eval_gij_dq1(q_zeta,qloc, q_zeta)*Y1 )   &
          -btz_sJ(i_mn,iGP)*2.0_wp*( hmap%eval_gij(    q_thet,qloc,Y1_zeta)        &
                                    +hmap%eval_gij(    q_zeta,qloc,Y1_thet)        &
                                    +hmap%eval_gij_dq1(q_thet,qloc, q_zeta)*Y1 ) 
      END DO !i_mn
    END DO !iGP
    F_GP(:)   = w_GP(:)*F_GP(:)    !*dthet_dzeta ... account below

    DO iGP=1,nGP
      F_ds_GP(iGP)=0.0_wp
      DO i_mn=1,mn_IP
        F_ds_GP(iGP) = F_ds_GP(iGP)    &
                       + dW(i_mn,iGP)*J_h(i_mn,iGP)*dX2_dthet(i_mn,iGP)*X1_base%f%base_IP(i_mn,iMode) !Y1
      END DO !i_mn
    END DO !iGP
    F_ds_GP(:)= w_GP(:)*F_ds_GP(:)

    iGP=1
    DO iElem=1,nElems
      ibase=X1_base%s%base_offset(iElem)
      F_X1(iBase:iBase+deg,iMode) = F_X1(iBase:iBase+deg,iMode) &
                                    + MATMUL(F_GP(   iGP:iGP+degGP),X1_base%s%base_GP(   0:degGP,0:deg,iElem)) &
                                    + MATMUL(F_ds_GP(iGP:iGP+degGP),X1_base%s%base_ds_GP(0:degGP,0:deg,iElem)) 
      iGP=iGP+(degGP+1)
    END DO !iElem
  END DO !iMode
  F_X1(:,:)=F_X1(:,:)*dthet_dzeta !scale with constants
  END ASSOCIATE !F_X1

  ASSOCIATE(F_X2=>F_MHD3D%X2)
  nBase = X2_base%s%nBase 
  modes = X2_base%f%modes
  deg   = X2_base%s%deg
  F_X2(:,:)=0.0_wp
  Y2_thet(1:3)=(/0.0_wp,0.0_wp,0.0_wp/) ! test (0,Y2_dthet,0)
  Y2_zeta(1:3)=(/0.0_wp,0.0_wp,1.0_wp/) ! test (0,Y2_dzeta,1)
  DO iMode=1,modes
    DO iGP=1,nGP
      F_GP(iGP)=0.0_wp
      DO i_mn=1,mn_IP
        qloc(1:3)   =(/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
        q_thet(1:2) =(/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP)/)
        q_zeta(1:2) =(/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP)/)
        Y2           =X2_base%f%base_IP(i_mn,iMode)
        Y2_thet(2)   =X2_base%f%base_dthet_IP(i_mn,iMode)
        Y2_zeta(2)   =X2_base%f%base_dzeta_IP(i_mn,iMode)
        F_GP(iGP) = F_GP(iGP) &
          +dW(    i_mn,iGP)*(  J_p(i_mn,iGP)*hmap%eval_Jh_dq2(qloc(:))*Y2          &
                             + J_h(i_mn,iGP)*dX1_ds(i_mn,iGP)         *Y2_thet(2) )&
          -btt_sJ(i_mn,iGP)*( 2.0_wp*hmap%eval_gij(    q_thet,qloc,Y2_thet)        &
                                    +hmap%eval_gij_dq2(q_thet,qloc, q_thet)*Y2 )   &
          -bzz_sJ(i_mn,iGP)*( 2.0_wp*hmap%eval_gij(    q_zeta,qloc,Y2_zeta)        & 
                                    +hmap%eval_gij_dq2(q_zeta,qloc, q_zeta)*Y2 )   &
          -btz_sJ(i_mn,iGP)*2.0_wp*( hmap%eval_gij(    q_thet,qloc,Y2_zeta)        &
                                    +hmap%eval_gij(    q_zeta,qloc,Y2_thet)        &
                                    +hmap%eval_gij_dq2(q_thet,qloc, q_zeta)*Y2 ) 
      END DO !i_mn
    END DO !iGP
    F_GP(:)   = w_GP(:)*F_GP(:)    !*dthet_dzeta ... account below

    DO iGP=1,nGP
      F_ds_GP(iGP)=0.0_wp
      DO i_mn=1,mn_IP
        F_ds_GP(iGP) = F_ds_GP(iGP)    &
                       - dW(i_mn,iGP)*J_h(i_mn,iGP)*dX1_dthet(i_mn,iGP)*X2_base%f%base_IP(i_mn,iMode) !Y2
      END DO !i_mn
    END DO !iGP
    F_ds_GP(:)= w_GP(:)*F_ds_GP(:)

    iGP=1
    DO iElem=1,nElems
      ibase=X2_base%s%base_offset(iElem)
      F_X2(iBase:iBase+deg,iMode) = F_X2(iBase:iBase+deg,iMode) &
                                    + MATMUL(F_GP(   iGP:iGP+degGP),X2_base%s%base_GP(     0:degGP,0:deg,iElem))      &
                                    + MATMUL(F_ds_GP(iGP:iGP+degGP),X2_base%s%base_ds_GP(  0:degGP,0:deg,iElem)) 
      iGP=iGP+(degGP+1)
    END DO !iElem
  END DO !iMode
  F_X2(:,:)=F_X2(:,:)*dthet_dzeta !scale with constants
  END ASSOCIATE !F_X2


  ASSOCIATE(F_LA=>F_MHD3D%LA)

  nBase = LA_base%s%nBase
  modes = LA_base%f%modes 
  deg   = LA_base%s%deg
 
  F_LA(:,:)=0.0_wp
  DO iMode=1,modes
    DO iGP=1,nGP
      F_GP(iGP)=0.0_wp
      DO i_mn=1,mn_IP
        F_GP(iGP)=F_GP(iGP)+ sJ_bcov_thet(i_mn,iGP)*LA_base%f%base_dzeta_IP(i_mn,iMode) &
                           - sJ_bcov_zeta(i_mn,iGP)*LA_base%f%base_dthet_IP(i_mn,iMode)  
      END DO !i_mn
    END DO !iGP
    F_GP(:)= (PhiPrime_GP(:)**2*w_GP(:)*F_GP(:))  !*2*dthet_dzeta ... account below
    iGP=1
    DO iElem=1,nElems
      ibase=LA_base%s%base_offset(iElem)
      F_LA(iBase:iBase+deg,iMode) = F_LA(iBase:iBase+deg,iMode) &
                                    + MATMUL(F_GP(iGP:iGP+degGP),LA_base%s%base_GP(0:degGP,0:deg,iElem))
      iGP=iGP+(degGP+1)
    END DO !iElem

  END DO !iMode
  F_LA(:,:)=F_LA(:,:)*(2.0_wp*s2mu_0*dthet_dzeta) !scale with constants


  END ASSOCIATE !F_LA

  SWRITE(UNIT_stdOut,'(A,3E21.11)')'Norm of force |X1|^2,|X2|^2,|LA|^2: ',F_MHD3D%norm_2()
  SWRITE(UNIT_stdOut,'(A,E21.11)')'... DONE: '
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE EvalForce

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeMHD3D_EvalFunc() 
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  nElems=-1
  degGP =-1
  nGP   =-1
  mn_IP =-1
  dthet_dzeta=0.0_wp
  SDEALLOCATE(s_GP         )
  SDEALLOCATE(w_GP         )
  SDEALLOCATE(zeta_IP      )
  SDEALLOCATE(pres_GP      )
  SDEALLOCATE(iota_GP      )
  SDEALLOCATE(PhiPrime_GP  )
  SDEALLOCATE(J_h          )
  SDEALLOCATE(J_p          )
  SDEALLOCATE(detJ         )
  SDEALLOCATE(X1_IP_GP     )
  SDEALLOCATE(X2_IP_GP     )
  SDEALLOCATE(dX1_ds       )
  SDEALLOCATE(dX2_ds       )
  SDEALLOCATE(dX1_dthet    )
  SDEALLOCATE(dX2_dthet    )
  SDEALLOCATE(dLA_dthet    )
  SDEALLOCATE(dX1_dzeta    )
  SDEALLOCATE(dX2_dzeta    )
  SDEALLOCATE(dLA_dzeta    )
  SDEALLOCATE(b_thet       )
  SDEALLOCATE(b_zeta       )
  SDEALLOCATE(sJ_bcov_thet )
  SDEALLOCATE(sJ_bcov_zeta )
!  SDEALLOCATE(g_tt         )
!  SDEALLOCATE(g_tz         )
!  SDEALLOCATE(g_zz         )

END SUBROUTINE FinalizeMHD3D_EvalFunc

END MODULE MOD_MHD3D_EvalFunc
