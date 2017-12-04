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
REAL(wp),ALLOCATABLE :: g_tt(     :,:)  !! metric tensor g_(theta,theta)
REAL(wp),ALLOCATABLE :: g_tz(     :,:)  !! metric tensor g_(theta,zeta )=g_(zeta,theta)
REAL(wp),ALLOCATABLE :: g_zz(     :,:)  !! metric tensor g_(zeta ,zeta )


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

  ASSOCIATE( nGP   => X1_base%s%nGP, & 
             mn_IP => X1_base%f%mn_IP  ) !same for all variables

  ALLOCATE(pres_GP(         nGP) )
  ALLOCATE(iota_GP(         nGP) )
  ALLOCATE(PhiPrime_GP(     nGP) )
  ALLOCATE(J_h(       mn_IP,nGP) )
  ALLOCATE(J_p(       mn_IP,nGP) )
  ALLOCATE(detJ(      mn_IP,nGP) )
  ALLOCATE(X1_IP_GP(  mn_IP,nGP) )
  ALLOCATE(X2_IP_GP(  mn_IP,nGP) )
  ALLOCATE(dX1_ds(    mn_IP,nGP) )
  ALLOCATE(dX2_ds(    mn_IP,nGP) )
  ALLOCATE(dX1_dthet( mn_IP,nGP) )
  ALLOCATE(dX2_dthet( mn_IP,nGP) )
  ALLOCATE(dLA_dthet( mn_IP,nGP) )
  ALLOCATE(dX1_dzeta( mn_IP,nGP) )
  ALLOCATE(dX2_dzeta( mn_IP,nGP) )
  ALLOCATE(dLA_dzeta( mn_IP,nGP) )
  ALLOCATE(b_thet(    mn_IP,nGP) )
  ALLOCATE(b_zeta(    mn_IP,nGP) )
  ALLOCATE(g_tt(      mn_IP,nGP) )
  ALLOCATE(g_tz(      mn_IP,nGP) )
  ALLOCATE(g_zz(      mn_IP,nGP) )
 
  END ASSOCIATE !mn_IP,nGP

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
  

END SUBROUTINE InitializeMHD3D_evalFunc

!===================================================================================================================================
!> Evaluate auxiliary variables at input state
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
  INTEGER   :: iGP,i_mn,nGP,mn_IP,IP_GP(2)
  REAL(wp)  :: qloc(3),q_thet(3),q_zeta(3)
  REAL(wp)  :: s_GP(   X1_base%s%nGP  )
  REAL(wp)  :: zeta_IP(X1_base%f%mn_IP)
!===================================================================================================================================
  nGP     = X1_base%s%nGP   ! same for all basis
  mn_IP   = X1_base%f%mn_IP ! same for all basis
  zeta_IP = X1_base%f%x_IP(2,:)
  s_GP    = X1_base%s%s_GP(:)
  X1_IP_GP   = X1_base%evalDOF(      0   ,Uin%X1)
  X2_IP_GP   = X2_base%evalDOF(      0   ,Uin%X2)
  dX1_ds    = X1_base%evalDOF(DERIV_S   ,Uin%X1)
  dX2_ds    = X2_base%evalDOF(DERIV_S   ,Uin%X2)
  dX1_dthet = X1_base%evalDOF(DERIV_THET,Uin%X1)
  dX2_dthet = X2_base%evalDOF(DERIV_THET,Uin%X2)
  dLA_dthet = LA_base%evalDOF(DERIV_THET,Uin%LA)
  dX1_dzeta = X1_base%evalDOF(DERIV_ZETA,Uin%X1)
  dX2_dzeta = X2_base%evalDOF(DERIV_ZETA,Uin%X2)
  dLA_dzeta = LA_base%evalDOF(DERIV_ZETA,Uin%LA)
  DO iGP=1,nGP
    iota_GP(    iGP) = Eval_iota(    s_GP(iGP))
    pres_GP(    iGP) = Eval_pres(    s_GP(iGP))
    PhiPrime_GP(iGP) = Eval_PhiPrime(s_GP(iGP))
  END DO !iGP
  q_thet(3)=0.0_wp
  q_zeta(3)=1.0_wp
  DO iGP=1,nGP
    DO i_mn=1,mn_IP
      J_p(  i_mn,iGP) = ( dX1_ds(i_mn,iGP)*dX2_dthet(i_mn,iGP) &
                         -dX2_ds(i_mn,iGP)*dX1_dthet(i_mn,iGP) )

      b_thet(i_mn,iGP) = (iota_GP(iGP)- dLA_dzeta(i_mn,iGP))    !b_theta
      b_zeta(i_mn,iGP) = (1.0_wp      + dLA_dthet(i_mn,iGP))    !b_zeta

      qloc(  1:3) = (/ X1_IP_GP( i_mn,iGP), X2_IP_GP( i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:2) = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP)/)
      q_zeta(1:2) = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP)/)

      J_h( i_mn,iGP) = hmap%eval_Jh(qloc)
      detJ(i_mn,iGP) = J_p(i_mn,iGP)*J_h(i_mn,iGP)
      g_tt(i_mn,iGP) = hmap%eval_gij(q_thet,qloc,q_thet)   !g_theta,theta
      g_tz(i_mn,iGP) = hmap%eval_gij(q_thet,qloc,q_zeta)   !g_theta,zeta =g_zeta,theta
      g_zz(i_mn,iGP) = hmap%eval_gij(q_zeta,qloc,q_zeta)   !g_zeta,zeta
    END DO !i_IP
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

END SUBROUTINE EvalAux

!===================================================================================================================================
!> Evaluate 3D MHD energy
!! NOTE: that auxiliary variables should be called before !!
!!
!===================================================================================================================================
FUNCTION EvalEnergy(Uin,callEvalAux) RESULT(W_MHD3D)
! MODULES
USE MOD_MHD3D_Vars, ONLY: X1_base,s2mu_0,sgammM1
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin   !! input solution 
  LOGICAL               , INTENT(IN   ) :: callEvalAux !! call auxiliary variable computation on current Uin
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)  :: w_MHD3D
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iGP,i_mn,nGP,mn_IP
  REAL(wp)  :: dthet_dzeta
  REAL(wp)  :: w_GP(     X1_base%s%nGP)    !! gauss weight
  REAL(wp)  :: Wmag_GP(  X1_base%s%nGP)    !! magnetic energy at gauss points 
                                           !! = 1/(dtheta*dzeta) * ( int [1/detJ * b_alpha*g_alpha,beta*b_beta]_iGP,dtheta,dzeta)
  REAL(wp)  :: v_prime_GP(X1_base%s%nGP)   !! =  1/(dtheta*dzeta) *( int detJ|_iGP ,dtheta dzeta)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'COMPUTE ENERGY...'
  IF(callEvalAux) CALL EvalAux(Uin)
  nGP    = X1_base%s%nGP   ! same for all basis
  mn_IP  = X1_base%f%mn_IP ! same for all basis
  w_GP   = X1_base%s%w_GP(:)
  dthet_dzeta  =X1_base%f%d_thet*X1_base%f%d_zeta
  Wmag_GP=0.0_wp
  v_prime_GP=0.0_wp
  DO iGP=1,nGP
    DO i_mn=1,mn_IP
      v_prime_GP(iGP) = v_prime_GP(iGP)+ detJ(i_mn,iGP)
      Wmag_GP(iGP)    = Wmag_GP(iGP) +  & 
                         ( b_thet(i_mn,iGP)*(        b_thet(i_mn,iGP)*g_tt(i_mn,iGP)     &
                                             +2.0_wp*b_zeta(i_mn,iGP)*g_tz(i_mn,iGP) )   &
                          +b_zeta(i_mn,iGP)*(        b_zeta(i_mn,iGP)*g_zz(i_mn,iGP))  )  / detJ(i_mn,iGP)
    END DO !i_mn
  END DO !iGP
  W_MHD3D= dthet_dzeta* (   s2mu_0 *SUM(PhiPrime_GP(:)**2*Wmag_GP(:)*w_GP(:)) &
                          + sgammM1*SUM(pres_GP(:)*v_prime_GP(:)    *w_GP(:)) )
  
  SWRITE(UNIT_stdOut,'(A,E21.11)')'... DONE: ',W_MHD3D
  SWRITE(UNIT_stdOut,fmt_sep)
END FUNCTION EvalEnergy

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
  SDEALLOCATE(pres_GP     )
  SDEALLOCATE(iota_GP     )
  SDEALLOCATE(PhiPrime_GP )
  SDEALLOCATE(J_h         )
  SDEALLOCATE(J_p         )
  SDEALLOCATE(detJ        )
  SDEALLOCATE(X1_IP_GP    )
  SDEALLOCATE(X2_IP_GP    )
  SDEALLOCATE(dX1_ds      )
  SDEALLOCATE(dX2_ds      )
  SDEALLOCATE(dX1_dthet   )
  SDEALLOCATE(dX2_dthet   )
  SDEALLOCATE(dLA_dthet   )
  SDEALLOCATE(dX1_dzeta   )
  SDEALLOCATE(dX2_dzeta   )
  SDEALLOCATE(dLA_dzeta   )
  SDEALLOCATE(b_thet      )
  SDEALLOCATE(b_zeta      )
  SDEALLOCATE(g_tt        )
  SDEALLOCATE(g_tz        )
  SDEALLOCATE(g_zz        )

END SUBROUTINE FinalizeMHD3D_EvalFunc

END MODULE MOD_MHD3D_EvalFunc
