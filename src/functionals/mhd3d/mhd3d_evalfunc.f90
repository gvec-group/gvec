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
REAL(wp),ALLOCATABLE :: X1_GPIP(:,:)    !! evaluation of X1
REAL(wp),ALLOCATABLE :: X2_GPIP(:,:)    !! evaluation of X2
REAL(wp),ALLOCATABLE :: J_h(:,:)        !! Jacobian of the mapping h (X1,X2,zeta) -->(x,y,z) (global Jacobian: J_h*J_p) 
REAL(wp),ALLOCATABLE :: J_p(:,:)        !! Jacobian of poloidal mapping: dX1_ds*dX2_dtheta - dX2_ds*dX1_theta
REAL(wp),ALLOCATABLE :: dX1_ds(:,:)     !! radial derivative of X1
REAL(wp),ALLOCATABLE :: dX2_ds(:,:)     !! radial derivative of X2
REAL(wp),ALLOCATABLE :: dX1_dthet(:,:)  !! theta  derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dthet(:,:)  !! theta  derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dthet(:,:)  !! theta  derivative of lambda
REAL(wp),ALLOCATABLE :: dX1_dzeta(:,:)  !! zeta   derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dzeta(:,:)  !! zeta   derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dzeta(:,:)  !! zeta   derivative of lambda
REAL(wp),ALLOCATABLE :: b_a (  :,:,:)   !! b=(iota-dlamba_dzeta,1+dlambda_dtheta), normalized contravariant magnetic field 
                                        !! components at all points, size(1:2,1:mn_IP,1:nGP)
REAL(wp),ALLOCATABLE :: g_ab(:,:,:,:)   !! metric tensor g_(alpha,beta), alpha/beta=theta and zeta, size(1:2,1:2,1:mn_IP,1:nGP)


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
  ALLOCATE(X1_GPIP(   mn_IP,nGP) )
  ALLOCATE(X2_GPIP(   mn_IP,nGP) )
  ALLOCATE(dX1_ds(    mn_IP,nGP) )
  ALLOCATE(dX2_ds(    mn_IP,nGP) )
  ALLOCATE(dX1_dthet( mn_IP,nGP) )
  ALLOCATE(dX2_dthet( mn_IP,nGP) )
  ALLOCATE(dLA_dthet( mn_IP,nGP) )
  ALLOCATE(dX1_dzeta( mn_IP,nGP) )
  ALLOCATE(dX2_dzeta( mn_IP,nGP) )
  ALLOCATE(dLA_dzeta( mn_IP,nGP) )
  ALLOCATE(b_a (    2,mn_IP,nGP) )
  ALLOCATE(g_ab(  2,2,mn_IP,nGP) )
 
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
USE MOD_MHD3D_Profiles,ONLY: Eval_iota,Eval_pres,Eval_phiPrime
USE MOD_MHD3D_vars,ONLY: X1_base,X2_base,LA_base,hmap
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin   !! input solution 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iGP,i_mn,nGP,mn_IP
  REAL(wp)  :: qloc(3),q_thet(3),q_zeta(3)
  REAL(wp)  :: s_GP(X1_base%s%nGP),zeta_IP(X1_base%f%mn_IP)
!===================================================================================================================================
  nGP    = X1_base%s%nGP   ! same for all basis
  mn_IP  = X1_base%f%mn_IP ! same for all basis
  zeta_IP= X1_base%f%x_IP(2,:)
  s_GP   = X1_base%s%s_GP(:)
  X1_GPIP   = X1_base%evalDOF(      0   ,Uin%X1)
  X2_GPIP   = X2_base%evalDOF(      0   ,Uin%X2)
  dX1_ds    = X1_base%evalDOF(DERIV_S   ,Uin%X1)
  dX2_ds    = X2_base%evalDOF(DERIV_S   ,Uin%X2)
  dX1_dthet = X1_base%evalDOF(DERIV_THET,Uin%X1)
  dX2_dthet = X2_base%evalDOF(DERIV_THET,Uin%X2)
  dLA_dthet = LA_base%evalDOF(DERIV_THET,Uin%LA)
  dX1_dzeta = X1_base%evalDOF(DERIV_ZETA,Uin%X1)
  dX2_dzeta = X2_base%evalDOF(DERIV_ZETA,Uin%X2)
  dLA_dzeta = LA_base%evalDOF(DERIV_ZETA,Uin%LA)
  q_thet(3)=0.0_wp
  q_zeta(3)=1.0_wp
  DO iGP=1,nGP
    iota_GP(    iGP) = Eval_iota(    s_GP(iGP))
    pres_GP(    iGP) = Eval_pres(    s_GP(iGP))
    PhiPrime_GP(iGP) = Eval_PhiPrime(s_GP(iGP))
    DO i_mn=1,mn_IP
      J_p(  i_mn,iGP) = ( dX1_ds(   i_mn,iGP)*dX2_dthet(i_mn,iGP) &
                         -dX1_dthet(i_mn,iGP)*dX2_ds(   i_mn,iGP) )

      b_a(1,i_mn,iGP) = (iota_GP(iGP)- dLA_dzeta(i_mn,iGP))    !b_theta
      b_a(2,i_mn,iGP) = (1.0_wp      + dLA_dzeta(i_mn,iGP))    !b_zeta

      qloc(  1:3) = (/X1_GPIP(  i_mn,iGP),X2_GPIP(  i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:2) = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP)/)
      q_zeta(1:2) = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP)/)

      J_h(     i_mn,iGP) = hmap%eval_Jh(qloc)
      g_ab(1,1,i_mn,iGP) = hmap%eval_gij(q_thet,qloc,q_thet)   !g_theta,theta
      g_ab(2,1,i_mn,iGP) = hmap%eval_gij(q_zeta,qloc,q_thet)   !g_zeta,theta
      g_ab(1,2,i_mn,iGP) = g_ab(2,1,i_mn,iGP)                  
      g_ab(2,2,i_mn,iGP) = hmap%eval_gij(q_zeta,qloc,q_zeta)   !g_zeta,zeta
    END DO !i_IP
  END DO !iGP

END SUBROUTINE EvalAux

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
  SDEALLOCATE(X1_GPIP     )
  SDEALLOCATE(X2_GPIP     )
  SDEALLOCATE(dX1_ds      )
  SDEALLOCATE(dX2_ds      )
  SDEALLOCATE(dX1_dthet   )
  SDEALLOCATE(dX2_dthet   )
  SDEALLOCATE(dLA_dthet   )
  SDEALLOCATE(dX1_dzeta   )
  SDEALLOCATE(dX2_dzeta   )
  SDEALLOCATE(dLA_dzeta   )
  SDEALLOCATE(b_a         )
  SDEALLOCATE(g_ab        )

END SUBROUTINE FinalizeMHD3D_EvalFunc

END MODULE MOD_MHD3D_EvalFunc
