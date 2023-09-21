!===================================================================================================================================
! Copyright (C) 2017 - 2022  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2021 - 2022  Tiago Ribeiro
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
MODULE MODgvec_MHD3D_evalFunc
  ! MODULES
  USE MODgvec_Globals,            ONLY : wp,abort,UNIT_stdOut,fmt_sep,MPIRoot
  USE sll_m_spline_matrix,        ONLY : sll_c_spline_matrix !for precond
  USE sll_m_spline_matrix_banded, ONLY : sll_t_spline_matrix_banded
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC::InitializeMHD3D_EvalFunc, InitProfilesGP, & 
          EvalEnergy, EvalForce,CheckEvalForce, EvalAux, EvalTotals, FinalizeMHD3D_EvalFunc
  
  !evaluations at radial gauss points, size(base%s%nGP_str:base%s%nGP_end)
  REAL(wp),ALLOCATABLE :: pres_GP(:)      !! mass profile 
  REAL(wp),ALLOCATABLE :: chiPrime_GP(:)  !! s derivative of poloidal flux 
  REAL(wp),ALLOCATABLE :: phiPrime_GP(:)  !! s derivative of toroidal flux 
  REAL(wp),ALLOCATABLE :: phiPrime2_GP(:) !! s derivative of toroidal flux : |Phi'(s)|^2
  
  !evaluations at all integration points, size(1:base%f%mn_IP,base%s%nGP_str:base%s%nGP_end)                                       
  REAL(wp),ALLOCATABLE :: X1_IP_GP(:,:)   !! evaluation of X1
  REAL(wp),ALLOCATABLE :: X2_IP_GP(:,:)   !! evaluation of X2
  REAL(wp),ALLOCATABLE :: J_h(:,:)        !! Jacobian of the mapping h (X1,X2,zeta) -->(x,y,z) 
  REAL(wp),ALLOCATABLE :: J_p(:,:)        !! Jacobian of poloidal mapping: dX1_ds*dX2_dtheta - dX2_ds*dX1_theta
  REAL(wp),ALLOCATABLE :: sJ_h(:,:)       !! 1/J_h
  REAL(wp),ALLOCATABLE :: sJ_p(:,:)       !! 1/J_p
  REAL(wp),ALLOCATABLE :: detJ(:,:)       !! global Jacobian: detJ=sqrt(det g)=J_h*J_p ) 
  REAL(wp),ALLOCATABLE :: sdetJ(:,:)      !! 1/detJ
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
  REAL(wp),ALLOCATABLE :: bbcov_sJ(:,:)      !! (b^alpha*g_{alpha,beta}*b^beta)/(detJ)
  REAL(wp),ALLOCATABLE :: g_tt(     :,:)     !! metric tensor g_(theta,theta)
  REAL(wp),ALLOCATABLE :: g_tz(     :,:)     !! metric tensor g_(theta,zeta )=g_(zeta,theta)
  REAL(wp),ALLOCATABLE :: g_zz(     :,:)     !! metric tensor g_(zeta ,zeta )
  
  INTEGER                         :: nGP
  INTEGER                         :: nGP_str, nGP_end !< for MPI
  INTEGER                         :: mn_IP
  REAL(wp)                        :: dthet_dzeta
  REAL(wp),ALLOCATABLE            :: w_GP(:)
  !private module variables, set in init
  INTEGER                ,PRIVATE :: nElems
  INTEGER                         :: nElems_str,nElems_end !< for MPI
  INTEGER                ,PRIVATE :: degGP
  REAL(wp),ALLOCATABLE   ,PRIVATE :: s_GP(:),zeta_IP(:)

  !FOR PRECONDITIONER
  REAL(wp),CONTIGUOUS,POINTER :: DX1_tt(:), DX1_tz(:), DX1_zz(:), DX1(:), DX1_ss(:)   
  REAL(wp),CONTIGUOUS,POINTER :: DX2_tt(:), DX2_tz(:), DX2_zz(:), DX2(:), DX2_ss(:)  
  REAL(wp),CONTIGUOUS,POINTER :: DLA_tt(:), DLA_tz(:), DLA_zz(:)
  REAL(wp),CONTIGUOUS,POINTER :: D_buf(:,:)    !! 2d array container for all 1d array abpove (is the one allocated)

  CLASS(sll_c_spline_matrix),PRIVATE,ALLOCATABLE :: precond_X1(:)  !! container for preconditioner matrices
  CLASS(sll_c_spline_matrix),PRIVATE,ALLOCATABLE :: precond_X2(:)  !! container for preconditioner matrices
  CLASS(sll_c_spline_matrix),PRIVATE,ALLOCATABLE :: precond_LA(:)  !! container for preconditioner matrices
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitializeMHD3D_evalFunc() 
! MODULES
  USE MODgvec_MHD3D_Vars,ONLY:X1_base,X2_base,LA_base,PrecondType
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iMode
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D_EVALFUNC...'
  !same for all basis
  nElems  = X1_base%s%grid%nElems
  nElems_str = X1_base%s%grid%nElems_str
  nElems_end = X1_base%s%grid%nElems_end
  nGP     = X1_base%s%nGP
  nGP_str = X1_base%s%nGP_str  !< for MPI
  nGP_end = X1_base%s%nGP_end  !< for MPI
  degGP   = X1_base%s%degGP
  mn_IP   = X1_base%f%mn_IP
  dthet_dzeta  =X1_base%f%d_thet*X1_base%f%d_zeta
  ALLOCATE(s_GP(1:nGP),w_GP(1:nGP),zeta_IP(1:mn_IP)) 
  s_GP    = X1_base%s%s_GP(1:nGP)
  w_GP    = X1_base%s%w_GP(1:nGP)
  zeta_IP = X1_base%f%x_IP(2,:)

  ALLOCATE(pres_GP(     1:nGP) )
  ALLOCATE(chiPrime_GP( 1:nGP) )
  ALLOCATE(phiPrime_GP( 1:nGP) )
  ALLOCATE(phiPrime2_GP(1:nGP) )
  ALLOCATE(J_h(         mn_IP,nGP_str:nGP_end) )
  ALLOCATE(J_p(         mn_IP,nGP_str:nGP_end) )
  ALLOCATE(sJ_h(        mn_IP,nGP_str:nGP_end) )
  ALLOCATE(sJ_p(        mn_IP,nGP_str:nGP_end) )
  ALLOCATE(detJ(        mn_IP,nGP_str:nGP_end) )
  ALLOCATE(sdetJ(       mn_IP,nGP_str:nGP_end) )
  ALLOCATE(X1_IP_GP(    mn_IP,nGP_str:nGP_end) )
  ALLOCATE(X2_IP_GP(    mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dX1_ds(      mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dX2_ds(      mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dX1_dthet(   mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dX2_dthet(   mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dLA_dthet(   mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dX1_dzeta(   mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dX2_dzeta(   mn_IP,nGP_str:nGP_end) )
  ALLOCATE(dLA_dzeta(   mn_IP,nGP_str:nGP_end) )
  ALLOCATE(b_thet(      mn_IP,nGP_str:nGP_end) )
  ALLOCATE(b_zeta(      mn_IP,nGP_str:nGP_end) )
  ALLOCATE(sJ_bcov_thet(mn_IP,nGP_str:nGP_end) )
  ALLOCATE(sJ_bcov_zeta(mn_IP,nGP_str:nGP_end) )
  ALLOCATE(bbcov_sJ    (mn_IP,nGP_str:nGP_end) )
  ALLOCATE(g_tt(        mn_IP,nGP_str:nGP_end) )
  ALLOCATE(g_tz(        mn_IP,nGP_str:nGP_end) )
  ALLOCATE(g_zz(        mn_IP,nGP_str:nGP_end) )
 
  IF(PrecondType.GT.0)THEN
    !WHEN CHANGED TO ALLGATHERV COMM IN BUILDPRECOND, THIS ALLOCATE WILL BE THE SAME.
    ! POINTERS HELP TO GATHER ALL DATA IN ONE ARRAY (buf)
    ALLOCATE(D_buf(1:nGP,13))
    ! this is where the "pointer allocation" occurs
    DX1_tt(1:nGP) => D_buf(1:nGP,1)
    DX1_tz(1:nGP) => D_buf(1:nGP,2)
    DX1_zz(1:nGP) => D_buf(1:nGP,3)
    DX1(1:nGP)    => D_buf(1:nGP,4)
    DX1_ss(1:nGP) => D_buf(1:nGP,5)
    
    DX2_tt(1:nGP) => D_buf(1:nGP,6)
    DX2_tz(1:nGP) => D_buf(1:nGP,7)
    DX2_zz(1:nGP) => D_buf(1:nGP,8)
    DX2(1:nGP)    => D_buf(1:nGP,9)
    DX2_ss(1:nGP) => D_buf(1:nGP,10)
    
    DLA_tt(1:nGP) => D_buf(1:nGP,11)
    DLA_tz(1:nGP) => D_buf(1:nGP,12)
    DLA_zz(1:nGP) => D_buf(1:nGP,13)

    !distribute the preconditioner per mode over the MPI tasks (modes_str:modes_end)
    ALLOCATE(sll_t_spline_matrix_banded :: precond_X1( X1_Base%f%modes_str:X1_base%f%modes_end))
    SELECT TYPE(precond_X1); TYPE IS(sll_t_spline_matrix_banded)
      DO iMode=X1_Base%f%modes_str,X1_base%f%modes_end
        CALL precond_X1(iMode)%init(X1_Base%s%nBase,X1_Base%s%deg,X1_Base%s%deg)
      END DO !iMode
    END SELECT !TYPE
    ALLOCATE(sll_t_spline_matrix_banded :: precond_X2( X2_Base%f%modes_str:X2_base%f%modes_end))
    SELECT TYPE(precond_X2); TYPE IS(sll_t_spline_matrix_banded)
      DO iMode=X2_Base%f%modes_str,X2_base%f%modes_end
        CALL precond_X2(iMode)%init(X2_Base%s%nBase,X2_Base%s%deg,X2_Base%s%deg)
      END DO !iMode
    END SELECT !TYPE
    ALLOCATE(sll_t_spline_matrix_banded :: precond_LA( LA_Base%f%modes_str:LA_base%f%modes_end))
    SELECT TYPE(precond_LA); TYPE IS(sll_t_spline_matrix_banded)
      DO iMode=LA_Base%f%modes_str,LA_base%f%modes_end
        CALL precond_LA(iMode)%init(LA_Base%s%nBase,LA_Base%s%deg,LA_Base%s%deg)
      END DO !iMode
    END SELECT !TYPE
  END IF !PrecondType>0

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
  

END SUBROUTINE InitializeMHD3D_evalFunc

!===================================================================================================================================
!> Initialise Profiles at GP!!!
!!
!===================================================================================================================================
SUBROUTINE InitProfilesGP()
! MODULES
  USE MODgvec_MPI             , ONLY: par_Bcast
  USE MODgvec_MHD3D_Profiles  , ONLY: Eval_pres,Eval_chiPrime,Eval_phiPrime
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                     :: iGP
!===================================================================================================================================
IF(MPIroot)THEN
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iGP) FIRSTPRIVATE(nGP)
  DO iGP=1,nGP
    chiPrime_GP( iGP) = Eval_chiPrime(s_GP(iGP))
    pres_GP(     iGP) = Eval_pres(    s_GP(iGP))
    PhiPrime_GP( iGP) = Eval_PhiPrime(s_GP(iGP))
    PhiPrime2_GP(iGP) = PhiPrime_GP(       iGP)**2
  END DO !iGP
!$OMP END PARALLEL DO
END  IF !MPIroot

CALL par_Bcast(chiPrime_GP,0)
CALL par_Bcast(pres_GP,0)
CALL par_Bcast(PhiPrime_GP,0)
CALL par_Bcast(PhiPrime2_GP,0)

END SUBROUTINE InitProfilesGP


!===================================================================================================================================
!> Evaluate auxiliary variables at input state, writes onto module variables!!!
!!
!===================================================================================================================================
SUBROUTINE EvalAux(Uin,JacCheck)
! MODULES
  USE MODgvec_MPI             , ONLY: par_AllReduce
  USE MODgvec_Globals         , ONLY: n_warnings_occured,myRank
  USE MODgvec_MHD3D_vars      , ONLY: X1_base,X2_base,LA_base,hmap
  USE MODgvec_sol_var_MHD3D   , ONLY: t_sol_var_MHD3D
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin      !! input solution 
  INTEGER               , INTENT(INOUT) :: JacCheck !! if 1 on input: abort if detJ<0. 
                                                    !! if 2 on input, no abort, unchanged if detJ>0 ,return -1 if detJ<=0
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iGP,i_mn,IP_GP(2)
  REAL(wp)  :: qloc(3),q_thet(3),q_zeta(3),min_detJ
!===================================================================================================================================
  IF(JacCheck.EQ.-1) THEN
      CALL abort(__STAMP__, &
          'You already called EvalAux, with a Jacobian smaller that  1.0e-12!!!' )
  END IF

  __PERFON('EvalAux')

  __PERFON('EvalDOF_1')
  !2D data: interpolation points x gauss-points
  CALL X1_base%evalDOF((/0,0/)         , Uin%X1, X1_IP_GP  )
  CALL X1_base%evalDOF((/DERIV_S,0/)   , Uin%X1, dX1_ds    )
  CALL X1_base%evalDOF((/0,DERIV_THET/), Uin%X1, dX1_dthet )
  CALL X2_base%evalDOF((/0,0/)         , Uin%X2, X2_IP_GP  )
  CALL X2_base%evalDOF((/DERIV_S,0/)   , Uin%X2, dX2_ds    )
  CALL X2_base%evalDOF((/0,DERIV_THET/), Uin%X2, dX2_dthet )
  __PERFOFF('EvalDOF_1')

  __PERFON('loop_1')
  min_detJ =HUGE(1.0_wp)
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(iGP,i_mn,qloc)  &
!$OMP   FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)  &
!$OMP   REDUCTION(min:min_detJ) &
!$OMP   SHARED(J_p,J_h,detJ,dX1_ds,dX2_dthet,dX2_ds,dX1_dthet,X1_IP_GP,X2_IP_GP,zeta_IP,hmap)
  DO iGP=nGP_str,nGP_end
    DO i_mn=1,mn_IP
      J_p(  i_mn,iGP) = ( dX1_ds(i_mn,iGP)*dX2_dthet(i_mn,iGP) &
                         -dX2_ds(i_mn,iGP)*dX1_dthet(i_mn,iGP) )

      qloc(  1:3) = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)

      J_h( i_mn,iGP) = hmap%eval_Jh(qloc)
      detJ(i_mn,iGP) = J_p(i_mn,iGP)*J_h(i_mn,iGP)
      min_detJ = MIN(min_detJ,detJ(i_mn,iGP))
    END DO !i_mn
  END DO !iGP
!$OMP END PARALLEL DO


  __PERFOFF('loop_1')

  !check Jacobian
  IF(min_detJ .LT.1.0e-12) THEN
    SELECT CASE(JacCheck)
    CASE(1)
      n_warnings_occured=n_warnings_occured+1
      IP_GP= MINLOC(detJ(:,nGP_str:nGP_end))
      WRITE(UNIT_stdOut,'(4X,A8,I8,4(A,E11.3))')'WARNING ',n_warnings_occured, &
           &                                       ' : min(J)= ',MINVAL(detJ(:,nGP_str:nGP_end)),' at s= ',s_GP(IP_GP(2)), &
           &                                                             ' theta= ',X1_base%f%x_IP(1,IP_GP(1)), &
           &                                                              ' zeta= ',X1_base%f%x_IP(2,IP_GP(1)) 
      IP_GP= MAXLOC(detJ(:,:))
      WRITE(UNIT_stdOut,'(4X,16X,4(A,E11.3))')'     ...max(J)= ',MAXVAL(detJ(:,nGP_str:nGP_end)),' at s= ',s_GP(IP_GP(2)), &
           &                                                             ' theta= ',X1_base%f%x_IP(1,IP_GP(1)), &
           &                                                              ' zeta= ',X1_base%f%x_IP(2,IP_GP(1)) 
      CALL abort(__STAMP__, &
           'EvalAux: Jacobian smaller that  1.0e-12 !!!', IntInfo=myRank )
    CASE(2) !quiet check, give back
      JacCheck=-1
    END SELECT
  ELSE
    JacCheck=1 !set to default for safety (abort if detJ<0)
  END IF
  CALL par_AllReduce(JacCheck,'MIN')
  IF(JacCheck.EQ.-1) THEN
    __PERFOFF('EvalAux')
    RETURN
  END IF

  __PERFON('EvalDOF_2')
  !2D data: interpolation points x gauss-points
  CALL X1_base%evalDOF((/0,DERIV_ZETA/),Uin%X1,dX1_dzeta)
  CALL X2_base%evalDOF((/0,DERIV_ZETA/),Uin%X2,dX2_dzeta)
  CALL LA_base%evalDOF((/0,DERIV_THET/),Uin%LA,dLA_dthet)
  CALL LA_base%evalDOF((/0,DERIV_ZETA/),Uin%LA,dLA_dzeta)
  __PERFOFF('EvalDOF_2')


  __PERFON('loop_2')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(iGP,i_mn,qloc,q_thet,q_zeta)  &
!$OMP   FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)  &
!$OMP   SHARED(b_thet,b_zeta,g_tt,g_tz,g_zz,sJ_p,sJ_h,sdetJ,sJ_bcov_thet,sJ_bcov_zeta,bbcov_sJ, &
!$OMP          J_p,J_h,hmap,dX1_dzeta,dX2_dzeta,dX1_dthet,dX2_dthet,chiPrime_GP,phiPrime_GP,dLA_dzeta,    &
!$OMP          dLA_dthet,X1_IP_GP,X2_IP_GP,zeta_IP)
  DO iGP=nGP_str,nGP_end
    DO i_mn=1,mn_IP
      b_thet(i_mn,iGP) = (chiPrime_GP(iGP)- phiPrime_GP(iGP)*dLA_dzeta(i_mn,iGP))    !b_theta
      b_zeta(i_mn,iGP) = phiPrime_GP(iGP)*(1.0_wp          + dLA_dthet(i_mn,iGP))    !b_zeta

      qloc(  1:3) = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:2) = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP)/) !dq(1:2)/dtheta
      q_thet(3)   = 0.0_wp                                      !dq(3)/dtheta
      q_zeta(1:2) = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP)/) !dq(1:2)/dzeta
      q_zeta(3)   = 1.0_wp                                      !dq(3)/zeta

      g_tt(i_mn,iGP)         = hmap%eval_gij(q_thet,qloc,q_thet)   !g_theta,theta
      g_tz(i_mn,iGP)         = hmap%eval_gij(q_thet,qloc,q_zeta)   !g_theta,zeta =g_zeta,theta
      g_zz(i_mn,iGP)         = hmap%eval_gij(q_zeta,qloc,q_zeta)   !g_zeta,zeta
      sJ_p(i_mn,iGP)         = 1.0_wp/J_p(i_mn,iGP)
      sJ_h(i_mn,iGP)         = 1.0_wp/J_h(i_mn,iGP)
      sdetJ(i_mn,iGP)        =  sJ_p(i_mn,iGP)*  sJ_h(i_mn,iGP)
      sJ_bcov_thet(i_mn,iGP) = (g_tt(i_mn,iGP)*b_thet(i_mn,iGP) + g_tz(i_mn,iGP)*b_zeta(i_mn,iGP))*sdetJ(i_mn,iGP)
      sJ_bcov_zeta(i_mn,iGP) = (g_tz(i_mn,iGP)*b_thet(i_mn,iGP) + g_zz(i_mn,iGP)*b_zeta(i_mn,iGP))*sdetJ(i_mn,iGP)
      bbcov_sJ(    i_mn,iGP) =  b_thet( i_mn,iGP)*sJ_bcov_thet(i_mn,iGP) &
                               +b_zeta( i_mn,iGP)*sJ_bcov_zeta(i_mn,iGP)
    END DO !i_mn
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('loop_2')

  __PERFOFF('EvalAux')

END SUBROUTINE EvalAux

!===================================================================================================================================
!> Evaluate total volume and average surface
!!
!===================================================================================================================================
SUBROUTINE EvalTotals(Uin,vol,surfAvg)
! MODULES
  USE MODgvec_globals      , ONLY: TWOPI
  USE MODgvec_MPI          , ONLY: par_Reduce
  USE MODgvec_sol_var_MHD3D, ONLY: t_sol_var_MHD3D
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN ) :: Uin  !! input solution 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)              , INTENT(OUT) :: vol      !! total integral of the volume
  REAL(wp)              , INTENT(OUT) :: surfAvg  !! average polodial surface 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER             :: iGP,i_mn,JacCheck
!===================================================================================================================================
  JacCheck=2
  CALL EvalAux(Uin,JacCheck)
  IF(JacCheck.EQ.-1) THEN
      CALL abort(__STAMP__, &
          ' detJ<0 in EvalAux, called from EvalTotals!!!' )
  END IF
  vol=0.0_wp
  surfAvg=0.0_wp
!$OMP PARALLEL DO       &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   REDUCTION(+:vol,surfAvg) PRIVATE(iGP,i_mn) FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)&
!$OMP   SHARED(J_h,J_p,w_GP)
  DO iGP=nGP_str,nGP_end
    DO i_mn=1,mn_IP
      vol    =vol    +ABS(J_h(i_mn,iGP)*J_p(i_mn,iGP))*w_GP(iGP)
      surfAvg=surfAvg+ABS(J_p(i_mn,iGP))*w_GP(iGP)
    END DO
  END DO 
!$OMP END PARALLEL DO
  CALL par_Reduce(vol,'SUM',0)
  CALL par_Reduce(surfAvg,'SUM',0)
  vol     = dthet_dzeta *vol
  surfAvg = dthet_dzeta *surfAvg /TWOPI

END SUBROUTINE EvalTotals

!===================================================================================================================================
!> Evaluate 3D MHD energy
!! NOTE: set callEvalaux >0 if not called before for the same Uin !!
!!
!===================================================================================================================================
FUNCTION EvalEnergy(Uin,callEvalAux,JacCheck) RESULT(W_MHD3D)
! MODULES
  USE MODgvec_MPI          , ONLY: par_AllReduce
  USE MODgvec_MHD3D_Vars   , ONLY: mu_0,sgammM1
  USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin         !! input solution 
  LOGICAL               , INTENT(IN   ) :: callEvalAux !! set True if evalAux was not called on Uin 
  INTEGER               , INTENT(INOUT) :: JacCheck !! if 1 on input: abort if detJ<0. 
                                                    !! if 2 on input, no abort, unchanged if detJ>0 ,return -1 if detJ<=0
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                              :: W_MHD3D     !! total integral of MHD3D energy 
                                                       !! W_MHD3D= int ( p/(gamma-1) + 1/(2mu_0) |B|^2) detJ ds dtheta dzeta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iGP,i_mn
  REAL(wp) :: Wmag_GP    !! magnetic energy at gauss points 
                         !! = 1/(dtheta*dzeta) * ( int [1/detJ * b_alpha*g_{alpha,beta}*b_beta]_iGP dtheta dzeta )
  REAL(wp) :: Vprime_GP  !! =  1/(dtheta*dzeta) *( int detJ|_iGP ,dtheta dzeta)
  REAL(wp) :: Wmag,Wpres
!===================================================================================================================================
  !WRITE(UNIT_stdOut,'(4X,A,I4)')'COMPUTE ENERGY... on rank:',myRank
  __PERFON('EvalEnergy')

  IF(callEvalAux) THEN
    CALL EvalAux(Uin,JacCheck)
    IF(JacCheck.EQ.-1) THEN
      W_MHD3D=1.0e30_wp
      WRITE(UNIT_stdOut,'(A)')'... detJ<0 in EvalAux '
      __PERFOFF('EvalEnergy')
      RETURN !accept detJ<0
    END IF
  ELSE
    IF(JacCheck.EQ.-1) THEN
        CALL abort(__STAMP__, &
            'You seem to have called EvalAux before, with a Jacobian smaller that  1.0e-12!!!' )
    END IF
  END IF

  Wmag = 0.0_wp
  Wpres= 0.0_wp
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC)  DEFAULT(NONE)  &
!$OMP   PRIVATE(iGP,i_mn,Wmag_GP,Vprime_GP)   &
!$OMP   FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)  &
!$OMP   REDUCTION(+:Wmag,Wpres)          &
!$OMP   SHARED(bbcov_sJ,detJ,pres_GP,w_GP)
  DO iGP=nGP_str,nGP_end
    Wmag_GP=0.0_wp
!$OMP SIMD REDUCTION(+:Wmag_GP)
    DO i_mn=1,mn_IP
      Wmag_GP   = Wmag_GP+ bbcov_sJ(i_mn,iGP)
    END DO
    Vprime_GP = 0.0_wp
!$OMP SIMD REDUCTION(+:Vprime_GP)
    DO i_mn=1,mn_IP
      Vprime_GP = Vprime_GP+ detJ(i_mn,iGP)
    END DO
    Wmag  = Wmag  + Wmag_GP*w_GP(iGP)
    Wpres = Wpres + Vprime_GP*pres_GP(iGP)*w_GP(iGP)
  END DO !iGP
!$OMP END PARALLEL DO 

  W_MHD3D= dthet_dzeta* (  0.5_wp      *Wmag + mu_0*sgammM1*Wpres)
  __PERFON('reduce_W_MHD3D')
  CALL par_AllReduce(W_MHD3D,'SUM')
  __PERFOFF('reduce_W_MHD3D')

   __PERFOFF('EvalEnergy')

!  SWRITE(UNIT_stdOut,'(A,E21.11)')'... DONE: ',W_MHD3D
END FUNCTION EvalEnergy

!===================================================================================================================================
!> Evaluate the variation of the Energy = Force... F_j=-(D_U W(U))_j = -DW(u_h)*testfunc_j
!! NOTE: set callEvalaux TRUE if not called before for the same Uin !!
!!
!===================================================================================================================================
SUBROUTINE EvalForce(Uin,callEvalAux,JacCheck,F_MHD3D,noBC)
! MODULES
  USE MODgvec_Globals,       ONLY : nRanks,myRank
  USE MODgvec_MPI,           ONLY : par_IReduce,par_IBcast,par_Wait,req1,req2,req3,par_Barrier,par_BCast
  USE MODgvec_MHD3D_Vars,    ONLY : X1_base,X2_base,LA_base,hmap,mu_0,PrecondType
  USE MODgvec_MHD3D_Vars,    ONLY : X1_BC_type,X2_BC_type,LA_BC_type
  USE MODgvec_sol_var_MHD3D, ONLY : t_sol_var_MHD3D
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin         !! input solution
  LOGICAL               , INTENT(IN   ) :: callEvalAux !! set True if evalAux was not called on Uin
  INTEGER               , INTENT(INOUT) :: JacCheck !! if 1 on input: abort if detJ<0.
                                                    !! if 2 on input, no abort, unchanged if detJ>0 ,return -1 if detJ<=0
  LOGICAL, OPTIONAL     , INTENT(IN)    :: noBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: F_MHD3D     !! variation of the energy projected onto the basis functions of Uin
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: ibase,nBase,iMode,modes,iGP,i_mn,Deg,iElem,modes_str,modes_end,iRank,offset_modes(0:nRanks)
  REAL(wp)  :: qloc(3),q_thet(3),q_zeta(3),w_GP_IP,p_mu_0
  REAL(wp)  :: hmap_g_t1,hmap_g_z1,hmap_Jh_dq1,hmap_g_tt_dq1,hmap_g_tz_dq1,hmap_g_zz_dq1
  REAL(wp)  :: hmap_g_t2,hmap_g_z2,hmap_Jh_dq2,hmap_g_tt_dq2,hmap_g_tz_dq2,hmap_g_zz_dq2
  REAL(wp)  ::    F_X1_GP_IP(nGP_str:nGP_end,1:X1_base%f%modes)
  REAL(wp)  ::  F_X1ds_GP_IP(nGP_str:nGP_end,1:X1_base%f%modes)
  REAL(wp)  ::    F_X2_GP_IP(nGP_str:nGP_end,1:X2_base%f%modes)
  REAL(wp)  ::  F_X2ds_GP_IP(nGP_str:nGP_end,1:X2_base%f%modes)
  REAL(wp)  ::    F_LA_GP_IP(nGP_str:nGP_end,1:LA_base%f%modes)
  REAL(wp)  ::    dW(1:mn_IP,nGP_str:nGP_end)        != p+1/2*B^2=p(s)+|Phi'(s)|^2 (b^alpha *g_{alpha,beta} *b^beta)/(2 *detJ^2)
  REAL(wp),DIMENSION(1:mn_IP,nGP_str:nGP_end)  :: btt_sJ,btz_sJ,bzz_sJ,  & != b^theta*b^theta/detJ, b^theta*b^zeta/detJ,b^zeta*b^zeta/detJ 
                                                  coefY,coefY_thet,coefY_zeta,coefY_s
!===================================================================================================================================
!  SWRITE(UNIT_stdOut,'(4X,A)',ADVANCE='NO')'COMPUTE FORCE...'
#if MPIDEBUG==1                                               
  WRITE(UNIT_stdOut,'(4X,A,I4)')'COMPUTE FORCE...',myRank
  CALL par_Barrier(beforeScreenOut="DEBUG ENTER FORCE")
#endif
  __PERFON('EvalForce')
  IF(callEvalAux) THEN
    CALL EvalAux(Uin,JacCheck)
  END IF
  IF(JacCheck.EQ.-1) THEN
    CALL abort(__STAMP__, &
         'negative Jacobian was found when you call EvalAux before!!!')
  END IF

  __PERFON('buildPrecond')
  IF(PrecondType.GT.0) CALL BuildPrecond()
  __PERFOFF('buildPrecond')

  !additional auxiliary variables for X1 and X2 force
  __PERFON('loop_prepare')
!$OMP PARALLEL DO    &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iGP,i_mn,p_mu_0) FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)
  DO iGP=nGP_str,nGP_end
    p_mu_0=mu_0*pres_GP(iGP)
    DO i_mn=1,mn_IP
      dW(    i_mn,iGP)=  0.5_wp*bbcov_sJ(i_mn,iGP)                 *sdetJ(i_mn,iGP) + p_mu_0 !=1/(2)*B^2+p
      btt_sJ(i_mn,iGP)=  0.5_wp*b_thet(  i_mn,iGP)*b_thet(i_mn,iGP)*sdetJ(i_mn,iGP)
      btz_sJ(i_mn,iGP)=  0.5_wp*b_thet(  i_mn,iGP)*b_zeta(i_mn,iGP)*sdetJ(i_mn,iGP)
      bzz_sJ(i_mn,iGP)=  0.5_wp*b_zeta(  i_mn,iGP)*b_zeta(i_mn,iGP)*sdetJ(i_mn,iGP)
    END DO !i_mn
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('loop_prepare')

  nBase = X1_Base%s%nBase
  modes = X1_Base%f%modes
  modes_str = X1_Base%f%modes_str
  modes_end = X1_Base%f%modes_end
  offset_modes = X1_Base%f%offset_modes
  deg   = X1_base%s%deg

  __PERFON('EvalForce_modes1')
  __PERFON('loop_prep_coefs')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(iGP,i_mn,qloc,q_thet,q_zeta,                                                &
!$OMP           hmap_Jh_dq1,hmap_g_t1,hmap_g_tt_dq1,hmap_g_z1,hmap_g_zz_dq1,hmap_g_tz_dq1)  &
!$OMP   FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)                                                 &
!$OMP   SHARED(X1_IP_GP,X2_IP_GP,zeta_IP,dX1_dthet,dX2_dthet,dX1_dzeta,dX2_dzeta, &
!$OMP          hmap,dW,J_h,J_p,dX2_ds,btt_sJ,bzz_sJ,btz_sJ,                                 &
!$OMP          coefY,coefY_thet,coefY_zeta,coefY_s)
  DO iGP=nGP_str,nGP_end
    DO i_mn=1,mn_IP
      qloc(1:3)      = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:3)    = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP),0.0_wp/)
      q_zeta(1:3)    = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP),1.0_wp/)
      !Y1tilde=(1,0,0)
      hmap_Jh_dq1   = hmap%eval_Jh_dq1(        qloc        ) !~Y1
      hmap_g_t1     = hmap%eval_gij(    q_thet,qloc,(/1.0_wp,0.0_wp,0.0_wp/)) !~Y1_thet
      hmap_g_z1     = hmap%eval_gij(    q_zeta,qloc,(/1.0_wp,0.0_wp,0.0_wp/)) !~Y1_zeta
      hmap_g_tt_dq1 = hmap%eval_gij_dq1(q_thet,qloc, q_thet) !~Y1
      hmap_g_tz_dq1 = hmap%eval_gij_dq1(q_thet,qloc, q_zeta) !~Y1
      hmap_g_zz_dq1 = hmap%eval_gij_dq1(q_zeta,qloc, q_zeta) !~Y1

! ADDED TO F_X1_GP_IP(iGP,iMode):
!                         -dW(    i_mn,iGP)*J_h(i_mn,iGP)*dX2_ds(     i_mn,iGP)*Y1_thet   & ![deltaJ]_Y1
!                         +dW(    i_mn,iGP)*J_p(i_mn,iGP)*hmap_Jh_dq1(i_mn,iGP)*Y1        & ![deltaJ]_Y1
!                         -btt_sJ(i_mn,iGP)*2.0_wp*hmap_g_t1(         i_mn,iGP)*Y1_thet   & ![delta g_tt]_Y1
!                         -btt_sJ(i_mn,iGP)*       hmap_g_tt_dq1(     i_mn,iGP)*Y1        & ![delta g_tt]_Y1
!                         -bzz_sJ(i_mn,iGP)*2.0_wp*hmap_g_z1(         i_mn,iGP)*Y1_zeta   & ![delta g_zz]_Y1
!                         -bzz_sJ(i_mn,iGP)*       hmap_g_zz_dq1(     i_mn,iGP)*Y1        & ![delta g_zz]_Y1
!                         -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_t1(         i_mn,iGP)*Y1_zeta   & !2*[delta g_tz]_y1
!                         -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_z1(         i_mn,iGP)*Y1_thet   & !2*[delta g_tz]_y1
!                         -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_tz_dq1(     i_mn,iGP)*Y1        & !2*[delta g_tz]_y1

      coefY     (i_mn,iGP)=( dW(    i_mn,iGP)*J_p(i_mn,iGP)*hmap_Jh_dq1    & ![deltaJ]_Y1
                            -btt_sJ(i_mn,iGP)*       hmap_g_tt_dq1         & ![delta g_tt]_Y1
                            -bzz_sJ(i_mn,iGP)*       hmap_g_zz_dq1         & ![delta g_zz]_Y1
                            -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_tz_dq1      )  !2*[delta g_tz]_y1

      coefY_thet(i_mn,iGP)=(-dW(i_mn,iGP)*J_h(i_mn,iGP)*dX2_ds(i_mn,iGP)   & ![deltaJ]_Y1
                            +2.0_wp*(-btt_sJ(i_mn,iGP)*hmap_g_t1           & ![delta g_tt]_Y1
                                     -btz_sJ(i_mn,iGP)*hmap_g_z1  )      )   !2*[delta g_tz]_y1

      coefY_zeta(i_mn,iGP)=( 2.0_wp*(-bzz_sJ(i_mn,iGP)*hmap_g_z1           & ![delta g_zz]_Y1
                                     -btz_sJ(i_mn,iGP)*hmap_g_t1  )      )   !2*[delta g_tz]_y1

      coefY_s   (i_mn,iGP)=(dW(i_mn,iGP)*J_h(i_mn,iGP)*dX2_dthet( i_mn,iGP))
    END DO !i_mn
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('loop_prep_coefs')

  __PERFON('fbase')
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iGP,w_GP_IP) FIRSTPRIVATE(nGP_str,nGP_end)
  DO iGP=nGP_str,nGP_end
    w_GP_IP=w_GP(iGP)*dthet_dzeta
    CALL X1_base%f%projectIPtoDOF(.FALSE.,w_GP_IP,         0,coefY(     :,iGP),F_X1_GP_IP(  iGP,:))
    CALL X1_base%f%projectIPtoDOF(.TRUE. ,w_GP_IP,DERIV_THET,coefY_thet(:,iGP),F_X1_GP_IP(  iGP,:))!d/dthet
    CALL X1_base%f%projectIPtoDOF(.TRUE. ,w_GP_IP,DERIV_ZETA,coefY_zeta(:,iGP),F_X1_GP_IP(  iGP,:))!d/dzeta
    CALL X1_base%f%projectIPtoDOF(.FALSE.,w_GP_IP,         0,coefY_s(   :,iGP),F_X1ds_GP_IP(iGP,:))
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('fbase')

  __PERFON('sbase')
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,iElem,iGP,iBase) FIRSTPRIVATE(modes,nElems_str,nElems_end)
  DO iMode=1,modes
    F_MHD3D%X1(:,iMode)=0.0_wp
    DO iElem=nElems_str,nElems_end 
      iGP=(iElem-1)*(degGP+1)+1
      iBase=X1_base%s%base_offset(iElem)
      F_MHD3D%X1(iBase:iBase+deg,iMode) = F_MHD3D%X1(iBase:iBase+deg,iMode) &
                                    + MATMUL(F_X1_GP_IP(  iGP:iGP+degGP,iMode),X1_base%s%base_GP(   0:degGP,0:deg,iElem)) &
                                    + MATMUL(F_X1ds_GP_IP(iGP:iGP+degGP,iMode),X1_base%s%base_ds_GP(0:degGP,0:deg,iElem))
    END DO !iElem
  END DO !iMode
!$OMP END PARALLEL DO
  __PERFOFF('sbase')

#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG BEFORE FIRST REDUCE")
#endif
  !. add up all the pieces of X1 calculated by the different MPI tasks
  __PERFON('reduce_solution_X1')
  !!!CALL par_AllReduce(F_MHD3D%X1,'SUM') !<< possible alternative
  DO iRank=0,nRanks-1 !<<<<
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      !!!CALL par_Reduce(F_MHD3D%X1(:,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank) !<<<< possible alternative
      CALL par_IReduce(F_MHD3D%X1(:,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank,req1(iRank)) !<<<< I-reduce different mode ranges to different ranks
  END DO
  __PERFOFF('reduce_solution_X1')

  __PERFOFF('EvalForce_modes1')

  __PERFON('EvalForce_modes2')
  __PERFON('loop_prep_coefs')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(iGP,i_mn,qloc,q_thet,q_zeta, &
!$OMP           hmap_Jh_dq2,hmap_g_t2,hmap_g_tt_dq2,hmap_g_z2,hmap_g_zz_dq2,hmap_g_tz_dq2) &
!$OMP   FIRSTPRIVATE(nGP_str,nGP_end,mn_IP)  &
!$OMP   SHARED(X1_IP_GP,X2_IP_GP,zeta_IP,dX1_dthet,dX1_dzeta,dX2_dthet,dX2_dzeta,&
!$OMP          hmap,dW,J_h,J_p,dX2_ds,btt_sJ,bzz_sJ,btz_sJ,dX1_ds,  &
!$OMP          coefY,coefY_thet,coefY_zeta,coefY_s)
  DO iGP=nGP_str,nGP_end
    DO i_mn=1,mn_IP

      qloc(1:3)      = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:3)    = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP),0.0_wp/)
      q_zeta(1:3)    = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP),1.0_wp/)
      !Y2tilde=(0,1,0)
      hmap_Jh_dq2    = hmap%eval_Jh_dq2(        qloc        ) !~Y2
      hmap_g_t2      = hmap%eval_gij(    q_thet,qloc,(/0.0_wp,1.0_wp,0.0_wp/)) !~Y2_thet
      hmap_g_z2      = hmap%eval_gij(    q_zeta,qloc,(/0.0_wp,1.0_wp,0.0_wp/)) !~Y2_zeta
      hmap_g_tt_dq2  = hmap%eval_gij_dq2(q_thet,qloc, q_thet) !~Y2
      hmap_g_tz_dq2  = hmap%eval_gij_dq2(q_thet,qloc, q_zeta) !~Y2
      hmap_g_zz_dq2  = hmap%eval_gij_dq2(q_zeta,qloc, q_zeta) !~Y2

! ADDED TO F_X2_GP_IP(iGP,iMode):
!                         +dW(    i_mn,iGP)*J_h(i_mn,iGP)*dX1_ds(i_mn,iGP)*Y2_thet  & ! [deltaJ]_Y2
!                         +dW(    i_mn,iGP)*J_p(i_mn,iGP)*hmap_Jh_dq2     *Y2       & ! [deltaJ]_Y2
!                         -btt_sJ(i_mn,iGP)*2.0_wp*hmap_g_t2              *Y2_thet  & ! [delta g_tt]_Y2
!                         -btt_sJ(i_mn,iGP)*       hmap_g_tt_dq2          *Y2       & ! [delta g_tt]_Y2
!                         -bzz_sJ(i_mn,iGP)*2.0_wp*hmap_g_z2              *Y2_zeta  & ! [delta g_zz]_Y2
!                         -bzz_sJ(i_mn,iGP)*       hmap_g_zz_dq2          *Y2       & ! [delta g_zz]_Y2
!                         -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_t2              *Y2_zeta  & ! 2*[delta g_tz]_Y1
!                         -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_z2              *Y2_thet  & ! 2*[delta g_tz]_Y1
!                         -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_tz_dq2          *Y2       & ! 2*[delta g_tz]_Y1

      coefY     (i_mn,iGP)=( dW(i_mn,iGP)*J_p(i_mn,iGP)*hmap_Jh_dq2       & ! [deltaJ]_Y2
                            -btt_sJ(i_mn,iGP)*       hmap_g_tt_dq2        & ! [delta g_tt]_Y2
                            -bzz_sJ(i_mn,iGP)*       hmap_g_zz_dq2        & ! [delta g_zz]_Y2
                            -btz_sJ(i_mn,iGP)*2.0_wp*hmap_g_tz_dq2      )   ! 2*[delta g_tz]_Y1

      coefY_thet(i_mn,iGP)=( dW(i_mn,iGP)*J_h(i_mn,iGP)*dX1_ds(i_mn,iGP)  & ! [deltaJ]_Y2
                            +2.0_wp*(-btt_sJ(i_mn,iGP)*hmap_g_t2          & ! [delta g_tt]_Y2
                                     -btz_sJ(i_mn,iGP)*hmap_g_z2 )      )   ! 2*[delta g_tz]_Y1

      coefY_zeta(i_mn,iGP)=( 2.0_wp*(-bzz_sJ(i_mn,iGP)*hmap_g_z2          & ! [delta g_zz]_Y2
                                     -btz_sJ(i_mn,iGP)*hmap_g_t2 )      )   ! 2*[delta g_tz]_Y1

      coefY_s   (i_mn,iGP)=(-dW(i_mn,iGP)*J_h(i_mn,iGP)*dX1_dthet( i_mn,iGP))
    END DO !i_mn
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('loop_prep_coefs')

  nBase = X2_base%s%nBase
  modes = X2_base%f%modes
  modes_str = X2_base%f%modes_str
  modes_end = X2_base%f%modes_end
  offset_modes = X2_Base%f%offset_modes
  deg   = X2_base%s%deg

  __PERFON('fbase')
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iGP,w_GP_IP) FIRSTPRIVATE(nGP_str,nGP_end) 
  DO iGP=nGP_str,nGP_end
    w_GP_IP = w_GP(iGP)*dthet_dzeta
    CALL X2_base%f%projectIPtoDOF(.FALSE.,w_GP_IP,         0,coefY(     :,iGP),F_X2_GP_IP(  iGP,:))
    CALL X2_base%f%projectIPtoDOF(.TRUE. ,w_GP_IP,DERIV_THET,coefY_thet(:,iGP),F_X2_GP_IP(  iGP,:))!d/dthet
    CALL X2_base%f%projectIPtoDOF(.TRUE. ,w_GP_IP,DERIV_ZETA,coefY_zeta(:,iGP),F_X2_GP_IP(  iGP,:))!d/dzeta
    CALL X2_base%f%projectIPtoDOF(.FALSE.,w_GP_IP,         0,coefY_s(   :,iGP),F_X2ds_GP_IP(iGP,:))
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('fbase')
  __PERFON('sbase')
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,iElem,iGP,iBase) FIRSTPRIVATE(modes,nElems_str,nElems_end)
  DO iMode=1,modes
    F_MHD3D%X2(:,iMode)=0.0_wp
    DO iElem=nElems_str,nElems_end
      iGP=(iElem-1)*(degGP+1)+1
      ibase=X2_base%s%base_offset(iElem)
      F_MHD3D%X2(iBase:iBase+deg,iMode) = F_MHD3D%X2(iBase:iBase+deg,iMode) &
                                    + MATMUL(F_X2_GP_IP(  iGP:iGP+degGP,iMode),X2_base%s%base_GP(   0:degGP,0:deg,iElem)) & 
                                    + MATMUL(F_X2ds_GP_IP(iGP:iGP+degGP,iMode),X2_base%s%base_ds_GP(0:degGP,0:deg,iElem)) 
    END DO !iElem
  END DO !iMode
!$OMP END PARALLEL DO 
  __PERFOFF('sbase')

#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG BEFORE X2 REDUCE")
#endif
  __PERFON('reduce_solution_X2')
  !. add up all the pieces of X2 calculated by the different MPI tasks
  !!!CALL par_AllReduce(F_MHD3D%X2,'SUM') !<< possible alternative
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      !!!CALL par_Reduce(F_MHD3D%X2(:,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank) !<< possible alternative
      CALL par_IReduce(F_MHD3D%X2(:,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank,req2(iRank)) !<<<< I-reduce different mode ranges to different ranks
  END DO
  __PERFOFF('reduce_solution_X2')

  __PERFOFF('EvalForce_modes2')


  __PERFON('EvalForce_modes3')

  nBase = LA_base%s%nBase
  modes = LA_base%f%modes
  modes_str = LA_base%f%modes_str
  modes_end = LA_base%f%modes_end
  offset_modes = LA_Base%f%offset_modes
  deg   = LA_base%s%deg

  __PERFON('fbase')
!   coefY_zeta(i_mn,iGP)= w_GP_IP*sJ_bcov_thet(i_mn,iGP)
!   coefY_thet(i_mn,iGP)=-w_GP_IP*sJ_bcov_zeta(i_mn,iGP)
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iGP,w_GP_IP) FIRSTPRIVATE(nGP_str,nGP_end)
  DO iGP=nGP_str,nGP_end
    w_GP_IP=PhiPrime_GP(iGP)*w_GP(iGP)*dthet_dzeta
    CALL LA_base%f%projectIPtoDOF(.FALSE., w_GP_IP,DERIV_ZETA, sJ_bcov_thet(:,iGP),F_LA_GP_IP(iGP,:)) !d/dzeta
    CALL LA_base%f%projectIPtoDOF(.TRUE. ,-w_GP_IP,DERIV_THET, sJ_bcov_zeta(:,iGP),F_LA_GP_IP(iGP,:)) !d/dthet
  END DO !iGP
!$OMP END PARALLEL DO
  __PERFOFF('fbase')
  __PERFON('sbase')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,iElem,iGP,iBase) FIRSTPRIVATE(modes,nElems_str,nElems_end)
  DO iMode=1,modes
    F_MHD3D%LA(:,iMode)=0.0_wp
    DO iElem=nElems_str,nElems_end
      iGP=(iElem-1)*(degGP+1)+1
      ibase=LA_base%s%base_offset(iElem)
      F_MHD3D%LA(iBase:iBase+deg,iMode) = F_MHD3D%LA(iBase:iBase+deg,iMode) &
                                    + MATMUL(F_LA_GP_IP(iGP:iGP+degGP,iMode),LA_base%s%base_GP(0:degGP,0:deg,iElem))
    END DO !iElem
  END DO !iMode
!$OMP END PARALLEL DO
  __PERFOFF('sbase')

#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG BEFORE LA REDUCE")
#endif
  __PERFON('reduce_solution_LA')
  !. Add up all the pieces of LA calculated by the different MPI tasks
  !!!CALL par_AllReduce(F_MHD3D%LA,'SUM') !<< possible alternative
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      !!!CALL par_Reduce(F_MHD3D%LA(:,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank) !<< possible alternative
      CALL par_IReduce(F_MHD3D%LA(:,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank,req3(iRank)) !<<<< I-reduce different mode ranges to different ranks
  END DO
  __PERFOFF('reduce_solution_LA')

  __PERFOFF('EvalForce_modes3')


  __PERFON('EvalForce_modes1')
  modes_str = X1_base%f%modes_str
  modes_end = X1_base%f%modes_end
  offset_modes = X1_Base%f%offset_modes

  __PERFON('reduce_solution_X1')
  CALL par_Wait(req1(0:nRanks-1))
#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG AFTER FINISH REDUCE")
#endif
  __PERFOFF('reduce_solution_X1')

  __PERFON('apply_precond')
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(F_MHD3D,X1_BC_type,X1_base) PRIVATE(iMode) FIRSTPRIVATE(modes_str,modes_end)
  DO iMode=modes_str,modes_end
    CALL X1_base%s%applyBCtoRHS(F_MHD3D%X1(:,iMode),X1_BC_type(:,iMode))
  END DO !iMode
!$OMP END PARALLEL DO

  IF(PrecondType.GT.0)THEN
    SELECT TYPE(precond_X1); TYPE IS(sll_t_spline_matrix_banded)
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode)  FIRSTPRIVATE(modes_str,modes_end,nBase) &
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(F_MHD3D,precond_X1)
#else
!$OMP   DEFAULT(SHARED)
#endif
      DO iMode=modes_str,modes_end !<<<<
        CALL ApplyPrecond(nBase,precond_X1(iMode),F_MHD3D%X1(:,iMode))
      END DO !iMode
!$OMP END PARALLEL DO
    END SELECT !TYPE(precond_X1)
  END IF !PrecondType.GT.0
  __PERFOFF('apply_precond')
  IF(PrecondType.LE.0)THEN
    !apply strong BC to X1 if no Precond
    CALL ApplyBC_Fstrong(1,F_MHD3D)  
  END IF 

  __PERFON('Bcast_solution_X1')
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      CALL par_Bcast(F_MHD3D%X1(:,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank) !<<<< broadcast different mode ranges to different ranks
      !CALL par_IBcast(F_MHD3D%X1(:,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank,req1(iRank)) !<<<< broadcast different mode ranges to different ranks
  END DO
  __PERFOFF('Bcast_solution_X1')

  __PERFOFF('EvalForce_modes1')

  __PERFON('EvalForce_modes2')

  modes_str = X2_base%f%modes_str
  modes_end = X2_base%f%modes_end
  offset_modes = X2_Base%f%offset_modes

  __PERFON('reduce_solution_X2')
  CALL par_Wait(req2(0:nRanks-1))
#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG AFTER FINISH REDUCE X2")
#endif
  __PERFOFF('reduce_solution_X2')

  __PERFON('apply_precond')
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(X2_base,F_MHD3D,X2_BC_type) PRIVATE(iMode)  FIRSTPRIVATE(modes_str,modes_end)
  DO iMode=modes_str,modes_end
    CALL X2_base%s%applyBCtoRHS(F_MHD3D%X2(:,iMode),X2_BC_type(:,iMode))
  END DO !iMode
!$OMP END PARALLEL DO

  IF(PrecondType.GT.0)THEN
    SELECT TYPE(precond_X2); TYPE IS(sll_t_spline_matrix_banded)
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode)  FIRSTPRIVATE(modes_str,modes_end,nBase) &
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(precond_X2,F_MHD3D)
#else
!$OMP   DEFAULT(SHARED)
#endif
      DO iMode=modes_str,modes_end
        CALL ApplyPrecond(nBase,precond_X2(iMode),F_MHD3D%X2(:,iMode))
      END DO !iMode
!$OMP END PARALLEL DO
    END SELECT !TYPE(precond_X2)
  END IF !PrecondType.GT.0
  __PERFOFF('apply_precond')
  IF(PrecondType.LE.0)THEN
    !apply strong BC to X2 if no Precond
    CALL ApplyBC_Fstrong(2,F_MHD3D)  
  END IF 

  __PERFON('Bcast_solution_X2')
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      CALL par_Bcast(F_MHD3D%X2(:,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank) !<<<< reduce different mode ranges to different ranks
      !CALL par_IBcast(F_MHD3D%X2(:,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank,req2(iRank)) !<<<< reduce different mode ranges to different ranks
  END DO
  __PERFOFF('Bcast_solution_X2')

  __PERFOFF('EvalForce_modes2')


  __PERFON('EvalForce_modes3')

  modes_str = LA_base%f%modes_str
  modes_end = LA_base%f%modes_end
  offset_modes = LA_Base%f%offset_modes

  __PERFON('reduce_solution_LA')
  CALL par_Wait(req3(0:nRanks-1))
#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG AFTER FINISH REDUCE LA")
#endif
  __PERFOFF('reduce_solution_LA')

  __PERFON('apply_precond')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(LA_base,F_MHD3D,LA_BC_type) PRIVATE(iMode)  FIRSTPRIVATE(modes_str,modes_end)
  DO iMode=modes_str,modes_end
    CALL LA_base%s%applyBCtoRHS(F_MHD3D%LA(:,iMode),LA_BC_type(:,iMode))
  END DO !iMode
!$OMP END PARALLEL DO

  IF(PrecondType.GT.0)THEN
    SELECT TYPE(precond_LA); TYPE IS(sll_t_spline_matrix_banded)
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode) FIRSTPRIVATE(modes_str,modes_end,nBase)&
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(precond_LA,F_MHD3D)
#else
!$OMP   DEFAULT(SHARED)
#endif
      DO iMode=modes_str,modes_end
        CALL ApplyPrecond(nBase,precond_LA(iMode),F_MHD3D%LA(:,iMode))
      END DO !iMode
!$OMP END PARALLEL DO
    END SELECT !TYPE(precond_LA)
  ELSE
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode) FIRSTPRIVATE(modes_str,modes_end)&
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(LA_base,F_MHD3D)
#else
!$OMP   DEFAULT(SHARED)
#endif
    DO iMode=modes_str,modes_end
      CALL LA_base%s%mass%solve_inplace(1,F_MHD3D%LA(:,iMode))
    END DO !iMode
!$OMP END PARALLEL DO
  END IF !PrecondType.GT.0
  __PERFOFF('apply_precond')
  IF(PrecondType.LE.0)THEN
    !apply strong BC to LA if no Precond
    CALL ApplyBC_Fstrong(3,F_MHD3D)  
  END IF 

#if MPIDEBUG==1
  WRITE(*,*)'DEBUG',myRank, offset_modes(myRank),offset_modes(myRank+1)
#endif
  __PERFON('Bcast_solution_LA')
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
!      CALL par_Bcast(F_MHD3D%LA(:,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank) !<< possible alternative
      CALL par_IBcast(F_MHD3D%LA(:,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank,req3(iRank)) !<<<< reduce different mode ranges to different ranks
  END DO

!  CALL par_Wait(req1(0:nRanks-1))
!  CALL par_Wait(req2(0:nRanks-1))
#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG BEFORE FINISH LA BCAST")
#endif
  CALL par_Wait(req3(0:nRanks-1))
#if MPIDEBUG==1
  CALL par_Barrier(beforeScreenOut="DEBUG AFTER FINISH BCASTS")
#endif
  __PERFOFF('Bcast_solution_LA')

  __PERFOFF('EvalForce_modes3')

  IF(PRESENT(noBC))THEN
    IF(noBC)THEN
      __PERFOFF('EvalForce')
      RETURN !DEBUG
    END IF
  END IF

  __PERFOFF('EvalForce')

!  SWRITE(UNIT_stdOut,'(A,3E21.11)')'... DONE: Norm of force |X1|,|X2|,|LA|: ',SQRT(F_MHD3D%norm_2())
END SUBROUTINE EvalForce


!===================================================================================================================================
!> Applies strong boundary condition to force DOF
!!
!===================================================================================================================================
SUBROUTINE ApplyBC_Fstrong(whichVar,F_MHD3D) 
! MODULES
  USE MODgvec_MHD3D_Vars,ONLY:X1_base,X2_base,LA_base
  USE MODgvec_MHD3D_Vars,ONLY:X1_BC_Type,X2_BC_Type,LA_BC_type
  USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: whichVar !! =1: X1, =2: X2,  =3: LA
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: F_MHD3D     !! variation of the energy projected onto the basis functions of Uin 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iMode
  REAL(wp)  :: BC_val(2)
!===================================================================================================================================
  !apply strong boundary conditions
  
  BC_val =(/      0.0_wp,      0.0_wp/)

  SELECT CASE(whichVar)  
  CASE(1)  !X1 BC
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(X1_base,F_MHD3D,X1_BC_type,BC_val) PRIVATE(iMode)
    DO imode=X1_base%f%modes_str,X1_base%f%modes_end
      CALL X1_base%s%applyBCtoDOF(F_MHD3D%X1(:,iMode),X1_BC_type(:,iMode),BC_val)
    END DO 
!$OMP END PARALLEL DO 
  
  CASE(2)  !X2 BC
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(X2_base,F_MHD3D,X2_BC_type,BC_val) PRIVATE(iMode)
    DO imode=X2_base%f%modes_str,X2_base%f%modes_end
      CALL X2_base%s%applyBCtoDOF(F_MHD3D%X2(:,iMode),X2_BC_type(:,iMode),BC_val)
    END DO 
!$OMP END PARALLEL DO 
  
  CASE(3)  !LA BC
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(LA_base,F_MHD3D,LA_BC_type,BC_val) PRIVATE(iMode)
    DO imode=LA_base%f%modes_str,LA_base%f%modes_end
      CALL LA_base%s%applyBCtoDOF(F_MHD3D%LA(:,iMode),LA_BC_type(:,iMode),BC_val)
    END DO 
!$OMP END PARALLEL DO 
  END SELECT !whichVar

END SUBROUTINE ApplyBC_Fstrong


!===================================================================================================================================
!> Build preconditioner matrices for X1,X2,LA and factorize, for all modes
!! the matrix is only radially dependent, and has the form
!! K_ij = int(s,0,1) d/ds sbase_i(s) <D_ss>(s) d/ds sbase_j(s)
!!                   + sbase_i(s) (<S>(s) + |Phi'(s)|^2 (-m^2 <D_tt>(s) - n^2 <D_zz>(s) ) ) sbase_j(s)
!! where < > denote an average over the angular coordinates
!!
!===================================================================================================================================
SUBROUTINE BuildPrecond() 
! MODULES
  USE MODgvec_MPI,        ONLY : par_AllReduce
  USE MODgvec_MHD3D_Vars, ONLY : X1_base,X2_base,LA_base,hmap
  USE MODgvec_MHD3D_Vars, ONLY : X1_BC_Type,X2_BC_Type,LA_BC_type
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                     :: ibase,nBase,iMode,modes_str,modes_end,iGP,i_mn,Deg,iElem,i,j
  INTEGER                     :: nD,tBC
  REAL(wp)                    :: qloc(3),q_thet(3),q_zeta(3),smn_IP,smn_IP_w_GP,norm_mn
  REAL(wp),DIMENSION(1:mn_IP) :: G11, G21, G31, G22, G32, dJh_dq1, dJh_dq2, bt_sJ, bz_sJ, &
       &                         b_dX1_tz,b_dX2_tz,gtt_dq1,gtz_dq1,gzz_dq1,gtt_dq2,gtz_dq2,gzz_dq2
!  REAL(wp),DIMENSION(nGP_str:nGP_end)   :: DX1_tt, DX1_tz, DX1_zz, DX1, DX1_ss
!  REAL(wp),DIMENSION(nGP_str:nGP_end)   :: DX2_tt, DX2_tz, DX2_zz, DX2, DX2_ss
!  REAL(wp),DIMENSION(nGP_str:nGP_end)   :: DLA_tt, DLA_tz, DLA_zz

  REAL(wp),ALLOCATABLE        :: D_mn(:),P_BCaxis(:,:), P_BCedge(:,:) !only needed on MPIroot
!===================================================================================================================================
!  WRITE(*,*)'BUILD PRECONDITIONER MATRICES'
  __PERFON('loop_1')
  smn_IP=1.0_wp/REAL(mn_IP,wp)
  
  !WHEN COMM CHANGED TO GATHER, ZEROING NOT NEEDED ANYMORE
  D_buf=0.0_wp
  
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED)   &
!$OMP   PRIVATE(iGP,i_mn,qloc,q_thet,q_zeta,G11,G21,G31,G22,G32,dJh_dq1,dJh_dq2,bt_sJ,bz_sJ,&
!$OMP           b_dX1_tz,b_dX2_tz,gtt_dq1,gtz_dq1,gzz_dq1,gtt_dq2,gtz_dq2,gzz_dq2,smn_IP_w_GP) &
!$OMP   FIRSTPRIVATE(smn_IP,mn_IP,nGP_str,nGP_end)
  loop_nGP: DO iGP=nGP_str,nGP_end  !<<<<
    !dont forget to average
    !additional variables
    smn_IP_w_GP=smn_IP*w_GP(iGP) !include gauss weight here!
    loop_mn_IP: DO i_mn=1,mn_IP
      qloc(1:3)     = (/ X1_IP_GP(i_mn,iGP), X2_IP_GP(i_mn,iGP),zeta_IP(i_mn)/)
      q_thet(1:3)   = (/dX1_dthet(i_mn,iGP),dX2_dthet(i_mn,iGP),0.0_wp/)
      q_zeta(1:3)   = (/dX1_dzeta(i_mn,iGP),dX2_dzeta(i_mn,iGP),1.0_wp/)
      gtt_dq1(i_mn) = hmap%eval_gij_dq1(q_thet,qloc, q_thet)
      gtz_dq1(i_mn) = hmap%eval_gij_dq1(q_thet,qloc, q_zeta)
      gzz_dq1(i_mn) = hmap%eval_gij_dq1(q_zeta,qloc, q_zeta)
      gtt_dq2(i_mn) = hmap%eval_gij_dq2(q_thet,qloc, q_thet)
      gtz_dq2(i_mn) = hmap%eval_gij_dq2(q_thet,qloc, q_zeta)
      gzz_dq2(i_mn) = hmap%eval_gij_dq2(q_zeta,qloc, q_zeta)
      G11(    i_mn) = hmap%eval_gij((/1.0_wp,0.0_wp,0.0_wp/),qloc,(/1.0_wp,0.0_wp,0.0_wp/))
      G21(    i_mn) = hmap%eval_gij((/0.0_wp,1.0_wp,0.0_wp/),qloc,(/1.0_wp,0.0_wp,0.0_wp/))
      G31(    i_mn) = hmap%eval_gij((/0.0_wp,0.0_wp,1.0_wp/),qloc,(/1.0_wp,0.0_wp,0.0_wp/))
     !G12=G21
      G22(    i_mn) = hmap%eval_gij((/0.,1.,0./),qloc,(/0.,1.,0./))
      G32(    i_mn) = hmap%eval_gij((/0.,0.,1./),qloc,(/0.,1.,0./))
      dJh_dq1(i_mn) = hmap%eval_Jh_dq1(qloc)
      dJh_dq2(i_mn) = hmap%eval_Jh_dq2(qloc)
      bt_sJ(  i_mn) = b_thet(i_mn,iGP)*sdetJ(i_mn,iGP)
      bz_sJ(  i_mn) = b_zeta(i_mn,iGP)*sdetJ(i_mn,iGP)
      b_dX1_tz(i_mn)= b_thet(i_mn,iGP)*dX1_dthet(i_mn,iGP)+b_zeta(i_mn,iGP)*dX1_dzeta(i_mn,iGP)
      b_dX2_tz(i_mn)= b_thet(i_mn,iGP)*dX2_dthet(i_mn,iGP)+b_zeta(i_mn,iGP)*dX2_dzeta(i_mn,iGP)
    END DO loop_mn_IP !i_mn
    !averaged quantities
    !X1
    DX1_ss(iGP) =smn_IP_w_GP*SUM(bbcov_sJ(:,iGP)*(   sJ_p(:,iGP)*dX2_dthet(:,iGP) )**2 )
    DX1(   iGP) =smn_IP_w_GP*SUM((sJ_h(:,iGP)*dJh_dq1(:))*(bbcov_sJ(:,iGP)*(  sJ_h(:,iGP)*dJh_dq1(:)) &
                                  -( bt_sJ(:)*(b_thet(:,iGP)*gtt_dq1(:)+2.0*b_zeta(:,iGP)*gtz_dq1(:))    &
                                    +bz_sJ(:)*                              b_zeta(:,iGP)*gzz_dq1(:))  ) )
    DX1_tt(iGP) =smn_IP_w_GP*SUM(bbcov_sJ(:,iGP)*(   sJ_p(:,iGP)*dX2_ds(   :,iGP) )**2  &
                                +bt_sJ(:)*( (2.0_wp*(sJ_p(:,iGP)*dX2_ds(   :,iGP) )      &
                                            *( b_dX1_tz(:    )*G11(:)   &
                                              +b_dX2_tz(:    )*G21(:)   &
                                              +b_zeta(  :,iGP)*G31(:))) &
                                           +b_thet(:,iGP)*G11(:))                                                   )
    DX1_tz(iGP) =smn_IP_w_GP*SUM(bz_sJ(:)*( (2.0_wp*(sJ_p(:,iGP)*dX2_ds(   :,iGP) )      &
                                            *( b_dX1_tz(:    )*G11(:)   &
                                              +b_dX2_tz(:    )*G21(:)   &
                                              +b_zeta(  :,iGP)*G31(:))) &
                                           +b_thet(:,iGP)*2.0*G11(:))                                               )
    DX1_zz(iGP) =smn_IP_w_GP*SUM(b_zeta(:,iGP)*bz_sJ(:)*G11(:))
    !X2
    DX2_ss(iGP) =smn_IP_w_GP*SUM(bbcov_sJ(:,iGP)*(   sJ_p(:,iGP)*dX1_dthet(:,iGP) )**2 )
    DX2(   iGP) =smn_IP_w_GP*SUM((sJ_h(:,iGP)*dJh_dq2(:))*(bbcov_sJ(:,iGP)*(  sJ_h(:,iGP)*dJh_dq2(:)) &
                                  -( bt_sJ(:)*(b_thet(:,iGP)*gtt_dq2(:)+2.0*b_zeta(:,iGP)*gtz_dq2(:))    &
                                    +bz_sJ(:)*                              b_zeta(:,iGP)*gzz_dq2(:))  ) )
    DX2_tt(iGP) =smn_IP_w_GP*SUM(bbcov_sJ(:,iGP)*(   sJ_p(:,iGP)*dX1_ds(   :,iGP) )**2  &
                                +bt_sJ(:)*(-(2.0_wp*(sJ_p(:,iGP)*dX1_ds(   :,iGP) )      &
                                            *( b_dX1_tz(:    )*G21(:)   &
                                              +b_dX2_tz(:    )*G22(:)   &
                                              +b_zeta(  :,iGP)*G32(:))) &
                                           +b_thet(:,iGP)*G22(:))                                                   )
    DX2_tz(iGP) =smn_IP_w_GP*SUM(bz_sJ(:)*(-(2.0_wp*(sJ_p(:,iGP)*dX1_ds(   :,iGP) )      &
                                            *( b_dX1_tz(:    )*G21(:)   &
                                              +b_dX2_tz(:    )*G22(:)   &
                                              +b_zeta(  :,iGP)*G32(:))) &
                                           +b_thet(:,iGP)*2.0*G22(:))                  )
    DX2_zz(iGP) =smn_IP_w_GP*SUM(b_zeta(:,iGP)*bz_sJ(:)*G22(:))
    !LA
    DLA_tt(iGP) =         smn_IP_w_GP*phiPrime2_GP(iGP)*SUM(g_zz(:,iGP)*sdetJ(:,iGP))
    DLA_tz(iGP) = -2.0_wp*smn_IP_w_GP*phiPrime2_GP(iGP)*SUM(g_tz(:,iGP)*sdetJ(:,iGP))
    DLA_zz(iGP) =         smn_IP_w_GP*phiPrime2_GP(iGP)*SUM(g_tt(:,iGP)*sdetJ(:,iGP))
  END DO loop_nGP !iGP
!$OMP END PARALLEL DO
  __PERFOFF('loop_1')

  __PERFON('par_Reduce_D_buf')
  !gather all D** (already stored contiguously in D_buf)
  !THIS SHOULD BE A ALLGATHERV of (nGP_str:nGP_end,:)<=>(1:nGP,:) , NOT A ALLREDUCE (BUT CHEAP ANYWAYS)!
  CALL par_AllReduce(D_buf,'SUM')
  __PERFOFF('par_Reduce_D_buf')
  
  ALLOCATE(D_mn(1:nGP))
  SELECT TYPE(precond_X1); TYPE IS(sll_t_spline_matrix_banded)
    nBase = X1_Base%s%nBase
    modes_str = X1_Base%f%modes_str
    modes_end = X1_Base%f%modes_end
    deg   = X1_base%s%deg
    ALLOCATE(P_BCaxis(1:deg+1,1:2*deg+1),P_BCedge(nBase-deg:nBase,nBase-2*deg:nBase))

      !CHECK =0
    IF(MPIroot)THEN
      IF(SUM(ABS(DX1_ss(:))).LT.REAL(nGP,wp)*1.0E-10)  &
         WRITE(*,*)'WARNING: very small DX1_ss: m,n,SUM(|DX1_ss|)= ',SUM(ABS(DX1_ss(:)))
    END IF

    __PERFON('modes_loop_1')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC)   &
!$OMP   PRIVATE(iMode,iGP,D_mn,iElem,i,j,iBase,tBC,nD,norm_mn) &
!$OMP   FIRSTPRIVATE(P_BCaxis,P_BCedge) & 
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(X1_base,precond_X1,DX1,DX1_tt,DX1_tz,DX1_zz,DX1_ss,X1_BC_Type, &
!$OMP          modes_str,modes_end,nElems,deg,degGP,nBase)
#else
!$OMP   DEFAULT(SHARED)
#endif
    DO iMode=modes_str,modes_end !<<<
      norm_mn=1.0_wp/X1_base%f%snorm_base(iMode)
      CALL precond_X1(iMode)%reset() !set all values to zero
      D_mn(:)=        (DX1+       (X1_Base%f%Xmn(1,iMode)**2)*DX1_tt(:)   &
                          -PRODUCT(X1_Base%f%Xmn(:,iMode))   *DX1_tz(:)   &  !correct sign of theta,zeta derivative!
                          +       (X1_Base%f%Xmn(2,iMode)**2)*DX1_zz(:) ) 
      iGP=1
      DO iElem=1,nElems
        iBase=X1_base%s%base_offset(iElem)
        DO i=0,deg
          DO j=0,deg
            CALL precond_X1(iMode)%add_element(iBase+i,iBase+j,                   &
                                   (SUM( X1_base%s%base_ds_GP(0:degGP,i,iElem)    &
                                        *DX1_ss(iGP:iGP+degGP)                    &
                                        *X1_base%s%base_ds_GP(0:degGP,j,iElem)    &
                                       + X1_base%s%base_GP(0:degGP,i,iElem)       &
                                        *D_mn(iGP:iGP+degGP)                      &
                                        *X1_base%s%base_GP(0:degGP,j,iElem)       &
                                   )*norm_mn)  )
          END DO !j=0,deg
        END DO !i=0,deg
        iGP=iGP+(degGP+1)
      END DO !iElem=1,nElems
      !ACCOUNT FOR BOUNDARY CONDITIONS!
      P_BCaxis=0.0_wp; P_BCedge=0.0_wp
      tBC = X1_BC_Type(BC_AXIS,iMode)
      nD  = X1_base%s%nDOF_BC(tBC) 
      IF(nD.GT.0)THEN
        !save 1:deg rows
        DO i=1,deg+1; DO j=1,MIN(deg+i,nBase)
          P_BCaxis(i,j)=precond_X1(iMode)%get_element(i,j) 
        END DO; END DO !j,i
        P_BCaxis(:,1:MIN(2*deg+1,nBase))     =MATMUL(X1_base%s%R_axis(:,:,tBC),P_BCaxis(:,1:MIN(2*deg+1,nBase))) !also sets rows 1:nD =0
        P_BCaxis(1:nD,1:deg+1) =X1_base%s%A_axis(1:nD,:,tBC)
        DO i=1,deg+1; DO j=1,MIN(deg+i,nBase)
          CALL Precond_X1(iMode)%set_element( i,j,P_BCaxis(i,j))
        END DO; END DO !j,i
      END IF !nDOF_BCaxis>0
      tBC = X1_BC_Type(BC_EDGE,iMode)
      nD  = X1_base%s%nDOF_BC(tBC) 
      IF(nD.GT.0)THEN
        !save nBase-deg:nBase rows
        DO i=nBase-deg,nBase; DO j=MAX(1,i-deg),nBase
          P_BCedge(i,j)=precond_X1(iMode)%get_element(i,j) 
        END DO; END DO !j,i
        P_BCedge(:,MAX(1,nBase-2*deg):nBase) =MATMUL(X1_base%s%R_edge(:,:,tBC),P_BCedge(:,MAX(1,nBase-2*deg):nBase)) !also sets rows nBase-nD+1:nBase =0 
        P_BCedge(nBase-nD+1:nBase,nBase-deg:nBase)=X1_base%s%A_edge(nBase-nD+1:nBase,:,tBC)
        DO i=nBase-deg,nBase; DO j=MAX(1,i-deg),nBase
          CALL Precond_X1(iMode)%set_element( i,j,P_BCedge(i,j) )
        END DO; END DO !j,i
      END IF !nDOF_BCedge>0
    END DO !iMode
!$OMP END PARALLEL DO
    __PERFOFF('modes_loop_1')

    DEALLOCATE(P_BCaxis,P_BCedge)
  END SELECT !TYPE X1

  SELECT TYPE(precond_X2); TYPE IS(sll_t_spline_matrix_banded)
    nBase = X2_Base%s%nBase 
    modes_str = X2_Base%f%modes_str
    modes_end = X2_Base%f%modes_end
    deg   = X2_base%s%deg
    ALLOCATE(P_BCaxis(1:deg+1,1:2*deg+1),P_BCedge(nBase-deg:nBase,nBase-2*deg:nBase))
  
    IF(MPIroot)THEN
      !CHECK =0
      IF(SUM(ABS(DX2_ss(:))).LT.REAL(nGP,wp)*1.0E-10)  &
           WRITE(*,*)'WARNING: very small DX2_ss: m,n,SUM(|DX2_ss|)= ',SUM(ABS(DX2_ss(:)))
    END IF
  
    __PERFON('modes_loop_2')
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC)  &
!$OMP   PRIVATE(iMode,iGP,D_mn,iElem,i,j,iBase,tBC,nD,norm_mn) &
!$OMP   FIRSTPRIVATE(P_BCaxis,P_BCedge)  &
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(X2_base,precond_X2,DX2,DX2_tt,DX2_tz,DX2_zz,DX2_ss,X2_BC_Type, &
!$OMP          modes_str,modes_end,nElems,deg,degGP,nBase) 
#else
!$OMP   DEFAULT(SHARED)
#endif
    DO iMode=modes_str,modes_end !<<<
      norm_mn=1.0_wp/X2_base%f%snorm_base(iMode)
      CALL precond_X2(iMode)%reset() !set all values to zero
      D_mn(:)=        (DX2+       (X2_Base%f%Xmn(1,iMode)**2)*DX2_tt(:)   &
                          -PRODUCT(X2_Base%f%Xmn(:,iMode))   *DX2_tz(:)   & !correct sign of theta,zeta derivative!
                          +       (X2_Base%f%Xmn(2,iMode)**2)*DX2_zz(:) ) 
      iGP=1
      DO iElem=1,nElems
        iBase=X2_base%s%base_offset(iElem)
        DO i=0,deg
          DO j=0,deg
            CALL precond_X2(iMode)%add_element(iBase+i,iBase+j,                   &
                                   (SUM( X2_base%s%base_ds_GP(0:degGP,i,iElem)    &
                                        *DX2_ss(iGP:iGP+degGP)                    &
                                        *X2_base%s%base_ds_GP(0:degGP,j,iElem)    &
                                       + X2_base%s%base_GP(0:degGP,i,iElem)       &
                                        *D_mn(iGP:iGP+degGP)                      &
                                        *X2_base%s%base_GP(0:degGP,j,iElem)       &
                                   )*norm_mn)  )
          END DO !j=0,deg
        END DO !i=0,deg
        iGP=iGP+(degGP+1)
      END DO !iElem=1,nElems
      !ACCOUNT FOR BOUNDARY CONDITIONS!
      P_BCaxis=0.0_wp; P_BCedge=0.0_wp
      tBC = X2_BC_Type(BC_AXIS,iMode)
      nD  = X2_base%s%nDOF_BC(tBC) 
      IF(nD.GT.0)THEN
        !save 1:deg rows
        DO i=1,deg+1; DO j=1,MIN(deg+i,nBase)
          P_BCaxis(i,j)=precond_X2(iMode)%get_element(i,j) 
        END DO; END DO !j,i
        P_BCaxis(:,1:MIN(2*deg+1,nBase))=MATMUL(X2_base%s%R_axis(:,:,tBC),P_BCaxis(:,1:MIN(2*deg+1,nBase))) !also sets rows 1:nD =0
        P_BCaxis(1:nD,1:deg+1) =X2_base%s%A_axis(1:nD,:,tBC)
        DO i=1,deg+1; DO j=1,MIN(deg+i,nBase)
          CALL Precond_X2(iMode)%set_element( i,j,P_BCaxis(i,j))
        END DO; END DO !j,i
      END IF !nDOF_BCaxis>0
      tBC = X2_BC_Type(BC_EDGE,iMode)
      nD  = X2_base%s%nDOF_BC(tBC) 
      IF(nD.GT.0)THEN
        !save nBase-deg:nBase rows
        DO i=nBase-deg,nBase; DO j=MAX(1,i-deg),nBase
          P_BCedge(i,j)=precond_X2(iMode)%get_element(i,j) 
        END DO; END DO !j,i
        P_BCedge(:,MAX(1,nBase-2*deg):nBase) =MATMUL(X2_base%s%R_edge(:,:,tBC),P_BCedge(:,MAX(1,nBase-2*deg):nBase)) !also sets rows nBase-nD+1:nBase =0 
        P_BCedge(nBase-nD+1:nBase,nBase-deg:nBase)=X2_base%s%A_edge(nBase-nD+1:nBase,:,tBC)
        DO i=nBase-deg,nBase; DO j=MAX(1,i-deg),nBase
          CALL Precond_X2(iMode)%set_element( i,j,P_BCedge(i,j) )
        END DO; END DO !j,i
      END IF !nDOF_BCedge>0
    END DO !iMode
!$OMP END PARALLEL DO
    __PERFOFF('modes_loop_2')
  
    DEALLOCATE(P_BCaxis,P_BCedge)
  END SELECT !TYPE X2
 
  SELECT TYPE(precond_LA); TYPE IS(sll_t_spline_matrix_banded)
    nBase = LA_Base%s%nBase 
    modes_str = LA_Base%f%modes_str
    modes_end = LA_Base%f%modes_end
    deg   = LA_base%s%deg
    ALLOCATE(P_BCaxis(1:deg+1,1:2*deg+1),P_BCedge(nBase-deg:nBase,nBase-2*deg:nBase))

    __PERFON('modes_loop_3')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC)   &
!$OMP   PRIVATE(iMode,iGP,D_mn,iElem,i,j,iBase,tBC,nD,norm_mn) &
!$OMP   FIRSTPRIVATE(P_BCaxis,P_BCedge) & 
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(LA_base,precond_LA,DLA_tt,DLA_tz,DLA_zz,LA_BC_Type, &
!$OMP          modes_str,modes_end,nElems,deg,degGP,nBase,nGP) 
#else
!$OMP   DEFAULT(SHARED)       
#endif       
    DO iMode=modes_str,modes_end !<<<
      norm_mn=1.0_wp/LA_base%f%snorm_base(iMode)
      CALL precond_LA(iMode)%reset() !set all values to zero
      IF(LA_base%f%zero_odd_even(iMode) .NE. MN_ZERO) THEN !MN_ZERO should not exist
        D_mn(:)=(        (        (LA_Base%f%Xmn(1,iMode)**2)*DLA_tt(:) &
                          -PRODUCT(LA_Base%f%Xmn(:,iMode))   *DLA_tz(:) & !correct sign of theta,zeta derivative!
                          +       (LA_Base%f%Xmn(2,iMode)**2)*DLA_zz(:) ))*norm_mn
        !CHECK =0
        IF(SUM(ABS(D_mn(:))).LT.REAL(nGP,wp)*1.0E-10) WRITE(*,*)'WARNING: small DLA: m,n,SUM(|DLA_mn|)= ', &
             LA_Base%f%Xmn(1:2,iMode),SUM(D_mn(:))
        iGP=1
        DO iElem=1,nElems
          iBase=LA_base%s%base_offset(iElem)
          DO i=0,deg
            DO j=0,deg
              CALL precond_LA(iMode)%add_element(iBase+i,iBase+j,               &
                                     (SUM( LA_base%s%base_GP(0:degGP,i,iElem)   &
                                          *D_mn(iGP:iGP+degGP)                  &
                                          *LA_base%s%base_GP(0:degGP,j,iElem))) )
            END DO !j=0,deg
          END DO !i=0,deg
          iGP=iGP+(degGP+1)
        END DO !iElem=1,nElems
      ELSE !safety, unit matrix, if MN_ZERO exists
        DO iBase=1,nBase
          CALL precond_LA(iMode)%set_element(iBase,iBase,1.0_wp)
        END DO
      END IF !MN_ZERO does not exists
      !ACCOUNT FOR BOUNDARY CONDITIONS!
      P_BCaxis=0.0_wp; P_BCedge=0.0_wp
      tBC = LA_BC_Type(BC_AXIS,iMode)
      nD  = LA_base%s%nDOF_BC(tBC) 
      IF(nD.GT.0)THEN
        !save 1:deg rows
        DO i=1,deg+1; DO j=1,MIN(deg+i,nBase)
          P_BCaxis(i,j)=precond_LA(iMode)%get_element(i,j) 
        END DO; END DO !j,i
        P_BCaxis(:,1:MIN(2*deg+1,nBase)) =MATMUL(LA_base%s%R_axis(:,:,tBC),P_BCaxis(:,1:MIN(2*deg+1,nBase))) !also sets rows 1:nD =0
        P_BCaxis(1:nD,1:deg+1) =LA_base%s%A_axis(1:nD,:,tBC)
        DO i=1,deg+1; DO j=1,MIN(deg+i,nBase)
          CALL Precond_LA(iMode)%set_element( i,j,P_BCaxis(i,j))
        END DO; END DO !j,i
      END IF !nDOF_BCaxis>0
      tBC = LA_BC_Type(BC_EDGE,iMode)
      nD  = LA_base%s%nDOF_BC(tBC) 
      IF(nD.GT.0)THEN
        !save nBase-deg:nBase rows
        DO i=nBase-deg,nBase; DO j=MAX(1,i-deg),nBase
          P_BCedge(i,j)=precond_LA(iMode)%get_element(i,j) 
        END DO; END DO !j,i
        P_BCedge(:,MAX(1,nBase-2*deg):nBase) =MATMUL(LA_base%s%R_edge(:,:,tBC),P_BCedge(:,MAX(1,nBase-2*deg):nBase)) !also sets rows nBase-nD+1:nBase =0 
        P_BCedge(nBase-nD+1:nBase,nBase-deg:nBase)=LA_base%s%A_edge(nBase-nD+1:nBase,:,tBC)
        DO i=nBase-deg,nBase; DO j=MAX(1,i-deg),nBase
          CALL Precond_LA(iMode)%set_element( i,j,P_BCedge(i,j) )
        END DO; END DO !j,i
      END IF !nDOF_BCedge>0
    END DO !iMode
!$OMP END PARALLEL DO
    __PERFOFF('modes_loop_3')

    DEALLOCATE(P_BCaxis,P_BCedge)
  END SELECT !TYPE LA

  DEALLOCATE(D_mn)

  !  WRITE(*,*)'    ---> FACTORIZE PRECONDITIONER MATRICES ...'

  SELECT TYPE(precond_X1); TYPE IS(sll_t_spline_matrix_banded)
    __PERFON('modes_loop_4')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode) &
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(precond_X1,X1_Base)
#else
!$OMP   DEFAULT(SHARED)       
#endif       
    DO iMode=X1_Base%f%modes_str,X1_Base%f%modes_end !<<<
      CALL precond_X1(iMode)%factorize()
    END DO !iMode
!$OMP END PARALLEL DO 
    __PERFOFF('modes_loop_4')
  END SELECT !TYPE X1

  SELECT TYPE(precond_X2); TYPE IS(sll_t_spline_matrix_banded)
    __PERFON('modes_loop_5')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode) &
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(precond_X2,X2_Base)
#else
!$OMP   DEFAULT(SHARED)       
#endif       
    DO iMode=X2_base%f%modes_str,X2_Base%f%modes_end !<<<
      CALL precond_X2(iMode)%factorize()
    END DO !iMode
!$OMP END PARALLEL DO 
    __PERFOFF('modes_loop_5')
  END SELECT !TYPE X2

  SELECT TYPE(precond_LA); TYPE IS(sll_t_spline_matrix_banded)
    __PERFON('modes_loop_6')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) PRIVATE(iMode) &
#ifdef __INTEL_COMPILER
!$OMP   DEFAULT(NONE) SHARED(precond_LA,LA_Base)
#else
!$OMP   DEFAULT(SHARED)       
#endif       
    DO iMode=LA_base%f%modes_str,LA_Base%f%modes_end !<<<
      CALL precond_LA(iMode)%factorize()
    END DO !iMode
!$OMP END PARALLEL DO 
    __PERFOFF('modes_loop_6')
  END SELECT !TYPE LA
    

END SUBROUTINE BuildPrecond


!===================================================================================================================================
!> Apply preconditioner matrix for single mode of one variable
!!
!===================================================================================================================================
SUBROUTINE ApplyPrecond(nBase,precond,F_inout)
! MODULES
  USE MODgvec_base,   ONLY: t_base
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER                         ,INTENT(IN   ) ::  nBase   !! length of inout vector
  TYPE(sll_t_spline_matrix_banded),INTENT(IN   ) ::  precond !! preconditioner matrix (factorized!)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),    INTENT(INOUT) :: F_inout(1:nBase)    !! apply preconditioner on this force vector
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL precond%solve_inplace(1,F_inout(1:nBase))
END SUBROUTINE ApplyPrecond


!===================================================================================================================================
!> check force with finite difference 
!!
!===================================================================================================================================
SUBROUTINE checkEvalForce(Uin,fileID)
! MODULES
  USE MODgvec_Globals       , ONLY: testlevel
  USE MODgvec_MHD3D_Vars    , ONLY: X1_base,X2_base,LA_base,PrecondType
  USE MODgvec_MHD3D_visu    , ONLY: WriteDataMN_visu
  USE MODgvec_sol_var_MHD3D , ONLY: t_sol_var_MHD3D
  USE MODgvec_Output_Vars   , ONLY: outputLevel
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin   !! input solution 
  INTEGER               , INTENT(IN   ) :: FileID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                               :: iBase,nBase,iMode,modes,JacCheck,PrecondTypeTmp
  REAL(wp)                              :: W_MHD3D_in,eps_glob,eps,pTg
  CLASS(t_sol_var_MHD3D),ALLOCATABLE    :: Ucopy
  CLASS(t_sol_var_MHD3D),ALLOCATABLE    :: Utest
  CLASS(t_sol_var_MHD3D),ALLOCATABLE    :: Ftest 
  CLASS(t_sol_var_MHD3D),ALLOCATABLE    :: Feval 
  CHARACTER(LEN=60)                     :: fname
!===================================================================================================================================
  IF(testLevel.EQ.-1) RETURN
  PrecondTypeTmp=PrecondType
  ALLOCATE(t_sol_var_MHD3D :: Utest)
  ALLOCATE(t_sol_var_MHD3D :: Ucopy)
  ALLOCATE(t_sol_var_MHD3D :: Ftest)
  ALLOCATE(t_sol_var_MHD3D :: Feval)

  CALL Ucopy%copy(Uin)
  CALL Utest%copy(Ucopy) !initialize
  CALL Ftest%copy(Ucopy) !initialize
  CALL Ftest%set_to(0.0_wp)
  CALL Feval%copy(Ftest)

  JacCheck=1 !abort if detJ<0 
  W_MHD3D_in=EvalEnergy(Utest,.TRUE.,JacCheck)
  eps_glob=1.0e-10_wp*W_MHD3D_in

  nBase = X1_Base%s%nBase 
  modes = X1_Base%f%modes
  SWRITE(*,*)'Test X1',nBase,modes
  DO iBase=1,nBase
    eps=eps_glob
    DO iMode=1,modes
      Utest%X1(iBase,iMode)= Ucopy%X1(iBase,iMode)+eps
      Utest%W_MHD3D        = EvalEnergy(Utest,.TRUE.,JacCheck)
      Utest%X1(iBase,iMode)= Ucopy%X1(iBase,iMode)

      Ftest%X1(iBase,iMode)= -(Utest%W_MHD3D-W_MHD3D_in)/eps
    END DO
  END DO 

  nBase = X2_Base%s%nBase 
  modes = X2_Base%f%modes
  SWRITE(*,*)'Test X2',nBase,modes
  DO iBase=1,nBase
    eps=eps_glob
    DO iMode=1,modes
      Utest%X2(iBase,iMode)= Ucopy%X2(iBase,iMode)+eps
      Utest%W_MHD3D        = EvalEnergy(Utest,.TRUE.,JacCheck)
      Utest%X2(iBase,iMode)= Ucopy%X2(iBase,iMode)

      Ftest%X2(iBase,iMode)= -(Utest%W_MHD3D-W_MHD3D_in)/eps
    END DO
  END DO 

  nBase = LA_Base%s%nBase 
  modes = LA_Base%f%modes
  SWRITE(*,*)'Test LA',nBase,modes
  DO iBase=1,nBase
    eps=eps_glob 
    DO iMode=1,modes
      Utest%LA(iBase,iMode)= Ucopy%LA(iBase,iMode)+eps
      Utest%W_MHD3D        = EvalEnergy(Utest,.TRUE.,JacCheck)
      Utest%LA(iBase,iMode)= Ucopy%LA(iBase,iMode)

      Ftest%LA(iBase,iMode)= -(Utest%W_MHD3D-W_MHD3D_in)/eps
    END DO
  END DO 

  SWRITE(UNIT_stdOut,'(A,3E21.11)')'Norm of test force without BC |X1|,|X2|,|LA|: ',SQRT(Ftest%norm_2())

  PrecondType=-1
  CALL EvalForce(Ucopy,.TRUE.,JacCheck,Feval,noBC=.FALSE.)
  SWRITE(UNIT_stdOut,'(A,3E21.11)')'Norm of eval force with BC |X1|,|X2|,|LA|: ',SQRT(Feval%norm_2())
  CALL EvalForce(Ucopy,.TRUE.,JacCheck,Feval,noBC=.TRUE.)
  SWRITE(UNIT_stdOut,'(A,3E21.11)')'Norm of eval force without BC |X1|,|X2|,|LA|: ',SQRT(Feval%norm_2())

  IF((testlevel.GE.2).AND.MPIroot)THEN
  WRITE(*,*)'-----------------------'
  modes = X1_Base%f%modes
  DO iMode=1,modes
    WRITE(*,'(A,2I4,A,E11.3)')'X1 mode',X1_base%f%Xmn(:,iMode),', maxUMode= ',&
                                            MAXVAL(ABS(Ucopy%X1(:,iMode)))
    WRITE(*,'(4X,2(A,E11.3))')    'minfeval= ',MINVAL((Feval%X1(:,iMode))), & 
                                ', maxfeval= ',MAXVAL((Feval%X1(:,iMode)))
    WRITE(*,'(4X,3(A,E11.3))')    'minftest= ',MINVAL((Ftest%X1(:,iMode))), & 
                                ', maxftest= ',MAXVAL((Ftest%X1(:,iMode))), & 
                     ', maxdiff force= ' , MAXVAL(ABS( Ftest%X1(:,iMode) &
                                                      -Feval%X1(:,iMode)))
  END DO
  WRITE(*,*)'-----------------------'
  modes = X2_Base%f%modes
  DO iMode=1,modes
    WRITE(*,'(A,2I4,A,E11.3)')'X2 mode',X2_base%f%Xmn(:,iMode),', maxUMode= ',&
                                            MAXVAL(ABS(Ucopy%X2(:,iMode)))
    WRITE(*,'(4X,2(A,E11.3))')    'minfeval= ',MINVAL((Feval%X2(:,iMode))), & 
                                ', maxfeval= ',MAXVAL((Feval%X2(:,iMode)))
    WRITE(*,'(4X,3(A,E11.3))')    'minftest= ',MINVAL((Ftest%X2(:,iMode))), & 
                                ', maxftest= ',MAXVAL((Ftest%X2(:,iMode))), & 
                     ', maxdiff force= ' , MAXVAL(ABS( Ftest%X2(:,iMode) &
                                                      -Feval%X2(:,iMode)))
  END DO
  WRITE(*,*)'-----------------------'
  modes = LA_Base%f%modes
  DO iMode=1,modes
    WRITE(*,'(A,2I4,A,E11.3)')'LA mode',LA_base%f%Xmn(:,iMode),', maxUMode= ',&
                                            MAXVAL(ABS(Ucopy%LA(:,iMode)))
    WRITE(*,'(4X,2(A,E11.3))')    'minfeval= ',MINVAL((Feval%LA(:,iMode))), & 
                                ', maxfeval= ',MAXVAL((Feval%LA(:,iMode)))
    WRITE(*,'(4X,3(A,E11.3))')    'minftest= ',MINVAL((Ftest%LA(:,iMode))), & 
                                ', maxftest= ',MAXVAL((Ftest%LA(:,iMode))), & 
                     ', maxdiff force= ' , MAXVAL(ABS( Ftest%LA(:,iMode) &
                                                      -Feval%LA(:,iMode)))
  END DO
  END IF !testlevel >=2
  
  IF(testlevel.GE.3.AND.MPIroot)THEN
    ASSOCIATE(np1d=>2*(X1_base%s%deg+3) )
    WRITE(fname,'(A,"_",I4.4,"_",I8.8)')'Ftest_X1_cos_',outputLevel,fileID
    CALL writeDataMN_visu(np1d,fname,'X1_cos_',0,X1_base,Ftest%X1)
    WRITE(fname,'(A,"_",I4.4,"_",I8.8)')'Ftest_X2_sin_',outputLevel,fileID
    CALL writeDataMN_visu(np1d,fname,'X2_sin_',0,X2_base,Ftest%X2)
    WRITE(fname,'(A,"_",I4.4,"_",I8.8)')'Ftest_LA_sin_',outputLevel,fileID
    CALL writeDataMN_visu(np1d,fname,'LA_sin_',0,LA_base,Ftest%LA)
    WRITE(fname,'(A,"_",I4.4,"_",I8.8)')'Feval_X1_cos_',outputLevel,fileID
    CALL writeDataMN_visu(np1d,fname,'X1_cos_',0,X1_base,Feval%X1)
    WRITE(fname,'(A,"_",I4.4,"_",I8.8)')'Feval_X2_sin_',outputLevel,fileID
    CALL writeDataMN_visu(np1d,fname,'X2_sin_',0,X2_base,Feval%X2)
    WRITE(fname,'(A,"_",I4.4,"_",I8.8)')'Feval_LA_sin_',outputLevel,fileID
    CALL writeDataMN_visu(np1d,fname,'LA_sin_',0,LA_base,Feval%LA)
    END ASSOCIATE !np1d
  END IF !testlevel >=3

  !check preconditioner positivity: preconditioned force^T * force should be <0
  PrecondType=-1
  CALL EvalForce(Ucopy,.TRUE.,JacCheck,Feval)
  PrecondType= 1
  CALL EvalForce(Ucopy,.TRUE.,JacCheck,Ftest)
  IF(MPIroot)THEN
    pTg= - SUM(Ftest%q*Feval%q)
    IF (pTg.GT.0.0_wp) THEN
      WRITE(*,'(A,E11.3)')'WARING: PRECONDITIONER is not positive symmetric',ptg
    ELSE
      WRITE(*,'(A,E11.3)')'OK: PRECONDITIONER is positive symmetric',ptg
    END IF
  END IF
  
  PrecondType=PrecondTypeTmp !save back 

  CALL Utest%free()
  CALL Ucopy%free()
  CALL Ftest%free()
  CALL Feval%free()


  DEALLOCATE(Utest)
  DEALLOCATE(Ucopy)
  DEALLOCATE(Ftest)
  DEALLOCATE(Feval)

END SUBROUTINE checkEvalForce


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeMHD3D_EvalFunc()
! MODULES
  USE MODgvec_MHD3D_Vars,ONLY:X1_base,X2_base,LA_base,PrecondType
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iMode
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
  SDEALLOCATE(chiPrime_GP  )
  SDEALLOCATE(phiPrime_GP  )
  SDEALLOCATE(phiPrime2_GP )
  SDEALLOCATE(J_h          )
  SDEALLOCATE(J_p          )
  SDEALLOCATE(sJ_h         )
  SDEALLOCATE(sJ_p         )
  SDEALLOCATE(detJ         )
  SDEALLOCATE(sdetJ        )
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
  SDEALLOCATE(bbcov_sJ     )
  SDEALLOCATE(g_tt         )
  SDEALLOCATE(g_tz         )
  SDEALLOCATE(g_zz         )

  IF(PrecondType.GT.0)THEN
    NULLIFY(DX1_tt); NULLIFY(DX1_tz); NULLIFY(DX1_zz); NULLIFY(DX1); NULLIFY(DX1_ss)
    NULLIFY(DX2_tt); NULLIFY(DX2_tz); NULLIFY(DX2_zz); NULLIFY(DX2); NULLIFY(DX2_ss)
    NULLIFY(DLA_tt); NULLIFY(DLA_tz); NULLIFY(DLA_zz)         
    DEALLOCATE(D_buf)
    IF(MPIroot)THEN
      SELECT TYPE(precond_X1); TYPE IS(sll_t_spline_matrix_banded)
        DO iMode=1,X1_Base%f%modes
          CALL precond_X1(iMode)%free()
        END DO !iMode
      END SELECT !TYPE
      DEALLOCATE(precond_X1)
     
      SELECT TYPE(precond_X2); TYPE IS(sll_t_spline_matrix_banded)
        DO iMode=1,X2_base%f%modes
          CALL precond_X2(iMode)%free()
        END DO !iMode
      END SELECT !TYPE
      DEALLOCATE(precond_X2)
     
      SELECT TYPE(precond_LA); TYPE IS(sll_t_spline_matrix_banded)
        DO iMode=1,LA_base%f%modes
          CALL precond_LA(iMode)%free()
        END DO !iMode
      END SELECT !TYPE
      DEALLOCATE(precond_LA)
    END IF !MPIroot
  END IF !PrecondType>0

END SUBROUTINE FinalizeMHD3D_EvalFunc

END MODULE MODgvec_MHD3D_EvalFunc
