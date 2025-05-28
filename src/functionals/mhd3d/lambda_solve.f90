!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **lambda_solve**
!!
!! CONTAINS routine to solve for lambda at a specific flux surface (for example for the boudnary condition at the last flux surface)
!!
!===================================================================================================================================
MODULE MODgvec_lambda_solve
! MODULES
USE MODgvec_Globals, ONLY:wp,UNIT_StdOut,abort
IMPLICIT NONE
PUBLIC

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Solve for lambda on one given flux surface (spos_in), using weak form of J^s=0: d/dzeta(B_theta)-d/dtheta(B_zeta)=0
!!
!! Note that the mapping defined by  X1 and X2 must be fully initialized, since derivatives in s must be taken!
!!
!===================================================================================================================================
SUBROUTINE Lambda_solve(spos_in,hmap_in,hmap_xv,X1_base_in,X2_base_in,LA_fbase_in,X1_in,X2_in,LA_s,phiPrime_s,chiPrime_s)
! MODULES
  USE MODgvec_Globals,       ONLY:n_warnings_occured
  USE MODgvec_base          ,ONLY: t_base
  USE MODgvec_fbase         ,ONLY: t_fbase
  USE MODgvec_hmap          ,ONLY: PP_T_HMAP,PP_T_HMAP_AUXVAR
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base),INTENT(IN)  :: X1_base_in,X2_base_in           !< base classes belong to solution X1_in,X2_in
  TYPE(t_fbase),INTENT(IN) :: LA_fbase_in                     !< base class belong to solution LA_s
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP), INTENT(IN) :: hmap_in
  TYPE(PP_T_HMAP_AUXVAR), INTENT(IN) :: hmap_xv(X1_base_in%f%mn_IP)  !< auxiliary variables for hmap, must be pre-computed
#else
  CLASS(PP_T_HMAP), INTENT(IN) :: hmap_in
  CLASS(PP_T_HMAP_AUXVAR), INTENT(IN) :: hmap_xv(X1_base_in%f%mn_IP)  !< auxiliary variables for hmap, must be pre-computed
#endif
  REAL(wp)     , INTENT(IN) :: spos_in                  !! s position to evaluate lambda
  REAL(wp)     , INTENT(IN) :: X1_in(1:X1_base_in%s%nBase,1:X1_base_in%f%modes) !! U%X1 variable, is reshaped to 2D at input
  REAL(wp)     , INTENT(IN) :: X2_in(1:X2_base_in%s%nBase,1:X2_base_in%f%modes) !! U%X2 variable, is reshaped to 2D at input
  REAL(wp)     , INTENT(IN) :: phiPrime_s !! toroidal flux derivative phi' at the position s
  REAL(wp)     , INTENT(IN) :: chiPrime_s !! poloidal flux derivative chi' at the position s
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)     , INTENT(  OUT) :: LA_s(1:LA_fbase_in%modes) !! lambda at spos
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                               :: iMode,i_mn,mn_IP
  REAL(wp)                              :: spos,minJ
  REAL(wp),DIMENSION(1:X1_base_in%f%modes) :: X1_s,X1_ds !! X1 solution at spos
  REAL(wp),DIMENSION(1:X2_base_in%f%modes) :: X2_s,X2_ds !! X1 solution at spos
  REAL(wp),DIMENSION(1:X1_base_in%f%mn_IP) :: X1_s_IP,dX1ds,dX1dthet,dX1dzeta, & !mn_IP should be same for all!
                                              X2_s_IP,dX2ds,dX2dthet,dX2dzeta, &
                                              detJ,gam_tt,gam_tz,gam_zz

!===================================================================================================================================
  __PERFON('lambda_solve')

  spos=MIN(1.0_wp-1.0e-12_wp,MAX(1.0e-04_wp,spos_in))
  mn_IP = X1_base_in%f%mn_IP
  IF(X2_base_in%f%mn_IP.NE.mn_IP) CALL abort(__STAMP__,&
                                             'X2 mn_IP /= X1 mn_IP')
  IF(LA_fbase_in%mn_IP .NE.mn_IP)  CALL abort(__STAMP__,&
                                             'LA mn_IP /= X1 mn_IP')


!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) &
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iMode)   &
!$OMP   SHARED(spos,X1_s,X1_ds,X1_base_in,X1_in)
  DO iMode=1,X1_base_in%f%modes
    X1_s( iMode)  = X1_base_in%s%evalDOF_s(spos,      0,X1_in(:,iMode))
    X1_ds(iMode)  = X1_base_in%s%evalDOF_s(spos,DERIV_S,X1_in(:,iMode))
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) &
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iMode)   &
!$OMP   SHARED(spos,X2_s,X2_ds,X2_base_in,X2_in)
  DO iMode=1,X2_base_in%f%modes
    X2_s( iMode)  = X2_base_in%s%evalDOF_s(spos,      0,X2_in(:,iMode))
    X2_ds(iMode)  = X2_base_in%s%evalDOF_s(spos,DERIV_S,X2_in(:,iMode))
  END DO
!$OMP END PARALLEL DO

  X1_s_IP  = X1_base_in%f%evalDOF_IP(         0,X1_s )
  dX1ds    = X1_base_in%f%evalDOF_IP(         0,X1_ds)
  dX1dthet = X1_base_in%f%evalDOF_IP(DERIV_THET,X1_s )
  dX1dzeta = X1_base_in%f%evalDOF_IP(DERIV_ZETA,X1_s )

  X2_s_IP  = X2_base_in%f%evalDOF_IP(         0,X2_s )
  dX2ds    = X2_base_in%f%evalDOF_IP(         0,X2_ds)
  dX2dthet = X2_base_in%f%evalDOF_IP(DERIV_THET,X2_s )
  dX2dzeta = X2_base_in%f%evalDOF_IP(DERIV_ZETA,X2_s )



!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) &
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)  &
!$OMP   SHARED(mn_IP,X1_s_IP,X2_s_IP,dX1ds,dX2ds,dX1dthet,dX2dthet,detJ,hmap_in,hmap_xv)
  DO i_mn=1,mn_IP
    detJ(i_mn)= (dX1ds(i_mn)*dX2dthet(i_mn)-dX1dthet(i_mn)*dX2ds(i_mn)) &
               *hmap_in%eval_Jh_aux(X1_s_IP(i_mn),X2_s_IP(i_mn),hmap_xv(i_mn)) !J_p*J_h
  END DO !i_mn
!$OMP END PARALLEL DO

  minJ=MINVAL(detJ)
  IF(minJ.LT.1.0e-12) THEN
    n_warnings_occured=n_warnings_occured+1
    i_mn= MINLOC(detJ(:),1)
    WRITE(UNIT_stdOut,'(4X,A8,I8,4(A,E11.3))')'WARNING ',n_warnings_occured, &
                                                 ' : min(J)= ',MINVAL(detJ),' at s= ',spos, &
                                                                       ' theta= ',X1_base_in%f%x_IP(1,i_mn), &
                                                                        ' zeta= ',X1_base_in%f%x_IP(2,i_mn)
    i_mn= MAXLOC(detJ(:),1)
    WRITE(UNIT_stdOut,'(4X,16X,4(A,E11.3))')'     ...max(J)= ',MAXVAL(detJ),' at s= ',spos, &
                                                                       ' theta= ',X1_base_in%f%x_IP(1,i_mn), &
                                                                        ' zeta= ',X1_base_in%f%x_IP(2,i_mn)
!    CALL abort(__STAMP__, &
!        'Lambda_solve: Jacobian smaller that  1.0e-12!!!' )
  END IF

!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) &
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)  &
!$OMP   SHARED(mn_IP,gam_tt,gam_tz,gam_zz,detJ,hmap_in,hmap_xv,X1_s_IP,X2_s_IP,dX1dthet,dX2dthet,dX1dzeta,dX2dzeta)
  DO i_mn=1,mn_IP
    gam_tt(i_mn) = hmap_in%eval_gij_aux(dX1dthet(i_mn),dX2dthet(i_mn),0.0_wp, &
                                         X1_s_IP(i_mn), X2_s_IP(i_mn),        &
                                        dX1dthet(i_mn),dX2dthet(i_mn),0.0_wp, &
                                        hmap_xv(i_mn)) / detJ(i_mn)
    gam_tz(i_mn) = hmap_in%eval_gij_aux(dX1dthet(i_mn),dX2dthet(i_mn),0.0_wp, &
                                         X1_s_IP(i_mn), X2_s_IP(i_mn),        &
                                        dX1dzeta(i_mn),dX2dzeta(i_mn),1.0_wp, &
                                        hmap_xv(i_mn)) / detJ(i_mn)
    gam_zz(i_mn) = hmap_in%eval_gij_aux(dX1dzeta(i_mn),dX2dzeta(i_mn),1.0_wp, &
                                         X1_s_IP(i_mn), X2_s_IP(i_mn),        &
                                        dX1dzeta(i_mn),dX2dzeta(i_mn),1.0_wp, &
                                        hmap_xv(i_mn)) / detJ(i_mn)
  END DO !i_mn
!$OMP END PARALLEL DO

  CALL lambda_setup_and_solve(LA_fbase_in,phiPrime_s,ChiPrime_s,gam_tt,gam_tz,gam_zz,LA_s)

  __PERFOFF('lambda_solve')

END SUBROUTINE Lambda_solve

SUBROUTINE Lambda_setup_and_solve(LA_fbase_in,phiPrime_s,ChiPrime_s,gam_tt,gam_tz,gam_zz,LA_s)
  ! MODULES
    USE MODgvec_LinAlg,ONLY: SOLVE
    USE MODgvec_fbase ,ONLY: t_fbase
    IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    TYPE(t_fbase),INTENT(IN)        :: LA_fbase_in           !< base classes belong to solution U_in
    REAL(wp),INTENT(IN)              :: phiPrime_s,ChiPrime_s   !! toroidal and poloidal flux s derivatives at s_pos
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP), INTENT(IN) :: gam_tt  !! g_tt/J evaluated on IP points
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP), INTENT(IN) :: gam_tz  !! g_tz/J evaluated on IP points
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP), INTENT(IN) :: gam_zz  !! g_zz/J evaluated on IP points
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)     , INTENT(  OUT) :: LA_s(1:LA_fbase_in%modes) !! lambda at spos
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER                               :: iMode,jMode,mn_IP,LA_modes
    REAL(wp)                              :: Amat(1:LA_fbase_in%modes,1:LA_fbase_in%modes)
    REAL(wp),DIMENSION(1:LA_fbase_in%modes) :: RHS,sAdiag
  !  REAL(wp)                              :: gam_ta_da,gam_za_da
    REAL(wp)                              :: sum_gam_ta_da,sum_gam_za_da
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP,1:LA_fbase_in%modes) :: gam_ta_da,gam_za_da
  !===================================================================================================================================
    __PERFON('setup_1')
    LA_modes=LA_fbase_in%modes
    mn_IP   =LA_fbase_in%mn_IP

  !$OMP PARALLEL DO        &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iMode)  &
  !$OMP   SHARED(mn_IP,LA_modes,LA_fbase_in,sAdiag)
    !estimate of 1/Adiag for preconditioning (=1 for m=n=0)
    DO iMode=1,LA_modes
    sAdiag(iMode)=1.0_wp/(MAX(1.0_wp,REAL((LA_fbase_in%Xmn(1,iMode))**2+LA_fbase_in%Xmn(2,iMode)**2 ,wp) )*REAL(mn_IP,wp))
        !sAdiag(iMode)=1.0_wp
    END DO !iMode
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO        &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iMode)        &
  !$OMP   SHARED(gam_ta_da,gam_za_da,LA_modes,gam_tt,gam_tz,gam_zz,LA_fbase_in)
    DO iMode=1,LA_modes
      gam_ta_da(:,iMode)=gam_tz(:)*LA_fbase_in%base_dthet_IP(:,iMode) - gam_tt(:)*LA_fbase_in%base_dzeta_IP(:,iMode)
      gam_za_da(:,iMode)=gam_zz(:)*LA_fbase_in%base_dthet_IP(:,iMode) - gam_tz(:)*LA_fbase_in%base_dzeta_IP(:,iMode)
    END DO!iMode
  !$OMP END PARALLEL DO

  __PERFOFF('setup_1')
  __PERFON('setup_Amat')

!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(jMode)        &
!$OMP   SHARED(Amat,gam_ta_da,gam_za_da,LA_fbase_in,LA_modes,PhiPrime_s,sAdiag)
    DO jMode=1,LA_modes
      !m=n=0 should not be in lambda, but check
      IF (LA_fbase_in%zero_odd_even(jMode).NE.MN_ZERO) THEN
        CALL LA_fbase_in%projectIPtoDOF(.FALSE., PhiPrime_s,DERIV_ZETA,gam_ta_da(:,jMode),Amat(:,jMode))
        CALL LA_fbase_in%projectIPtoDOF(.TRUE. ,-PhiPrime_s,DERIV_THET,gam_za_da(:,jMode),Amat(:,jMode))
        Amat(:,jMode) = Amat(:,jMode) *sAdiag(:)
      ELSE
        Amat(:    ,jMode)=0.0_wp
        Amat(jMode,jMode)=1.0_wp
      END IF
    END DO!jMode
!$OMP END PARALLEL DO

  __PERFOFF('setup_Amat')
  __PERFON('setup_rhs')

  !$OMP PARALLEL DO        &
  !$OMP   SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(iMode,sum_gam_ta_da,sum_gam_za_da)        &
  !$OMP   SHARED(RHS,gam_ta_da,gam_za_da,LA_fbase_in,LA_modes,chiPrime_s,PhiPrime_s,sAdiag)
    DO iMode=1,LA_modes
      !m=n=0 should not be in lambda, but check
      IF (LA_fbase_in%zero_odd_even(iMode).NE.MN_ZERO) THEN
        sum_gam_ta_da=SUM(gam_ta_da(:,iMode))
        sum_gam_za_da=SUM(gam_za_da(:,iMode))
        ! 1/J( iota (g_thet,zeta dsigma_dthet - g_thet,thet dsigma_dzeta )
        !          +(g_zeta,zeta dsigma_dthet - g_zeta,thet dsigma_dzeta ) )
        RHS(iMode) = (chiPrime_s*sum_gam_ta_da +phiPrime_s*sum_gam_za_da) *sAdiag(iMode)
      ELSE
        RHS(iMode) = 0.0_wp
      END IF
    END DO!iMode
  !$OMP END PARALLEL DO

    __PERFOFF('setup_rhs')
    __PERFON('solve')
    LA_s=SOLVE(Amat,RHS)
    __PERFOFF('solve')
END SUBROUTINE Lambda_setup_and_solve

END MODULE MODgvec_lambda_solve
