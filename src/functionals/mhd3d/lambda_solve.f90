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
SUBROUTINE Lambda_solve(spos_in,hmap_in,X1_base_in,X2_base_in,LA_fbase_in,X1_in,X2_in,LA_s,phiPrime_s,chiPrime_s)
! MODULES
  USE MODgvec_Globals,       ONLY:n_warnings_occured
  USE MODgvec_base          ,ONLY: t_base
  USE MODgvec_fbase         ,ONLY: t_fbase
  USE MODgvec_hmap          ,ONLY: c_hmap
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base),INTENT(IN)  :: X1_base_in,X2_base_in           !< base classes belong to solution X1_in,X2_in
  CLASS(t_fbase),INTENT(IN) :: LA_fbase_in                     !< base class belong to solution LA_s
  CLASS(c_hmap), INTENT(INOUT) :: hmap_in                         
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
  REAL(wp)                              :: spos,Jh,minJ,qloc(3),dqdthet(3),dqdzeta(3)
  REAL(wp),DIMENSION(1:X1_base_in%f%modes) :: X1_s,X1_ds !! X1 solution at spos 
  REAL(wp),DIMENSION(1:X2_base_in%f%modes) :: X2_s,X2_ds !! X1 solution at spos 
  REAL(wp),DIMENSION(1:X1_base_in%f%mn_IP) :: X1_s_IP,dX1ds,dX1dthet,dX1dzeta, & !mn_IP should be same for all!
                                              X2_s_IP,dX2ds,dX2dthet,dX2dzeta, &
                                              detJ,gam_tt,gam_tz,gam_zz,zeta_IP
!===================================================================================================================================
  __PERFON('lambda_solve')
  spos=MIN(1.0_wp-1.0e-12_wp,MAX(1.0e-04,spos_in))
  mn_IP = X1_base_in%f%mn_IP
  IF(X2_base_in%f%mn_IP.NE.mn_IP) STOP 'X2 mn_IP /= X1 mn_IP'
  IF(LA_fbase_in%mn_IP .NE.mn_IP)  STOP 'LA mn_IP /= X1 mn_IP'
  zeta_IP  = X1_base_in%f%x_IP(2,:) 

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
!$OMP   PRIVATE(i_mn,qloc,Jh)  &
!$OMP   SHARED(mn_IP,X1_s_IP,X2_s_IP,zeta_IP,dX1ds,dX2ds,dX1dthet,dX2dthet,hmap_in,detJ)
  DO i_mn=1,mn_IP
    qloc(1:3) = (/X1_s_IP(i_mn) , X2_s_IP(i_mn) , zeta_IP(i_mn)/)
    Jh=hmap_in%eval_Jh( qloc ) !X1,X2,zeta
    detJ(i_mn)=(dX1ds(i_mn)*dX2dthet(i_mn)-dX1dthet(i_mn)*dX2ds(i_mn))*Jh !J_p*J_h
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
!$OMP   PRIVATE(i_mn,qloc,dqdthet,dqdzeta)  &
!$OMP   SHARED(mn_IP,X1_s_IP,X2_s_IP,zeta_IP,dX1dthet,dX2dthet,dX1dzeta,dX2dzeta,hmap_in,gam_tt,gam_tz,gam_zz,detJ)
  !account for 1/J here
  DO i_mn=1,mn_IP
    qloc   (1:3) =(/  X1_s_IP(i_mn), X2_s_IP(i_mn), zeta_IP(i_mn)/) 
    dqdthet(1:3) =(/ dX1dthet(i_mn),dX2dthet(i_mn), 0.0_wp       /) 
    dqdzeta(1:3) =(/ dX1dzeta(i_mn),dX2dzeta(i_mn), 1.0_wp       /) 

    gam_tt(i_mn) = (hmap_in%eval_gij(dqdthet,qloc,dqdthet))/detJ(i_mn)
    gam_tz(i_mn) = (hmap_in%eval_gij(dqdthet,qloc,dqdzeta))/detJ(i_mn)
    gam_zz(i_mn) = (hmap_in%eval_gij(dqdzeta,qloc,dqdzeta))/detJ(i_mn)
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
    CLASS(t_fbase),INTENT(IN)        :: LA_fbase_in           !< base classes belong to solution U_in
    REAL(wp),INTENT(IN)              :: phiPrime_s,ChiPrime_s   !! toroidal and poloidal flux s derivatives at s_pos
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP), INTENT(IN) :: gam_tt  !! g_tt/J evaluated on IP points
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP), INTENT(IN) :: gam_tz  !! g_tz/J evaluated on IP points
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP), INTENT(IN) :: gam_zz  !! g_zz/J evaluated on IP points
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)     , INTENT(  OUT) :: LA_s(1:LA_fbase_in%modes) !! lambda at spos 
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER                               :: iMode,jMode,i_mn,mn_IP,LA_modes
    REAL(wp)                              :: Amat(1:LA_fbase_in%modes,1:LA_fbase_in%modes)
    REAL(wp),DIMENSION(1:LA_fbase_in%modes) :: RHS,sAdiag
  !  REAL(wp)                              :: gam_ta_da,gam_za_da
    REAL(wp)                              :: sum_gam_ta_da,sum_gam_za_da,sum_mn
    REAL(wp),DIMENSION(1:LA_fbase_in%mn_IP,1:LA_fbase_in%modes) :: gam_ta_da,gam_za_da
  !===================================================================================================================================
    __PERFON('setup')
    LA_modes=LA_fbase_in%modes
    mn_IP   =LA_fbase_in%mn_IP

  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) & 
  !$OMP   DEFAULT(NONE)    &
  !$OMP   PRIVATE(iMode)  &
  !$OMP   SHARED(mn_IP,LA_modes,LA_fbase_in,sAdiag)
    !estimate of 1/Adiag for preconditioning (=1 for m=n=0)
    DO iMode=1,LA_modes
        ASSOCIATE(nfp=> LA_fbase_in%nfp,m=>LA_fbase_in%Xmn(1,iMode),n=>LA_fbase_in%Xmn(2,iMode))
        sAdiag(iMode)=1.0_wp/(MAX(1.0_wp,REAL((nfp*m)**2+n**2 ,wp) )*REAL(mn_IP,wp))
        !sAdiag(iMode)=1.0_wp 
        END ASSOCIATE
    END DO !iMode
  !$OMP END PARALLEL DO 

    ASSOCIATE(sigma_dthet => LA_fbase_in%base_dthet_IP, &
              sigma_dzeta => LA_fbase_in%base_dzeta_IP  )
  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) & 
  !$OMP   DEFAULT(NONE)    &
  !$OMP   PRIVATE(iMode,i_mn)        &
  !$OMP   SHARED(gam_ta_da,gam_za_da,mn_IP,LA_modes,gam_tt,gam_tz,gam_zz,LA_fbase_in)
    DO iMode=1,LA_modes
      DO i_mn=1,mn_IP
        gam_ta_da(i_mn,iMode)=gam_tz(i_mn)*sigma_dthet(i_mn,iMode) - gam_tt(i_mn)*sigma_dzeta(i_mn,iMode)
        gam_za_da(i_mn,iMode)=gam_zz(i_mn)*sigma_dthet(i_mn,iMode) - gam_tz(i_mn)*sigma_dzeta(i_mn,iMode)
      END DO !i_mn
    END DO!iMode
  !$OMP END PARALLEL DO 

  
  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) & 
  !$OMP   DEFAULT(NONE)    &
  !$OMP   PRIVATE(iMode,jMode,i_mn,sum_mn)        &
  !$OMP   SHARED(Amat,gam_ta_da,gam_za_da,LA_fbase_in,mn_IP,LA_modes,PhiPrime_s,sAdiag) 
    DO jMode=1,LA_modes
      !m=n=0 should not be in lambda, but check
      IF (LA_fbase_in%zero_odd_even(jMode).NE.MN_ZERO) THEN 
        DO iMode=1,LA_modes
          sum_mn=0.
  !$OMP SIMD REDUCTION(+:sum_mn)
          DO i_mn=1,mn_IP
            sum_mn = sum_mn + ( gam_ta_da(i_mn,jMode)*sigma_dzeta(i_mn,iMode) &
                               -gam_za_da(i_mn,jMode)*sigma_dthet(i_mn,iMode))
          END DO !i_mn
          ! 1/J ( (g_thet,zeta dsigma_dthet -g_thet,thet dsigma_dzeta ) dlambdaSIN_dzeta
          !      -(g_zeta,zeta dsigma_dthet -g_zeta,thet dsigma_dzeta ) dlambdaSIN_dthet)
          Amat(iMode,jMode) = sum_mn *PhiPrime_s*sAdiag(iMode)
        END DO !iMode
      ELSE
        Amat(:    ,jMode)=0.0_wp
        Amat(jMode,jMode)=1.0_wp
      END IF
    END DO!jMode
  !$OMP END PARALLEL DO 
  
    END ASSOCIATE !sigma_dthet,sigma_dzeta

  !$OMP PARALLEL DO        &  
  !$OMP   SCHEDULE(STATIC) & 
  !$OMP   DEFAULT(NONE)    &
  !$OMP   PRIVATE(iMode,i_mn,sum_gam_ta_da,sum_gam_za_da)        &
  !$OMP   SHARED(RHS,gam_ta_da,gam_za_da,LA_fbase_in,mn_IP,LA_modes,chiPrime_s,PhiPrime_s,sAdiag) 
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
  
    __PERFOFF('setup')
    __PERFON('solve')
    LA_s=SOLVE(Amat,RHS)  
    __PERFOFF('solve')
END SUBROUTINE Lambda_setup_and_solve

END MODULE MODgvec_lambda_solve
