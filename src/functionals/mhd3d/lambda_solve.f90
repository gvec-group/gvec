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
MODULE MOD_lambda_solve
! MODULES
USE MOD_Globals, ONLY:wp,UNIT_StdOut,abort
IMPLICIT NONE
PUBLIC

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Solve for lambda on one given flux surface
!!
!! Note that the mapping defined by  X1 and X2 must be fully initialized, since derivatives in s must be taken!
!!
!===================================================================================================================================
SUBROUTINE Lambda_solve(spos,iota_s,X1_in,X2_in,LA_s) 
! MODULES
USE MOD_sol_var_MHD3D, ONLY: t_sol_var_MHD3D
USE MOD_LinAlg,     ONLY: SOLVE
USE MOD_MHD3D_Vars, ONLY: hmap,X1_base,X2_base,LA_base
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp)     , INTENT(IN   ) :: spos                    !! s position to evaluate lambda
REAL(wp)     , INTENT(IN   ) :: iota_s                  !! iota at s_pos
REAL(wp)     , INTENT(IN   ) :: X1_in(1:X1_base%s%nBase,1:X1_base%f%modes) !! U%X1 variable, is reshaped to 2D at input
REAL(wp)     , INTENT(IN   ) :: X2_in(1:X2_base%s%nBase,1:X2_base%f%modes) !! U%X2 variable, is reshaped to 2D at input 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)     , INTENT(  OUT) :: LA_s(1:LA_base%f%modes) !! lambda at spos 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: iMode,jMode,i_mn,mn_IP,LA_modes
REAL(wp)                              :: Jh,minJ,gta_da,gza_da,qloc(3),dqdthet(3),dqdzeta(3)
REAL(wp),DIMENSION(1:X1_base%f%modes) :: X1_s,X1_ds !! X1 solution at spos 
REAL(wp),DIMENSION(1:X2_base%f%modes) :: X2_s,X2_ds !! X1 solution at spos 
REAL(wp),DIMENSION(1:X1_base%f%mn_IP) :: X1_s_IP,dX1ds,dX1dthet,dX1dzeta, & !mn_IP should be same for all!
                                         X2_s_IP,dX2ds,dX2dthet,dX2dzeta, &
                                         detJ,g_tt,g_tz,g_zz,zeta_IP
REAL(wp)                              :: Amat(1:LA_base%f%modes,1:LA_base%f%modes)
REAL(wp),DIMENSION(1:LA_base%f%modes) :: RHS,sAdiag
!===================================================================================================================================
  mn_IP = X1_base%f%mn_IP
  IF(X2_base%f%mn_IP.NE.mn_IP) STOP'X2 mn_IP /= X1 mn_IP'
  IF(LA_base%f%mn_IP.NE.mn_IP) STOP'LA mn_IP /= X1 mn_IP'
  zeta_IP  = X1_base%f%x_IP(2,:) 

  DO iMode=1,X1_base%f%modes
    X1_s( iMode)  = X1_base%s%evalDOF_s(spos,      0,X1_in(:,iMode))
    X1_ds(iMode)  = X1_base%s%evalDOF_s(spos,DERIV_S,X1_in(:,iMode))
  END DO
  DO iMode=1,X2_base%f%modes
    X2_s( iMode)  = X2_base%s%evalDOF_s(spos,      0,X2_in(:,iMode))
    X2_ds(iMode)  = X2_base%s%evalDOF_s(spos,DERIV_S,X2_in(:,iMode))
  END DO
  X1_s_IP  = X1_base%f%evalDOF_IP(         0,X1_s )
  dX1ds    = X1_base%f%evalDOF_IP(         0,X1_ds)
  dX1dthet = X1_base%f%evalDOF_IP(DERIV_THET,X1_s )
  dX1dzeta = X1_base%f%evalDOF_IP(DERIV_ZETA,X1_s )

  X2_s_IP  = X2_base%f%evalDOF_IP(         0,X2_s )
  dX2ds    = X2_base%f%evalDOF_IP(         0,X2_ds)
  dX2dthet = X2_base%f%evalDOF_IP(DERIV_THET,X2_s )
  dX2dzeta = X2_base%f%evalDOF_IP(DERIV_ZETA,X2_s )


  minJ=1.0E+12_wp
  DO i_mn=1,mn_IP
    qloc(1:3) = (/X1_s_IP(i_mn) , X2_s_IP(i_mn) , zeta_IP(i_mn)/)
    Jh=hmap%eval_Jh( qloc ) !X1,X2,zeta
    detJ(i_mn)=(dX1ds(i_mn)*dX2dthet(i_mn)-dX1dthet(i_mn)*dX2ds(i_mn))*Jh !J_p*J_h
  END DO !i_mn
  IF(MINVAL(detJ) .LT.1.0e-12) THEN
    i_mn= MINLOC(detJ(:),1)
    WRITE(UNIT_stdOut,'(4X,4(A,E11.3))')'WARNING min(J)= ',MINVAL(detJ),' at s= ',spos, &
                                                                       ' theta= ',X1_base%f%x_IP(1,i_mn), &
                                                                        ' zeta= ',X1_base%f%x_IP(2,i_mn) 
    i_mn= MAXLOC(detJ(:),1)
    WRITE(UNIT_stdOut,'(4X,4(A,E11.3))')'     ...max(J)= ',MAXVAL(detJ),' at s= ',spos, &
                                                                       ' theta= ',X1_base%f%x_IP(1,i_mn), &
                                                                        ' zeta= ',X1_base%f%x_IP(2,i_mn) 
!    CALL abort(__STAMP__, &
!        'Lambda_solve: Jacobian smaller that  1.0e-12!!!' )
  END IF
  !account for 1/J here
  DO i_mn=1,mn_IP
    qloc   (1:3) =(/  X1_s_IP(i_mn), X2_s_IP(i_mn), zeta_IP(i_mn)/) 
    dqdthet(1:3) =(/ dX1dthet(i_mn),dX2dthet(i_mn), 0.0_wp       /) 
    dqdzeta(1:3) =(/ dX1dzeta(i_mn),dX2dzeta(i_mn), 1.0_wp       /) 

    g_tt(i_mn) = (hmap%eval_gij(dqdthet,qloc,dqdthet))/detJ(i_mn)
    g_tz(i_mn) = (hmap%eval_gij(dqdthet,qloc,dqdzeta))/detJ(i_mn)
    g_zz(i_mn) = (hmap%eval_gij(dqdzeta,qloc,dqdzeta))/detJ(i_mn)
  END DO !i_mn
  
  LA_modes=LA_base%f%modes
  !estimate of 1/Adiag
  DO iMode=1,LA_modes
    ASSOCIATE(nfp=> LA_base%f%nfp,m=>LA_base%f%Xmn(1,iMode),n=>LA_base%f%Xmn(2,iMode))
    sAdiag(iMode)=1.0_wp/(MAX(1.0_wp,REAL((nfp*m)**2+n**2 ,wp) )*REAL(mn_IP,wp))
    END ASSOCIATE
  END DO
  Amat(:,:)=0.0_wp
  RHS(:)   =0.0_wp
  ASSOCIATE(sigma_dthet => LA_base%f%base_dthet_IP, &
            sigma_dzeta => LA_base%f%base_dzeta_IP  )
  DO jMode=1,LA_modes
    !m=n=0 should not be in lambda, but check
    IF (LA_base%f%zero_odd_even(jMode).NE.MN_ZERO) THEN 
      DO i_mn=1,mn_IP
        gta_da=g_tz(i_mn)*sigma_dthet(i_mn,jMode) - g_tt(i_mn)*sigma_dzeta(i_mn,jMode)
        gza_da=g_zz(i_mn)*sigma_dthet(i_mn,jMode) - g_tz(i_mn)*sigma_dzeta(i_mn,jMode)
        DO iMode=1,LA_modes
          ! 1/J ( (g_thet,zeta dsigma_dthet -g_thet,thet dsigma_dzeta ) dlambdaSIN_dzeta
          !      -(g_zeta,zeta dsigma_dthet -g_zeta,thet dsigma_dzeta ) dlambdaSIN_dthet)
          Amat(iMode,jMode) = Amat(iMode,jMode) +&
                              ( gta_da*sigma_dzeta(i_mn,iMode) &
                               -gza_da*sigma_dthet(i_mn,iMode)) *sAdiag(iMode)
        END DO !iMode
        ! 1/J( iota (g_thet,zeta dsigma_dthet - g_thet,thet dsigma_dzeta )
        !          +(g_zeta,zeta dsigma_dthet - g_zeta,thet dsigma_dzeta ) )
        RHS(jMode)      =   RHS(jMode)+ (iota_s*gta_da +gza_da) *sAdiag(jMode)
      END DO !i_mn
    ELSE
      Amat(jMode,jMode)=1.0_wp
      RHS(       jMode)=0.0_wp
    END IF
  END DO!jMode
  END ASSOCIATE !sigma_dthet,sigma_dzeta
 
  LA_s=SOLVE(Amat,RHS)  

END SUBROUTINE Lambda_solve


END MODULE MOD_lambda_solve
