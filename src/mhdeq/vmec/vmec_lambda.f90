!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>

! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
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

!===================================================================================================================================
!>
!!# Module **VMEC lambda**
!!
!! recompute lambda at the full flux surfaces via an elliptic solution on 
!! each flux surface
!! integration over angles done wih trapezoidal rule, cosine and sine are
!! sampled at theta/zeta = [0,...,k/(2pi),...,(np-1)/(2pi)]
!!
!===================================================================================================================================
MODULE MOD_VMEC_Lambda
  USE MOD_Globals,ONLY:wp
  IMPLICIT NONE
  PUBLIC

! ----------------------------------------------------------------------
INTERFACE RecomputeLambda
  MODULE PROCEDURE RecomputeLambda
END INTERFACE

CONTAINS


!=====================================================================
!> Recompute lambda, given the VMEC R,Z data and its recontructed splines 
!! lambda coefficients lmns and lmnc will be overwritten!
!!
!===================================================================================================================================
SUBROUTINE RecomputeLambda(np_m,np_n)
!=====================================================================
USE MOD_LinAlg,     ONLY: SOLVE
USE MOD_VMEC_vars,  ONLY: xmabs,rho,Rmnc_Spl,Rmns_Spl,Zmnc_Spl,Zmns_Spl
USE MOD_VMEC_Readin,ONLY: Rmnc,Rmns,Zmnc,Zmns,lmnc,lmns,iotaf,lasym,xm,xn
USE MOD_VMEC_Readin,ONLY: nfp,nFluxVMEC,mn_mode,twoPi
USE SPLINE1_MOD,    ONLY: SPLINE1_EVAL
  IMPLICIT NONE
!---------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER          :: np_m !< number of points for integration in theta 
  INTEGER          :: np_n !< number of points for integration in zeta
!---------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                   :: iFlux,iMode,jMode,mn0Mode
  INTEGER                                   :: i_m,i_n,i_mn,np_mn,iGuess
  REAL(wp),DIMENSION(mn_mode,nFluxVMEC)     :: dRmnc_drho,dRmns_drho,dZmnc_drho,dZmns_drho
  REAL(wp),DIMENSION(1:np_m*np_n,mn_mode)   :: cosMN,sinMN
  REAL(wp),DIMENSION(1:np_m*np_n)           :: R,dRdrho,dRdtheta,dRdzeta
  REAL(wp),DIMENSION(1:np_m*np_n)           :: dZdrho,dZdtheta,dZdzeta
  REAL(wp),DIMENSION(1:np_m*np_n)           :: J,g_tt,g_tz,g_zz 
  REAL(wp),ALLOCATABLE                      :: Amat(:,:),RHS(:),lambda(:),sAdiag(:)
  REAL(wp)                                  :: rho_p,rhom,drhom,splOut(3) !for weighted spline interpolation
!=====================================================================
  WRITE(*,'(4X,A)')'RECOMPUTE LAMBDA ON FULL GRID...'
  WRITE(*,'(4X,A,2I8)')' # int. points in m and n : ',np_m,np_n
  IF(lasym)THEN !cosine & sine
    ALLOCATE(Amat(2*mn_mode,2*mn_mode),RHS(2*mn_mode),lambda(2*mn_mode))
  ELSE !only sine
    ALLOCATE(Amat(mn_mode,mn_mode),RHS(mn_mode),lambda(mn_mode))
  END IF
  !estimate of 1/Adiag
  ALLOCATE(sAdiag(1:mn_mode)) !sAdiag( mn_mode+iMode)=sAdiag(iMode)
  DO iMode=1,mn_mode
    sAdiag(iMode)=1.0_wp/(MAX(1.0_wp,((REAL(nfp,wp)*xm(iMode))**2+(xn(iMode))**2) )*REAL(np_m*np_n,wp))
  END DO

  DO iMode=1,mn_mode
    IF((xm(iMode).EQ.0).AND.(xn(iMode).EQ.0)) mn0Mode=iMode
  END DO
  np_mn=np_m*np_n 
  DO i_n=1,np_n; DO i_m=1,np_m
    i_mn=1+(i_m-1) + (i_n-1)*np_m
    cosMN(i_mn,:)=COS(twoPi*( xm(:)*(REAL(i_m-1,wp)/REAL(np_m,wp)) & 
                             -xn(:)*(REAL(i_n-1,wp)/REAL(nfp*np_n,wp))))
    sinMN(i_mn,:)=SIN(twoPi*( xm(:)*(REAL(i_m-1,wp)/REAL(np_m,wp)) &
                             -xn(:)*(REAL(i_n-1,wp)/REAL(nfp*np_n,wp))))
  END DO; END DO !i_m,i_n
  cosMN(:,mn0Mode)=1.0_wp
  sinMN(:,mn0Mode)=0.0_wp

  !compute radial derivatives of each R,Z mode, using splines of R and Z in rho direction
  DO iFlux=2,nFluxVMEC
    rho_p=MIN(1.0_wp,rho(iFlux)) !rho=sqrt(psi_norm)
    DO iMode=1,mn_mode
      SELECT CASE(xmabs(iMode))
      CASE(0)
        rhom=1.0_wp
        drhom=0.0_wp
      CASE(1)
        rhom=rho_p
        drhom=1.0_wp
      CASE(2)
        rhom=rho_p*rho_p
        drhom=2.0_wp*rho_p
      CASE DEFAULT
        rhom=rho_p**xmabs(iMode)
        drhom=REAL(xmabs(iMode) , wp)*rho_p**(xmabs(iMode)-1)
      END SELECT
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmnc_Spl(:,:,iMode),iGuess,splout) 
      dRmnc_drho(iMode,iFlux)=(rhom*splout(2)+splout(1)*drhom)
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmns_Spl(:,:,iMode),iGuess,splout) 
      dZmns_drho(iMode,iFlux)=(rhom*splout(2)+splout(1)*drhom)
      IF(lasym)THEN
        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmns_Spl(:,:,iMode),iGuess,splout) 
        dRmns_drho(iMode,iFlux)=(rhom*splout(2)+splout(1)*drhom)
        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmnc_Spl(:,:,iMode),iGuess,splout) 
        dZmnc_drho(iMode,iFlux)= (rhom*splout(2)+splout(1)*drhom)
      END IF !lasym
    END DO!iMode
  END DO!iFlux

  !compute metrics on flux surface and Amat + RHS
  DO iFlux=2,nFluxVMEC
    WRITE(0,'(I6,A4,I6,A27,A1)',ADVANCE='NO')iFlux, ' of ',nFluxVMEC,' flux surfaces evaluated...',ACHAR(13)
    !sample derivatives on mesh
    IF(.not.lasym)THEN
      !$OMP PARALLEL 
      DO i_mn=1,np_mn
        R       (i_mn) = 0.0_wp
        dRdrho  (i_mn) = 0.0_wp
        dRdtheta(i_mn) = 0.0_wp
        dRdzeta (i_mn) = 0.0_wp
        dZdrho  (i_mn) = 0.0_wp
        dZdtheta(i_mn) = 0.0_wp
        dZdzeta (i_mn) = 0.0_wp
        DO iMode=1,mn_mode
          R       (i_mn) =       R(i_mn)           +Rmnc(iMode,iFlux) *cosMN(i_mn,iMode)
          dRdrho  (i_mn) =  dRdrho(i_mn)+     dRmnc_drho(iMode,iFlux) *cosMN(i_mn,iMode)
          dRdtheta(i_mn) =dRdtheta(i_mn)-(xm(iMode)*Rmnc(iMode,iFlux))*sinMN(i_mn,iMode)
          dRdzeta (i_mn) = dRdzeta(i_mn)+(xn(iMode)*Rmnc(iMode,iFlux))*sinMN(i_mn,iMode)
      
          dZdrho  (i_mn) =  dZdrho(i_mn)+     dZmns_drho(iMode,iFlux) *sinMN(i_mn,iMode)
          dZdtheta(i_mn) =dZdtheta(i_mn)+(xm(iMode)*Zmns(iMode,iFlux))*cosMN(i_mn,iMode)
          dZdzeta (i_mn) = dZdzeta(i_mn)-(xn(iMode)*Zmns(iMode,iFlux))*cosMN(i_mn,iMode)
        END DO!iMode
      END DO !i_mn
      !$OMP END PARALLEL
    ELSE !lasym
      !$OMP PARALLEL 
      DO i_mn=1,np_mn
        R       (i_mn) = 0.0_wp
        dRdrho  (i_mn) = 0.0_wp
        dRdtheta(i_mn) = 0.0_wp
        dRdzeta (i_mn) = 0.0_wp
        dZdrho  (i_mn) = 0.0_wp
        dZdtheta(i_mn) = 0.0_wp
        dZdzeta (i_mn) = 0.0_wp
        DO iMode=1,mn_mode
          R       (i_mn) =       R(i_mn)         +  Rmnc(iMode,iFlux) *cosMN(i_mn,iMode) &
                                                 +  Rmns(iMode,iFlux) *sinMN(i_mn,iMode)
          dRdrho  (i_mn) =  dRdrho(i_mn)    + dRmnc_drho(iMode,iFlux) *cosMN(i_mn,iMode) &
                                            + dRmns_drho(iMode,iFlux) *sinMN(i_mn,iMode)
          dRdtheta(i_mn) =dRdtheta(i_mn)-(xm(iMode)*Rmnc(iMode,iFlux))*sinMN(i_mn,iMode) &
                                        +(xm(iMode)*Rmns(iMode,iFlux))*cosMN(i_mn,iMode)
          dRdzeta (i_mn) = dRdzeta(i_mn)+(xn(iMode)*Rmnc(iMode,iFlux))*sinMN(i_mn,iMode) &
                                        -(xn(iMode)*Rmns(iMode,iFlux))*cosMN(i_mn,iMode)
      
          dZdrho  (i_mn) =  dZdrho(i_mn)    + dZmns_drho(iMode,iFlux) *sinMN(i_mn,iMode) &
                                            + dZmnc_drho(iMode,iFlux) *cosMN(i_mn,iMode) 
          dZdtheta(i_mn) =dZdtheta(i_mn)+(xm(iMode)*Zmns(iMode,iFlux))*cosMN(i_mn,iMode) &
                                        -(xm(iMode)*Zmnc(iMode,iFlux))*sinMN(i_mn,iMode)
          dZdzeta (i_mn) = dZdzeta(i_mn)-(xn(iMode)*Zmns(iMode,iFlux))*cosMN(i_mn,iMode) &
                                        +(xn(iMode)*Zmnc(iMode,iFlux))*sinMN(i_mn,iMode)
        END DO!iMode
      END DO !i_mn
      !$OMP END PARALLEL 
    END IF !lasym
    !J=sqrtG=R*(dRdtheta*dZds-dRds*dZdtheta) =  1/(2*rho_p)*R*(dRdtheta*dZdrho-dRdrho*dZdtheta) 
    !1/(2rho_p) is constant on flux, not needed
    !$OMP PARALLEL
    DO i_mn=1,np_mn
      J(i_mn)=R(i_mn)*(dRdtheta(i_mn)*dZdrho(i_mn)-dRdrho(i_mn)*dZdtheta(i_mn)) 
    END DO !i_mn
    !$OMP END PARALLEL
    IF(MINVAL(ABS(J)).LT.1.0e-12) STOP 'Jacobian near zero! J< 1.0e-12'
    !account for 1/J here
    !$OMP PARALLEL
    DO i_mn=1,np_mn
      g_tt(i_mn) =(dRdtheta(i_mn)*dRdtheta(i_mn) + dZdtheta(i_mn)*dZdtheta(i_mn) )/J(i_mn)
      g_tz(i_mn) =(dRdtheta(i_mn)*dRdzeta (i_mn) + dZdtheta(i_mn)*dZdzeta (i_mn) )/J(i_mn) !=g_zt
      g_zz(i_mn) =(dRdzeta (i_mn)*dRdzeta (i_mn) + dZdzeta (i_mn)*dZdzeta (i_mn) + R(i_mn)*R(i_mn))/J(i_mn)
    END DO !i_mn
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    DO jMode=1,mn_mode
      DO iMode=1,mn_mode
        Amat(iMode,jMode)=0.0_wp
        DO i_mn=1,np_mn
          ! 1/J ( (g_thet,zeta  dlambdaSIN_dzeta -g_zeta,zeta  dlambdaSIN_dthet) dsigmaSIN_dthet
          !      -(g_thet,thet  dlambdaSIN_dzeta -g_zeta,thet  dlambdaSIN_dthet) dsigmaSIN_dzeta)
          Amat(iMode,jMode) = Amat(iMode,jMode) +&
!                            ( ( g_tz(i_mn)*(-xn(iMode)) -g_zz(i_mn)*( xm(iMode)))*( xm(jMode))    &  
!                             -( g_tt(i_mn)*(-xn(iMode)) -g_tz(i_mn)*( xm(iMode)))*(-xn(jMode)) )  &  
!optimized
                            (  g_tz(i_mn)*(((-xn(iMode))*( xm(jMode))) +            (( xm(iMode))*(-xn(jMode))))   &  
                             - g_tt(i_mn)* ((-xn(iMode))*(-xn(jMode))) - g_zz(i_mn)*(( xm(iMode))*( xm(jMode))) )  &  
                                     *cosMN(i_mn,iMode)*cosMN(i_mn,jMode)*sAdiag(iMode)
        END DO !i_mn
      END DO !iMode
      RHS(jMode)=0.0_wp
      DO i_mn=1,np_mn
        ! 1/J( (iota g_thet,zeta  + g_zeta,zeta )  dsigmaSIN_dthet
        !     -(iota g_thet,thet  + g_zeta,thet )  dsigmaSIN_dzeta
        RHS(jMode)      =   RHS(jMode)+&
!                           ( ( g_tz(i_mn)*iotaf(iFlux) +g_zz(i_mn)  )*( xm(jMode))                &  
!                            -( g_tt(i_mn)*iotaf(iFlux) +g_tz(i_mn)  )*(-xn(jMode)))               &  
!optimized
                            ( g_tz(i_mn)*(iotaf(iFlux)*( xm(jMode))  -            (-xn(jMode)))   &  
                             -g_tt(i_mn)*(iotaf(iFlux)*(-xn(jMode))) + g_zz(i_mn)*( xm(jMode)) )  & 
                                                *cosMN(i_mn,jMode) *sAdiag(jMode)
      END DO !i_mn
    END DO!jMode
    !$OMP END PARALLEL
    Amat(mn0mode,:)      =0.0_wp
    Amat(mn0mode,mn0mode)=1.0_wp
    RHS(         mn0mode)=0.0_wp
    IF(lasym)THEN
      !lambda(1:mn_mode)=lmns(:,iFlux)
      !lambda(mn_mode+1:2*mn_mode)=lmnc(:,iFlux)
      !$OMP PARALLEL
      DO jMode=1,mn_mode
        DO iMode=1,mn_mode
          Amat(mn_mode+iMode,jMode)=0.0_wp 
          Amat(iMode,mn_mode+jMode)=0.0_wp
          Amat(mn_mode+iMode,mn_mode+jMode)=0.0_wp
          DO i_mn=1,np_mn
            ! 1/J ( (g_thet,thet dlambdaCOS_dzeta-g_thet,zeta  dlambdaCOS_dthet) dsigmaSIN_dthet
            !      -(g_zeta,thet dlambdaCOS_dzeta-g_zeta,zeta  dlambdaCOS_dthet) dsigmaSIN_dzeta)
            Amat(mn_mode+iMode,jMode)= Amat(mn_mode+iMode,jMode)+&
!                            ( ( g_tz(i_mn)*( xn(iMode)) -g_zz(i_mn)*(-xm(iMode)))*( xm(jMode))    &  
!                             -( g_tt(i_mn)*( xn(iMode)) -g_tz(i_mn)*(-xm(iMode)))*(-xn(jMode)) )  &  
!optimized
                            (  g_tz(i_mn)*(( xn(iMode))*( xm(jMode))  +            (-xm(iMode))*(-xn(jMode)))    &  
                             - g_tt(i_mn)*(( xn(iMode))*(-xn(jMode))) -g_zz(i_mn)*((-xm(iMode))*( xm(jMode))) )  &  
                                    *sinMN(i_mn,iMode)*cosMN(i_mn,jMode) *sAdiag(iMode)

            ! 1/J ( (g_thet,thet dlambdaSIN_dzeta-g_thet,zeta  dlambdaSIN_dthet) dsigmaCOS_dthet
            !      -(g_zeta,thet dlambdaSIN_dzeta-g_zeta,zeta  dlambdaSIN_dthet) dsigmaCOS_dzeta)
            Amat(iMode,mn_mode+jMode)= Amat(iMode,mn_mode+jMode)+ &
!                            ( ( g_tz(i_mn)*(-xn(iMode)) -g_zz(i_mn)*( xm(iMode)))*(-xm(jMode))    &  
!                             -( g_tt(i_mn)*(-xn(iMode)) -g_tz(i_mn)*( xm(iMode)))*( xn(jMode)) )  &  
!optimized
                            (  g_tz(i_mn)*((-xn(iMode))*(-xm(jMode))  +            ( xm(iMode))*( xn(jMode)))    &  
                             - g_tt(i_mn)*((-xn(iMode))*( xn(jMode))) -g_zz(i_mn)*(( xm(iMode))*(-xm(jMode))))  &  
                                    *cosMN(i_mn,iMode)*sinMN(i_mn,jMode) *sAdiag(iMode)

            ! 1/J ( (g_thet,thet dlambdaCOS_dzeta-g_thet,zeta  dlambdaCOS_dthet) dsigmaCOS_dthet
            !      -(g_zeta,thet dlambdaCOS_dzeta-g_zeta,zeta  dlambdaCOS_dthet) dsigmaCOS_dzeta)
            Amat(mn_mode+iMode,mn_mode+jMode)= Amat(mn_mode+iMode,mn_mode+jMode)+&
!                            ( ( g_tz(i_mn)*( xn(iMode)) -g_zz(i_mn)*(-xm(iMode)))*(-xm(jMode))    &  
!                             -( g_tt(i_mn)*( xn(iMode)) -g_tz(i_mn)*(-xm(iMode)))*( xn(jMode)) )  &  
!optimized
                            (  g_tz(i_mn)*(( xn(iMode))*(-xm(jMode))  +            (-xm(iMode))*( xn(jMode)))    &  
                             - g_tt(i_mn)*(( xn(iMode))*( xn(jMode))) -g_zz(i_mn)*((-xm(iMode))*(-xm(jMode))) )  &  
                                    *sinMN(i_mn,iMode)*sinMN(i_mn,jMode) *sAdiag(iMode)
        
        
          END DO!i_mn
        END DO !iMode
        RHS(mn_mode+jMode) =0.0_wp 
        DO i_mn=1,np_mn
          ! 1/J( (iota g_thet,thet+ g_thet,zeta)  dsigmaCOS_dthet
          !     -(iota g_zeta,thet+ g_zeta,zeta)  dsigmaCOS_dzeta
          RHS(mn_mode+jMode) = RHS(mn_mode+jMode) + &
!                               ( ( g_tz(i_mn)*iotaf(iFlux) +g_zz(i_mn)  )*(-xm(jMode))               &  
!                                -( g_tt(i_mn)*iotaf(iFlux) +g_tz(i_mn)  )*( xn(jMode)))              &  
!optimized
                                (   g_tz(i_mn)*((iotaf(iFlux)*(-xm(jMode))) -           ( xn(jMode)))     &  
                                 -  g_tt(i_mn)* (iotaf(iFlux)*( xn(jMode))) +g_zz(i_mn)*(-xm(jMode))   )  &  
                                                       *sinMN(i_mn,jMode) *sAdiag(jMode)
       
        END DO !i_mn
      END DO!jMode
      !$OMP END PARALLEL
      Amat(mn_mode+mn0Mode,:              )=0.0_wp
      Amat(mn_mode+mn0Mode,mn_mode+mn0Mode)=1.0_wp
      RHS(                 mn_mode+mn0Mode)=0.0_wp
    END IF !lasym
   
    lambda=SOLVE(Amat,RHS)  
!    WRITE(*,'(A,I6,2E21.11)')'DEBUG',iFlux,MAXVAL(ABS(lambda(1:mn_mode)-lmns(:,iFlux)))
!    DO iMode=1,mn_mode
!      WRITE(*,'(A,I6,I6,3E21.11)')'DEBUG',NINT(xm(iMode)),NINT(xn(iMode)),lambda(iMode),lmns(iMode,iFlux),ABS(lambda(iMode)-lmns(iMode,iFlux))
!    END DO
!    WRITE(*,*)'***************************'

    !overwrite
    lmns(:,iFlux)=lambda(1:mn_mode)
    IF(lasym)THEN
!    WRITE(*,'(A,I6,2E21.11)')'DEBUG',iFlux,MAXVAL(ABS(lambda(mn_mode+1:2*mn_mode)-lmnc(:,iFlux)))
!    DO iMode=1,mn_mode
!      WRITE(*,'(A,I6,I6,3E21.11)')'DEBUG',NINT(xm(iMode)),NINT(xn(iMode)),lambda(mn_mode+iMode),lmnc(iMode,iFlux) &
!                            ,lambda(mn_mode+iMode)/(1.0e-13+lmnc(iMode,iFlux))
!    END DO
!    WRITE(*,*)'###########################'

      lmnc(:,iFlux)=lambda(mn_mode+1:2*mn_mode)
    END IF
  END DO !iFlux
  !on axis
  lmns(:,1)=0.0_wp
  IF(lasym)THEN
    lmnc(:,1)=0.0_wp
  END IF
!STOP
  DEALLOCATE(lambda,Amat,RHS,sAdiag)
  WRITE(*,'(A)')'                                                       '
  WRITE(*,'(4X,A)')'...DONE.'
END SUBROUTINE RecomputeLambda

END MODULE MOD_VMEC_lambda 
