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
  REAL(wp),DIMENSION(1:np_m*np_n,mn_mode)   :: dcosMN_dthet,dcosMN_dzeta
  REAL(wp),DIMENSION(1:np_m*np_n,mn_mode)   :: dsinMN_dthet,dsinMN_dzeta
  REAL(wp),DIMENSION(1:np_m*np_n)           :: R,dRdrho,dRdthet,dRdzeta
  REAL(wp),DIMENSION(1:np_m*np_n)           :: dZdrho,dZdthet,dZdzeta
  REAL(wp),DIMENSION(1:np_m*np_n)           :: J,g_tt,g_tz,g_zz 
  REAL(wp),ALLOCATABLE                      :: Amat(:,:),RHS(:),lambda(:),sAdiag(:)
  REAL(wp)                                  :: dthet,dzeta
  REAL(wp)                                  :: rho_p,rhom,drhom,splOut(3) !for weighted spline interpolation
  REAL(wp)                                  :: gta_dsinda,gza_dsinda
  REAL(wp)                                  :: gta_dcosda,gza_dcosda
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
  dthet=twoPi/REAL(np_m,wp)     ![0;2pi]
  dzeta=twoPi/REAL(nfp*np_m,wp) ![0;2pi/nfp]

  DO i_n=1,np_n; DO i_m=1,np_m
    i_mn=1+(i_m-1) + (i_n-1)*np_m
    cosMN(i_mn,:)=COS( xm(:)*(dthet*(REAL(i_m,wp)-0.5_wp)) & 
                      -xn(:)*(dzeta*(REAL(i_n,wp)-0.5_wp)))
    sinMN(i_mn,:)=SIN( xm(:)*(dthet*(REAL(i_m,wp)-0.5_wp)) &
                      -xn(:)*(dzeta*(REAL(i_n,wp)-0.5_wp)))
  END DO; END DO !i_m,i_n
  cosMN(:,mn0Mode)=1.0_wp
  sinMN(:,mn0Mode)=0.0_wp

  DO iMode=1,mn_mode
    DO i_mn=1,np_mn
      dcosMN_dthet(i_mn,iMode) = -xm(iMode)*sinMN(i_mn,iMode)
      dcosMN_dzeta(i_mn,iMode) =  xn(iMode)*sinMN(i_mn,iMode)
      dsinMN_dthet(i_mn,iMode) =  xm(iMode)*cosMN(i_mn,iMode)
      dsinMN_dzeta(i_mn,iMode) = -xn(iMode)*cosMN(i_mn,iMode)
    END DO !i_mn
  END DO !iMode

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
    DO i_mn=1,np_mn
      R      (i_mn) = 0.0_wp
      dRdrho (i_mn) = 0.0_wp
      dRdthet(i_mn) = 0.0_wp
      dRdzeta(i_mn) = 0.0_wp
      dZdrho (i_mn) = 0.0_wp
      dZdthet(i_mn) = 0.0_wp
      dZdzeta(i_mn) = 0.0_wp
      DO iMode=1,mn_mode
        R      (i_mn) =      R(i_mn)      +Rmnc(iMode,iFlux)* cosMN(      i_mn,iMode)
        dRdrho (i_mn) = dRdrho(i_mn)+dRmnc_drho(iMode,iFlux)* cosMN(      i_mn,iMode)
        dRdthet(i_mn) =dRdthet(i_mn)+      Rmnc(iMode,iFlux)*dcosMN_dthet(i_mn,iMode)
        dRdzeta(i_mn) =dRdzeta(i_mn)+      Rmnc(iMode,iFlux)*dcosMN_dzeta(i_mn,iMode)
    
        dZdrho (i_mn) = dZdrho(i_mn)+ dZmns_drho(iMode,iFlux)*sinMN(       i_mn,iMode)
        dZdthet(i_mn) =dZdthet(i_mn)+       Zmns(iMode,iFlux)*dsinMN_dthet(i_mn,iMode)
        dZdzeta(i_mn) =dZdzeta(i_mn)+       Zmns(iMode,iFlux)*dsinMN_dzeta(i_mn,iMode)
      END DO!iMode
    END DO !i_mn
    IF(lasym)THEN
      DO i_mn=1,np_mn
        DO iMode=1,mn_mode
          R      (i_mn) =      R(i_mn) +       Rmns(iMode,iFlux)* sinMN(      i_mn,iMode)
          dRdrho (i_mn) = dRdrho(i_mn) + dRmns_drho(iMode,iFlux)* sinMN(      i_mn,iMode)
          dRdthet(i_mn) =dRdthet(i_mn) +       Rmns(iMode,iFlux)*dsinMN_dthet(i_mn,iMode)
          dRdzeta(i_mn) =dRdzeta(i_mn) +       Rmns(iMode,iFlux)*dsinMN_dzeta(i_mn,iMode)
      
          dZdrho (i_mn) = dZdrho(i_mn) + dZmnc_drho(iMode,iFlux)* cosMN(      i_mn,iMode) 
          dZdthet(i_mn) =dZdthet(i_mn) +       Zmnc(iMode,iFlux)*dcosMN_dthet(i_mn,iMode)
          dZdzeta(i_mn) =dZdzeta(i_mn) +       Zmnc(iMode,iFlux)*dcosMN_dzeta(i_mn,iMode)
        END DO!iMode
      END DO !i_mn
    END IF !lasym
    !J=sqrtG=R*(dRdthet*dZds-dRds*dZdthet) =  1/(2*rho_p)*R*(dRdthet*dZdrho-dRdrho*dZdthet) 
    !1/(2rho_p) is constant on flux, not needed
    DO i_mn=1,np_mn
      J(i_mn)=R(i_mn)*(dRdthet(i_mn)*dZdrho(i_mn)-dRdrho(i_mn)*dZdthet(i_mn)) 
    END DO !i_mn
    IF(MINVAL(ABS(J)).LT.1.0e-12) STOP 'Jacobian near zero! J< 1.0e-12'
    !account for 1/J here
    DO i_mn=1,np_mn
      g_tt(i_mn) =(dRdthet(i_mn)*dRdthet(i_mn) + dZdthet(i_mn)*dZdthet(i_mn) )/J(i_mn)
      g_tz(i_mn) =(dRdthet(i_mn)*dRdzeta(i_mn) + dZdthet(i_mn)*dZdzeta(i_mn) )/J(i_mn) !=g_zt
      g_zz(i_mn) =(dRdzeta(i_mn)*dRdzeta(i_mn) + dZdzeta(i_mn)*dZdzeta(i_mn) + R(i_mn)*R(i_mn))/J(i_mn)
    END DO !i_mn
    
    Amat(:,:)=0.0_wp
    RHS(:)   =0.0_wp
    DO jMode=1,mn_mode
      DO i_mn=1,np_mn
        gta_dsinda=g_tz(i_mn)*dsinMN_dthet(i_mn,jMode) - g_tt(i_mn)*dsinMN_dzeta(i_mn,jMode)
        gza_dsinda=g_zz(i_mn)*dsinMN_dthet(i_mn,jMode) - g_tz(i_mn)*dsinMN_dzeta(i_mn,jMode)
        DO iMode=1,mn_mode
          ! 1/J ( (g_thet,zeta dsigmaSIN_dthet -g_thet,thet dsigmaSIN_dzeta ) dlambdaSIN_dzeta
          !      -(g_zeta,zeta dsigmaSIN_dthet -g_zeta,thet dsigmaSIN_dzeta ) dlambdaSIN_dthet)
          Amat(iMode,jMode) = Amat(iMode,jMode) +&
                              ( gta_dsinda*dsinMN_dzeta(i_mn,iMode) &
                               -gza_dsinda*dsinMN_dthet(i_mn,iMode)) *sAdiag(iMode)
        END DO !iMode
        ! 1/J( iota (g_thet,zeta dsigmaSIN_dthet - g_thet,thet dsigmaSIN_dzeta )
        !          +(g_zeta,zeta dsigmaSIN_dthet - g_zeta,thet dsigmaSIN_dzeta ) )
        RHS(jMode)      =   RHS(jMode)+ (iotaf(iFlux)*gta_dsinda +gza_dsinda) *sAdiag(jMode)
      END DO !i_mn
    END DO!jMode
    !m,n=0 Mode !=0
    Amat(mn0mode,:)      =0.0_wp
    Amat(mn0mode,mn0mode)=1.0_wp
    RHS(         mn0mode)=0.0_wp
    IF(lasym)THEN
      DO jMode=1,mn_mode
        DO i_mn=1,np_mn
          gta_dsinda=g_tz(i_mn)*dsinMN_dthet(i_mn,jMode) - g_tt(i_mn)*dsinMN_dzeta(i_mn,jMode)
          gza_dsinda=g_zz(i_mn)*dsinMN_dthet(i_mn,jMode) - g_tz(i_mn)*dsinMN_dzeta(i_mn,jMode)
          gta_dcosda=g_tz(i_mn)*dcosMN_dthet(i_mn,jMode) - g_tt(i_mn)*dcosMN_dzeta(i_mn,jMode)
          gza_dcosda=g_zz(i_mn)*dcosMN_dthet(i_mn,jMode) - g_tz(i_mn)*dcosMN_dzeta(i_mn,jMode)
          DO iMode=1,mn_mode
            ! 1/J ( (g_thet,thet dsigmaSIN_dthet -g_zeta,thet dsigmaSIN_dzeta ) dlambdaCOS_dzeta
            !      -(g_thet,zeta dsigmaSIN_dthet -g_zeta,zeta dsigmaSIN_dzeta ) dlambdaCOS_dthet)
            Amat(mn_mode+iMode,jMode)= Amat(mn_mode+iMode,jMode)+&
                                      ( gta_dsinda*dcosMN_dzeta(i_mn,iMode)  & 
                                       -gza_dsinda*dcosMN_dthet(i_mn,iMode)) *sAdiag(iMode)

            ! 1/J ( (g_thet,thet dsigmaCOS_dthet-g_zeta,thet dsigmaCOS_dzeta) dlambdaSIN_dzeta
            !      -(g_thet,zeta dsigmaCOS_dthet-g_zeta,zeta dsigmaCOS_dzeta) dlambdaSIN_dthet)
            Amat(iMode,mn_mode+jMode)= Amat(iMode,mn_mode+jMode)+ &
                                     ( gta_dcosda*dsinMN_dzeta(i_mn,iMode) &
                                      -gza_dcosda*dsinMN_dthet(i_mn,iMode)) *sAdiag(iMode)

            ! 1/J ( (g_thet,thet dsigmaCOS_dthet-g_zeta,thet dsigmaCOS_dzeta) dlambdaCOS_dzeta
            !      -(g_thet,zeta dsigmaCOS_dthet-g_zeta,zeta dsigmaCOS_dzeta) dlambdaCOS_dthet)
            Amat(mn_mode+iMode,mn_mode+jMode)= Amat(mn_mode+iMode,mn_mode+jMode)+&
                                     ( gta_dcosda*dcosMN_dzeta(i_mn,iMode)       & 
                                      -gza_dcosda*dcosMN_dthet(i_mn,iMode)) *sAdiag(iMode)
          END DO !iMode
          ! 1/J( iota (g_thet,thet dsigmaCOS_dthet - g_zeta,thet dsigmaCOS_dzeta)
          !          + g_thet,zeta dsigmaCOS_dthet - g_zeta,zeta dsigmaCOS_dzeta)
          RHS(mn_mode+jMode) = RHS(mn_mode+jMode) + (iotaf(iFlux)*gza_dcosda + gta_dcosda)*sAdiag(jMode)
        END DO !i_mn
      END DO!jMode
      !m,n=0 Mode !=0 
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

  DEALLOCATE(lambda,Amat,RHS,sAdiag)
  WRITE(*,'(A)')'                                                       '
  WRITE(*,'(4X,A)')'...DONE.'
END SUBROUTINE RecomputeLambda

END MODULE MOD_VMEC_lambda 
