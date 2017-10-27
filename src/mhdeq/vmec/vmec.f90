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
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **VMEC**
!!
!! Initializes variables to evaluate a VMEC dataset. 
!! In radial direction, a cubic spline is used to interpolate the data. 
!! Calls readin of VMEC "wout" datafile (netcdf format). 
!!
!===================================================================================================================================
MODULE MOD_VMEC
! MODULES
USE MOD_Globals,ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitVMEC 
  MODULE PROCEDURE InitVMEC 
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToVMEC 
!  MODULE PROCEDURE MapToVMEC 
!END INTERFACE

PUBLIC::InitVMEC
PUBLIC::MapToVMEC
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Initialize VMEC module
!!
!===================================================================================================================================
SUBROUTINE InitVMEC 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,abort
USE MOD_ReadInTools
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
USE MOD_VMEC_Vars
USE MOD_VMEC_lambda, ONLY:recomputeLambda
USE MOD_VMEC_Readin
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iMode,nyq,np_m,np_n
LOGICAL              :: useFilter,reLambda
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)')'  INIT VMEC INPUT ...'

  !VMEC "wout*.nc"  file
  VMECdataFile=GETSTR("VMECwoutfile")

  !use StraightFieldline mapping 
  useSFL=GETLOGICAL("VMEC_useSFL",".FALSE.")
 
  ! read VMEC 2000 output (netcdf)
  CALL ReadVmec(VMECdataFile)


  !gmnc not needed anymore
  !ALLOCATE(gmnc_half_nyq(1:nFluxVMEC,mn_mode_nyq))
  !gmnc_half_nyq     = gmnc
  ! half data is stored from 2:nFluxVMEC

  !toroidal flux from VMEC, now called PHI!
  ALLOCATE(Phi_prof(nFluxVMEC))
  Phi_prof = phi

  !normalized toroidal flux (=flux variable s [0;1] in VMEC)
  ALLOCATE(Phinorm_prof(nFluxVMEC))
  Phinorm_prof=(Phi_prof-Phi_prof(1))/(Phi_prof(nFluxVMEC)-Phi_prof(1))
  WRITE(UNIT_stdOut,'(4X,A,3F10.4)')'normalized toroidal flux of first three flux surfaces',Phinorm_prof(2:4)
  !poloidal flux from VMEC
  ALLOCATE(chi_prof(nFluxVMEC))
  chi_prof=chi
  WRITE(UNIT_stdOut,'(4X,A,3F10.4)')'min/max toroidal flux',MINVAL(phi),MAXVAL(phi)
  WRITE(UNIT_stdOut,'(4X,A,3F10.4)')'min/max poloidal flux',MINVAL(chi),MAXVAL(chi)

  WRITE(UNIT_stdOut,'(4X,A, I6)')'Total Number of mn-modes:',mn_mode
  WRITE(UNIT_stdOut,'(4X,A,3I6)')'Max Mode m,n,nfp: ',NINT(MAXVAL(xm)),NINT(MAXVAL(xn)),nfp
!  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes (Nyquist):',mn_mode_nyq
!  WRITE(UNIT_stdOut,*)'   Max Mode m,n: ',MAXVAL(xm_nyq),MAXVAL(xn_nyq)

  useFilter=.TRUE. !GETLOGICAL('VMECuseFilter','.TRUE.') !SHOULD BE ALWAYS TRUE...

  ALLOCATE(xmabs(mn_mode))
  DO iMode=1,mn_mode
    xmabs(iMode)=ABS(NINT(xm(iMode)))
    IF(useFilter)THEN
      IF(xmabs(iMode) > 3) THEN !Filtering for |m| > 3
        IF(MOD(xmabs(iMode),2) == 0) THEN
          xmabs(iMode)=2 !Even mode, remove rho**2
        ELSE
          xmabs(iMode)=3 !Odd mode, remove rho**3
        END IF
      END IF
    END IF !usefilter
  END DO !iMode=1,mn_mode

  !prepare Spline interpolation
  ALLOCATE(rho(1:nFluxVMEC))
  rho(:)=SQRT(Phinorm_prof(:))
  

  ALLOCATE(Rmnc_Spl(4,1:nFluxVMEC,mn_mode)) !first dim is for spline interpolation
  CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Rmnc,Rmnc_Spl)

  ALLOCATE(Zmns_Spl(4,1:nFluxVMEC,mn_mode))
  CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Zmns,Zmns_Spl)

  IF(lasym)THEN
    WRITE(*,'(4X,A)')'LASYM=TRUE : R,Z,lambda in cos and sin!'
    ALLOCATE(Rmns_Spl(4,1:nFluxVMEC,mn_mode)) 
    CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Rmns,Rmns_Spl)
    
    ALLOCATE(Zmnc_Spl(4,1:nFluxVMEC,mn_mode))
    CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Zmnc,Zmnc_Spl)
    
  END IF

  relambda=GETLOGICAL("VMEC_relambda",".TRUE.")
  IF(relambda) nyq=GETINT("VMEC_Lam_nyq","4")

  ALLOCATE(lmns_Spl(4,1:nFluxVMEC,mn_mode))
  IF(lasym) ALLOCATE(lmnc_Spl(4,1:nFluxVMEC,mn_mode))
  IF(reLambda)THEN
    np_m=1+nyq*2*(NINT(MAXVAL(ABS(xm)))/2) !m_points [0,2pi]
    np_n=1+nyq*2*(NINT(MAXVAL(ABS(xn)))/(2*nfp)) !n_points [0,2pi/nfs] ->  
    !recompute lambda on FULL GRID
    CALL RecomputeLambda(np_m,np_n) 
    CALL           FitSpline(mn_mode,nFluxVMEC,xmAbs,lmns,lmns_Spl)
    IF(lasym) CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,lmnc,lmnc_Spl)
  ELSE
    !lambda given on half grid
    CALL           FitSplineHalf(mn_mode,nFluxVMEC,xmAbs,lmns,lmns_Spl)
    IF(lasym) CALL FitSplineHalf(mn_mode,nFluxVMEC,xmAbs,lmnc,lmnc_Spl)
  END IF

  ALLOCATE(pres_spl(4,1:nFluxVMEC))
  pres_spl(1,:)=presf(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,pres_Spl(:,:), K_BC1=3, K_BCN=0)

  ALLOCATE(Phi_spl(4,1:nFluxVMEC))
  Phi_spl(1,:)=Phi_Prof(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,Phi_Spl(:,:), K_BC1=3, K_BCN=0)

  ALLOCATE(chi_spl(4,1:nFluxVMEC))
  chi_spl(1,:)=chi_Prof(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,chi_Spl(:,:), K_BC1=3, K_BCN=0)

  WRITE(*,'(4X,A,3F10.4)')'iota axis/middle/edge',iotaf(1),iotaf(nFluxVMEC/2),iotaf(nFluxVMEC)
  

  WRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE InitVMEC


!===================================================================================================================================
!> Fit disrete data along flux surfaces as spline for each fourier mode
!!
!===================================================================================================================================
SUBROUTINE FitSpline(modes,nFlux,mAbs,Xmn,Xmn_Spl)
! MODULES
USE MOD_Globals
USE MOD_VMEC_Vars,     ONLY: rho
USE SPLINE1_MOD,       ONLY: SPLINE1_FIT 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)  :: modes             !! number of modes
INTEGER, INTENT(IN)  :: nFlux             !! number of flux surfaces
INTEGER, INTENT(IN)  :: mabs(modes)       !! filtered m-mode value
REAL(wp), INTENT(IN) :: Xmn(modes,nFlux)  !! fourier coefficients at all flux surfaces 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp), INTENT(OUT):: Xmn_Spl(4,nFlux,modes)  !!  spline fitted fourier coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
!===================================================================================================================================
Xmn_Spl=0.0_wp
DO iMode=1,modes
  !scaling with rho^|m|
  DO iFlux=2,nFlux
    IF(mabs(iMode).EQ.0)THEN
      Xmn_Spl(1,iFlux,iMode)=Xmn(iMode,iFlux)
    ELSE
      Xmn_Spl(1,iFlux,iMode)=Xmn(iMode,iFlux) /(rho(iFlux)**mabs(iMode))
    END IF
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
  Xmn_Spl(1,1,iMode)=(Xmn_Spl(1,2,iMode)*rho(3)**2-Xmn_Spl(1,3,iMode)*rho(2)**2) /(rho(3)**2-rho(2)**2)
!  !Quadratic extrapolation to axis with dx'(rho=0)=0.0_wp
!  r1=rho(2)**2*rho(3)**4-rho(2)**4*rho(3)**2
!  r2=rho(2)**2*rho(4)**4-rho(2)**4*rho(4)**2
!  Xmn_Spl(1,1,iMode)= ( r1*(Xmn_Spl(1,2,iMode)*rho(4)**4-Xmn_Spl(1,4,iMode)*rho(2)**4) &
!                       -r2*(Xmn_Spl(1,2,iMode)*rho(3)**4-Xmn_Spl(1,3,iMode)*rho(2)**4)) &
!                     /( r1*(rho(4)**4-rho(2)**4)-r2*(rho(3)**4-rho(2)**4))
  CALL SPLINE1_FIT(nFlux,rho,Xmn_Spl(:,:,iMode), K_BC1=3, K_BCN=0)
END DO !iMode 

END SUBROUTINE FitSpline


!===================================================================================================================================
!> Fit disrete data along flux surfaces as spline for each fourier mode
!! input is given on the half mesh 2:nFluxVMEC 
!!
!===================================================================================================================================
SUBROUTINE FitSplineHalf(modes,nFlux,mabs,Xmn_half,Xmn_Spl)
! MODULES
USE MOD_Globals
USE MOD_VMEC_Vars,     ONLY: rho,Phinorm_prof
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
USE SPLINE1_MOD,       ONLY:SPLINE1_INTERP 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: modes                  !! number of modes
INTEGER, INTENT(IN) :: nFlux                  !! number of flux surfaces
INTEGER, INTENT(IN) :: mabs(modes)            !! filtered m-mode value
REAL(wp),INTENT(IN) :: Xmn_half(modes,nFlux)  !! fourier coefficients at all flux surfaces 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT):: Xmn_Spl(4,nFlux,modes) !!  spline fitted fourier coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
REAL(wp)          :: Xmn_half_Spl(4,nFlux+1)  ! spline fitted fourier coefficients 
REAL(wp)          :: rho_half(1:nFlux+1)
INTEGER           :: iFlag
CHARACTER(len=100):: message
!===================================================================================================================================
DO iFlux=1,nFlux-1
  rho_half(iFlux+1)=SQRT(0.5_wp*(Phinorm_prof(iFlux+1)+Phinorm_prof(iFlux))) !0.5*(rho(iFlux)+rho(iFlux+1))
END DO
!add end points
rho_half(1)=0.0_wp
rho_half(nFlux+1)=1.0_wp

DO iMode=1,modes
  !scaling with rho^|m|
  DO iFlux=2,nFlux
    IF(mabs(iMode).EQ.0)THEN
      Xmn_half_Spl(1,iFlux)=Xmn_half(iMode,iFlux)
    ELSE
      Xmn_half_Spl(1,iFlux)=Xmn_half(iMode,iFlux) /(rho_half(iFlux)**mabs(iMode))
    END IF
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
  Xmn_Half_Spl(1,1)=(Xmn_Half_Spl(1,2)*rho_half(3)**2-Xmn_Half_Spl(1,3)*rho_half(2)**2) /(rho_half(3)**2-rho_half(2)**2)
  !Extrapolate to Edge 
  Xmn_Half_Spl(1,nFlux+1)= ( Xmn_half_Spl(1,nFlux  )*(rho_half(nFlux+1)-rho_half(nFlux-1))     &
                            -Xmn_half_Spl(1,nFlux-1)*(rho_half(nFlux+1)-rho_half(nFlux  )) )   &
                               /(  rho_half(nFlux)   -rho_half(nFlux-1) )
  CALL SPLINE1_FIT(nFlux+1,rho_half,Xmn_half_Spl(:,:), K_BC1=3, K_BCN=0)
  iflag=0
  message=''
  CALL SPLINE1_INTERP((/1,0,0/),nFlux+1,rho_half,Xmn_half_Spl, &
                                nFlux  ,rho     ,Xmn_Spl(:,:,iMode),       &
                          iflag,message, K_BC1=3,K_BCN=0)
  !respline
  Xmn_Spl(2:4,:,iMode)=0.0_wp
  Xmn_Spl(1,1,iMode)  =(Xmn_Spl(1,2,iMode)*rho(3)**2-Xmn_Spl(1,3,iMode)*rho(2)**2) /(rho(3)**2-rho(2)**2)
  CALL SPLINE1_FIT(nFlux,rho,Xmn_Spl(:,:,iMode), K_BC1=3, K_BCN=0)
END DO !iMode 

END SUBROUTINE FitSplineHalf



!===================================================================================================================================
!> Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from VMEC data. 
!! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!!
!===================================================================================================================================
SUBROUTINE MapToVMEC(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars,    ONLY: nVarMHDEQ
USE MOD_MHDEQ_Vars,    ONLY: nRhoCoefs,RhoFluxVar,RhoCoefs
USE MOD_MHDEQ_Tools,   ONLY: Eval1DPoly
USE MOD_VMEC_Vars
USE MOD_VMEC_Readin,   ONLY: mu0,lasym,mn_mode,xm,xn,nFluxVMEC
USE MOD_Newton,        ONLY: NewtonRoot1D_FdF
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nTotal          !! total number of points
REAL(wp),INTENT(IN):: x_in(3,nTotal)  !! input coordinates represent a cylinder: 
INTEGER, INTENT(IN):: InputCoordSys   !!  0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
                                      !! 1: x_in(1:3) are (r,z,phi) coordinates r= [0;1], z= [0;1], phi=[0;1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT):: x_out(3,nTotal) !! mapped x,y,z coordinates with vmec data
REAL(wp),INTENT(OUT):: MHDEQdata(nVarMHDEQ,nTotal) !! vector of equilibrium variables, see definition in mhdeq_vars.f90
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,percent
INTEGER :: iMode
INTEGER :: iGuess=1
REAL(wp):: CosMN(mn_mode)          ! =cos(m*theta-n*zeta) for all modes
REAL(wp):: SinMN(mn_mode)          ! =sin(m*theta-n*zeta) for all modes
REAL(wp):: r_p                     ! radius in cylindrical coordinate system
REAL(wp):: Phinorm                 ! normalized TOROIDAL flux (=flux coordinate s [0,1]), use map r_p ~ sqrt(Phinorm)
                                   ! NOTE: Phi is called PHI in VMEC
REAL(wp):: chinorm                 ! normalized POLOIDAL flux [0,1]
REAL(wp):: theta                   ! poloidal angle [0,2pi]
REAL(wp):: zeta                    ! toroidal angle [0,2pi]
REAL(wp):: rho_p,rhom,drhom,splOut(3) !for weighted spline interpolation
REAL(wp):: dRdrho,dRdtheta,dRdzeta ! derivatives of R,Z, from splines
REAL(wp):: dZdrho,dZdtheta,dZdzeta ! 
REAL(wp):: lam,dldtheta,dldzeta    ! spline interpolation of lambda function and derivatives
REAL(wp):: sqrtGr                  ! (part of) mapping Jacobian between VMEC and cylinder coords. 
                                   ! (R,Z,phi) <->(rho,theta,zeta)
REAL(wp):: Phi_int,dPhi_drho_int   ! toroidal flux and derivative from interpolation at rho_p 
                                   ! NOTE: Phi is called PHI in VMEC
REAL(wp):: chi_int,dchi_drho_int   ! poloidal flux and derivative from interpolation at rho_p
REAL(wp):: iota_int                ! rotational transform, >0, iota=-chi'/Phi'
REAL(wp):: Br,Bz,Btor              ! mangetic field components in (R,Z,phi) system (phi=R*zeta)
                                   ! Br=dRdtheta*Btheta+dRdzeta*Bszeta 
                                   ! Bz=dZdtheta*Btheta+dZdzeta*Bzeta
                                   ! Bphi=R*Bzeta
REAL(wp):: Bcart(3)                ! magnetic field components in (X,Y,Z) system 
                                   ! Bx=Br*cos(torangle) - Btor*sin(torangle)
                                   ! By=Br*sin(torangle) + Btor*cos(torangle)
                                   ! Bz=Bz
REAL(wp):: Ar,Az,Ator,Acart(3)     ! R,Z,torangle and cartesian components of the magnetic vector potential
REAL(wp):: Density,Pressure        ! from density and pressure profiles
REAL(wp):: coszeta,sinzeta         ! cos(zeta),sin(zeta)
REAL(wp):: R,Z                     !
REAL(wp):: theta_star              ! for SFL
REAL(wp):: lam_rho(mn_mode,2)      ! for SFL
REAL(wp):: dldtheta_rho(mn_mode,2) ! for SFL
REAL(wp):: Rmin,Rmax,Zmin,Zmax 
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO VMEC DATA FROM ',TRIM(VMECdataFile),' ...'
Rmin=1.0E+12; Rmax=-1.0E12
Zmin=1.0E+12; Zmax=-1.0E12
percent=0
DO iNode=1,nTotal
  ! output of progress in %
  IF((nTotal.GT.10000).AND.(MOD(iNode,(nTotal/100)).EQ.0)) THEN
    percent=percent+1
    WRITE(0,'(I4,A23,A1)',ADVANCE='NO')percent, ' % of nodes evaluated...',ACHAR(13)
  END IF
  SELECT CASE(InputCoordSys)
  CASE(0)!x_in(1:3) = x,y,z of cylinder with r<1 and z=[0;1]
    r_p   = SQRT(x_in(1,iNode)**2+x_in(2,iNode)**2) 
    theta = ATAN2(x_in(2,iNode),x_in(1,iNode))
    zeta  = -twoPi*x_in(3,iNode) 
  CASE(1) !x_in(1:3) = (r,z,phi) with r= [0;1], z= [0;1], phi=[0;1] 
    r_p =  x_in(1,iNode) ! =r
    theta = twoPi*x_in(3,iNode) ! =2*pi*phi
    zeta  = twoPi*x_in(2,iNode) ! =2*pi*z
  END SELECT 

  ! Phinorm ~ r_p**2 , use scaling of radius to Phi toroidal flux evaluation variable
  ! rho_p = SQRT(Phinorm)=r_p

  rho_p=MIN(1.0_wp,MAX(r_p,1.0E-4)) ! ~ Phinorm [1.0-e08,1.]
  
  IF(useSFL)THEN ! straight-field line coordinates: theta is actually theta^*, find corresponding theta from 
              ! -theta^*+theta+lambda(s,theta,zeta) != 0
    theta_star=theta
    ! prepare evaluation of all modes of lambda at position rho
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
        drhom=REAL(xmabs(iMode),wp)*rho_p**(xmabs(iMode)-1)
      END SELECT
      ! lambda
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,lmns_Spl(:,:,iMode),iGuess,splout) 
      lam_rho(     iMode,1)= rhom*splout(1)
      dldtheta_rho(iMode,1)= rhom*splout(1)*xm(iMode) 
      IF(lasym)THEN
        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,lmnc_Spl(:,:,iMode),iGuess,splout) 
        lam_rho(     iMode,2)= rhom*splout(1)
        dldtheta_rho(iMode,2)= rhom*splout(1)*xm(iMode) 
      END IF
    END DO !iMode

    
    theta=NewtonRoot1D_FdF(1.0e-12,theta_star-Pi,theta_star+Pi,theta_star,theta_star,FRdFR)

  END IF
  CosMN(:)      = COS(    xm(:) * theta -     xn(:) * zeta)
  SinMN(:)      = SIN(    xm(:) * theta -     xn(:) * zeta) 
    ! not needed anymore (only for gmnc)
  ! CosMN_nyq(:)  = COS(xm_nyq(:) * theta - xn_nyq(:) * zeta)
  
  
  R         =0.0_wp
  Z         =0.0_wp
  dRdtheta  =0.0_wp
  dRdzeta   =0.0_wp
  dRdrho    =0.0_wp
  dZdtheta  =0.0_wp
  dZdzeta   =0.0_wp
  dZdrho    =0.0_wp
  lam       =0.0_wp
  dldtheta  =0.0_wp
  dldzeta   =0.0_wp
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
      drhom=REAL(xmabs(iMode),wp)*rho_p**(xmabs(iMode)-1)
    END SELECT
    CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmnc_Spl(:,:,iMode),iGuess,splout) 
    R    = R    + rhom*splout(1)*CosMN(iMode)
    ! dR/drho = sum_mn [ rho**m * dR_mn/drho + R_mn *d(rho**m)/drho ] * cos(m*theta-n*zeta)
    ! dR/dtheta (cos(m*theta-n*zeta))=-m*sin(m*theta-n*zeta)
    ! dR/dzeta  (cos(m*theta-n*zeta))=n*sin(m*theta-n*zeta)
    dRdrho   = dRdrho   + (rhom*splout(2)+splout(1)*drhom)*CosMN(iMode)
    dRdtheta = dRdtheta - rhom*splout(1)*xm(iMode)*SinMN(iMode) 
    dRdzeta  = dRdzeta  + rhom*splout(1)*xn(iMode)*SinMN(iMode)
    IF(lasym) THEN
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmns_Spl(:,:,iMode),iGuess,splout) 
      R    = R    + rhom*splout(1)*SinMN(iMode)
      dRdrho   = dRdrho   + (rhom*splout(2)+splout(1)*drhom)*SinMN(iMode)
      dRdtheta = dRdtheta + rhom*splout(1)*xm(iMode)*CosMN(iMode)
      dRdzeta  = dRdzeta  - rhom*splout(1)*xn(iMode)*CosMN(iMode)
    END IF
    CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmns_Spl(:,:,iMode),iGuess,splout) 
    Z    = Z    + rhom*splout(1)*SinMN(iMode)
    ! dZ/drho = sum_mn [ rho**m * dZ_mn/drho + Z_mn *d(rho**m)/drho ] * sin(m*theta-n*zeta)
    ! dZ/dtheta (sin(m*theta-n*zeta))=m*cos(m*theta-n*zeta)
    ! dZ/dzeta  (sin(m*theta-n*zeta))=-n*cos(m*theta-n*zeta)
    dZdrho   = dZdrho   + (rhom*splout(2)+splout(1)*drhom)*SinMN(iMode)
    dZdtheta = dZdtheta + rhom*splout(1)*xm(iMode)*CosMN(iMode)
    dZdzeta  = dZdzeta  - rhom*splout(1)*xn(iMode)*CosMN(iMode)
    IF(lasym) THEN
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmnc_Spl(:,:,iMode),iGuess,splout) 
      Z    = Z    + rhom*splout(1)*CosMN(iMode)
      dZdrho   = dZdrho   + (rhom*splout(2)+splout(1)*drhom)*CosMN(iMode)
      dZdtheta = dZdtheta - rhom*splout(1)*xm(iMode)*SinMN(iMode)
      dZdzeta  = dZdzeta  + rhom*splout(1)*xn(iMode)*SinMN(iMode)
    END IF
    ! lambda
    CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,lmns_Spl(:,:,iMode),iGuess,splout) 
    lam      = lam      + rhom*splout(1)*SinMN(iMode)
    ! derivatives of lambda is rho not needed
    ! dl/drho = sum_mn [ rho**m * dl_mn/drho + l_mn *d(rho**m)/drho ] * sin(m*theta-n*zeta)
    ! dldrho   = dldrho   + (rhom*splout(2)+splout(1)*drhom)*SinMN(iMode)
    ! dL/dtheta (sin(m*theta-n*zeta))=m*cos(m*theta-n*zeta)
    ! dl/dzeta  (sin(m*theta-n*zeta))=-n*cos(m*theta-n*zeta)
    dldtheta = dldtheta + rhom*splout(1)*CosMN(iMode)*xm(iMode) 
    dldzeta  = dldzeta  - rhom*splout(1)*CosMN(iMode)*xn(iMode)
    IF(lasym) THEN
      CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,lmnc_Spl(:,:,iMode),iGuess,splout) 
      lam      = lam      + rhom*splout(1)*CosMN(iMode)
      dldtheta = dldtheta - rhom*splout(1)*SinMN(iMode)*xm(iMode) 
      dldzeta  = dldzeta  + rhom*splout(1)*SinMN(iMode)*xn(iMode)
    END IF
  END DO !iMode=1,mn_mode 

!  sqrtG     =0.
!  DO iMode=1,mn_mode_nyq
!    IF(xmabs_nyq(iMode).EQ.0)THEN
!      rhom=1.
!      drhom=0.
!    ELSEIF(xmabs_nyq(iMode).EQ.1)THEN
!      rhom=rho_p
!      drhom=1.
!    ELSE
!      rhom=rho_p**xmabs_nyq(iMode)
!      drhom=xmabs_nyq(iMode)*rho_p**(xmabs_nyq(iMode)-1)
!    END IF
!    CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,gmnc_nyq_Spl(:,:,iMode),iGuess,splout) 
!    sqrtG    = sqrtG    + rhom*splout(1)*CosMN_nyq(iMode)
!  END DO !iMode=1,mn_mode_nyq 

  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,pres_Spl(:,:),iGuess,splout) 
  pressure=splout(1)
  !save way to compute dPhi_ds=dPhi_drho*drho/ds, drho/ds=1/(2*rho), where s=Phinorm
  CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Phi_Spl(:,:),iGuess,splout) 
  Phi_int=splout(1)
  dPhi_drho_int=splout(2)
!  dPhi_ds_int=dPhi_drho_int/(2.*rho_p)
!  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,dPhi_ds_Spl(:,:),iGuess,splout) 
!  IF(ABS(dPhi_ds_int-splout(1)).GT.1.0E-07) &
!    WRITE(*,*)'DEBUG,ABS(dPhi_ds-dPhi_drho/(2rho))>1.0-07',dPhi_ds_int,splout(1)

  CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,chi_Spl(:,:),iGuess,splout) 
  chi_int=splout(1)
  dchi_drho_int=splout(2)

!  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,iota_Spl(:,:),iGuess,splout) 
!  iota_int=splout(1)
!   iota should be >0, chi is growing radially, but Phi is decreasing radially 
  iota_int = dchi_drho_int/dPhi_drho_int 

  ! compute magnetic field, following Michael Kraus formulas: 
  ! B^s     = 0 
  ! B^theta = dPhi/ds*(iota-dlambda/dzeta)  : dPhi_ds*(iotaf -  (lVmnc) ) 
  ! B^zeta = dPhi/ds*(1+dlambda/dtheta)     : dPhi_ds*(1+ (lUmnc) )
  !   ...( lmns is overwritten to full mesh and then d/dtheta d/dzeta is applied )
  

  ! compute sqrtG=Jacobian: R*(dRdtheta*dZds-dRds*dZdtheta)
  ! sqrtG=R*(dRdtheta*dZds-dRds*dZdtheta) =  1/(2*rho_p)*R*(dRdtheta*dZdrho-dRdrho*dZdtheta) 
  ! sqrtG=R/(2*rho_p)*sqrtGr 
  sqrtGr= (dRdtheta*dZdrho-dRdrho*dZdtheta) 

! !CHECK WITH INTERPOLATION
! IF(ABS(1-2*rho_p/R*sqrtG/sqrtGr).GT.1.0E-02)  &
!    WRITE(*,'(A,E11.5,A,F11.5)')'rel.err. sqrtG: |1 - 2*rho_p/R*sqrtG/sqrtGr|>1.0E-02 ',1- 2*rho_p/R*sqrtG/sqrtGr, ' Phinorm= ' ,Phinorm_p

  ! contravariant components of B  ! 
  ! B^s    = 0
  ! B^theta = (dchi_ds_int - dPhi_ds_int*dldzeta) /sqrtG
  ! B^zeta  = dPhi_ds_int*(1.  + dldtheta)  /sqrtG
  ! 
  ! Br   =  dRdtheta*Btheta+dRdzeta*Bzeta
  ! Bz   =  dZdtheta*Btheta+dZdzeta*Bzeta
  ! Bphi =               dPhi_dzeta*Bzeta = R*Bzeta

  ! cylindrical components of B 
  Br   =  (  (dchi_drho_int-dPhi_drho_int*dldzeta)*dRdtheta &
                               + dPhi_drho_int*(1.0_wp  + dldtheta)*dRdzeta     )/(R*sqrtGr)
  Bz   =  (  (dchi_drho_int-dPhi_drho_int*dldzeta)*dZdtheta &
                               + dPhi_drho_int*(1.0_wp  + dldtheta)*dZdzeta     )/(R*sqrtGr)
  Btor =                         dPhi_drho_int*(1.0_wp  + dldtheta)              /sqrtGr    ! *R/R


  ! compute cylindrical components of A (R,Z, phi) ,  
  ! needs metric tensor to rho,theta,zeta (from inverse of Jacobian DR/Drho)
  ! coordinates directions (rho,theta,zeta)  expressed in (R,Z,phi)
  !
  ! detJ = dphi_dzeta*(dR_dtheta*dZ_drho - dR_drho*dZ_dtheta) , dphi_dzeta=R
  ! grad(rho  ) = 1/detJ * (/-R*dZ_dtheta, 
  !                           R*dR_dtheta, 
  !                          (dR_dzeta*dZ_dtheta-dR_dtheta*dZ_dzeta) /)
  ! grad(theta) = 1/detJ * (/ R*dZ_drho  , 
  !                          -R*dR_drho  ,
  !                          (dR_drho  *dZ_dzeta-dR_dzeta*dZ_drho  ) /)
  ! grad(zeta ) = (/0 , 0, 1/R /)
 
  ! the vector potential reads as
  ! A=Phi*grad(theta) + Phi*grad(lambda) - chi*grad(zeta) 
  !
  !    with  Phi*grad(lambda)=Phi*(lambda_rho*grad(rho)+lambda_theta*grad(theta)+lamdba_zeta*grad(zeta))
  !
  !A=Phi* [ (1+lambda_theta)*grad(theta) + lambda_rho*grad(rho) ] + (Phi*lambda_zeta-chi)*grad(zeta)


  ! Ar   = Phi_int *(dldrho*(-dZdtheta) + (1.+dldtheta)*( dZdrho)    )/(dRdtheta*dZdrho - dRdrho*dZdtheta) ! *R/R
  ! Az   = Phi_int *(dldrho*(+dRdtheta) + (1.+dldtheta)*(-dRdrho)    )/(dRdtheta*dZdrho - dRdrho*dZdtheta) ! *R/R
  ! Ator = Phi_int *(dldrho*(dRdzeta*dZdtheta-dRdtheta*dZdzeta)                                             & 
  !                                     + (1.+dldtheta)*(dRdrho*dZdzeta-dRdzeta*dZdrho))                    &
  !                                                                /(R*(dRdtheta*dZdrho - dRdrho*dZdtheta)) &
  !        + (Phi_int*dldzeta-chi_int)/R 

  ! OR also, much better, without lambda derivatives:
  ! A= Phi grad(theta) - lambda grad(Phi) - chi grad(zeta) , where grad(Phi) = dPhi_drho *grad(rho)

  Ar   = (-lam*dPhi_drho_int*(-dZdtheta)  + Phi_int *( dZdrho)    )/(sqrtGr) ! *R/R
  Az   = (-lam*dPhi_drho_int*( dRdtheta)  + Phi_int *(-dRdrho)    )/(sqrtGr) ! *R/R
  Ator = (-lam*dPhi_drho_int*(dRdzeta*dZdtheta-dRdtheta*dZdzeta)                                         &
                                          + Phi_int *(dRdrho*dZdzeta-dRdzeta*dZdrho) )                   &
                                                                /(R*(sqrtGr)) &
         -chi_int/R 
 
  ! convert to cartesian
  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  
  x_out(1,iNode)= R*coszeta
  x_out(2,iNode)= R*sinzeta
  x_out(3,iNode)= Z
  
  Bcart(1)= Br*coszeta-Btor*sinzeta
  Bcart(2)= Br*sinzeta+Btor*coszeta
  Bcart(3)= Bz

  Acart(1)= Ar*coszeta-Ator*sinzeta
  Acart(2)= Ar*sinzeta+Ator*coszeta
  Acart(3)= Az
 
  Rmin=MIN(Rmin,R);  Rmax=MAX(Rmax,R) 
  Zmin=MIN(Zmin,Z);  Zmax=MAX(Zmax,Z) 

  !NOTE: Phi is the toroidal flux, called PHI in VMEC
  Phinorm=(Phi_int-Phi_prof(1))/(Phi_prof(nFluxVMEC)-Phi_prof(1))
  ! chi is the poloidal flux 
  chinorm=(chi_int-chi_prof(1))/(chi_prof(nFluxVMEC)-chi_prof(1))

  SELECT CASE (RhoFluxVar)
  CASE(0) !use normalized toroidal flux 
    Density=Eval1DPoly(nRhoCoefs,RhoCoefs,Phinorm) 
  CASE(1) !use normalized poloidal flux 
    Density=Eval1DPoly(nRhoCoefs,RhoCoefs,chinorm) 
  CASE(2) !use sqrt of normalized toroidal flux 
    Density=Eval1DPoly(nRhoCoefs,RhoCoefs,sqrt(Phinorm))
  CASE(3) !use sqrt of normalized poloidal flux 
    Density=Eval1DPoly(nRhoCoefs,RhoCoefs,sqrt(chinorm))
  CASE DEFAULT
    STOP 'RhoFluxVar for vmec eq. not between 0 or 3'
  END SELECT

  MHDEQdata(  1,iNode)=Density
  MHDEQdata(  2,iNode)=pressure*mu0 !pressure transformed to mu0=1
  MHDEQdata( 3:5,iNode)=Bcart(:)
  MHDEQdata(   6,iNode)=chi_int !poloidal flux
  MHDEQdata(   7,iNode)=Phi_int !toroidal flux
  MHDEQdata(8:10,iNode)=Acart(:)

END DO !iNode=1,nTotal
WRITE(UNIT_stdOut,'(A,4(1X,E12.5))')'  Rmin/Zmax, Zmin/Zmax ', Rmin,Rmax,Zmin,Zmax 

WRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '

!for iteration on theta^*
CONTAINS 

  FUNCTION FRdFR(theta_iter)
    IMPLICIT NONE
    REAL:: theta_iter
    REAL:: FRdFR(2)
    CosMN(:)      = COS(    xm(:) * theta_iter -     xn(:) * zeta)
    SinMN(:)      = SIN(    xm(:) * theta_iter -     xn(:) * zeta) 
    lam        =0.0_wp
    dldtheta   =0.0_wp
    DO iMode=1,mn_mode
      lam      = lam      + lam_rho(iMode,1)*SinMN(iMode)
      dldtheta = dldtheta + dldtheta_rho(iMode,1)*CosMN(iMode)
    END DO !iMode=1,mn_mode 
    IF(lasym)THEN
      DO iMode=1,mn_mode
        lam      = lam      + lam_rho(iMode,2)*CosMN(iMode)
        dldtheta = dldtheta - dldtheta_rho(iMode,2)*SinMN(iMode)
      END DO !iMode=1,mn_mode 
    END IF
    FRdFR(1)=theta_iter+lam
    FRdFR(2)=1.0_wp+dldtheta
  END FUNCTION FRdFR

END SUBROUTINE MapToVmec 


END MODULE MOD_VMEC
