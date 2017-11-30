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
!!# Module **Soloviev**
!!
!! Evaluate Soloviev equilibrium
!!
!===================================================================================================================================
MODULE MOD_Solov
! MODULES
USE MOD_Globals,ONLY:UNIT_StdOut,wp
IMPLICIT NONE
PRIVATE

INTERFACE InitSolov 
  MODULE PROCEDURE InitSolov 
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToSolov 
!  MODULE PROCEDURE MapToSolov 
!END INTERFACE

INTERFACE FinalizeSolov 
  MODULE PROCEDURE FinalizeSolov 
END INTERFACE

PUBLIC::InitSolov
PUBLIC::MapToSolov
PUBLIC::FinalizeSolov
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Sololviev module
!!
!===================================================================================================================================
SUBROUTINE InitSolov 
! MODULES
USE MOD_ReadInTools
USE MOD_Solov_Vars
USE MOD_CCint,ONLY:initCCint
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'  INIT SOLOV INPUT ...'

setup =GETINT("Solov_setup",0)
IF(setup.EQ.0)THEN
  P_R0    =GETREAL("Solov_R0")
  P_eps   =GETREAL("Solov_eps")
  P_kappa =GETREAL("Solov_kappa")
  P_delta =GETREAL("Solov_delta")
  P_A     =GETREAL("Solov_A")
  P_B0    =GETREAL("Solov_B0")
  P_qaxis =GETREAL("Solov_qaxis")
  P_paxis =GETREAL("Solov_paxis")
ELSE
  !preselected setups (defaults are set, can still be modified in inifile)
  SELECT CASE(setup)
  CASE(1) !circular
    SWRITE(UNIT_stdout,*)'circular equilibrium setup:'
    P_R0    =GETREAL("Solov_R0"   ,10.0_wp)
    P_eps   =GETREAL("Solov_eps"  ,0.4_wp)
    P_kappa =GETREAL("Solov_kappa",1.0_wp)
    P_delta =GETREAL("Solov_delta",0.0_wp)
    P_A     =GETREAL("Solov_A"    ,0.955_wp)
    P_B0    =GETREAL("Solov_B0"   ,1.0_wp)
    P_qaxis =GETREAL("Solov_qaxis",1.89_wp)
    P_paxis =GETREAL("Solov_paxis",2.0e-03_wp)
  CASE(2) !parameter ITER-like
    SWRITE(UNIT_stdout,*)'Iter-like equilibrium setup:'
    P_R0    =GETREAL("Solov_R0"   ,6.2_wp)
    P_eps   =GETREAL("Solov_eps"  ,0.32_wp)
    P_kappa =GETREAL("Solov_kappa",1.7_wp)
    P_delta =GETREAL("Solov_delta",0.33_wp)
    P_A     =GETREAL("Solov_A"    ,-0.155_wp)
    P_B0    =GETREAL("Solov_B0"   ,1.0_wp)
    P_qaxis =GETREAL("Solov_qaxis",1.6_wp)
    P_paxis =GETREAL("Solov_paxis",8.0e-02_wp)
  CASE DEFAULT
    SWRITE(UNIT_stdout,*) "WARNING: This soloviev setup does not exist" ,setup
    STOP
  END SELECT
END IF !setup=0
IF(p_A    .GT.1.0_wp)  SWRITE(Unit_stdOut,*) 'WARNING, USING POSITIVE PRESSURE GRADIENT p_axis< PresEdge'
IF(p_eps  .GT.0.98_wp) SWRITE(Unit_stdOut,*) 'WARNING, EPSILON >0.98, makes no sense!'
IF(p_delta.GT.0.98_wp) SWRITE(Unit_stdOut,*) 'WARNING, DELTA >0.98, makes no sense!'
  
asin_delta=ASIN(P_delta)

!initialize clenshaw-curtis module for integration later in MapToSolov
CALL InitCCint()

!initialize Soloviev parameters
CALL InitSolovievEquilibrium()

SWRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE InitSolov


!===================================================================================================================================
!> Solve for psiCoefs depending on input data
!!
!===================================================================================================================================
SUBROUTINE InitSolovievEquilibrium 
! MODULES
USE MOD_Solov_Vars
USE MOD_PsiEval, ONLY:EvalPsi
USE MOD_PsiEval, ONLY:EvaldPsi
USE MOD_Newton, ONLY:NewtonMin1D
USE MOD_Newton, ONLY:NewtonMin2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)  :: tol,a(2),b(2)
REAL(wp)  :: qaxis_fact,mu0,dp_dpsin,qedge
REAL(wp)  :: F_edge
REAL(wp)  :: d2Psi_dx2_axis  !d^2/dx^2 Psi, second derivatives at axis
REAL(wp)  :: d2Psi_dy2_axis  !d^2/dy^2 Psi
!===================================================================================================================================

CALL SolveForPsiCoefs()
SWRITE(Unit_stdOut,'(4X,A,4(E22.15,1X))')"PsiCoefs(1:4)        : ", PsiCoefs(1:4)
SWRITE(Unit_stdOut,'(4X,A,3(E22.15,1X))')"PsiCoefs(5:7)        : ", PsiCoefs(5:7)

!                tolerance,
tol=1.0e-15
a(1)=1.0_wp-0.1_wp*p_eps !bound x- 
b(1)=1.0_wp+0.4_wp*p_eps !bound x+
a(2)=-0.1_wp*p_kappa*p_eps !bound y-
b(2)=-a(2)              !bound y+
             

!1D minimum find, since we know here y=0
xaxis(1)=1.0_wp
xaxis(2)=0.0_wp
psi_axis=NewtonMin1D(tol,a(1),b(1),xaxis(1),F1,dF1,ddF1 )
SWRITE(Unit_stdOut,'(A,2E22.15)')'   magentic axis x,y    : ',xaxis(:)

!2D minimum search
xaxis(1)=1.0_wp
xaxis(2)=0.0_wp
psi_axis=NewtonMin2D(tol,a,b,xaxis,F2,dF2,ddF2 )

psi_edge=EvalPsi(1.0_wp+p_eps,0.0_wp)
SWRITE(UNIT_stdOut,'(A,2E22.15)')'   magentic axis x,y    : ',xaxis(:)
SWRITE(UNIT_stdOut,'(A,E22.15)') '   psi at magentic axis : ',psi_axis
SWRITE(UNIT_StdOut,'(A,E22.15)') '   psi at domain edge   : ',psi_edge

!find out scaling of psi function
F_axis=p_B0*p_R0

! q_axis=F_axis/(R_axis*(d_rr(psireal)*d_zz(psireal))^1/2), psireal=psi*psi0, evaluated at axis...
!       =F_axis/psi_scale* 1/(R_axis*(d_rr(psi)*d_zz(psi))^1/2), from book Freidberg, ideal MHD, p.134
! R=x*R0,Z=y*R0
! d/dR=(d/dx)*1/R0, d/dZ=(d/dy)*1/R0

d2Psi_dx2_axis = EvaldPsi(1,2,xaxis(1),xaxis(2))
d2Psi_dy2_axis = EvaldPsi(2,2,xaxis(1),xaxis(2))
qaxis_fact=p_R0/(xaxis(1)*SQRT(d2Psi_dx2_axis*d2Psi_dy2_axis) )

psi_scale=qaxis_fact*F_axis/p_qaxis

SWRITE(UNIT_stdOut,'(A,E22.15)') '   psi_scale            : ',psi_scale
SWRITE(UNIT_StdOut,'(A,E22.15)') '   real psi on axis     : ',psi_scale*psi_axis

deltaF2 = 2*p_A*psi_scale**2*psi_axis/(p_R0**2)
F_edge=SQRT(F_axis**2+deltaF2)

SWRITE(UNIT_stdOut,'(A,E22.15)') '   F on axis            : ',F_axis
SWRITE(UNIT_StdOut,'(A,E22.15)') '   F on edge            : ',F_edge

!gradient of the pressure over normalized flux, p=p0+dpdpsin*psin 
mu0=1.
dp_dpsin=(1-p_A)*psi_scale**2*psi_axis/(mu0*p_R0**4)
PresEdge=p_paxis+dp_dpsin
SWRITE(UNIT_stdOut,'(A,E22.15)') '   pressure on axis     : ',p_paxis
SWRITE(UNIT_StdOut,'(A,E22.15)') '   pressure on edge     : ',PresEdge

IF(PresEdge.LT.0.0_wp) STOP 'Pressure negative in domain, increase paxis input'

CALL Eval_qedge(F_edge,qedge)

SWRITE(UNIT_stdOut,'(A,E22.15)') '   q-factor on axis     : ',p_qaxis
SWRITE(UNIT_StdOut,'(A,E22.15)') '   q-factor on edge     : ',qedge

CONTAINS
  !for 1d NewtonMin
  FUNCTION F1(x)
    IMPLICIT NONE
    REAL(wp):: x
    REAL(wp):: F1
    F1=EvalPsi(x,0.0_wp)
  END FUNCTION F1

  FUNCTION dF1(x)
    IMPLICIT NONE
    REAL(wp):: x
    REAL(wp):: dF1
    dF1=EvaldPsi(1,1,x,0.0_wp)
  END FUNCTION dF1

  FUNCTION ddF1(x)
    IMPLICIT NONE
    REAL(wp):: x
    REAL(wp):: ddF1
    ddF1=EvaldPsi(1,2,x,0.0_wp)
  END FUNCTION ddF1

  !for 2d NewtonMin
  FUNCTION F2(x)
    IMPLICIT NONE
    REAL(wp):: x(2)
    REAL(wp):: F2
    F2=EvalPsi(x(1),x(2))
  END FUNCTION F2

  FUNCTION dF2(x)
    IMPLICIT NONE
    REAL(wp):: x(2)
    REAL(wp):: dF2(2)
    dF2(1)=EvaldPsi(1,1,x(1),x(2))
    dF2(2)=EvaldPsi(2,1,x(1),x(2))
  END FUNCTION dF2

  FUNCTION ddF2(x)
    IMPLICIT NONE
    REAL(wp):: x(2)
    REAL(wp):: ddF2(2,2)
 
    ddF2(1,1)=EvaldPsi(1,2,x(1),x(2))
    ddF2(1,2)=EvaldPsi(3,1,x(1),x(2))
    ddF2(2,1)=ddF2(1,2)
    ddF2(2,2)=EvaldPsi(2,2,x(1),x(2))
  END FUNCTION ddF2

END SUBROUTINE InitSolovievEquilibrium 


!===================================================================================================================================
!> compute qfactor on edge, needs an circular integration over a flux surface, here given by boundary curve  
!!
!===================================================================================================================================
SUBROUTINE Eval_qedge(Fedge,qedge)
! MODULES
USE MOD_globals   , ONLY:Pi
USE MOD_Solov_Vars, ONLY:asin_delta,p_eps,p_kappa,xaxis,p_R0,psi_scale
USE MOD_PsiEval   , ONLY:EvalPsi,EvaldPsi
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)  :: Fedge
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: qedge
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: i,nmax
REAL(wp) :: x,y,rcoord,BthetaR
REAL(wp) :: tau,dtau
REAL(wp) :: psimax,psimin,psi_dx,psi_dy,intgr
!===================================================================================================================================
nmax=1000
dtau=2*pi/REAL(nmax) 
tau=0.0_wp
qedge=0.0_wp
DO i=1,nmax
  x = 1.0_wp + p_eps * COS(tau + asin_delta* SIN(tau))
  y = p_eps * p_kappa * SIN(tau)
  rcoord=SQRT((x-xaxis(1))**2+(y-xaxis(2))**2)
  psimin=MIN(psimin,EvalPsi(x,y))
  psimax=MAX(psimax,EvalPsi(x,y))
  psi_dx=EvaldPsi(1,1,x,y)
  psi_dy=EvaldPsi(2,1,x,y)
  !BR=-1/R*dpsireal_dZ, Bz=1/R*dpsireal_dR 
  ! Btheta=-BR*sin(theta)+BZ*cos(theta)
  ! Btheta*R= dpsireal_dZ*sin(theta)+dpsireal_dR*cos(theta), dpsiReal_dZ=psi_scale/R0*dpsi_dy,  sin=(y-y_a)/r, cos=(x-x_a)/r
  BthetaR=psi_scale/(p_R0*rcoord)*(psi_dy*(y-xaxis(2))+psi_dx*(x-xaxis(1))) !/R=(x*R0)
  
  intgr = rcoord/(x*BthetaR)  ! oint r/(R^2*Btheta) dtheta, r=rcoord*R0, R=x*R0
  qedge = qedge+intgr
  tau=tau+dtau
END DO !i
qedge=Fedge*qedge/REAL(nmax) !=Fedge/2*pi*qedge*dtau

END SUBROUTINE Eval_qedge


!===================================================================================================================================
!> Approximately flux-aligned coordinates, rr[0,1] ~sqrt(psi), theta [0,2pi]
!!
!===================================================================================================================================
FUNCTION ApproxFluxMap(rr,beta) RESULT(xyCoords)
! MODULES
USE MOD_Solov_Vars,ONLY:asin_delta,p_eps,p_kappa,xaxis
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN) :: rr     !! linear radial coordinate [0,1] ~ sqrt(psi)
REAL(wp), INTENT(IN) :: beta   !! poloidal angle (not atan(y,x)!) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)  :: xyCoords(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)  :: sinbeta
!===================================================================================================================================
sinbeta=SIN(beta)
xyCoords(1) =1.0_wp+(xaxis(1)-1.0_wp)*(1.0_wp-rr**2)+p_eps * rr*COS(beta + asin_delta*rr * sinbeta ) !edge_X-1, mit p_delta=rr*p_delta
xyCoords(2) = rr*p_eps * p_kappa * sinbeta
         
END FUNCTION ApproxFluxMap


!===================================================================================================================================
!> Approximately flux-aligned coordinates, rr[0,1] ~sqrt(psi), theta [0,2pi]
!!
!===================================================================================================================================
FUNCTION ApproxFluxMap_dr(rr,beta) RESULT(xyCoords_dr)
! MODULES
USE MOD_Solov_Vars,ONLY:asin_delta,p_eps,p_kappa,xaxis
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN) :: rr     !! linear radial coordinate [0,1] ~ sqrt(psi)
REAL(wp), INTENT(IN) :: beta   !! poloidal angle (not atan(y,x)!) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)  :: xyCoords_dr(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)  :: sinbeta
!===================================================================================================================================
sinbeta=SIN(beta)

xyCoords_dr(1) =-2.0_wp*rr*(xaxis(1)-1.0_wp)+p_eps * (COS(beta + asin_delta*rr * sinbeta ) &
                                             -rr*(SIN(beta + asin_delta*rr * sinbeta ))*asin_delta*sinbeta)
xyCoords_dr(2) = p_eps * p_kappa * sinbeta
         
END FUNCTION ApproxFluxMap_dr


!===================================================================================================================================
!> evaluate normalized flux [0,1], 0 on axis, 1 on edge
!!
!===================================================================================================================================
FUNCTION PsiToPsiNorm(psi) RESULT(PsiNorm)
! MODULES
USE MOD_Solov_Vars,ONLY:psi_axis,psi_edge
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN) :: psi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp) :: PsiNorm
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
PsiNorm=(Psi-psi_axis)/(psi_edge-psi_axis)
END FUNCTION PsiToPsiNorm


!===================================================================================================================================
!> evaluate normalized flux [0,1], 0 on axis, 1 on edge
!!
!===================================================================================================================================
FUNCTION PsiNormToPsi(psiNorm) RESULT(Psi)
! MODULES
USE MOD_Solov_Vars,ONLY:psi_axis,psi_edge
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN) :: psiNorm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp) :: Psi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Psi=PsiNorm*(psi_edge-psi_axis)+psi_axis
END FUNCTION PsiNormToPsi


!===================================================================================================================================
!> Builds up the linear system and solves for psiCoefs(0:7) 
!!
!===================================================================================================================================
SUBROUTINE SolveForPsiCoefs()
! MODULES
USE MOD_Solov_Vars
USE MOD_LinAlg,  ONLY:SOLVE
USE MOD_PsiEval, ONLY:EvalPsiVec,EvaldPsidxVec,Evald2PsidxVec
USE MOD_PsiEval, ONLY:EvaldPsidyVec,Evald2PsidyVec
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
REAL(wp)            :: SysMat(7,7),RHS(7)
REAL(wp)            :: xm,xc,xp,yp
REAL(wp)            :: N1,N2,N3
REAL(wp),DIMENSION(0:7) :: Psi_xmy0,dPsidx_xmy0,d2Psidy_xmy0
REAL(wp),DIMENSION(0:7) :: Psi_xpy0,dPsidx_xpy0,d2Psidy_xpy0
REAL(wp),DIMENSION(0:7) :: Psi_xcyp,dPsidy_xcyp,d2Psidx_xcyp,dPsidx_xcyp
!===================================================================================================================================
PsiCoefs(0)=1.0_wp

xm=1.0_wp-p_eps
xc=1.0_wp-p_delta*p_eps
xp=1.0_wp+p_eps
yp=p_kappa*p_eps
N1=  - (1.0_wp + asin_delta)**2 / (p_eps * p_kappa**2)
N2 = + (1.0_wp - asin_delta)**2 / (p_eps * p_kappa**2)
N3 = - p_kappa / (p_eps * COS(asin_delta)**2)

    Psi_xmy0=EvalPsiVec(xm,0.0_wp)
 dPsidx_xmy0=EvaldPsidxVec(xm,0.0_wp)
d2Psidy_xmy0=Evald2PsidyVec(xm,0.0_wp)

    Psi_xpy0=EvalPsiVec(xp,0.0_wp)
 dPsidx_xpy0=EvaldPsidxVec(xp,0.0_wp)
d2Psidy_xpy0=Evald2PsidyVec(xp,0.0_wp)

    Psi_xcyp=EvalPsiVec(xc,yp)
 dPsidx_xcyp=EvaldPsidxVec(xc,yp)
 dPsidy_xcyp=EvaldPsidyVec(xc,yp)
d2Psidx_xcyp=Evald2PsidxVec(xc,yp)

RHS(1) = -Psi_xpy0(0)
RHS(2) = -Psi_xmy0(0)
RHS(3) = -Psi_xcyp(0)
RHS(4) = -dPsidx_xcyp(0)
RHS(5) = -(d2Psidy_xpy0(0)+ N1*dPsidx_xpy0(0))
RHS(6) = -(d2Psidy_xmy0(0)+ N2*dPsidx_xmy0(0))
RHS(7) = -(d2Psidx_xcyp(0)+ N3*dPsidy_xcyp(0))

DO i=1,7
  sysmat(1,i) = Psi_xpy0(i)
  sysmat(2,i) = Psi_xmy0(i)
  sysmat(3,i) = Psi_xcyp(i)
  sysmat(4,i) = dPsidx_xcyp(i)
  sysmat(5,i) = d2Psidy_xpy0(i)+ N1*dPsidx_xpy0(i)
  sysmat(6,i) = d2Psidy_xmy0(i)+ N2*dPsidx_xmy0(i)
  sysmat(7,i) = d2Psidx_xcyp(i)+ N3*dPsidy_xcyp(i)
END DO !i=1,7

PsiCoefs(1:7)=SOLVE(sysmat,RHS)

END SUBROUTINE SolveForPsiCoefs


!===================================================================================================================================
!> Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from a tokamak soloviev equilibrium. 
!! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!! for a fixed theta, uses a newton method to find the exact location  in r of psi(x(r),y(r))-Psi0=0, psi0=psi(psinorm=r_p^2),
!! using an approximate map(R,Z)<->(rho,theta) of the soloviev equilibrium. 
!!
!===================================================================================================================================
SUBROUTINE MapToSolov(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars,  ONLY:nVarMHDEQ
USE MOD_MHDEQ_Vars,  ONLY: nRhoCoefs,RhoFluxVar,RhoCoefs
USE MOD_Newton,      ONLY:NewtonRoot1D
USE MOD_Solov_Vars,  ONLY:p_R0,p_kappa,p_paxis,PresEdge
USE MOD_Solov_Vars,  ONLY:F_axis,deltaF2,xaxis,psi_scale
USE MOD_PsiEval,     ONLY:EvalPsi,EvaldPsi
USE MOD_CCint,       ONLY:CCint
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,  INTENT(IN) :: nTotal         !! total number of points
REAL(wp), INTENT(IN) :: x_in(3,nTotal) !! input coordinates represent a cylinder: 
INTEGER,  INTENT(IN) :: InputCoordSys  !! =0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
                                       !! =1: x_in(1:3) are (r,z,phi) coordinates r= [0;1], z= [0;1], phi=[0;1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: x_out(3,nTotal) !! mapped x,y,z coordinates with vmec data
REAL(wp),INTENT(OUT) :: MHDEQdata(nVarMHDEQ,nTotal) !! MHD equilibrium data 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)  :: r_p    ! raduis in cylindrical coordinate system
REAL(wp)  :: theta  ! some poloidal angle [0,2pi] ,THIS IS NOT =atan((Z-Zaxis)/(R-Raxis))
REAL(wp)  :: zeta   ! toroidal angle [0,2pi]
INTEGER   :: iNode
INTEGER   :: percent
REAL(wp)  :: coszeta,sinzeta
REAL(wp)  :: psiVal,psiNorm 
REAL(wp)  :: tol,rNewton,xPos(2)
REAL(wp)  :: Density
REAL(wp)  :: R,dPsi_dx,dPsi_dy,Fval,alpha_sRho2
REAL(wp)  :: BR,BZ,Bphi,Bcart(3) 
REAL(wp)  :: AR,AZ,Aphi,Acart(3)
LOGICAL   :: converged
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO SOLOVIEV EQUILIBRIUM'
percent=0
tol=1.0E-13_wp
DO iNode=1,nTotal
  ! output of progress in %
  IF((nTotal.GT.10000).AND.(MOD(iNode,(nTotal/100)).EQ.0)) THEN
    percent=percent+1
    WRITE(0,'(I4,A23,A1)',ADVANCE='NO')percent, ' % of nodes evaluated...',ACHAR(13)
  END IF
  SELECT CASE(InputCoordSys)
  CASE(0)!x_in(1:3) = x,y,z of cylinder with r<1 and z=[0;1]
    r_p   = SQRT(x_in(1,iNode)**2+x_in(2,iNode)**2) 
    Theta = ATAN2(x_in(2,iNode),x_in(1,iNode))
    zeta  = -2.0_wp*Pi*x_in(3,iNode) 
  CASE(1) !x_in(1:3) = (r,z,phi) with r= [0;1], z= [0;1], phi=[0;1] 
    r_p =  x_in(1,iNode) !=r
    Theta =  2.0_wp*Pi*x_in(3,iNode) !=2*pi*phi
    zeta  = -2.0_wp*Pi*x_in(2,iNode) !=2*pi*z
  END SELECT 
  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  psiNorm=r_p**2 !r ~ sqrt(psiNorm), r^2 ~ psiNorm 
  psiVal = PsiNormToPsi(psiNorm)  

  theta = theta -0.1_wp*(p_kappa-1.0_wp)*SIN(2.0_wp*theta) !slight angle correction with ellipticity
  IF(r_p.LT.tol) THEN
    rNewton=r_p
  ELSE
    rNewton= NewtonRoot1D(tol,0.0_wp,1.1_wp,r_p,psiVal,FR1,dFR1)
  END IF
  xPos=ApproxFluxMap(rNewton,theta)
  IF(ABS(PsiVal-EvalPsi(xPos(1),xPos(2))).GT.10*tol) THEN
    !WRITE(*,*)psinorm,xPos,PsiVal,EvalPsi(xPos(1),xPos(2)),ABS(PsiVal-EvalPsi(xPos(1),xPos(2)))
    WRITE(*,*) 'Newton converged, but (PsiVal-Psi(x,y)) still > 10*tol'
    STOP
  END IF
  
  R=p_R0*xPos(1)
  !Z=p_R0*xPos(2)

  x_out(1,iNode)= R*COS(zeta)
  x_out(2,iNode)= R*SIN(zeta)
  x_out(3,iNode)= (p_R0*xPos(2))

  SELECT CASE (RhoFluxVar)
  CASE(1) !use normalized poloidal flux (here psi!!!)
    Density=Eval1DPoly(nRhoCoefs,RhoCoefs,psinorm) 
  CASE(3) !use sqrt of normalized poloidal flux (here psi!!!)
    Density=Eval1DPoly(nRhoCoefs,RhoCoefs,sqrt(psinorm))
  CASE DEFAULT
    STOP 'RhoFluxVar for Solov eq. can only depend poloidal flux(=1) or sqrt(norm.pol.flux) (=3)'
  END SELECT
  !BR=-1/R*dpsiReal_dZ = -1/(x*R0)*psi_scale*dpsi_dy*dy/dZ = -psiscale/(x*R0^2)*dpsi_dy
  !BZ= 1/R*dpsiReal_dR =  1/(x*R0)*psi_scale*dpsi_dx*dx/dR =  psiscale/(x*R0^2)*dpsi_dx
  !Bphi=F/R
  dpsi_dx=EvaldPsi(1,1,xPos(1),xPos(2)) 
  dpsi_dy=EvaldPsi(2,1,xPos(1),xPos(2)) 

  Fval = SQRT(F_axis**2+deltaF2*psiNorm) 

  BR= -psi_scale/(p_R0*R)*dpsi_dy
  BZ=  psi_scale/(p_R0*R)*dpsi_dx
  Bphi= Fval/R
  
  !AR=AR0+deltaAR, AZ=AZ0+deltaAZ
  ! curl(A0)=Faxis/R => AR0=0 AZ0 = -Faxis*log(R)
  !
  ! curl(deltaA) = deltaF(Psi)/R, deltaF=sqrt(Faxis^2-dF2*psinorm)-Faxis
  ! introduce a polar coordinate rho=sqrt((R-Raxis)^2+(Z-Zaxis)^2) and the 
  ! real poloidal angle phi=atan2((Z-Zaxis)/(R-Raxis)), and integrate along the line from axis to 
  ! current point
  ! alpha(rho,phi)= \int_0^rho deltaF(Psi)/R rhohat drhohat , Psi=Psi(x(rhohat,phi),y(rhohat,phi)), R=R0*x(rhohat,phi)
  ! dphi/dR = -(Z-Zaxis)/rho^2, dphi/dZ = (R-Raxis)/rho^2
  !
  ! => deltaAR =  -alpha * dphi/dR, deltaAZ=-alpha *dphi/dZ
  alpha_sRho2=CCint(tol,FI1,converged) ! *rho^2,  since integral normalized to [0,1]

  AR  =  alpha_sRho2*P_R0*(xPos(2)-xaxis(2))
  AZ  = -alpha_sRho2*P_R0*(xPos(1)-xaxis(1))
  Aphi= psi_scale*PsiVal/R

!DEBUG 1: poloidal Bfield
!  BR= -psi_scale/(p_R0*R)*dpsi_dy
!  BZ=  psi_scale/(p_R0*R)*dpsi_dx
!  Bphi= 0.
!  AR  = 0.
!  AZ  = 0. 
!  Aphi= psi_scale*PsiVal/R

!DEBUG 2: toroidal full F field, integrated
!  BR= 0.
!  BZ= 0.
!  Bphi= Fval/R
!  AR  =  alpha_sRho2*P_R0*(xPos(2)-xaxis(2))
!  AZ  = -alpha_sRho2*P_R0*(xPos(1)-xaxis(1))
!  Aphi= 0.

!DEBUG 3: toroidal  field analytical to constant F_axis
!  BR= 0.
!  BZ= 0.
!  Bphi= F_axis/R
!  AR  = 0.
!  AZ  = -F_axis*LOG(R) 
!  Aphi= 0.

!DEBUG 4: toroidal deltaF field (integrate only (Fval-Faxis))
!  BR= 0.
!  BZ= 0.
!  Bphi= (Fval-F_axis)/R
!  AR  =  alpha_sRho2*P_R0*(xPos(2)-xaxis(2))
!  AZ  = -alpha_sRho2*P_R0*(xPos(1)-xaxis(1))
!  Aphi= 0.

  Bcart(1)= BR*coszeta-Bphi*sinzeta
  Bcart(2)= BR*sinzeta+Bphi*coszeta
  Bcart(3)= BZ

  Acart(1)= AR*coszeta-Aphi*sinzeta
  Acart(2)= AR*sinzeta+Aphi*coszeta
  Acart(3)= AZ

  MHDEQdata(:,iNode)=0.0_wp
  MHDEQdata(  1,iNode)=Density
  MHDEQdata(  2,iNode)=p_paxis + (PresEdge-p_paxis)*psiNorm !linear pressure profile
  MHDEQdata( 3:5,iNode)=Bcart(:)
  MHDEQdata(   6,iNode)=psi_scale*PsiVal !polodial flux
  MHDEQdata(   7,iNode)=0.0_wp           !toroidal flux, not known
  MHDEQdata(8:10,iNode)=Acart(:)
END DO !iNode


SWRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '

CONTAINS
!for newton search
  FUNCTION FR1(r)
    !USE MOD_PsiEval,ONLY:EvalPsi
    !uses ApproxFluxMap function
    !uses current theta from subroutine
    IMPLICIT NONE
    !-----------------------------------------------------------------
    REAL(wp) :: r
    REAL(wp) :: FR1
    !local
    REAL(wp) :: xloc(2)
    !-----------------------------------------------------------------
    xloc=ApproxFluxMap(r,Theta) !theta from subroutine!
    FR1=EvalPsi(xloc(1),xloc(2))
  END FUNCTION FR1

  FUNCTION dFR1(r)
    !USE MOD_PsiEval,ONLY:EvaldPsi
    !uses ApproxFluxMap function
    !uses ApproxFluxMap_dr function
    !uses current theta from subroutine
    IMPLICIT NONE
    !-----------------------------------------------------------------
    REAL(wp) :: r
    REAL(wp) :: dFR1
    !local
    REAL(wp) :: xloc(2),dxloc_dr(2)
    !-----------------------------------------------------------------
    xloc=ApproxFluxMap(r,Theta) !theta from subroutine!
    dxloc_dr=ApproxFluxMap_dr(r,Theta) !theta from subroutine!
    !d/dr(psi(x(r),y(r))) = dpsi_dx(r)*dx_dr(r) +dpsi_dy*dy_dr(r)
    dFR1= EvaldPsi(1,1,xloc(1),xloc(2))*dxloc_dr(1) &
          +EvaldPsi(2,1,xloc(1),xloc(2))*dxloc_dr(2)
  END FUNCTION dFR1

  !Function for cc integration to compute poloidal vector potential
  FUNCTION FI1(nn,rhohat)
    !USE MOD_PsiEval,     ONLY:EvalPsi
    !uses PsiToPsiNorm(psi) 
    !uses current point position xPos(1:2) from subroutine
    !uses xaxis,F_axis,deltaF2,p_R0 from subroutine
    IMPLICIT NONE
    !-----------------------------------------------------------------
    INTEGER  :: nn
    REAL(wp) :: rhohat(1:nn) !integration interval, should be [0,1]
    REAL(wp) :: FI1(1:nn)
    !local
    INTEGER  :: i
    REAL(wp) :: xloc(2),Fpsi
    !-----------------------------------------------------------------
    !Integrand alpha = F(Psi(x,y))*rhohat/R(x) 
    DO i=1,nn
      xloc=rhohat(i)*(xPos(:)-xaxis(:))+xaxis(:) !position at rhohat
      Fpsi=SQRT(F_axis**2+deltaF2*PsiToPsiNorm(EvalPsi(xloc(1),xloc(2))))
      FI1(i)=Fpsi*rhohat(i)/(P_R0*xloc(1))  !xloc(1)is always >0
    END DO !i=1,nn
  END FUNCTION FI1


END SUBROUTINE MapToSolov 

!===================================================================================================================================
!> Finalize Sololviev module
!!
!===================================================================================================================================
SUBROUTINE FinalizeSolov
! MODULES
USE MOD_CCint,ONLY:FinalizeCCint
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL FinalizeCCint()


END SUBROUTINE FinalizeSolov


END MODULE MOD_Solov
