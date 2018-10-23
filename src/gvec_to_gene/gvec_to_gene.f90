!===================================================================================================================================
! Copyright (C) 2018  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2018  Maurice Maurer <maurice_maurer@gmx.de>
! Copyright (C) 2018  Alejandro Banon Navarro <abanonna@ipp.mpg.de>
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
!!# Module **gvec_to_gene**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_gene
! MODULES
USE MODgvec_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE init_gvec_to_gene
  MODULE PROCEDURE init_gvec_to_gene
END INTERFACE
!
INTERFACE gvec_to_gene_scalars
  MODULE PROCEDURE gvec_to_gene_scalars
END INTERFACE

INTERFACE gvec_to_gene_profile
  MODULE PROCEDURE gvec_to_gene_profile
END INTERFACE

INTERFACE gvec_to_gene_coords
  MODULE PROCEDURE gvec_to_gene_coords
END INTERFACE

INTERFACE gvec_to_gene_metrics
  MODULE PROCEDURE gvec_to_gene_metrics
END INTERFACE

INTERFACE finalize_gvec_to_gene
  MODULE PROCEDURE finalize_gvec_to_gene
END INTERFACE

PUBLIC::init_gvec_to_gene
PUBLIC::gvec_to_gene_scalars
PUBLIC::gvec_to_gene_profile
PUBLIC::gvec_to_gene_coords
PUBLIC::gvec_to_gene_metrics
PUBLIC::finalize_gvec_to_gene

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE init_gvec_to_gene(fileName) 
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_ReadState,ONLY: ReadState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileName !< name of GVEC file
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT EVAL GVEC ...'

  CALL ReadState(fileName)

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE init_gvec_to_gene


!===================================================================================================================================
!> Scalar variables of the equilibrium
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_scalars(Fa,minor_r,PhiPrime_edge,q_edge,n0_global)
! MODULES
USE MODgvec_globals,ONLY: TWOPI
USE MODgvec_ReadState_Vars,ONLY: a_minor,X1_base_r,profiles_1d
USE MODgvec_MHD3D_profiles,ONLY: Eval_PhiPrime,Eval_chiPrime
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: Fa                !! toroidal flux at the edge
REAL(wp),INTENT(OUT) :: minor_r           !! length scale, minor radius
REAL(wp),INTENT(OUT) :: phiPrime_edge     !! toroidal flux derivative dPhi/ds at the edge. phi'=chi'*q
REAL(wp),INTENT(OUT) :: q_edge            !! q-profile evaluated at the edge. q=phi'/chi'
INTEGER, INTENT(OUT) :: n0_global         !! number of field periods
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Fa = TWOPI*X1_base_r%s%evalDOF_s(1.0, 0,profiles_1d(:,1)) !phi(s=1)
minor_r   = a_minor
n0_global = X1_base_r%f%nfp
PhiPrime_edge = X1_base_r%s%evalDOF_s(1.0, DERIV_S,profiles_1d(:,1))
q_edge    = 1./(X1_base_r%s%evalDOF_s(1.0,       0,profiles_1d(:,3)) ) !q=1/iota

END SUBROUTINE gvec_to_gene_scalars


!===================================================================================================================================
!> Evaluate only s dependend variables 
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_profile(spos,q,q_prime,p,p_prime)
! MODULES
USE MODgvec_ReadState_Vars,ONLY: profiles_1d,X1_base_r
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: spos            !! radial position (sqrt(phi_norm)), phi_norm: normalized toroidal flux [0,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),OPTIONAL,INTENT(OUT) :: q                !! q=1/iota profile
REAL(wp),OPTIONAL,INTENT(OUT) :: q_prime          !! dq/ds=-(d/ds iota)/iota^2=-(d/ds iota)*q^2
REAL(wp),OPTIONAL,INTENT(OUT) :: p                !! p, pressure profile
REAL(wp),OPTIONAL,INTENT(OUT) :: p_prime          !! dp/ds, derivative of pressure profile 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(PRESENT(q)) q       = 1./(  X1_base_r%s%evalDOF_s(spos,       0,profiles_1d(:,3)) ) !q=1/iota
IF(PRESENT(q_prime)) q_prime = -q*q*(X1_base_r%s%evalDOF_s(spos, DERIV_S,profiles_1d(:,3)) ) !q'=-iota'/iota^2
IF(PRESENT(p)) p =      (X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,4)) ) !pressure
IF(PRESENT(p_prime)) p_prime =      (X1_base_r%s%evalDOF_s(spos, DERIV_S,profiles_1d(:,4)) ) !pressure'

END SUBROUTINE gvec_to_gene_profile


!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_coords(nthet,nzeta,spos_in,theta_star_in,zeta_in,theta_out,cart_coords)
! MODULES
USE MODgvec_ReadState_Vars
USE MODgvec_globals, ONLY: PI
USE MODgvec_Newton,  ONLY: NewtonRoot1D_FdF
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER              :: nthet          !! number of points in theta_star
INTEGER              :: nzeta          !! number of points in zeta
REAL(wp),INTENT( IN) :: spos_in        !! radial position (sqrt(phi_norm)), phi_norm: normalized toroidal flux [0,1]
REAL(wp),INTENT( IN) :: theta_star_in(nthet,nzeta)  !! thetaStar poloidal angle (straight field line angle PEST)
REAL(wp),INTENT( IN) :: zeta_in(      nthet,nzeta)  !! zeta toroidal angle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp), INTENT(OUT) :: theta_out(nthet,nzeta)
REAL(wp),INTENT(OUT) :: cart_coords(3,nthet,nzeta)  !! x,y,z cartesian coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iMode,ithet,izeta
REAL(wp)    :: theta_star,zeta
REAL(wp)    :: xp(2),qvec(3)
REAL(wp)    :: X1_s(   1:X1_base_r%f%modes)
REAL(wp)    :: X2_s(   1:X2_base_r%f%modes)
REAL(wp)    :: LA_s(   1:LA_base_r%f%modes)
REAL(wp)    :: X1_int,X2_int,spos
!===================================================================================================================================
spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,spos_in)) !for satefy reasons at the axis and edge
!interpolate first in s direction
DO iMode=1,X1_base_r%f%modes
  X1_s(iMode)      =X1_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode)) !R
END DO
DO iMode=1,X2_base_r%f%modes
  X2_s(iMode)      =X2_base_r%s%evalDOF_s(spos,      0,X2_r(:,iMode)) !Z
END DO
DO iMode=1,LA_base_r%f%modes
  LA_s(iMode)      =LA_base_r%s%evalDOF_s(spos,      0,LA_r(:,iMode)) !lambda
END DO

DO izeta=1,nzeta; DO ithet=1,nthet
  theta_star = theta_star_in(ithet,izeta) !theta_star depends on zeta!!
  zeta       = zeta_in(      ithet,izeta)
  !find angle theta from straight field line angle (PEST) theta_star=theta+lambda(s,theta,zeta) 
  ! 1D Newton uses derivative function FRdFR defined below... solves FR(1)-F0=0, FR(2)=dFR(1)/dtheta
  !                      tolerance , lower bound , upper bound , maxstep, start value , F0
  theta_out(ithet,izeta)=NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1*PI  ,theta_star   , theta_star,FRdFR)

  xp=(/theta_out(ithet,izeta),zeta/)

  X1_int      = X1_base_r%f%evalDOF_x(xp,0,X1_s)
  X2_int      = X2_base_r%f%evalDOF_x(xp,0,X2_s)

  qvec = (/X1_int,X2_int,zeta/)
  cart_coords(:,ithet,izeta)=hmap_r%eval(qvec)

END DO; END DO !ithet,izeta

!for iteration on theta^*
CONTAINS 

  FUNCTION FRdFR(theta_iter)
    !uses current zeta where newton is called, and LA_s from subroutine above
    IMPLICIT NONE
    REAL(wp) :: theta_iter
    REAL(wp) :: FRdFR(2) !output
    !--------------------------------------------------- 
    FRdFR(1)=theta_iter+LA_base_r%f%evalDOF_x((/theta_iter,zeta/),0,LA_s)  !theta_iter+lambda
    FRdFR(2)=1.0_wp+LA_base_r%f%evalDOF_x((/theta_iter,zeta/),DERIV_THET,LA_s) !1+dlambda/dtheta
  END FUNCTION FRdFR

END SUBROUTINE gvec_to_gene_coords

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_metrics(nthet,nzeta,spos_in,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
! MODULES
USE MODgvec_ReadState_Vars
USE MODgvec_globals, ONLY: PI,CROSS
USE MODgvec_Newton,  ONLY: NewtonRoot1D_FdF
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER              :: nthet          !! number of points in theta_star
INTEGER              :: nzeta          !! number of points in zeta
REAL(wp),INTENT( IN) :: spos_in        !! radial position (sqrt(phi_norm)), phi_norm: normalized toroidal flux [0,1]
REAL(wp),INTENT( IN) :: theta_star_in(nthet,nzeta)  !! thetaStar poloidal angle (straight field line angle PEST)
REAL(wp),INTENT( IN) :: zeta_in(      nthet,nzeta)  !! zeta toroidal angle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: grad_s(         3,nthet,nzeta)  !! gradient in cartesian space, of the radial coordinate
REAL(wp),INTENT(OUT) :: grad_theta_star(3,nthet,nzeta)  !! gradient in cartesian space, of the theta_star coordinate
REAL(wp),INTENT(OUT) :: grad_zeta(      3,nthet,nzeta)  !! gradient in cartesian space, of the zeta coordinate
REAL(wp),INTENT(OUT) :: Bfield(         3,nthet,nzeta)  !! magnetic field in cartesian space
REAL(wp),INTENT(OUT) :: grad_absB(      3,nthet,nzeta)  !! gradient in cartesian space, of the magnetic field magnitude |B|
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iMode,ithet,izeta
REAL(wp) :: iota_int,PhiPrime_int,theta_star,theta,zeta
REAL(wp) :: iota_int_eps,PhiPrime_int_eps
REAL(wp) :: xp(2),qvec(3)
REAL(wp),DIMENSION(1:X1_base_r%f%modes) :: X1_s,dX1ds_s,X1_s_eps,dX1ds_s_eps
REAL(wp),DIMENSION(1:X2_base_r%f%modes) :: X2_s,dX2ds_s,X2_s_eps,dX2ds_s_eps
REAL(wp),DIMENSION(1:LA_base_r%f%modes) :: LA_s,dLAds_s,LA_s_eps
REAL(wp) :: X1_int,dX1ds,dX1dthet,dX1dzeta
REAL(wp) :: X2_int,dX2ds,dX2dthet,dX2dzeta
REAL(wp) :: dLAds,dLAdthet,dLAdzeta
REAL(wp) :: e_s(3),e_thet(3),e_zeta(3),grad_thet(3)
REAL(wp) :: sqrtG,spos
REAL(wp) :: absB,absB_ds,absB_dthet,absB_dzeta
REAL(wp) :: eps=1.0e-08
!===================================================================================================================================
spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,spos_in)) !for satefy reasons at the axis and edge
!interpolate first in s direction

DO iMode=1,X1_base_r%f%modes
  X1_s(   iMode) =X1_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode))
  dX1ds_s(iMode) =X1_base_r%s%evalDOF_s(spos,DERIV_S,X1_r(:,iMode))
END DO
DO iMode=1,X2_base_r%f%modes
  X2_s(   iMode) =X2_base_r%s%evalDOF_s(spos,      0,X2_r(:,iMode))
  dX2ds_s(iMode) =X2_base_r%s%evalDOF_s(spos,DERIV_S,X2_r(:,iMode))
END DO
DO iMode=1,LA_base_r%f%modes
  LA_s(   iMode) =LA_base_r%s%evalDOF_s(spos,      0,LA_r(:,iMode))
  dLAds_s(iMode) =LA_base_r%s%evalDOF_s(spos,DERIV_S,LA_r(:,iMode))
END DO

iota_int     = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))
PhiPrime_int = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))

!for FD in s-direction
DO iMode=1,X1_base_r%f%modes
  X1_s_eps(   iMode)  =X1_base_r%s%evalDOF_s(spos+eps,      0,X1_r(:,iMode))
  dX1ds_s_eps(iMode)  =X1_base_r%s%evalDOF_s(spos+eps,DERIV_S,X1_r(:,iMode))
END DO
DO iMode=1,X2_base_r%f%modes
  X2_s_eps(   iMode)  =X2_base_r%s%evalDOF_s(spos+eps,      0,X2_r(:,iMode))
  dX2ds_s_eps(iMode)  =X2_base_r%s%evalDOF_s(spos+eps,DERIV_S,X2_r(:,iMode))
END DO
DO iMode=1,LA_base_r%f%modes
  LA_s_eps(   iMode)  =LA_base_r%s%evalDOF_s(spos+eps,      0,LA_r(:,iMode))
END DO

iota_int_eps     = X1_base_r%s%evalDOF_s(spos+eps, 0,profiles_1d(:,3))
PhiPrime_int_eps = X1_base_r%s%evalDOF_s(spos+eps, DERIV_S ,profiles_1d(:,1))

DO izeta=1,nzeta; DO ithet=1,nthet
  theta_star = theta_star_in(ithet,izeta) !theta_star depends on zeta!!
  zeta = zeta_in(ithet,izeta)
  !find angle theta from straight field line angle (PEST) theta_star=theta+lambda(s,theta,zeta) 
  ! 1D Newton uses derivative function FRdFR defined below... solves FR(1)-F0=0, FR(2)=dFR(1)/dtheta
  !                      tolerance , lower bound , upper bound , maxstep, start value , F0
  theta=NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1*PI  ,theta_star   , theta_star,FRdFR)

  xp=(/theta,zeta/)

  X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
  dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
  dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
  dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

  X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
  dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
  dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
  dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

  dLAds    =LA_base_r%f%evalDOF_x(xp,          0,dLAds_s)
  dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
  dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

  qvec=(/X1_int,X2_int,zeta/)
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,0.0_wp/))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 

  grad_s(:,ithet,izeta)   = CROSS(e_thet,e_zeta)/sqrtG
  grad_thet(:)            = CROSS(e_zeta,e_s   )/sqrtG
  grad_zeta(:,ithet,izeta)= CROSS(e_s   ,e_thet)/sqrtG
  grad_theta_star(:,ithet,izeta)=  grad_s(   :,ithet,izeta)*dLAds         &
                                  +grad_thet(:)            *(1.+dLAdthet) &
                                  +grad_zeta(:,ithet,izeta)*dLAdzeta

  Bfield(:,ithet,izeta)= (  e_thet(:)*(iota_int-dLAdzeta )  & 
                          + e_zeta(:)*(1.0_wp+dLAdthet   ) )*(PhiPrime_int/sqrtG)

  absB=SQRT(SUM((Bfield(:,ithet,izeta))**2))

  !-----------TO COMPUTE grad|B|, we do a finite difference in s,theta,zeta ----------

  !variation of |B| in s coordinate (using _eps variables evaluated at spos+eps, above) 
  xp=(/theta,zeta/)

  X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s_eps  )
  dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s_eps)
  dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s_eps  )
  dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s_eps  )

  X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s_eps  )
  dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s_eps)
  dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s_eps  )
  dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s_eps  )

  dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s_eps)
  dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s_eps)

  qvec   = (/X1_int,X2_int,zeta/)
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 

  absB_ds =(SQRT(SUM(((  e_thet(:)*(iota_int_eps-dLAdzeta )  & 
                       + e_zeta(:)*(1.0_wp      +dLAdthet )  )*(PhiPrime_int_eps/sqrtG))**2)) &
            -absB)

  !variation of |B| in theta 

  xp=(/theta+eps,zeta/)

  X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
  dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
  dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
  dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

  X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
  dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
  dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
  dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

  dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
  dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

  qvec   = (/X1_int,X2_int,zeta/)
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 

  absB_dthet =(SQRT(SUM(((   e_thet(:)*(iota_int-dLAdzeta )  & 
                           + e_zeta(:)*(1.0_wp  +dLAdthet ) )*(PhiPrime_int/sqrtG))**2)) &
               -absB)

  !variation of |B| in zeta 

  xp=(/theta,zeta+eps/)

  X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
  dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
  dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
  dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

  X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
  dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
  dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
  dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

  dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
  dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

  qvec=(/X1_int,X2_int,zeta/)
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 

  absB_dzeta =(SQRT(SUM(((  e_thet(:)*(iota_int-dLAdzeta )  & 
                          + e_zeta(:)*(1.0_wp  +dLAdthet ) )*(PhiPrime_int/sqrtG))**2)) &
               -absB)


  grad_absB(:,ithet,izeta)=( absB_ds   *grad_s(:,ithet,izeta)  &
                            +absB_dthet*grad_thet(:)           &
                            +absB_dzeta*grad_zeta(:,ithet,izeta) )/eps

END DO; END DO !ithet,izeta

!for iteration on theta^*
CONTAINS 

  FUNCTION FRdFR(theta_iter)
    !uses current zeta where newton is called, and LA_s from subroutine above
    IMPLICIT NONE
    REAL(wp) :: theta_iter
    REAL(wp) :: FRdFR(2) !output
    !--------------------------------------------------- 
    FRdFR(1)=theta_iter+LA_base_r%f%evalDOF_x((/theta_iter,zeta/),0,LA_s)  !theta_iter+lambda
    FRdFR(2)=1.0_wp+LA_base_r%f%evalDOF_x((/theta_iter,zeta/),DERIV_THET,LA_s) !1+dlambda/dtheta
  END FUNCTION FRdFR

END SUBROUTINE gvec_to_gene_metrics

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE finalize_gvec_to_gene 
! MODULES
USE MODgvec_readState, ONLY: finalize_readState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL Finalize_ReadState()

END SUBROUTINE finalize_gvec_to_gene

END MODULE MODgvec_gvec_to_gene
