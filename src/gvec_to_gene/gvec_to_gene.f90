!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
USE MODgvec_output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_base, ONLY: t_base
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
PUBLIC::WriteDataToVTK

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE init_gvec_to_gene(fileName,SFLcoord_in,factorSFL_in)
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_ReadState,ONLY: ReadState,eval_phiPrime_r,eval_iota_r
USE MODgvec_ReadState_Vars,ONLY: hmap_r,X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r
USE MODgvec_transform_sfl     ,ONLY: transform_sfl_new
USE MODgvec_gvec_to_gene_vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileName !< name of GVEC file
INTEGER,INTENT(IN),OPTIONAL :: SFLcoord_in !< type of straight field line coordinate (0: 'old way' PEST, 1: PEST, 2: BOOZER)
INTEGER,INTENT(IN),OPTIONAL :: factorSFL_in !< factor on fourier modes for SFL>0
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: mn_max(2)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT EVAL GVEC ...'

  IF(PRESENT(SFLcoord_in))THEN
    SFLcoord=SFLcoord_in
  ELSE
    SFLcoord=0
  END IF

  CALL ReadState(fileName)

  IF(SFLcoord.NE.0)THEN
    IF(PRESENT(factorSFL_in))THEN
      factorSFL=factorSFL_in
    ELSE
      factorSFL = 4  !!DEFAULT, added as optional input argument!
    END IF

    mn_max(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
    mn_max(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))

    !mn_max       = mn_max*factorSFL !*SFLfactor on modes of GVEC solution

    IF(mn_max(2).GT.0)THEN
      !IMPORTANT FIX: same maximal number of modes in both directions, since SFL coordinates are coupled (if not a tokamak)!
      mn_max       = MAXVAL(mn_max)*factorSFL !*SFLfactor on modes of GVEC solution
    END IF
    CALL transform_sfl_new(trafoSFL,mn_max,SFLcoord,X1_base_r%s%deg,X1_base_r%s%continuity, &
                           X1_base_r%s%degGP,X1_base_r%s%grid ,hmap_r,X1_base_r,X2_base_r,LA_base_r,eval_phiPrime_r,eval_iota_r)
    CALL trafoSFL%buildTransform(X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r)
  END IF


  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE init_gvec_to_gene


!===================================================================================================================================
!> Scalar variables of the equilibrium
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_scalars(Fa,minor_r,PhiPrime_edge,q_edge,n0_global,major_R)
! MODULES
USE MODgvec_globals,ONLY: TWOPI
USE MODgvec_ReadState_Vars,ONLY: a_minor,X1_base_r,sbase_prof,profiles_1d,r_major
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
REAL(wp),INTENT(OUT),OPTIONAL :: major_R  !! major radius
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Fa = TWOPI*sbase_prof%evalDOF_s(1.0, 0,profiles_1d(:,1)) !phi(s=1)
minor_r   = a_minor
n0_global = X1_base_r%f%nfp
PhiPrime_edge = sbase_prof%evalDOF_s(1.0, DERIV_S,profiles_1d(:,1))
q_edge    = 1./(sbase_prof%evalDOF_s(1.0,       0,profiles_1d(:,3)) ) !q=1/iota
IF(PRESENT(major_R)) major_R=r_major

END SUBROUTINE gvec_to_gene_scalars


!===================================================================================================================================
!> Evaluate only s dependend variables
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_profile(spos,q,q_prime,p,p_prime)
! MODULES
USE MODgvec_ReadState_Vars,ONLY: profiles_1d,sbase_prof
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
IF(PRESENT(q))       q       = 1./(  sbase_prof%evalDOF_s(spos,       0,profiles_1d(:,3)) ) !q=1/iota
IF(PRESENT(q_prime)) q_prime = -q*q*(sbase_prof%evalDOF_s(spos, DERIV_S,profiles_1d(:,3)) ) !q'=-iota'/iota^2
IF(PRESENT(p))             p =      (sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,4)) ) !pressure
IF(PRESENT(p_prime)) p_prime =      (sbase_prof%evalDOF_s(spos, DERIV_S,profiles_1d(:,4)) ) !pressure'

END SUBROUTINE gvec_to_gene_profile

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_coords(nthet,nzeta,spos_in,theta_star_in,zeta_in,theta_out,cart_coords)
! MODULES
USE MODgvec_gvec_to_gene_vars,ONLY:SFLcoord
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
!===================================================================================================================================
SELECT CASE(SFLcoord)
CASE(0) ! 'old way' of PEST
  CALL gvec_to_gene_coords_old(nthet,nzeta,spos_in,theta_star_in,zeta_in,theta_out,cart_coords)
CASE(1,2)  !SFL of PEST or BOOZER, fourier angles are already SFL angles!
  CALL gvec_to_gene_coords_sfl(nthet,nzeta,spos_in,theta_star_in,zeta_in,cart_coords)
  theta_out=0.  !'GVEC THETA ANGLE'!
CASE DEFAULT
  STOP 'SFLcoord not known'
END SELECT !SFLcoord

END SUBROUTINE gvec_to_gene_coords

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_coords_old(nthet,nzeta,spos_in,theta_star_in,zeta_in,theta_out,cart_coords)
! MODULES
USE MODgvec_ReadState_Vars
USE MODgvec_globals, ONLY: PI
USE MODgvec_Transform_SFL, ONLY: get_pest_newton
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
INTEGER     :: ithet,izeta
REAL(wp)    :: theta_star,zeta
REAL(wp)    :: xp(2),qvec(3)
REAL(wp)    :: X1_s(   1:X1_base_r%f%modes)
REAL(wp)    :: X2_s(   1:X2_base_r%f%modes)
REAL(wp)    :: LA_s(   1:LA_base_r%f%modes)
REAL(wp)    :: X1_int,X2_int,spos
!===================================================================================================================================
spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,spos_in)) !for satefy reasons at the axis and edge
!interpolate first in s direction
X1_s(:)      =X1_base_r%s%evalDOF2D_s(spos,X1_base_r%f%modes,      0,X1_r(:,:)) !R
X2_s(:)      =X2_base_r%s%evalDOF2D_s(spos,X2_base_r%f%modes,      0,X2_r(:,:)) !Z
LA_s(:)      =LA_base_r%s%evalDOF2D_s(spos, LA_base_r%f%modes,     0,LA_r(:,:)) !lambda

DO izeta=1,nzeta; DO ithet=1,nthet
  theta_star = theta_star_in(ithet,izeta) !theta_star depends on zeta!!
  zeta       = zeta_in(      ithet,izeta)
  !find angle theta from straight field line angle (PEST) theta_star=theta+lambda(s,theta,zeta)
  theta_out(ithet,izeta)=get_pest_newton(theta_star, zeta, LA_base_r%f, LA_s)

  xp=(/theta_out(ithet,izeta),zeta/)

  X1_int      = X1_base_r%f%evalDOF_x(xp,0,X1_s)
  X2_int      = X2_base_r%f%evalDOF_x(xp,0,X2_s)

  qvec = (/X1_int,X2_int,zeta/)
  cart_coords(:,ithet,izeta)=hmap_r%eval(qvec)

END DO; END DO !ithet,izeta
END SUBROUTINE gvec_to_gene_coords_old

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_coords_sfl(nthet,nzeta,spos_in,theta_star_in,zeta_in,cart_coords)
! MODULES
USE MODgvec_ReadState_Vars,ONLY: hmap_r
USE MODgvec_gvec_to_gene_vars,ONLY:SFLcoord
USE MODgvec_gvec_to_gene_vars,ONLY:trafoSFL
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
REAL(wp),INTENT(OUT) :: cart_coords(3,nthet,nzeta)  !! x,y,z cartesian coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: ithet,izeta
REAL(wp)    :: zeta,zetastar
REAL(wp)    :: xp(2),qvec(3)
REAL(wp)    :: X1sfl_s(   1:trafoSFL%X1sfl_base%f%modes)
REAL(wp)    :: X2sfl_s(   1:trafoSFL%X2sfl_base%f%modes)
REAL(wp),ALLOCATABLE :: GZsfl_s(:)
REAL(wp)    :: X1_int,X2_int,spos
!===================================================================================================================================
ASSOCIATE(X1sfl_base=>trafoSFL%X1sfl_base,X1sfl=>trafoSFL%X1sfl,&
          X2sfl_base=>trafoSFL%X2sfl_base,X2sfl=>trafoSFL%X2sfl,&
          GZsfl_base=>trafoSFL%GZsfl_base,GZsfl=>trafoSFL%GZsfl )
spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,spos_in)) !for satefy reasons at the axis and edge
X1sfl_s(:)      =X1sfl_base%s%evalDOF2D_s(spos,X1sfl_base%f%modes,      0,X1sfl(:,:)) !R
X2sfl_s(:)      =X2sfl_base%s%evalDOF2D_s(spos,X2sfl_base%f%modes,      0,X2sfl(:,:)) !Z
IF(SFLcoord.EQ.2)THEN !BOOZER
  ALLOCATE(GZsfl_s(1:trafoSFL%GZsfl_base%f%modes))
  GZsfl_s(:)      =GZsfl_base%s%evalDOF2D_s(spos,GZsfl_base%f%modes,      0,GZsfl(:,:))
END IF

DO izeta=1,nzeta; DO ithet=1,nthet
  zetastar=zeta_in( ithet,izeta) !=zetastar = zeta+Gsfl(thetastar,zetastar)
  xp=(/theta_star_in(ithet,izeta),zetastar/)

  X1_int      = X1sfl_base%f%evalDOF_x(xp,0,X1sfl_s)
  X2_int      = X2sfl_base%f%evalDOF_x(xp,0,X2sfl_s)
  IF(SFLcoord.EQ.2)THEN !BOOZER
    zeta = zetastar- GZsfl_base%f%evalDOF_x(xp,0,GZsfl_s)
  ELSE  !PEST / GVEC
    zeta = zetastar
  END IF

  qvec = (/X1_int,X2_int,zeta/)
  cart_coords(:,ithet,izeta)=hmap_r%eval(qvec)

END DO; END DO !ithet,izeta

END ASSOCIATE !X1sfl,X2sfl

END SUBROUTINE gvec_to_gene_coords_sfl

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_metrics(nthet,nzeta,spos_in,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
! MODULES
USE MODgvec_gvec_to_gene_vars,ONLY:SFLcoord
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
!===================================================================================================================================
SELECT CASE(SFLcoord)
CASE(0) ! 'old way' of PEST
  CALL gvec_to_gene_metrics_old(nthet,nzeta,spos_in,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
CASE(1,2)  !SFL of PEST or BOOZER, fourier angles are already SFL angles!
  CALL gvec_to_gene_metrics_sfl(nthet,nzeta,spos_in,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
CASE DEFAULT
  STOP 'SFLcoord not known'
END SELECT !SFLcoord
END SUBROUTINE gvec_to_gene_metrics

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_metrics_old(nthet,nzeta,spos_in,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
! MODULES
USE MODgvec_ReadState_Vars
USE MODgvec_globals, ONLY: PI,CROSS
USE MODgvec_Transform_SFL, ONLY: get_pest_newton
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
INTEGER  :: ithet,izeta
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

!! INTERMEDIATE FLAG FOR TESTING THE NEW implementation of the Finite Difference...
!!! #define PP_test_new_delta

#ifdef PP_test_new_delta
INTEGER  :: iElem
REAL(wp),PARAMETER :: eps_s    = 1.0e-4 !
REAL(wp),PARAMETER :: eps_thet = 1.0e-8 !theta,zeta
REAL(wp),PARAMETER :: eps_zeta = 1.0e-8 !theta,zeta
#else
REAL(wp) :: eps=1.0e-08
#endif
REAL(wp) :: delta_s,delta_thet,delta_zeta
REAL(wp) :: g_tt,g_tz,g_zz,g_tt_FD,g_tz_FD,g_zz_FD,Bt,Bz,Bt_FD,Bz_FD,q_thet(3),q_zeta(3)
!===================================================================================================================================
spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,spos_in)) !for satefy reasons at the axis and edge
!interpolate first in s direction

#ifdef PP_test_new_delta
iElem=sgrid_r%find_elem(spos)
delta_s = MERGE(-1,1,(iElem.EQ.sgrid_r%nElems))* eps_s*sgrid_r%ds(iElem)
#else
delta_s=eps
#endif

X1_s(   :)      = X1_base_r%s%evalDOF2D_s(spos        ,X1_base_r%f%modes,       0,X1_r(:,:))
dX1ds_s(:)      = X1_base_r%s%evalDOF2D_s(spos        ,X1_base_r%f%modes, DERIV_S,X1_r(:,:))
X1_s_eps(   :)  = X1_base_r%s%evalDOF2D_s(spos+delta_s,X1_base_r%f%modes,       0,X1_r(:,:))
dX1ds_s_eps(:)  = X1_base_r%s%evalDOF2D_s(spos+delta_s,X1_base_r%f%modes, DERIV_S,X1_r(:,:))

X2_s(   :)      = X2_base_r%s%evalDOF2D_s(spos        ,X2_base_r%f%modes,       0,X2_r(:,:))
dX2ds_s(:)      = X2_base_r%s%evalDOF2D_s(spos        ,X2_base_r%f%modes, DERIV_S,X2_r(:,:))
X2_s_eps(   :)  = X2_base_r%s%evalDOF2D_s(spos+delta_s,X2_base_r%f%modes,       0,X2_r(:,:))
dX2ds_s_eps(:)  = X2_base_r%s%evalDOF2D_s(spos+delta_s,X2_base_r%f%modes, DERIV_S,X2_r(:,:))

LA_s(   :)      = LA_base_r%s%evalDOF2D_s(spos        ,LA_base_r%f%modes,       0,LA_r(:,:))
dLAds_s(:)      = LA_base_r%s%evalDOF2D_s(spos        ,LA_base_r%f%modes, DERIV_S,LA_r(:,:))
LA_s_eps(   :)  = LA_base_r%s%evalDOF2D_s(spos+delta_s,LA_base_r%f%modes,       0,LA_r(:,:))

iota_int     = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,3))
PhiPrime_int = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))

iota_int_eps     = sbase_prof%evalDOF_s(spos+delta_s, 0,profiles_1d(:,3))
PhiPrime_int_eps = sbase_prof%evalDOF_s(spos+delta_s, DERIV_S ,profiles_1d(:,1))

DO izeta=1,nzeta; DO ithet=1,nthet
  theta_star = theta_star_in(ithet,izeta) !theta_star depends on zeta!!
  zeta = zeta_in(ithet,izeta)
  !find angle theta from straight field line angle (PEST) theta_star=theta+lambda(s,theta,zeta)
  theta = get_pest_newton(theta_star, zeta, LA_base_r%f, LA_s)

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

  !absB=SQRT(SUM((Bfield(:,ithet,izeta))**2))


  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp/)

  g_tt        = hmap_r%eval_gij(q_thet,qvec,q_thet)   !g_theta,theta
  g_tz        = hmap_r%eval_gij(q_thet,qvec,q_zeta)   !g_theta,zeta =g_zeta,theta
  g_zz        = hmap_r%eval_gij(q_zeta,qvec,q_zeta)   !g_zeta,zeta
  Bt          = (iota_int-dLAdzeta )*(PhiPrime_int/sqrtG)
  Bz          = (1.0_wp  +dLAdthet )*(PhiPrime_int/sqrtG)
  absB        = SQRT(Bt*(Bt*g_tt+2.0_wp*Bz*g_tz)+Bz*Bz*g_zz)

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

  !absB_ds =(SQRT(SUM(((  e_thet(:)*(iota_int_eps-dLAdzeta )  &
  !                     + e_zeta(:)*(1.0_wp      +dLAdthet )  )*(PhiPrime_int_eps/sqrtG))**2)) &
  !          -absB)

  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp/)

  g_tt_FD  = (hmap_r%eval_gij(q_thet,qvec,q_thet) - g_tt)/delta_s
  g_tz_FD  = (hmap_r%eval_gij(q_thet,qvec,q_zeta) - g_tz)/delta_s
  g_zz_FD  = (hmap_r%eval_gij(q_zeta,qvec,q_zeta) - g_zz)/delta_s

  Bt_FD    = ( (iota_int_eps-dLAdzeta )*(PhiPrime_int_eps/sqrtG) - Bt)/delta_s
  Bz_FD    = ( (1.0_wp      +dLAdthet )*(PhiPrime_int_eps/sqrtG) - Bz)/delta_s

  absB_ds = (2.0_wp*((Bt*g_tt+ Bz*g_tz)*Bt_FD +(Bt*g_tz + Bz*g_zz)*Bz_FD +Bt*Bz*g_tz_FD) +Bt*Bt*g_tt_FD + Bz*Bz*g_zz_FD)! / (2.0_wp*absB)


  !variation of |B| in theta

#ifdef PP_test_new_delta
  delta_thet = eps_thet*SQRT(SUM(grad_thet*grad_thet))
#else
  delta_thet = eps
#endif
  xp=(/theta+delta_thet,zeta/)

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
  !e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  !e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet)

  !absB_dthet =(SQRT(SUM(((   e_thet(:)*(iota_int-dLAdzeta )  &
  !                         + e_zeta(:)*(1.0_wp  +dLAdthet ) )*(PhiPrime_int/sqrtG))**2)) &
  !             -absB)
  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp/)

  g_tt_FD  = (hmap_r%eval_gij(q_thet,qvec,q_thet) - g_tt)/delta_thet
  g_tz_FD  = (hmap_r%eval_gij(q_thet,qvec,q_zeta) - g_tz)/delta_thet
  g_zz_FD  = (hmap_r%eval_gij(q_zeta,qvec,q_zeta) - g_zz)/delta_thet

  Bt_FD    = ( (iota_int-dLAdzeta )*(PhiPrime_int/sqrtG) - Bt)/delta_thet
  Bz_FD    = ( (1.0_wp  +dLAdthet )*(PhiPrime_int/sqrtG) - Bz)/delta_thet

  absB_dthet = (2.0_wp*((Bt*g_tt+ Bz*g_tz)*Bt_FD +(Bt*g_tz + Bz*g_zz)*Bz_FD +Bt*Bz*g_tz_FD) +Bt*Bt*g_tt_FD + Bz*Bz*g_zz_FD)! / (2.0_wp*absB)

  !variation of |B| in zeta

#ifdef PP_test_new_delta
  delta_zeta = eps_zeta*SQRT(SUM(grad_zeta*grad_zeta))
#else
  delta_zeta = eps
#endif
  xp=(/theta,zeta+delta_zeta/)

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

  qvec=(/X1_int,X2_int,xp(2)/) ! USE THE CORRECT PERTURBED ZETA POSITION!  xp(2)=zeta+delta_zeta !!
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
  sqrtG  = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet)

  !absB_dzeta =(SQRT(SUM(((  e_thet(:)*(iota_int-dLAdzeta )  &
  !                        + e_zeta(:)*(1.0_wp  +dLAdthet ) )*(PhiPrime_int/sqrtG))**2)) &
  !             -absB)
  !
  !grad_absB(:,ithet,izeta)=( absB_ds   *grad_s(:,ithet,izeta)/delta_s &
  !                          +absB_dthet*grad_thet(:)/delta_thet          &
  !                          +absB_dzeta*grad_zeta(:,ithet,izeta)/delta_zeta)

  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp/)

  g_tt_FD  = (hmap_r%eval_gij(q_thet,qvec,q_thet) - g_tt)/delta_zeta
  g_tz_FD  = (hmap_r%eval_gij(q_thet,qvec,q_zeta) - g_tz)/delta_zeta
  g_zz_FD  = (hmap_r%eval_gij(q_zeta,qvec,q_zeta) - g_zz)/delta_zeta

  Bt_FD    = ( (iota_int-dLAdzeta )*(PhiPrime_int/sqrtG) - Bt)/delta_zeta
  Bz_FD    = ( (1.0_wp  +dLAdthet )*(PhiPrime_int/sqrtG) - Bz)/delta_zeta

  absB_dzeta = (2.0_wp*((Bt*g_tt+ Bz*g_tz)*Bt_FD +(Bt*g_tz + Bz*g_zz)*Bz_FD +Bt*Bz*g_tz_FD) +Bt*Bt*g_tt_FD + Bz*Bz*g_zz_FD)! / (2.0_wp*absB)

  grad_absB(:,ithet,izeta)=( absB_ds   *grad_s(:,ithet,izeta)   &
                            +absB_dthet*grad_thet(:)            &
                            +absB_dzeta*grad_zeta(:,ithet,izeta))/(2.0_wp*absB)
END DO; END DO !ithet,izeta
END SUBROUTINE gvec_to_gene_metrics_old


!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_metrics_sfl(nthet,nzeta,spos_in,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
! MODULES
USE MODgvec_globals, ONLY: CROSS
USE MODgvec_ReadState_Vars,ONLY: hmap_r,sbase_prof,profiles_1d
USE MODgvec_gvec_to_gene_vars,ONLY:SFLcoord
USE MODgvec_gvec_to_gene_vars,ONLY:trafoSFL
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
INTEGER  :: ithet,izeta
REAL(wp) :: iota_int,PhiPrime_int,theta,zeta
REAL(wp) :: iota_int_eps,PhiPrime_int_eps
REAL(wp) :: xp(2),qvec(3)
REAL(wp),DIMENSION(1:trafoSFL%X1sfl_base%f%modes) :: X1_s,dX1ds_s,X1_s_eps,dX1ds_s_eps
REAL(wp),DIMENSION(1:trafoSFL%X2sfl_base%f%modes) :: X2_s,dX2ds_s,X2_s_eps,dX2ds_s_eps
!REAL(wp),DIMENSION(1:GZsfl_base%f%modes) :: GZ_s,dGZds_s,GZ_s_eps,dGZds_s_eps
REAL(wp),ALLOCATABLE :: GZ_s(:),dGZds_s(:),GZ_s_eps(:),dGZds_s_eps(:)
REAL(wp) :: X1_int,dX1ds,dX1dthet,dX1dzeta
REAL(wp) :: X2_int,dX2ds,dX2dthet,dX2dzeta
REAL(wp) :: GZ_int,dGZds,dGZdthet,dGZdzeta
REAL(wp) :: e_s(3),e_thet(3),e_zeta(3),grad_thet(3)
REAL(wp) :: sqrtG,spos
REAL(wp) :: absB,absB_ds,absB_dthet,absB_dzeta
#ifdef PP_test_new_delta
  LOGICAL  :: test_new_delta=.TRUE.  !!! False: old style FD using eps only, True: scaled delta for FD
  INTEGER  :: iElem
  REAL(wp),PARAMETER :: eps_s    = 1.0e-4 !
  REAL(wp),PARAMETER :: eps_thet = 1.0e-8 !theta,zeta
  REAL(wp),PARAMETER :: eps_zeta = 1.0e-8 !theta,zeta
#else
REAL(wp) :: eps=1.0e-08
#endif
REAL(wp) :: delta_s,delta_thet,delta_zeta
REAL(wp) :: g_tt,g_tz,g_zz,g_tt_FD,g_tz_FD,g_zz_FD,Bt,Bz,Bt_FD,Bz_FD,q_thet(3),q_zeta(3)
!===================================================================================================================================

spos=MAX(1.0e-08_wp,MIN(1.0_wp-1.0e-12_wp,spos_in)) !for satefy reasons at the axis and edge

!!!! new delta in s-direction, scaled with the size of the element, change direction of FD in the last element!
#ifdef PP_test_new_delta
iElem=sgrid_sfl%find_elem(spos)
delta_s = MERGE(-1,1,(iElem.EQ.sgrid_sfl%nElems))* eps_s*sgrid_sfl%ds(iElem)
#else
delta_s=eps
#endif

!interpolate first in s direction
ASSOCIATE(X1sfl_base=>trafoSFL%X1sfl_base,X1sfl=>trafoSFL%X1sfl,&
          X2sfl_base=>trafoSFL%X2sfl_base,X2sfl=>trafoSFL%X2sfl,&
          GZsfl_base=>trafoSFL%GZsfl_base,GZsfl=>trafoSFL%GZsfl )

X1_s(   :)      = X1sfl_base%s%evalDOF2D_s(spos        ,X1sfl_base%f%modes,       0,X1sfl(:,:))
dX1ds_s(:)      = X1sfl_base%s%evalDOF2D_s(spos        ,X1sfl_base%f%modes, DERIV_S,X1sfl(:,:))
X1_s_eps(   :)  = X1sfl_base%s%evalDOF2D_s(spos+delta_s,X1sfl_base%f%modes,       0,X1sfl(:,:))
dX1ds_s_eps(:)  = X1sfl_base%s%evalDOF2D_s(spos+delta_s,X1sfl_base%f%modes, DERIV_S,X1sfl(:,:))

X2_s(   :)      = X2sfl_base%s%evalDOF2D_s(spos        ,X2sfl_base%f%modes,       0,X2sfl(:,:))
dX2ds_s(:)      = X2sfl_base%s%evalDOF2D_s(spos        ,X2sfl_base%f%modes, DERIV_S,X2sfl(:,:))
X2_s_eps(   :)  = X2sfl_base%s%evalDOF2D_s(spos+delta_s,X2sfl_base%f%modes,       0,X2sfl(:,:))
dX2ds_s_eps(:)  = X2sfl_base%s%evalDOF2D_s(spos+delta_s,X2sfl_base%f%modes, DERIV_S,X2sfl(:,:))

IF(SFLcoord.EQ.2)THEN !BOOZER
  ALLOCATE(GZ_s(       1:GZsfl_base%f%modes))
  ALLOCATE(dGZds_s(    1:GZsfl_base%f%modes))
  ALLOCATE(GZ_s_eps(   1:GZsfl_base%f%modes))
  ALLOCATE(dGZds_s_eps(1:GZsfl_base%f%modes))
  GZ_s(   :)      = GZsfl_base%s%evalDOF2D_s(spos        ,GZsfl_base%f%modes,       0,GZsfl(:,:))
  dGZds_s(:)      = GZsfl_base%s%evalDOF2D_s(spos        ,GZsfl_base%f%modes, DERIV_S,GZsfl(:,:))
  GZ_s_eps(   :)  = GZsfl_base%s%evalDOF2D_s(spos+delta_s,GZsfl_base%f%modes,       0,GZsfl(:,:))
  dGZds_s_eps(:)  = GZsfl_base%s%evalDOF2D_s(spos+delta_s,GZsfl_base%f%modes, DERIV_S,GZsfl(:,:))
END IF

iota_int     = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,3))
PhiPrime_int = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))

iota_int_eps     = sbase_prof%evalDOF_s(spos+delta_s, 0,profiles_1d(:,3))
PhiPrime_int_eps = sbase_prof%evalDOF_s(spos+delta_s, DERIV_S ,profiles_1d(:,1))

!FOR SFLcoord=1
GZ_int   =0.
dGZds    =0.
dGZdthet =0.
dGZdzeta =0.

DO izeta=1,nzeta; DO ithet=1,nthet
  theta = theta_star_in(ithet,izeta) !theta_star is the SFL angle
  zeta = zeta_in(ithet,izeta)

  xp=(/theta,zeta/)

  X1_int  =X1sfl_base%f%evalDOF_x(xp,          0, X1_s  )
  dX1ds   =X1sfl_base%f%evalDOF_x(xp,          0,dX1ds_s)
  dX1dthet=X1sfl_base%f%evalDOF_x(xp, DERIV_THET, X1_s  )
  dX1dzeta=X1sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

  X2_int  =X2sfl_base%f%evalDOF_x(xp,          0, X2_s  )
  dX2ds   =X2sfl_base%f%evalDOF_x(xp,          0,dX2ds_s)
  dX2dthet=X2sfl_base%f%evalDOF_x(xp, DERIV_THET, X2_s  )
  dX2dzeta=X2sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

  IF(SFLcoord.EQ.2)THEN !BOOZER
    GZ_int   =GZsfl_base%f%evalDOF_x(xp,          0, GZ_s)
    dGZds    =GZsfl_base%f%evalDOF_x(xp,          0,dGZds_s)
    dGZdthet =GZsfl_base%f%evalDOF_x(xp, DERIV_THET, GZ_s)
    dGZdzeta =GZsfl_base%f%evalDOF_x(xp, DERIV_ZETA, GZ_s)
  END IF

  qvec=(/X1_int,X2_int,zeta-GZ_int/)
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,      -dGZds/))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,      -dGZdthet/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/))
  sqrtG  = SUM(e_s*CROSS(e_thet,e_zeta))

  grad_s(:,ithet,izeta)   = CROSS(e_thet,e_zeta)/sqrtG
  grad_thet(:)            = CROSS(e_zeta,e_s   )/sqrtG
  grad_zeta(:,ithet,izeta)= CROSS(e_s   ,e_thet)/sqrtG

  !!!! theta=theta_sfl!
  grad_theta_star(:,ithet,izeta)= grad_thet(:)

  Bfield(:,ithet,izeta)= ( e_thet(:)*iota_int + e_zeta(:) )*(PhiPrime_int/sqrtG)

  !absB=SQRT(SUM((Bfield(:,ithet,izeta))**2))

  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp-dGZdthet/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/)

  g_tt        = hmap_r%eval_gij(q_thet,qvec,q_thet)   !g_theta,theta
  g_tz        = hmap_r%eval_gij(q_thet,qvec,q_zeta)   !g_theta,zeta =g_zeta,theta
  g_zz        = hmap_r%eval_gij(q_zeta,qvec,q_zeta)   !g_zeta,zeta
  Bt          = iota_int*(PhiPrime_int/sqrtG)
  Bz          = 1.0_wp  *(PhiPrime_int/sqrtG)
  absB        = SQRT(Bt*(Bt*g_tt+2.0_wp*Bz*g_tz)+Bz*Bz*g_zz)


  !-----------TO COMPUTE grad|B|, we do a finite difference in s,theta,zeta ----------

  !variation of |B| in s coordinate (using _eps variables evaluated at spos+eps, above)
  xp=(/theta,zeta/)

  X1_int  =X1sfl_base%f%evalDOF_x(xp,          0, X1_s_eps  )
  dX1ds   =X1sfl_base%f%evalDOF_x(xp,          0,dX1ds_s_eps)
  dX1dthet=X1sfl_base%f%evalDOF_x(xp, DERIV_THET, X1_s_eps  )
  dX1dzeta=X1sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X1_s_eps  )

  X2_int  =X2sfl_base%f%evalDOF_x(xp,          0, X2_s_eps  )
  dX2ds   =X2sfl_base%f%evalDOF_x(xp,          0,dX2ds_s_eps)
  dX2dthet=X2sfl_base%f%evalDOF_x(xp, DERIV_THET, X2_s_eps  )
  dX2dzeta=X2sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X2_s_eps  )

  IF(SFLcoord.EQ.2)THEN !BOOZER
    GZ_int   =GZsfl_base%f%evalDOF_x(xp,          0, GZ_s_eps)
    dGZds    =GZsfl_base%f%evalDOF_x(xp,          0,dGZds_s_eps)
    dGZdthet =GZsfl_base%f%evalDOF_x(xp, DERIV_THET, GZ_s_eps)
    dGZdzeta =GZsfl_base%f%evalDOF_x(xp, DERIV_ZETA, GZ_s_eps)
  END IF

  qvec=(/X1_int,X2_int,zeta-GZ_int/)
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,      -dGZds/))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,      -dGZdthet/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/))
  sqrtG = SUM(e_s*CROSS(e_thet,e_zeta))

  !absB_ds =(SQRT(SUM(((e_thet(:)*iota_int_eps + e_zeta(:))*(PhiPrime_int_eps/sqrtG))**2)) &
  !          -absB)

  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp-dGZdthet/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/)

  g_tt_FD  = (hmap_r%eval_gij(q_thet,qvec,q_thet) - g_tt)/delta_s
  g_tz_FD  = (hmap_r%eval_gij(q_thet,qvec,q_zeta) - g_tz)/delta_s
  g_zz_FD  = (hmap_r%eval_gij(q_zeta,qvec,q_zeta) - g_zz)/delta_s

  Bt_FD    = (iota_int_eps*(PhiPrime_int_eps/sqrtG) - Bt)/delta_s
  Bz_FD    = (1.0_wp      *(PhiPrime_int_eps/sqrtG) - Bz)/delta_s

  absB_ds = (2.0_wp*((Bt*g_tt+ Bz*g_tz)*Bt_FD +(Bt*g_tz + Bz*g_zz)*Bz_FD +Bt*Bz*g_tz_FD) +Bt*Bt*g_tt_FD + Bz*Bz*g_zz_FD)! / (2.0_wp*absB)

  !variation of |B| in theta

#ifdef PP_test_new_delta
  delta_thet = eps_thet*SQRT(SUM(grad_thet*grad_thet))
#else
  delta_thet = eps
#endif
  xp=(/theta+delta_thet,zeta/)

  X1_int  =X1sfl_base%f%evalDOF_x(xp,          0, X1_s  )
  dX1ds   =X1sfl_base%f%evalDOF_x(xp,          0,dX1ds_s)
  dX1dthet=X1sfl_base%f%evalDOF_x(xp, DERIV_THET, X1_s  )
  dX1dzeta=X1sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

  X2_int  =X2sfl_base%f%evalDOF_x(xp,          0, X2_s  )
  dX2ds   =X2sfl_base%f%evalDOF_x(xp,          0,dX2ds_s)
  dX2dthet=X2sfl_base%f%evalDOF_x(xp, DERIV_THET, X2_s  )
  dX2dzeta=X2sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

  IF(SFLcoord.EQ.2)THEN
    GZ_int   =GZsfl_base%f%evalDOF_x(xp,          0, GZ_s)
    dGZds    =GZsfl_base%f%evalDOF_x(xp,          0,dGZds_s)
    dGZdthet =GZsfl_base%f%evalDOF_x(xp, DERIV_THET, GZ_s)
    dGZdzeta =GZsfl_base%f%evalDOF_x(xp, DERIV_ZETA, GZ_s)
  END IF

  qvec=(/X1_int,X2_int,zeta-GZ_int/)
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,      -dGZds/))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,      -dGZdthet/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/))
  sqrtG  = SUM(e_s*CROSS(e_thet,e_zeta))

  !absB_dthet =(SQRT(SUM(((e_thet(:)*iota_int + e_zeta(:))*(PhiPrime_int/sqrtG))**2)) &
  !             -absB)

  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp-dGZdthet/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/)

  g_tt_FD  = (hmap_r%eval_gij(q_thet,qvec,q_thet) - g_tt)/delta_thet
  g_tz_FD  = (hmap_r%eval_gij(q_thet,qvec,q_zeta) - g_tz)/delta_thet
  g_zz_FD  = (hmap_r%eval_gij(q_zeta,qvec,q_zeta) - g_zz)/delta_thet

  Bt_FD    = ( iota_int*(PhiPrime_int/sqrtG) - Bt)/delta_thet
  Bz_FD    = ( 1.0_wp  *(PhiPrime_int/sqrtG) - Bz)/delta_thet

  absB_dthet = (2.0_wp*((Bt*g_tt+ Bz*g_tz)*Bt_FD +(Bt*g_tz + Bz*g_zz)*Bz_FD +Bt*Bz*g_tz_FD) +Bt*Bt*g_tt_FD + Bz*Bz*g_zz_FD)! / (2.0_wp*absB)

  !variation of |B| in zeta

#ifdef PP_test_new_delta
  delta_zeta = eps_zeta*SQRT(SUM(grad_zeta*grad_zeta))
#else
  delta_zeta = eps
#endif
  xp=(/theta,zeta+delta_zeta/)

  X1_int  =X1sfl_base%f%evalDOF_x(xp,          0, X1_s  )
  dX1ds   =X1sfl_base%f%evalDOF_x(xp,          0,dX1ds_s)
  dX1dthet=X1sfl_base%f%evalDOF_x(xp, DERIV_THET, X1_s  )
  dX1dzeta=X1sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

  X2_int  =X2sfl_base%f%evalDOF_x(xp,          0, X2_s  )
  dX2ds   =X2sfl_base%f%evalDOF_x(xp,          0,dX2ds_s)
  dX2dthet=X2sfl_base%f%evalDOF_x(xp, DERIV_THET, X2_s  )
  dX2dzeta=X2sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

  IF(SFLcoord.EQ.2)THEN
    GZ_int   =GZsfl_base%f%evalDOF_x(xp,          0, GZ_s)
    dGZds    =GZsfl_base%f%evalDOF_x(xp,          0,dGZds_s)
    dGZdthet =GZsfl_base%f%evalDOF_x(xp, DERIV_THET, GZ_s)
    dGZdzeta =GZsfl_base%f%evalDOF_x(xp, DERIV_ZETA, GZ_s)
  END IF

  qvec=(/X1_int,X2_int,xp(2)-GZ_int/) ! USE THE CORRECT PERTURBED ZETA POSITION!  xp(2)=zeta+delta_zeta !!
  e_s    = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,      -dGZds/))
  e_thet = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,      -dGZdthet/))
  e_zeta = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/))
  sqrtG  = SUM(e_s*CROSS(e_thet,e_zeta))

  !absB_dzeta =(SQRT(SUM(((e_thet(:)*iota_int + e_zeta(:))*(PhiPrime_int/sqrtG))**2)) &
  !             -absB)


  !grad_absB(:,ithet,izeta)=( absB_ds   *grad_s(:,ithet,izeta)/delta_s  &
  !                          +absB_dthet*grad_thet(:)/delta_thet           &
  !                          +absB_dzeta*grad_zeta(:,ithet,izeta)/delta_zeta )

  q_thet(1:3) = (/dX1dthet,dX2dthet,0.0_wp-dGZdthet/)
  q_zeta(1:3) = (/dX1dzeta,dX2dzeta,1.0_wp-dGZdzeta/)

  g_tt_FD  = (hmap_r%eval_gij(q_thet,qvec,q_thet) - g_tt)/delta_zeta
  g_tz_FD  = (hmap_r%eval_gij(q_thet,qvec,q_zeta) - g_tz)/delta_zeta
  g_zz_FD  = (hmap_r%eval_gij(q_zeta,qvec,q_zeta) - g_zz)/delta_zeta

  Bt_FD    = ( iota_int*(PhiPrime_int/sqrtG) - Bt)/delta_zeta
  Bz_FD    = ( 1.0_wp  *(PhiPrime_int/sqrtG) - Bz)/delta_zeta

  absB_dzeta = (2.0_wp*((Bt*g_tt+ Bz*g_tz)*Bt_FD +(Bt*g_tz + Bz*g_zz)*Bz_FD +Bt*Bz*g_tz_FD) +Bt*Bt*g_tt_FD + Bz*Bz*g_zz_FD)! / (2.0_wp*absB)

  grad_absB(:,ithet,izeta)=( absB_ds   *grad_s(:,ithet,izeta)   &
                            +absB_dthet*grad_thet(:)            &
                            +absB_dzeta*grad_zeta(:,ithet,izeta))/(2.0_wp*absB)
END DO; END DO !ithet,izeta

IF(SFLcoord.EQ.2)THEN !BOOZER
  DEALLOCATE(GZ_s,dGZds_s,GZ_s_eps,dGZds_s_eps)
END IF
END ASSOCIATE !X1sfl_base,X1sfl,X2sfl_base,X2sfl,GZsfl_base,GZsfl

END SUBROUTINE gvec_to_gene_metrics_sfl

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE finalize_gvec_to_gene
! MODULES
USE MODgvec_readState, ONLY: finalize_readState
USE MODgvec_gvec_to_gene_vars, ONLY: SFLcoord,trafoSFL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL Finalize_ReadState()
  IF(SFLcoord.NE.0) THEN
    CALL trafoSFL%free()
    DEALLOCATE(trafoSFL)
  END IF

END SUBROUTINE finalize_gvec_to_gene

END MODULE MODgvec_gvec_to_gene
