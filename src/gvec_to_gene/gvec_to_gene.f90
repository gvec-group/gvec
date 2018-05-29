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
USE MODgvec_gvec_to_gene_Vars
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
!> read an input solution and initialize U(0) (X1,X2,LA) of size X1/X2/LA_base , from an ascii .dat file 
!! if size of grid/X1/X2/LA  not equal X1/X2/X3_base
!! interpolate readin solution to the current base of Uin
!!
!===================================================================================================================================
SUBROUTINE ReadState(fileString)
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MODgvec_gvec_to_gene_Vars
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base, base_new
USE MODgvec_fbase,  ONLY: sin_cos_map 
USE MODgvec_hmap,  ONLY: hmap_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=*)    , INTENT(IN   ) :: fileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: fileID_r,OutputLevel_r
  INTEGER              :: ioUnit,iMode,is,nElems_r,grid_type_r,nfp_r,degGP_r,mn_nyq_r(2),which_hmap_r 
  INTEGER              :: X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  INTEGER              :: X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  INTEGER              :: LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r
  INTEGER,ALLOCATABLE  :: X1_mn_r(:,:),X2_mn_r(:,:),LA_mn_r(:,:)
  REAL(wp),ALLOCATABLE :: sp_r(:),profiles_IP(:,:)
  INTEGER              :: X1_mn_max_r(2),X2_mn_max_r(2),LA_mn_max_r(2)
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)')'   READ SOLUTION VARIABLE FROM FILE    "'//TRIM(FileString)//'" ...'

  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit         ,&
     FILE     = TRIM(FileString) ,&
     STATUS   = 'OLD'            ,&
     ACTION   = 'READ'           ,&
     ACCESS   = 'SEQUENTIAL' ) 

  READ(ioUnit,*) !## MHD3D Solution file
  READ(ioUnit,*) outputLevel_r,fileID_r
  READ(ioUnit,*) !## grid: nElems, grid_type
  READ(ioUnit,*) nElems_r,grid_type_r
  ALLOCATE(sp_r(0:nElems_r))

  READ(ioUnit,*) !## grid: sp(0:nElems)
  READ(ioUnit,*)sp_r(:)
  READ(ioUnit,*) !## global: nfp, degGP, mn_nyq
  READ(ioUnit,*) nfp_r, degGP_r,mn_nyq_r,which_hmap_r
  READ(ioUnit,*) !## X1_base: 
  READ(ioUnit,*) X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  READ(ioUnit,*) !## X2_base:                 
  READ(ioUnit,*) X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  READ(ioUnit,*) !## LA_base:                 
  READ(ioUnit,*) LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r

  ALLOCATE(X1_r(1:X1_nbase_r,1:X1_modes_r))
  ALLOCATE(X2_r(1:X2_nbase_r,1:X2_modes_r))
  ALLOCATE(LA_r(1:LA_nbase_r,1:LA_modes_r))

  ALLOCATE(profiles_IP(1:X1_nbase_r,5))
  ALLOCATE(X1_mn_r(2,1:X1_modes_r))
  ALLOCATE(X2_mn_r(2,1:X2_modes_r))
  ALLOCATE(LA_mn_r(2,1:LA_modes_r))
  READ(ioUnit,*) !## X1: 
  DO iMode=1,X1_modes_r
    READ(ioUnit,*)X1_mn_r(:,iMode),X1_r(:,iMode)
  END DO
  READ(ioUnit,*) !## X2: 
  DO iMode=1,X2_modes_r
    READ(ioUnit,*)X2_mn_r(:,iMode),X2_r(:,iMode)
  END DO
  READ(ioUnit,*) !## LA: 
  DO iMode=1,LA_modes_r
    READ(ioUnit,*)LA_mn_r(:,iMode),LA_r(:,iMode)
  END DO
  READ(ioUnit,*) !## profiles at X1_base IP points : spos,phi,chi,iota,pressure 
  DO is=1,X1_nbase_r
    READ(ioUnit,*)profiles_IP(is,:)
  END DO
  !!TODO!! READ a_minor
  READ(ioUnit,*) !## a_minor,r_major,volume 
  READ(ioUnit,*)a_minor,r_major,volume

  CLOSE(ioUnit)

  CALL hmap_new(hmap_r,which_hmap_r)
  ! check if input has changed:

  CALL sgrid_r%init(nElems_r,grid_type_r)

  !needed to build base of restart file
  X1_mn_max_r = (/MAXVAL(X1_mn_r(1,:)),MAXVAL(X1_mn_r(2,:))/nfp_r/)
  X2_mn_max_r = (/MAXVAL(X2_mn_r(1,:)),MAXVAL(X2_mn_r(2,:))/nfp_r/)
  LA_mn_max_r = (/MAXVAL(LA_mn_r(1,:)),MAXVAL(LA_mn_r(2,:))/nfp_r/)

  CALL base_new(X1_base_r,X1_deg_r,X1_cont_r,sgrid_r,degGP_r,X1_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(X1_sin_cos_r),(X1_excl_mn_zero_r.EQ.1))
  CALL base_new(X2_base_r,X2_deg_r,X2_cont_r,sgrid_r,degGP_r,X2_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(X2_sin_cos_r),(X2_excl_mn_zero_r.EQ.1))
  CALL base_new(LA_base_r,LA_deg_r,LA_cont_r,sgrid_r,degGP_r,LA_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(LA_sin_cos_r),(LA_excl_mn_zero_r.EQ.1))

  ALLOCATE(profiles_1d(1:X1_nbase_r,4))
  !convert to spline DOF
  profiles_1d(:,1) =X1_base_r%s%initDOF( profiles_IP(:,1+1) ) !phi
  profiles_1d(:,2) =X1_base_r%s%initDOF( profiles_IP(:,1+2) ) !chi
  profiles_1d(:,3) =X1_base_r%s%initDOF( profiles_IP(:,1+3) ) !iota
  profiles_1d(:,4) =X1_base_r%s%initDOF( profiles_IP(:,1+4) ) !pressure

  DEALLOCATE(sp_r,profiles_IP)
  DEALLOCATE(X1_mn_r)
  DEALLOCATE(X2_mn_r)
  DEALLOCATE(LA_mn_r)


  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE ReadState


!===================================================================================================================================
!> Scalar variables of the equilibrium
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_scalars(Fa,minor_r,n0_global)
! MODULES
USE MODgvec_globals,ONLY: TWOPI
USE MODgvec_gvec_to_gene_Vars,ONLY: a_minor,X1_base_r,profiles_1d
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: Fa                !! toroidal flux at the edge
REAL(wp),INTENT(OUT) :: minor_r           !! length scale, minor radius
INTEGER,INTENT(OUT) :: n0_global         !! number of field periods
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Fa = TWOPI*X1_base_r%s%evalDOF_s(1.0, 0,profiles_1d(:,1)) !phi(s=1)
minor_r=a_minor
n0_global = X1_base_r%f%nfp

END SUBROUTINE gvec_to_gene_scalars


!===================================================================================================================================
!> Evaluate only s dependend variables 
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_profile(spos,q,q_prime,p_prime)
! MODULES
USE MODgvec_globals,ONLY: TWOPI
USE MODgvec_gvec_to_gene_Vars,ONLY: a_minor,profiles_1d,X1_base_r
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: spos              !! radial position (sqrt(phi_norm)), phi_norm: normalized toroidal flux [0,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),OPTIONAL,INTENT(OUT) :: q                !! q=1/iota profile
REAL(wp),OPTIONAL,INTENT(OUT) :: q_prime          !! dq/ds=-(d/ds iota)/iota^2=-(d/ds iota)*q^2
REAL(wp),OPTIONAL,INTENT(OUT) :: p_prime          !! dp/ds, derivative of pressure profile 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
q       = 1./(  X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3)) ) !q=1/iota
q_prime = -q*q*(X1_base_r%s%evalDOF_s(spos, DERIV_S,profiles_1d(:,3)) ) !q'=-iota'/iota^2
p_prime =      (X1_base_r%s%evalDOF_s(spos, DERIV_S,profiles_1d(:,4)) ) !pressure'

END SUBROUTINE gvec_to_gene_profile


!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_coords(nthet,nzeta,spos,theta_star_in,zeta_in,cart_coords)
! MODULES
USE MODgvec_gvec_to_gene_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER              :: nthet          !! number of points in theta_star
INTEGER              :: nzeta          !! number of points in zeta
REAL(wp),INTENT( IN) :: spos           !! radial position (sqrt(phi_norm)), phi_norm: normalized toroidal flux [0,1]
REAL(wp),INTENT( IN) :: theta_star_in(nthet,nzeta)  !! thetaStar poloidal angle
REAL(wp),INTENT( IN) :: zeta_in(      nthet,nzeta)  !! zeta toroidal angle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: cart_coords(3,nthet,nzeta)  !! x,y,z cartesian coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iMode,ithet,izeta
REAL(wp)    :: iota_int,theta_star,theta,zeta
REAL(wp)    :: xp(2),qvec(3)
REAL(wp)    :: X1_s(   1:X1_base_r%f%modes)
REAL(wp)    :: X2_s(   1:X2_base_r%f%modes)
REAL(wp)    :: LA_s(   1:LA_base_r%f%modes)
REAL(wp)    :: X1_int,X2_int
!===================================================================================================================================
!interpolate first in s direction
iota_int = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))

DO iMode=1,X1_base_r%f%modes
  X1_s(iMode)      =X1_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode))
END DO
DO iMode=1,X2_base_r%f%modes
  X2_s(iMode)      =X2_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode))
END DO
DO iMode=1,LA_base_r%f%modes
  LA_s(iMode)      =LA_base_r%s%evalDOF_s(spos,      0,LA_r(:,iMode))
END DO

DO izeta=1,nzeta; DO ithet=1,nthet
  theta_star = theta_star_in(ithet,izeta) !theta_star depends on zeta!!
  zeta = zeta_in(ithet,izeta)
  !find angle theta from straight field line angle theta_star=theta+lambda(s,theta,zeta) 
  theta = theta_star !TODO!!!

  xp=(/theta,zeta/)

  X1_int      =X1_base_r%f%evalDOF_x(xp,0,X1_s)
  X2_int      =X2_base_r%f%evalDOF_x(xp,0,X2_s)

  qvec=(/X1_int,X2_int,zeta/)
  cart_coords(:,ithet,izeta)=hmap_r%eval(qvec)

END DO; END DO !ithet,izeta

END SUBROUTINE gvec_to_gene_coords

!===================================================================================================================================
!> Evaluate gvec state at a list of theta,zeta positions and a fixed s position
!!
!===================================================================================================================================
SUBROUTINE gvec_to_gene_metrics(nthet,nzeta,spos,theta_star_in,zeta_in,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
! MODULES
USE MODgvec_gvec_to_gene_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER          :: nthet          !! number of points in theta_star
INTEGER          :: nzeta          !! number of points in zeta
REAL(wp),INTENT( IN) :: spos           !! radial position (sqrt(phi_norm)), phi_norm: normalized toroidal flux [0,1]
REAL(wp),INTENT( IN) :: theta_star_in(nthet,nzeta)  !! thetaStar poloidal angle
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
REAL(wp) :: iota_int,theta_star,theta,zeta
REAL(wp) :: xp(2),qvec(3)
REAL(wp) :: X1_s(   1:X1_base_r%f%modes)
REAL(wp) :: dX1ds_s(   1:X1_base_r%f%modes)
REAL(wp) :: X2_s(   1:X2_base_r%f%modes)
REAL(wp) :: dX2ds_s(   1:X2_base_r%f%modes)
REAL(wp) :: LA_s(   1:LA_base_r%f%modes)
REAL(wp) :: X1_int,X2_int
REAL(wp) :: dX1ds_int,dX2ds_int
!===================================================================================================================================
!interpolate first in s direction
iota_int = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))

DO iMode=1,X1_base_r%f%modes
  X1_s(   iMode)      =X1_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode))
  dX1ds_s(iMode)      =X1_base_r%s%evalDOF_s(spos,DERIV_S,X1_r(:,iMode))
END DO
DO iMode=1,X2_base_r%f%modes
  X2_s(   iMode)      =X2_base_r%s%evalDOF_s(spos,      0,X2_r(:,iMode))
  dX2ds_s(iMode)      =X2_base_r%s%evalDOF_s(spos,DERIV_S,X2_r(:,iMode))
END DO
DO iMode=1,LA_base_r%f%modes
  LA_s(   iMode)      =LA_base_r%s%evalDOF_s(spos,      0,LA_r(:,iMode))
END DO

DO izeta=1,nzeta; DO ithet=1,nthet
  theta_star = theta_star_in(ithet,izeta) !theta_star depends on zeta!!
  zeta = zeta_in(ithet,izeta)
  !find angle theta from straight field line angle theta_star=theta+lambda(s,theta,zeta) 
  theta = theta_star !TODO!!!

  xp=(/theta,zeta/)

  X1_int      =X1_base_r%f%evalDOF_x(xp, 0, X1_s  )
  dX1ds_int   =X1_base_r%f%evalDOF_x(xp, 0,dX1ds_s)
  X2_int      =X2_base_r%f%evalDOF_x(xp, 0, X2_s   )
  dX2ds_int   =X2_base_r%f%evalDOF_x(xp, 0,dX2ds_s)

  qvec=(/X1_int,X2_int,zeta/)

END DO; END DO !ithet,izeta

END SUBROUTINE gvec_to_gene_metrics

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE finalize_gvec_to_gene 
! MODULES
USE MODgvec_gvec_to_gene_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  DEALLOCATE(X1_r)
  DEALLOCATE(X2_r)
  DEALLOCATE(LA_r)
  DEALLOCATE(profiles_1d)
  CALL sgrid_r%free()
  CALL hmap_r%free()
  CALL X1_base_r%free()
  CALL X2_base_r%free()
  CALL LA_base_r%free()
  DEALLOCATE(hmap_r)
  DEALLOCATE(X1_base_r)
  DEALLOCATE(X2_base_r)
  DEALLOCATE(LA_base_r)

END SUBROUTINE finalize_gvec_to_gene

END MODULE MODgvec_gvec_to_gene
