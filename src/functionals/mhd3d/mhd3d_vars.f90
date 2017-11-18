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
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_MHD3D_Vars
! MODULES
USE MOD_Globals,ONLY: PI,wp,Unit_stdOut,abort
USE MOD_sgrid,  ONLY: t_sgrid
USE MOD_base,   ONLY: t_base
USE MOD_Sol_Var,ONLY: t_sol_var
USE MOD_c_hmap, ONLY: c_hmap
IMPLICIT NONE
PUBLIC


!-----------------------------------------------------------------------------------------------------------------------------------
! Globally used variables
!-----------------------------------------------------------------------------------------------------------------------------------
CLASS(t_base),ALLOCATABLE :: X1base            !! container for base of variable X1
CLASS(t_base),ALLOCATABLE :: X2base            !! container for base of variable X2
CLASS(t_base),ALLOCATABLE :: LAbase            !! container for base of variable lambda 
                             
TYPE(t_sgrid)    :: sgrid             !! only one grid up to now

TYPE(t_sol_var),ALLOCATABLE  :: U(:)  !! solutions at levels (k-1),(k),(k+1)
TYPE(t_sol_var)              :: dUdt  !! solution update

INTEGER          :: X1X2_BC(2)        !! BC axis (0) and edge (1)   for variables X1 and X2
INTEGER          :: LA_BC(2)          !! BC axis (0) and edge (1)   for variable lambda
INTEGER          :: nDOF_X1           !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_X2           !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_LA           !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 


!===================================================================================================================================
! locally used in functional evaluation only

! input parameters for functional
INTEGER              :: NFP           !! number of field periods
REAL(wp)             :: Phi_edge      !! toroidal flux at the last flux surface of the domain
INTEGER              :: n_mass_coefs  !! number of polynomial coeffients for mass profile
INTEGER              :: n_iota_coefs  !! number of polynomial coeffients for iota profile
REAL(wp),ALLOCATABLE :: mass_coefs(:) !! polynomial coefficients of the mass profile
REAL(wp),ALLOCATABLE :: iota_coefs(:) !! polynomial coefficients of the iota profile
!constants
REAL(wp)             :: mu_0          !! permeability
REAL(wp)             :: gamm          !! isentropic exponent

! --- store some variables to help/speed up evaluation of the functional
CLASS(c_hmap),ALLOCATABLE :: hmap    !! type containing subroutines for evaluating the map h (Omega_p x S^1) --> Omega

!evaluations at radial gauss points,  size(1:base%s%nGP)                                       
REAL(wp),ALLOCATABLE :: mass_GP(:)    !! mass profile 
REAL(wp),ALLOCATABLE :: iota_GP(:)    !! iota profile 
REAL(wp),ALLOCATABLE :: PhiPrime_GP(:)!! s derivative of toroidal flux : Phi'(s)
REAL(wp),ALLOCATABLE :: Vprime_GP(:)  !! s derivative of volume: V'(s)
!evaluations at all integration points,  size(1:base%f%mn_IP,1:base%s%nGP)                                       
REAL(wp),ALLOCATABLE :: J_h(:,:)      !! Jacobian of the mapping h (X1,X2,zeta) -->(x,y,z) (global Jacobian: J_h*J_p) 
REAL(wp),ALLOCATABLE :: J_p(:,:)      !! Jacobian of poloidal mapping: dX1_ds*dX2_dtheta - dX2_ds*dX1_theta
REAL(wp),ALLOCATABLE :: dX1_ds(:,:)   !! radial derivative of X1
REAL(wp),ALLOCATABLE :: dX2_ds(:,:)   !! radial derivative of X2
REAL(wp),ALLOCATABLE :: dX1_dthet(:,:)!! theta  derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dthet(:,:)!! theta  derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dthet(:,:)!! theta  derivative of lambda
REAL(wp),ALLOCATABLE :: dX1_dzeta(:,:)!! zeta   derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dzeta(:,:)!! zeta   derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dzeta(:,:)!! zeta   derivative of lambda
REAL(wp),ALLOCATABLE :: b_a (  :,:,:) !! b=(iota-dlamba_dzeta,1+dlambda_dtheta), normalized contravariant magnetic field components
                                      !! at all points, size(1:2,1:mn_IP,1:nGP)
REAL(wp),ALLOCATABLE :: g_ab(:,:,:,:) !! metric tensor g_(alpha,beta), alpha/beta=theta and zeta, size(1:2,1:2,1:mn_IP,1:nGP)

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type hmap with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_new( sf, which_hmap)
! MODULES
USE MOD_hmap_RZ,   ONLY:t_hmap_RZ
USE MOD_hmap_knot, ONLY:t_hmap_knot
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: which_hmap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_hmap), ALLOCATABLE,INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(which_hmap)
  CASE(1)
    ALLOCATE(t_hmap_RZ :: sf)
  CASE(2)
    ALLOCATE(t_hmap_knot :: sf)
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "this hmap choice does not exist  !")
  END SELECT 
  sf%which_hmap=which_hmap
  CALL sf%init()

END SUBROUTINE hmap_new


END MODULE MOD_MHD3D_Vars

