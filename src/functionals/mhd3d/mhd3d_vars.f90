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
USE MOD_fbase,  ONLY: t_fbase
USE MOD_Sol_Var_MHD3D,ONLY: t_sol_var_MHD3D
USE MOD_c_hmap, ONLY: c_hmap
IMPLICIT NONE
PUBLIC


!-----------------------------------------------------------------------------------------------------------------------------------
! derived type variables
!-----------------------------------------------------------------------------------------------------------------------------------
CLASS(t_base),  ALLOCATABLE :: X1_base   !! container for base of variable X1
CLASS(t_base),  ALLOCATABLE :: X2_base   !! container for base of variable X2
CLASS(t_base),  ALLOCATABLE :: LA_base   !! container for base of variable lambda 
                            
CLASS(t_fbase), ALLOCATABLE :: X1_b_base !! container for base of boundaries o X1
CLASS(t_fbase), ALLOCATABLE :: X2_b_base !! container for base of boundaries o X1
CLASS(t_fbase), ALLOCATABLE :: LA_b_base !! container for base of boundaries o LA
                             
                             
TYPE(t_sgrid)               :: sgrid     !! only one grid up to now
                                         
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: U(:)      !! solutions at levels (k-1),(k),(k+1)
TYPE(t_sol_var_MHD3D)             :: dUdt      !! solution update
INTEGER                     :: nDOF_X1   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER                     :: nDOF_X2   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER                     :: nDOF_LA   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
                                         
                                         
CLASS(c_hmap),  ALLOCATABLE :: hmap      !! type containing subroutines for evaluating the map h (Omega_p x S^1) --> Omega

!===================================================================================================================================
! input parameters for functional
INTEGER              :: which_init      !! select initialization. 0: only using input parameter, 1: using a VMEC equilibrium
INTEGER              :: NFP             !! number of field periods
REAL(wp)             :: Phi_edge        !! toroidal flux at the last flux surface of the domain
INTEGER              :: n_mass_coefs    !! number of polynomial coeffients for mass profile
INTEGER              :: n_iota_coefs    !! number of polynomial coeffients for iota profile
REAL(wp),ALLOCATABLE :: mass_coefs(:)   !! polynomial coefficients of the mass profile
REAL(wp),ALLOCATABLE :: iota_coefs(:)   !! polynomial coefficients of the iota profile
!constants
REAL(wp)             :: mu_0            !! permeability
REAL(wp)             :: gamm            !! isentropic exponent

INTEGER              :: X1X2_BC(2)      !! BC axis (0) and edge (1)   for variables X1 and X2 (default(0,1))
INTEGER              :: LA_BC(2)        !! BC axis (0) and edge (1)   for variable lambda     (default(0,0))
                                        
REAL(wp),ALLOCATABLE :: X1_b(:)         !! fourier modes of the boundary for X1
REAL(wp),ALLOCATABLE :: X2_b(:)         !! fourier modes of the boundary for X2
REAL(wp),ALLOCATABLE :: LA_b(:)         !! fourier modes of the boundary for LA

!evaluations at radial gauss points, size(1:base%s%nGP)                                       
REAL(wp),ALLOCATABLE :: mass_GP(:)      !! mass profile 
REAL(wp),ALLOCATABLE :: iota_GP(:)      !! iota profile 
REAL(wp),ALLOCATABLE :: PhiPrime_GP(:)  !! s derivative of toroidal flux : Phi'(s)
REAL(wp),ALLOCATABLE :: Vprime_GP(:)    !! s derivative of volume: V'(s)

!evaluations at all integration points, size(1:base%f%mn_IP,1:base%s%nGP)                                       
REAL(wp),ALLOCATABLE :: J_h(:,:)        !! Jacobian of the mapping h (X1,X2,zeta) -->(x,y,z) (global Jacobian: J_h*J_p) 
REAL(wp),ALLOCATABLE :: J_p(:,:)        !! Jacobian of poloidal mapping: dX1_ds*dX2_dtheta - dX2_ds*dX1_theta
REAL(wp),ALLOCATABLE :: dX1_ds(:,:)     !! radial derivative of X1
REAL(wp),ALLOCATABLE :: dX2_ds(:,:)     !! radial derivative of X2
REAL(wp),ALLOCATABLE :: dX1_dthet(:,:)  !! theta  derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dthet(:,:)  !! theta  derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dthet(:,:)  !! theta  derivative of lambda
REAL(wp),ALLOCATABLE :: dX1_dzeta(:,:)  !! zeta   derivative of X1
REAL(wp),ALLOCATABLE :: dX2_dzeta(:,:)  !! zeta   derivative of X2
REAL(wp),ALLOCATABLE :: dLA_dzeta(:,:)  !! zeta   derivative of lambda
REAL(wp),ALLOCATABLE :: b_a (  :,:,:)   !! b=(iota-dlamba_dzeta,1+dlambda_dtheta), normalized contravariant magnetic field 
                                        !! components at all points, size(1:2,1:mn_IP,1:nGP)
REAL(wp),ALLOCATABLE :: g_ab(:,:,:,:)   !! metric tensor g_(alpha,beta), alpha/beta=theta and zeta, size(1:2,1:2,1:mn_IP,1:nGP)

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
USE MOD_hmap_cyl, ONLY:t_hmap_cyl
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
    ALLOCATE(t_hmap_cyl  :: sf)
  CASE(10)
    ALLOCATE(t_hmap_knot :: sf)
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "this hmap choice does not exist  !")
  END SELECT 
  sf%which_hmap=which_hmap
  CALL sf%init()

END SUBROUTINE hmap_new


END MODULE MOD_MHD3D_Vars

