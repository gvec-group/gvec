!===================================================================================================================================
! Copyright (C) 2017 - 2024  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2021 - 2022  Tiago Ribeiro
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
!!# Module ** splProfile **
!!
!! Defines a 1-D profile in rho^2 via B-Spline knots and coefficients
!===================================================================================================================================
MODULE MODgvec_splProfile
! MODULES
USE MODgvec_Globals ,ONLY: wp,Unit_stdOut,abort,MPIRoot
USE sll_m_bsplines  ,ONLY: sll_s_bsplines_new, sll_c_bsplines

IMPLICIT NONE

PUBLIC 

TYPE :: t_splProfile
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER               :: n_knots !! number of knots, including repeated edge knots
  INTEGER               :: n_coefs !! number of basis coefficients
  REAL(wp), ALLOCATABLE :: knots(:)   !! knot values, includinng edge knots
  REAL(wp), ALLOCATABLE :: coefs(:)   !! basis coefficients
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER                           :: deg = 0!! degree of the B-splines
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl
  
  CONTAINS

  !PROCEDURE :: init             => splProfile_init
  !PROCEDURE :: copy             => splProfile_copy
  !PROCEDURE :: compare          => splProfile_compare
  PROCEDURE :: eval_at_phi_norm        => splProfile_eval_at_phi_norm
  PROCEDURE :: eval_prime_at_phi_norm  => splProfile_eval_prime_at_phi_norm
  FINAL     :: splProfile_free
END TYPE t_splProfile

interface t_splProfile
  module procedure :: splProfile_init
end interface t_splProfile

!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize the type splProfile
!!
!===================================================================================================================================
FUNCTION splProfile_init(knots, n_knots, coefs, n_coefs) RESULT(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: n_knots !! number of knots
  INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
  REAL,    INTENT(IN) :: knots(n_knots)  !! knots of the B-Spline with repeated start and end points
  REAL,    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_splProfile)  :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
  __PERFON("splProfile_init")
  sf%deg   = COUNT(knots==knots(1))-1 ! multiplicity of the first knot determines the degree
  sf%n_knots = n_knots
  sf%n_coefs = n_coefs
  ALLOCATE(sf%knots(1:n_knots))
  sf%knots = knots
  ALLOCATE(sf%coefs(1:n_coefs))
  sf%coefs = coefs
  CALL sll_s_bsplines_new(sf%bspl, sf%deg, .FALSE., & 
                          sf%knots(1),sf%knots(n_knots),&
                          sf%deg , & 
                          sf%knots(sf%deg+1:n_knots-sf%deg)) ! remove repeated edge knots
  __PERFOFF("splProfile_init")

END FUNCTION splProfile_init

!===================================================================================================================================
!> evaluate the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_splProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= s^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: basis_values(sf%deg+1) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
  CALL sf%bspl%eval_basis(phi_norm,basis_values,first_non_zero_bspl)
  profile_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*basis_values(:sf%deg+1))
END FUNCTION splProfile_eval_at_phi_norm

!===================================================================================================================================
!> evaluate the first derivative of the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_prime_at_phi_norm( sf, phi_norm ) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_splProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (s=phi/phi_edge= rho^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: deriv_values(sf%deg+1) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
  CALL sf%bspl%eval_deriv(phi_norm,deriv_values,first_non_zero_bspl)
  profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*deriv_values(:sf%deg+1))
END FUNCTION splProfile_eval_prime_at_phi_norm

!===================================================================================================================================
!> finalize the type splProfile
!!
!===================================================================================================================================
SUBROUTINE splProfile_free(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_splProfile), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!=================================================================================================================================== 
  SDEALLOCATE(sf%knots)
  SDEALLOCATE(sf%coefs)
  IF (ALLOCATED(sf%bspl)) CALL sf%bspl%free()
  SDEALLOCATE(sf%bspl)
  sf%deg        =-1
  sf%n_coefs    =-1
  sf%n_knots    =-1
END SUBROUTINE splProfile_free


END MODULE MODgvec_splProfile

