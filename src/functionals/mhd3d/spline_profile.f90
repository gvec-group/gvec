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

PRIVATE
PUBLIC t_splProfile,splProfile_new

TYPE :: t_splProfile
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER               :: n_knots !! number of knots, including repeated edge knots
  INTEGER               :: n_coefs !! number of basis coefficients
  REAL(wp), ALLOCATABLE :: knots(:)   !! knot values, includinng edge knots
  REAL(wp), ALLOCATABLE :: coefs(:)   !! basis coefficients
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL                           :: initialized
  INTEGER                           :: deg !! degree of the B-splines
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl
  
  CONTAINS

  PROCEDURE :: init             => splProfile_init
  PROCEDURE :: free             => splProfile_free
  !PROCEDURE :: copy             => splProfile_copy
  !PROCEDURE :: compare          => splProfile_compare
  PROCEDURE :: eval_at_phi_norm        => splProfile_eval_at_phi_norm
  PROCEDURE :: eval_prime_at_phi_norm  => splProfile_eval_prime_at_phi_norm
END TYPE t_splProfile



!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> allocate the type splProfile 
!!
!===================================================================================================================================
SUBROUTINE splProfile_new( sf, knots, n_knots, coefs, n_coefs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER , INTENT(IN) :: n_knots !! number of knots
  INTEGER , INTENT(IN) :: n_coefs !! number of coefficients
  REAL(wp), INTENT(IN) :: knots(n_knots)  !! knots of the B-Spline with repeated edge knots
  REAL(wp), INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_splProfile), ALLOCATABLE,INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(sf)
  __PERFON("splProfile_new")
  CALL sf%init(knots, n_knots, coefs, n_coefs)

  __PERFOFF("splProfile_new")
END SUBROUTINE splProfile_new

!===================================================================================================================================
!> initialize the type splProfile maximum mode numbers, number of integration points, type of basis (sin/cos or sin and cos) 
!!
!===================================================================================================================================
SUBROUTINE splProfile_init( sf, knots, n_knots, coefs, n_coefs )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: n_knots !! number of knots
  INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
  REAL,    INTENT(IN) :: knots(n_knots)  !! knots of the B-Spline with repeated start and end points
  REAL,    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_splProfile), INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
  IF(sf%initialized) THEN
    CALL abort(__STAMP__, &
    "Trying to reinit splProfile!") 
  END IF
  sf%deg   = COUNT(knots==knots(1))-1
  sf%n_knots = n_knots
  sf%n_coefs = n_coefs
  CALL splProfile_alloc(sf)
  sf%knots = knots
  sf%coefs = coefs
  CALL sll_s_bsplines_new(sf%bspl, sf%deg, .FALSE., & 
                          sf%knots(1),sf%knots(n_knots),&
                          sf%deg , & 
                          sf%knots(sf%deg+1:n_knots-sf%deg)) ! remove repeated edge knots
  
  sf%initialized = .TRUE.
  
END SUBROUTINE splProfile_init


!===================================================================================================================================
!> allocate all variables in splProfile
!!
!===================================================================================================================================
SUBROUTINE splProfile_alloc( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_splProfile), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ASSOCIATE(&
              n_knots   => sf%n_knots     &
            , n_coefs   => sf%n_coefs     &
            )
      ALLOCATE(sf%knots(1:n_knots))
      ALLOCATE(sf%coefs(1:n_coefs))
  END ASSOCIATE !m_nyq,n_nyq,modes
END SUBROUTINE splProfile_alloc

!===================================================================================================================================
!> evaluate the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
IMPLICIT NONE
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
IMPLICIT NONE
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
SUBROUTINE splProfile_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_splProfile), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN
  !allocatables  
  SDEALLOCATE(sf%knots)
  SDEALLOCATE(sf%coefs)
  CALL sf%bspl%free()
  SDEALLOCATE(sf%bspl)
  
  sf%deg        =-1
  sf%n_coefs    =-1
  sf%n_knots    =-1
  sf%initialized=.FALSE.

END SUBROUTINE splProfile_free


END MODULE MODgvec_splProfile

