!===================================================================================================================================
! Copyright (C) 2024 Robert Koeberl <robert.koeberl@ipp.mpg.de>
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
!!# Module ** rProfile **
!!
!! Defines a 1-D profile in rho^2 via B-Spline knots and coefficients
!===================================================================================================================================
MODULE MODgvec_rProfile_bspl
! MODULES
USE MODgvec_Globals ,ONLY: wp
USE MODgvec_rProfile_base, ONLY: c_rProfile
USE sll_m_bsplines  ,ONLY: sll_s_bsplines_new, sll_c_bsplines
IMPLICIT NONE

PUBLIC

TYPE, EXTENDS(c_rProfile) :: t_rProfile_bspl
  
    INTEGER               :: n_knots !! number of knots, including repeated edge knots
    REAL(wp), ALLOCATABLE :: knots(:)   !! knot values, includinng edge knots
    CLASS(sll_c_bsplines),ALLOCATABLE :: bspl !! b-spline class
    
  CONTAINS
  
  PROCEDURE :: eval_at_rho2        => bsplProfile_eval_at_rho2
  FINAL :: bsplProfile_free
  
END TYPE t_rProfile_bspl

INTERFACE t_rProfile_bspl
    MODULE PROCEDURE bsplProfile_new
END INTERFACE t_rProfile_bspl

CONTAINS

!===================================================================================================================================
!> initialize the rProfile of type bspline
!!
!===================================================================================================================================
FUNCTION bsplProfile_new(knots, n_knots, coefs, n_coefs) RESULT(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: n_knots !! number of knots
    INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
    REAL(wp),    INTENT(IN) :: knots(n_knots)  !! knots of the B-Spline with repeated start and end points
    REAL(wp),    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    TYPE(t_rProfile_bspl) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
    sf%deg   = COUNT(knots==knots(1))-1 ! multiplicity of the first knot determines the degree
    sf%n_knots = n_knots
    sf%n_coefs = n_coefs
    sf%knots = knots
    sf%coefs = coefs
    IF (sf%deg>0) THEN
      CALL sll_s_bsplines_new(sf%bspl, sf%deg, .FALSE., & 
                              sf%knots(1),sf%knots(n_knots),&
                              size(sf%knots(sf%deg+1:n_knots-sf%deg))-1 , & ! number of knots handed to the library
                              sf%knots(sf%deg+1:n_knots-sf%deg)) ! remove repeated edge knots
    END IF

END FUNCTION bsplProfile_new

!===================================================================================================================================
!> evaluate the n-th derivative of the bsplProfile at position s
!!
!===================================================================================================================================
FUNCTION bsplProfile_eval_at_rho2( sf, rho2, deriv ) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile_bspl), INTENT(IN)  :: sf !! self
  REAL(wp)              , INTENT(IN)  :: rho2 !! evaluation point in the toroidal flux coordinate (rho2=phi/phi_edge= spos^2)
  INTEGER , OPTIONAL    , INTENT(IN)    :: deriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp), ALLOCATABLE            :: deriv_values(:,:) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
  INTEGER                          :: deriv_case
!===================================================================================================================================
    IF (PRESENT(deriv)) THEN
      deriv_case = deriv
    ELSE
      deriv_case = 0
    END IF

    ALLOCATE(deriv_values(deriv_case+1,sf%deg+1))
    IF (sf%deg>0) THEN
      CALL sf%bspl%eval_basis_and_n_derivs(rho2,deriv_case,deriv_values,first_non_zero_bspl)
      profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*deriv_values(deriv_case+1,:sf%deg+1))
    ELSE 
      profile_prime_value = 0.0_wp
    END IF
    SDEALLOCATE(deriv_values)
END FUNCTION bsplProfile_eval_at_rho2

!===================================================================================================================================
!> finalize the type rProfile
!!
!===================================================================================================================================
SUBROUTINE bsplProfile_free(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_rProfile_bspl), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!=================================================================================================================================== 
  IF (ALLOCATED(sf%bspl)) CALL sf%bspl%free()
END SUBROUTINE bsplProfile_free

END MODULE MODgvec_rProfile_bspl
