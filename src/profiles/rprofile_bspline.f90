!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** rProfile **
!!
!! Defines a 1-D profile in rho^2 via B-Spline knots and coefficients
!===================================================================================================================================
MODULE MODgvec_rProfile_bspl
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp, abort
  USE MODgvec_rProfile_base, ONLY: c_rProfile
  USE sll_m_bsplines  ,ONLY: sll_s_bsplines_new, sll_c_bsplines
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_rProfile) :: t_rProfile_bspl

    INTEGER               :: n_knots !! number of knots, including repeated edge knots
    !INTEGER               :: n_coefs !! number of B-Spline coefficients, part of abstract type
    INTEGER               :: deg = 0
    REAL(wp), ALLOCATABLE :: knots(:)   !! knot values, includinng edge knots
    !REAL(wp), ALLOCATABLE :: coefs(:)   !! B-Spline coefficients, part of abstract type
    CLASS(sll_c_bsplines),ALLOCATABLE :: bspl !! b-spline class

    CONTAINS
    PROCEDURE :: eval_at_rho2        => bsplProfile_eval_at_rho2
    PROCEDURE :: antiderivative      => bsplProfile_antiderivative
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
  FUNCTION bsplProfile_new(knots,  coefs) RESULT(sf)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp),    INTENT(IN) :: knots(:)  !! knots of the B-Spline with repeated start and end points
    REAL(wp),    INTENT(IN) :: coefs(:)  !! B-Spline coefficients
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    TYPE(t_rProfile_bspl) :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: n_knots, n_coefs
  !===================================================================================================================================
    n_knots=SIZE(knots)
    n_coefs=SIZE(coefs)
    sf%deg   = COUNT((ABS(knots-knots(1)).LE.1E-12))-1 ! multiplicity of the first knot determines the degree
    IF(COUNT((ABS(knots-knots(n_knots)).LE.1E-12)).NE.sf%deg+1) THEN
      CALL abort(__STAMP__, &
                 "The Bspline knot sequence need the same multiplicity at the beginning and end (=degree+1).")
    END IF
    IF (n_coefs .NE. n_knots-sf%deg-1) THEN
      CALL abort(__STAMP__, &
                 "Number of Bspline coeffcients must be number of knots - (degree+1)!")
    END IF
    sf%n_knots = n_knots
    sf%n_coefs = n_coefs
    ALLOCATE(sf%knots(1:n_knots), sf%coefs(1:n_coefs))
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
    REAL(wp)              , INTENT(IN)  :: rho2 !! evaluation point in the toroidal flux coordinate (rho2=phi/phi_edge= rhopos^2)
    INTEGER , OPTIONAL    , INTENT(IN)  :: deriv !! derivative of bspline(rho^2) in rho^2
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)                         :: profile_prime_value
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp)           :: base_eval(0:sf%deg,0:sf%deg) !! value and derivatives of the (deg+1) B-splines that contribute at s_pos
    INTEGER            :: first_non_zero_bspl !! index offset for the coefficients
    INTEGER             :: deriv_case
  !===================================================================================================================================
    IF (PRESENT(deriv)) THEN
      deriv_case = deriv
    ELSE
      deriv_case = 0
    END IF
    IF(deriv_case.LE.sf%deg)THEN
      CALL sf%bspl%eval_basis_and_n_derivs(rho2,deriv_case,base_eval(0:deriv_case,:),first_non_zero_bspl)
      profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*base_eval(deriv_case,:))
    ELSE
      profile_prime_value = 0.0_wp
    END IF
  END FUNCTION bsplProfile_eval_at_rho2

  !===================================================================================================================================
  !> get the exact spline antiderivative, with respect to rho2
  !! the knotspan is increased by an extra multiplicity on both ends, and the new coefficients are computed as
  !! beta(i) = beta(i-1) + alpha(i)*(t(i+degree+1)-t(i))/(degree+1)
  !! From deBoor, "A practical guide to Splines", p.128
  !!
  !===================================================================================================================================
  FUNCTION bsplProfile_antiderivative(sf) RESULT(antideriv)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_rProfile_bspl), INTENT(IN)  :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(c_rProfile),ALLOCATABLE :: antideriv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp) :: coefs(sf%n_coefs+1), knots(sf%n_knots+2), intermid
    INTEGER :: i, n_coefs, n_knots, deg
  !===================================================================================================================================
    coefs = 0.0_wp
    knots = -42.0_wp

    n_coefs = sf%n_coefs+1
    n_knots = sf%n_knots+2
    deg = sf%deg+1
    ! increase multiplicity at the edges
    knots(2:n_knots-1) = sf%knots
    knots(1) = sf%knots(1)
    knots(n_knots) = sf%knots(sf%n_knots)

    DO i=1,sf%n_coefs
      intermid = sf%coefs(i)*(sf%knots(i+deg)-sf%knots(i))/deg
      coefs(i+1) = coefs(i) + intermid
    END DO
    antideriv = t_rProfile_bspl(knots,  coefs)
  END FUNCTION bsplProfile_antiderivative

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
    SDEALLOCATE(sf%knots)
    SDEALLOCATE(sf%coefs)
  END SUBROUTINE bsplProfile_free

END MODULE MODgvec_rProfile_bspl
