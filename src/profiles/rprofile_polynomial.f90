!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** polyProfile **
!!
!! Defines a 1-D profile in rho^2 via a power polynomial.
!===================================================================================================================================
MODULE MODgvec_rProfile_poly
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_rProfile_base, ONLY: c_rProfile, poly_derivative_prefactor
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_rProfile) :: t_rProfile_poly
    !INTEGER               :: n_coefs !! number of polynomial coefficients, part of abstract type
    INTEGER               :: deg = 0
    !REAL(wp), ALLOCATABLE :: coefs(:)   !! polynomial coefficients, part of abstract type
    CONTAINS
    PROCEDURE :: eval_at_rho2        => polyProfile_eval_at_rho2
    PROCEDURE :: antiderivative      => polyProfile_antiderivative
    FINAL     :: polyProfile_free

  END TYPE t_rProfile_poly

  INTERFACE t_rProfile_poly
      MODULE PROCEDURE polyProfile_new
  END INTERFACE t_rProfile_poly

  CONTAINS

  FUNCTION polyProfile_new(coefs) RESULT(sf)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL,    INTENT(IN) :: coefs(:)  !! B-Spline coefficients
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    TYPE(t_rProfile_poly) :: sf
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: n_coefs
  !===================================================================================================================================
    n_coefs=SIZE(coefs)
    sf%deg   = n_coefs-1
    sf%n_coefs = n_coefs
    ALLOCATE(sf%coefs(1:n_coefs))
    sf%coefs = coefs
  END FUNCTION polyProfile_new

  !===================================================================================================================================
  !> evaluate the n-th derivative of a power polynomial
  !!
  !===================================================================================================================================
  FUNCTION polyProfile_eval_at_rho2(sf, rho2, deriv) RESULT(profile_prime_value)
  ! MODULES
    USE MODgvec_Globals, ONLY: Eval1DPoly,Eval1DPoly_deriv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(IN)  :: sf !! self
    REAL(wp)              , INTENT(IN)  :: rho2 !! evaluation point in the toroidal flux coordinate (rho2=phi/phi_edge= spos^2)
    INTEGER , OPTIONAL    , INTENT(IN)  :: deriv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)                         :: profile_prime_value
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp)                         :: prefactor
    INTEGER                          :: d
    INTEGER                          :: deriv_case
  !===================================================================================================================================
    IF (PRESENT(deriv)) THEN
        deriv_case = deriv
    ELSE
        deriv_case = 0
    END IF

    IF (deriv_case>sf%deg) THEN
        profile_prime_value = 0.0_wp
    ELSE IF (deriv_case==0) THEN
        profile_prime_value = EVAL1DPOLY(sf%n_coefs, sf%coefs, rho2)
    ELSE IF (deriv_case==1) THEN
        profile_prime_value = EVAL1DPOLY_deriv(sf%n_coefs, sf%coefs, rho2)
    ELSE
        profile_prime_value = 0.0_wp
        DO d=sf%deg+1, deriv_case+1,-1
            prefactor=poly_derivative_prefactor(d-1,deriv_case)
            profile_prime_value = profile_prime_value*rho2+prefactor*sf%coefs(d)
        END DO
    END IF
  END FUNCTION polyProfile_eval_at_rho2

  !===================================================================================================================================
  !> get the exact polynomial antiderivative, with respect to rho2
  !!
  !===================================================================================================================================
  FUNCTION polyProfile_antiderivative(sf) RESULT(antideriv)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(IN)  :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(c_rProfile),ALLOCATABLE :: antideriv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp) :: coefs(sf%n_coefs+1)
    INTEGER :: i
  !===================================================================================================================================
    coefs = 0.0_wp
    DO i=1,sf%n_coefs
      coefs(i+1) = sf%coefs(i)/REAL(i,wp)
    END DO
    antideriv = t_rProfile_poly(coefs)
  END FUNCTION polyProfile_antiderivative

  !===================================================================================================================================
  !> finalize the type rProfile
  !!
  !===================================================================================================================================
  SUBROUTINE polyProfile_free(sf)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    TYPE(t_rProfile_poly), INTENT(INOUT)  :: sf !! self
  !===================================================================================================================================
    SDEALLOCATE(sf%coefs)
  END SUBROUTINE polyProfile_free

END MODULE MODgvec_rProfile_poly
