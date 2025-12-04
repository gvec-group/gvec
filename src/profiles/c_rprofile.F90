!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** rProfile **
!!
!! Abstract class for radial profiles
!===================================================================================================================================
MODULE MODgvec_rProfile_base
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp, abort
  IMPLICIT NONE

  PUBLIC

  TYPE, ABSTRACT :: c_rProfile
      INTEGER               :: n_coefs
      REAL(wp), ALLOCATABLE :: coefs(:)

      contains

      PROCEDURE(i_fun_eval_at_rho2), DEFERRED :: eval_at_rho2
      PROCEDURE(i_fun_antiderivative), DEFERRED :: antiderivative

      PROCEDURE :: eval_at_rho => rProfile_eval_at_rho
      ! hard coded derivatives with respect to rho=sqrt(phi/phi_edge)
      PROCEDURE, PRIVATE :: rProfile_drho2
      PROCEDURE, PRIVATE :: rProfile_drho3
      PROCEDURE, PRIVATE :: rProfile_drho4

  end type c_rProfile

  ABSTRACT INTERFACE

      FUNCTION i_fun_eval_at_rho2( sf, rho2, deriv ) RESULT(profile_value)
          IMPORT c_rProfile
          IMPORT wp
          CLASS(c_rProfile), INTENT(IN)  :: sf
          REAL(wp)         , INTENT(IN)  :: rho2
          INTEGER, OPTIONAL, INTENT(IN)  :: deriv
          REAL(wp)                       :: profile_value
      END FUNCTION i_fun_eval_at_rho2

      FUNCTION i_fun_antiderivative(sf) RESULT(antideriv)
          IMPORT c_rProfile
          IMPORT wp
          CLASS(c_rProfile), INTENT(IN)   :: sf
          CLASS(c_rProfile), ALLOCATABLE  :: antideriv
      END FUNCTION i_fun_antiderivative
  END INTERFACE

  CONTAINS
  !===================================================================================================================================
  !> calculate the prefactor for the d-th coefficient of the n-th derivative of a polynomial
  !!
  !===================================================================================================================================
  PURE FUNCTION poly_derivative_prefactor(D,deriv) RESULT(prefactor)
    INTEGER, INTENT(IN) :: D,deriv
    INTEGER :: i
    REAL(wp) :: prefactor
    !===================================================================================================================================
    prefactor = 1.0_wp
    DO i=D-deriv+1,D
        prefactor = prefactor*i
    END DO
  END FUNCTION poly_derivative_prefactor

  !===================================================================================================================================
  !> evaluate the n-th derivative of (rho^2) with respect to rho ~sqrt(magnetic flux).
  !!
  !===================================================================================================================================
  PURE FUNCTION rho2_derivative(rho,deriv) RESULT(rho2_prime)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(IN) :: rho                !! rho position rho ~sqrt(magnetic flux)
    INTEGER, INTENT(IN)  :: deriv  !! derivative in rho
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp) :: rho2_prime     !!n-th derivative of rho2 with respect to rho ~sqrt(magnetic flux).
  !===================================================================================================================================
    IF (deriv>2) THEN
      rho2_prime = 0.0_wp
  ELSE
      rho2_prime = poly_derivative_prefactor(2,deriv)*rho**(2-deriv)
  END IF
  END FUNCTION rho2_derivative

  !===================================================================================================================================
  !> evaluate the 2nd derivative of a radial profile with respect to rho ~sqrt(magnetic flux).
  !!
  !===================================================================================================================================
  FUNCTION rProfile_drho2(sf, rho) RESULT(derivative)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(c_rProfile)  :: sf !! self
    REAL(wp), INTENT(IN) :: rho                !! rho position rho ~sqrt(magnetic flux)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp) :: derivative   !! 2nd derivative of a radial profile with respect to rho ~sqrt(magnetic flux).
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp) :: rho2
  !===================================================================================================================================
    rho2 = rho2_derivative(rho,deriv=0)
    ! d^2/dx^2 f(g(x)) = [g'(x)]^2 * f''(g(x))+g''(x) * f'(g(x))
    derivative =  rho2_derivative(rho,deriv=1)**2*sf%eval_at_rho2(rho2, deriv=2) &
                + rho2_derivative(rho,deriv=2)*   sf%eval_at_rho2(rho2, deriv=1)
  END FUNCTION rProfile_drho2

  !===================================================================================================================================
  !> evaluate the 3rd derivative of a radial profile with respect to rho ~sqrt(magnetic flux).
  !!
  !===================================================================================================================================
  FUNCTION rProfile_drho3(sf, rho) RESULT(derivative)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(c_rProfile)  :: sf !! self
    REAL(wp), INTENT(IN) :: rho                !! rho position rho ~sqrt(magnetic flux)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp) :: derivative          !! 3rd derivative of a radial profile with respect to rho ~sqrt(magnetic flux)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp) :: rho2
  !===================================================================================================================================
    rho2 = rho2_derivative(rho,deriv=0)
    ! d^3/dx^3 f(g(x)) = 3g'(x) * f''(g(x)) * g''(x) + [g'(x)]^3*f'''(g(x)) + f'(x) * g'''(x)
    derivative = 3*rho2_derivative(rho,deriv=1)*   sf%eval_at_rho2(rho2, deriv=2)*rho2_derivative(rho,deriv=2) &
                +  rho2_derivative(rho,deriv=1)**3*sf%eval_at_rho2(rho2, deriv=3) &
                +  rho2_derivative(rho,deriv=3)*   sf%eval_at_rho2(rho2, deriv=1)
  END FUNCTION rProfile_drho3

  !===================================================================================================================================
  !> evaluate the 4th derivative of a radial profile with respect to rho ~sqrt(magnetic flux)
  !!
  !===================================================================================================================================
  FUNCTION rProfile_drho4(sf, rho) RESULT(derivative)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(c_rProfile) :: sf !! self
    REAL(wp), INTENT(IN) :: rho                !! rho position rho ~sqrt(magnetic flux)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp) :: derivative   !! 4th derivative of a radial profile with respect to rho ~sqrt(magnetic flux).
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp) :: rho2
  !===================================================================================================================================
    rho2 = rho2_derivative(rho,deriv=0)
    ! d^4/dx^4 f(g(x)) = f''''(g(x))g'(x)**4
    !                   + 6f'''(g(x))g''(x)g'(x)^2
    !                   + 3f''(g(x))g''(x)^2
    !                   + 4f''(g(x))g'''(x)g'(x)
    !                   + f'(g(x))g''''(x)
    derivative =    sf%eval_at_rho2(rho2, deriv=4)*rho2_derivative(rho,deriv=1)**4 &
                + 6*sf%eval_at_rho2(rho2, deriv=3)*rho2_derivative(rho,deriv=2)*rho2_derivative(rho,deriv=1)**2 &
                + 3*sf%eval_at_rho2(rho2, deriv=2)*rho2_derivative(rho,deriv=2)**2 &
                + 4*sf%eval_at_rho2(rho2, deriv=2)*rho2_derivative(rho,deriv=3)*rho2_derivative(rho,deriv=1) &
                +   sf%eval_at_rho2(rho2, deriv=1)*rho2_derivative(rho,deriv=4)
  END FUNCTION rProfile_drho4

  !===================================================================================================================================
  !> evaluate the n-th derivative of a radial profile with respect to rho ~sqrt(magnetic flux).
  !! NOTE: n has to be in [0,4] due to an explicit implementation of the product rule.
  !===================================================================================================================================
  FUNCTION rProfile_eval_at_rho(sf, rho, deriv) RESULT(derivative)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(c_rProfile)  :: sf !! self
    REAL(wp), INTENT(IN) :: rho                !! rho position rho ~sqrt(magnetic flux)
    INTEGER,  OPTIONAL, INTENT(IN)   :: deriv  !! derivative in rho
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp) :: derivative  !! n-th derivative of a radial profile with respect to rho ~sqrt(magnetic flux).
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: deriv_case
    REAL(wp) :: rho2
  !===================================================================================================================================
    IF (PRESENT(deriv)) THEN
        deriv_case = deriv
    ELSE
        deriv_case = 0
    END IF

    rho2 = rho2_derivative(rho,deriv=0)
    SELECT CASE(deriv_case)
    CASE(0)
        derivative = sf%eval_at_rho2(rho2, deriv=0)
    CASE(1)
        derivative = sf%eval_at_rho2(rho2,deriv=1)*rho2_derivative(rho,deriv=1)
    CASE(2)
        derivative = sf%rProfile_drho2(rho)
    CASE(3)
        derivative = sf%rProfile_drho3(rho)
    CASE(4)
        derivative = sf%rProfile_drho4(rho)
    CASE DEFAULT
        CALL abort(__STAMP__,&
            "error in rprofile: derivatives higher than 4 with respect to rho=sqrt(phi/phi_edge) are not implemented!")
    END SELECT
  END FUNCTION rProfile_eval_at_rho

END MODULE MODgvec_rProfile_base
