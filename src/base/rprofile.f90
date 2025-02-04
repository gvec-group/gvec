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
!!# Module ** rProfile **
!!
!! Defines a 1-D profile in rho^2 via B-Spline knots and coefficients
!===================================================================================================================================
MODULE MODgvec_rProfile
! MODULES
USE MODgvec_Globals ,ONLY: wp,Unit_stdOut,abort,MPIRoot
USE sll_m_bsplines  ,ONLY: sll_s_bsplines_new, sll_c_bsplines

IMPLICIT NONE

PUBLIC 

TYPE :: t_rProfile
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER               :: n_knots !! number of knots, including repeated edge knots
  INTEGER               :: n_coefs !! number of basis coefficients
  REAL(wp), ALLOCATABLE :: knots(:)   !! knot values, includinng edge knots
  REAL(wp), ALLOCATABLE :: coefs(:)   !! basis coefficients
  INTEGER               :: profile_type
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER                           :: deg = 0!! degree of the B-splines
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl !! b-spline class
  
  CONTAINS
  ! general procedures
  PROCEDURE :: init_spline             => splProfile_init
  PROCEDURE :: eval_at_phi_norm        => rProfile_eval_at_phi_norm
  PROCEDURE :: eval_prime_at_phi_norm  => rProfile_eval_prime_at_phi_norm
  PROCEDURE :: eval_derivative         => rProfile_eval_derivative
  PROCEDURE :: eval_rho_derivative     => rProfile_eval_rho_derivative

  ! b-spline specific evaluation routines
  PROCEDURE :: splProfile_eval_at_phi_norm
  PROCEDURE :: splProfile_eval_prime_at_phi_norm
  PROCEDURE :: splProfile_eval_derivative
  
  ! polynomial specififc routine
  PROCEDURE :: polyProfile_eval_derivative

  ! hard coded derivatives with respect to rho=sqrt(phi/phi_edge)
  PROCEDURE :: rProfile_drho2
  PROCEDURE :: rProfile_drho3
  PROCEDURE :: rProfile_drho4

  FINAL     :: rProfile_free
END TYPE t_rProfile

interface t_rProfile
  module procedure :: rProfile_init
end interface t_rProfile

!===================================================================================================================================

CONTAINS

FUNCTION rProfile_init(profile_type, coefs, n_coefs, knots, n_knots) RESULT(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER             :: profile_type !! defines the type of profile: 0=power polynomial, 1=B-Splines
  REAL,    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients
  INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
! OPTIONAL
  REAL,    INTENT(IN), OPTIONAL :: knots(:)  !! knots of the B-Spline with repeated start and end points
  INTEGER, INTENT(IN), OPTIONAL :: n_knots !! number of knots
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_rProfile) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF ((.NOT. PRESENT(knots)) .AND. (profile_type==1)) THEN
    WRITE(UNIT_stdOut,'(A)')'WARNING: Profile type "b_spline" specified without providing B-spline knots! Interpreting profile coefficients as for profile type "power_poly".'
    profile_type = 0
  END IF

  sf%profile_type = profile_type

  SELECT CASE(profile_type)
  CASE(0)
    sf%deg   = n_coefs-1 ! multiplicity of the first knot determines the degree
    sf%n_coefs = n_coefs
    sf%n_knots = -1
    ALLOCATE(sf%coefs(1:n_coefs))
    sf%coefs = coefs
  CASE(1)
    CALL sf%init_spline(knots, n_knots, coefs, n_coefs)
  END SELECT
END FUNCTION rProfile_init

!===================================================================================================================================
!> evaluate the rProfile at position phi_norm
!!
!===================================================================================================================================
FUNCTION rProfile_eval_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
  USE MODgvec_Globals    ,ONLY: EVAL1DPOLY
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE (sf%profile_type)
  CASE(0)
    profile_value = EVAL1DPOLY(sf%n_coefs,sf%coefs,phi_norm)
  CASE(1)
    profile_value = sf%splProfile_eval_at_phi_norm(phi_norm)
  END SELECT
END FUNCTION rProfile_eval_at_phi_norm

!===================================================================================================================================
!> evaluate the first derivative of the rProfile at position phi_norm
!!
!===================================================================================================================================
FUNCTION rProfile_eval_prime_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
  USE MODgvec_Globals    ,ONLY: EVAL1DPOLY_DERIV
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE (sf%profile_type)
  CASE(0)
    profile_value = EVAL1DPOLY_DERIV(sf%n_coefs,sf%coefs,phi_norm)
  CASE(1)
    profile_value = sf%splProfile_eval_prime_at_phi_norm(phi_norm)
  END SELECT
END FUNCTION rProfile_eval_prime_at_phi_norm

!===================================================================================================================================
!> evaluate the n-th derivative of the rProfile at position phi_norm
!!
!===================================================================================================================================
FUNCTION rProfile_eval_derivative(sf, phi_norm, n ) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
  INTEGER, OPTIONAL    :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE (sf%profile_type)
  CASE(0)
    profile_prime_value = sf%polyProfile_eval_derivative(phi_norm, n)
  CASE(1)
    profile_prime_value = sf%splProfile_eval_derivative(phi_norm, n)
  END SELECT
END FUNCTION rProfile_eval_derivative

!===================================================================================================================================
!> initialize the rProfile of type bspline
!!
!===================================================================================================================================
SUBROUTINE splProfile_init(sf, knots, n_knots, coefs, n_coefs)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(INOUT)  :: sf !! self
  INTEGER, INTENT(IN) :: n_knots !! number of knots
  INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
  REAL,    INTENT(IN) :: knots(n_knots)  !! knots of the B-Spline with repeated start and end points
  REAL,    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                :: nth_deriv_at_axis
  INTEGER             :: i
!===================================================================================================================================
  __PERFON("splProfile_init")
  sf%deg   = COUNT(knots==knots(1))-1 ! multiplicity of the first knot determines the degree
  sf%n_knots = n_knots
  sf%n_coefs = n_coefs
  ALLOCATE(sf%knots(1:n_knots))
  sf%knots = knots
  ALLOCATE(sf%coefs(1:n_coefs))
  sf%coefs = coefs
  IF (sf%deg>0) THEN
    CALL sll_s_bsplines_new(sf%bspl, sf%deg, .FALSE., & 
                            sf%knots(1),sf%knots(n_knots),&
                            size(sf%knots(sf%deg+1:n_knots-sf%deg))-1 , & ! number of knots handed to the library
                            sf%knots(sf%deg+1:n_knots-sf%deg)) ! remove repeated edge knots
  END IF
  __PERFOFF("splProfile_init")

END SUBROUTINE splProfile_init

!===================================================================================================================================
!> evaluate the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: basis_values(sf%deg+1) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
  IF (sf%deg > 0) THEN
    CALL sf%bspl%eval_basis(phi_norm,basis_values,first_non_zero_bspl)
    profile_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*basis_values(:sf%deg+1))
  ELSE
    profile_value = sf%coefs(1)
  END IF
END FUNCTION splProfile_eval_at_phi_norm

!===================================================================================================================================
!> evaluate the first derivative of the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_prime_at_phi_norm( sf, phi_norm ) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: deriv_values(sf%deg+1) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
  IF (sf%deg>0) THEN
    CALL sf%bspl%eval_deriv(phi_norm,deriv_values,first_non_zero_bspl)
    profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*deriv_values(:sf%deg+1))
    !profile_prime_value = sf%eval_derivative(phi_norm, n=1)
  ELSE 
    profile_prime_value = 0.0_wp
  END IF
END FUNCTION splProfile_eval_prime_at_phi_norm

!===================================================================================================================================
!> evaluate the n-th derivative of the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_derivative( sf, phi_norm, n ) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
  INTEGER, OPTIONAL    :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp), ALLOCATABLE            :: deriv_values(:,:) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
    IF (.NOT.PRESENT(n)) THEN
      n = 0
    END IF

    ALLOCATE(deriv_values(n+1,sf%deg+1))
    IF (sf%deg>0) THEN
      CALL sf%bspl%eval_basis_and_n_derivs(phi_norm,n,deriv_values,first_non_zero_bspl)
      profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*deriv_values(n+1,:sf%deg+1))
    ELSE 
      profile_prime_value = 0.0_wp
    END IF
    SDEALLOCATE(deriv_values)
END FUNCTION splProfile_eval_derivative

!===================================================================================================================================
!> calculate the prefactor for the d-th coefficient of the n-th derivative of a polynomial
!!
!===================================================================================================================================
FUNCTION poly_derivative_prefactor(D,n) RESULT(prefactor)
  INTEGER, INTENT(IN) :: D,n
  INTEGER :: i
  REAL(wp) :: prefactor
  prefactor = 1.0_wp
  DO i=D-n+1,D
    prefactor = prefactor*i
  END DO
END FUNCTION poly_derivative_prefactor

!===================================================================================================================================
!> evaluate the n-th derivative of phi_norm with respect to rho=sqrt(phi/phi_edge)
!!
!===================================================================================================================================
FUNCTION phi_norm_derivative(spos,n) RESULT(phi_norm_prime)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN) :: spos
  INTEGER, INTENT(IN)  :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) :: phi_norm_prime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF (n>2) THEN
    phi_norm_prime = 0.0_wp
  ELSE
    phi_norm_prime = poly_derivative_prefactor(2,n)*spos**(2-n)
  END IF
END FUNCTION phi_norm_derivative

!===================================================================================================================================
!> evaluate the n-th derivative of a power polynomial
!!
!===================================================================================================================================
FUNCTION polyProfile_eval_derivative(sf, phi_norm, n) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
  INTEGER            , INTENT(IN)  :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: prefactor
  INTEGER                          :: d
!===================================================================================================================================
  IF (n>sf%deg) THEN
    profile_prime_value = 0.0_wp
  ELSE 
    profile_prime_value = 0.0_wp
    DO d=n,sf%deg
      prefactor=poly_derivative_prefactor(d,n)
      profile_prime_value = profile_prime_value +prefactor*sf%coefs(d+1)*(phi_norm**(d-n))
    END DO
  END IF
END FUNCTION polyProfile_eval_derivative

!===================================================================================================================================
!> evaluate the n-th derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! NOTE: n has to be in [0,4] due to a lazy implementation of the product rule.
!===================================================================================================================================
FUNCTION rProfile_eval_rho_derivative(sf, spos, n) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp), INTENT(IN) :: spos
  INTEGER,  INTENT(IN) :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) :: derivative
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm = phi_norm_derivative(spos,n=0)
  SELECT CASE(n)
  CASE(0)
    derivative = sf%eval_at_phi_norm(phi_norm)
  CASE(1)
    derivative = sf%eval_prime_at_phi_norm(phi_norm)*phi_norm_derivative(spos,n=1)
  CASE(2)
    derivative = sf%rProfile_drho2(spos)
  CASE(3)
    derivative = sf%rProfile_drho3(spos)
  CASE(4)
    derivative = sf%rProfile_drho4(spos)
  CASE DEFAULT
    CALL abort(__STAMP__,&
        "error in rprofile: derivatives higher than 4 with respect to rho=sqrt(phi/phi_edge) are not implemented!")
  END SELECT
END FUNCTION rProfile_eval_rho_derivative

! The following functions hardcode the first iterations of Faà di Brunos formula.
! If at any point phi_norm is not spos**2 one would need to adapt phi_norm_derivative and these should still be valid.

!===================================================================================================================================
!> evaluate the 2nd derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! 
!===================================================================================================================================
FUNCTION rProfile_drho2(sf, spos) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp), INTENT(IN) :: spos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) :: derivative
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm = phi_norm_derivative(spos,n=0)
  ! d^2/dx^2 f(g(x)) = [g'(x)]^2 * f''(g(x))+g''(x) * f'(g(x))
  derivative = phi_norm_derivative(spos,n=1)**2*sf%eval_derivative(phi_norm, n=2) &
             + phi_norm_derivative(spos,n=2)*sf%eval_derivative(phi_norm, n=1)
END FUNCTION rProfile_drho2

!===================================================================================================================================
!> evaluate the 3rd derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! 
!===================================================================================================================================
FUNCTION rProfile_drho3(sf, spos) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp), INTENT(IN) :: spos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) :: derivative
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm = phi_norm_derivative(spos,n=0)
  ! d^3/dx^3 f(g(x)) = 3g'(x) * f''(g(x)) * g''(x) + [g'(x)]^3*f'''(g(x)) + f'(x) * g'''(x)
  derivative = 3*phi_norm_derivative(spos,n=1)*sf%eval_derivative(phi_norm, n=2)*phi_norm_derivative(spos,n=2) &
             + phi_norm_derivative(spos,n=1)**3*sf%eval_derivative(phi_norm, n=3) &
             + phi_norm_derivative(spos,n=3)*sf%eval_derivative(phi_norm, n=1)
END FUNCTION rProfile_drho3

!===================================================================================================================================
!> evaluate the 4th derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! 
!===================================================================================================================================
FUNCTION rProfile_drho4(sf, spos) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_rProfile), INTENT(IN)  :: sf !! self
  REAL(wp), INTENT(IN) :: spos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) :: derivative
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm = phi_norm_derivative(spos,n=0)
  ! d^4/dx^4 f(g(x)) = f''''(g(x))g'(x)**4 
  !                   + 6f'''(g(x))g''(x)g'(x)^2
  !                   + 3f''(g(x))g''(x)^2 
  !                   + 4f''(g(x))g'''(x)g'(x)
  !                   + f'(g(x))g''''(x)
  derivative = sf%eval_derivative(phi_norm, n=4)*phi_norm_derivative(spos,n=1)**4 &
             + 6*sf%eval_derivative(phi_norm, n=3)*phi_norm_derivative(spos,n=2)*phi_norm_derivative(spos,n=1)**2 &
             + 3*sf%eval_derivative(phi_norm, n=2)*phi_norm_derivative(spos,n=2)**2 &
             + 4*sf%eval_derivative(phi_norm, n=2)*phi_norm_derivative(spos,n=3)*phi_norm_derivative(spos,n=1) &
             + sf%eval_derivative(phi_norm, n=1)*phi_norm_derivative(spos,n=4)
END FUNCTION rProfile_drho4

!===================================================================================================================================
!> finalize the type rProfile
!!
!===================================================================================================================================
SUBROUTINE rProfile_free(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_rProfile), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!=================================================================================================================================== 
  SDEALLOCATE(sf%knots)
  SDEALLOCATE(sf%coefs)
  IF (ALLOCATED(sf%bspl)) CALL sf%bspl%free()
  SDEALLOCATE(sf%bspl)
  sf%deg          =-1
  sf%n_coefs      =-1
  sf%n_knots      =-1
  sf%profile_type =-1
END SUBROUTINE rProfile_free


END MODULE MODgvec_rProfile

