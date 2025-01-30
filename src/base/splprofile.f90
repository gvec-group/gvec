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

  PROCEDURE :: eval_at_rho        => splProfile_eval_at_rho
  PROCEDURE :: eval_prime_at_rho  => splProfile_eval_prime_at_rho
  PROCEDURE :: eval_derivative    => splProfile_eval_derivative
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
! LOCAL VARIABLES
  REAL                :: nth_deriv_at_axis
  INTEGER             :: i
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
  IF (sf%deg>0) THEN
    CALL sll_s_bsplines_new(sf%bspl, sf%deg, .FALSE., & 
                            sf%knots(1),sf%knots(n_knots),&
                            size(sf%knots(sf%deg+1:n_knots-sf%deg))-1 , & ! number of knots handed to the library
                            sf%knots(sf%deg+1:n_knots-sf%deg)) ! remove repeated edge knots
    ! Test if the first odd derivatives are zero at the magnetic axis to avoid discontinuous profiles
    DO i = 1,sf%deg-1,2
      nth_deriv_at_axis = sf%eval_derivative(0.0_wp,i)
      IF (ABS(nth_deriv_at_axis) > 1E-12) THEN
        WRITE(UNIT_stdOut,'(6X,A,I4,A,3E11.4)')'WARNING: n-th odd derivative at axis of spline profile is not zero. n= ',i," derivative value:",nth_deriv_at_axis
      END IF
    END DO
  END IF
  __PERFOFF("splProfile_init")

END FUNCTION splProfile_init

!===================================================================================================================================
!> evaluate the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_at_rho( sf, rho ) RESULT(profile_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_splProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: rho !! evaluation point in the toroidal flux coordinate (rho=sqrt(phi/phi_edge)= spos)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: basis_values(sf%deg+1) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
  IF (sf%deg > 0) THEN
    CALL sf%bspl%eval_basis(rho,basis_values,first_non_zero_bspl)
    profile_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*basis_values(:sf%deg+1))
    !profile_value = sf%eval_derivative(rho, n=0)
  ELSE
    profile_value = sf%coefs(1)
  END IF
END FUNCTION splProfile_eval_at_rho

!===================================================================================================================================
!> evaluate the first derivative of the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_prime_at_rho( sf, rho ) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_splProfile), INTENT(IN)  :: sf !! self
  REAL(wp)           , INTENT(IN)  :: rho !! evaluation point in the toroidal flux coordinate (spos=sqrt(phi/phi_edge) = rho)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                         :: deriv_values(sf%deg+1) !! values of the (deg+1) B-splines that contribute at s_pos
  INTEGER                          :: first_non_zero_bspl !! index offset for the coefficients
!===================================================================================================================================
  IF (sf%deg>0) THEN
    CALL sf%bspl%eval_deriv(rho,deriv_values,first_non_zero_bspl)
    profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*deriv_values(:sf%deg+1))
    !profile_prime_value = sf%eval_derivative(rho, n=1)
  ELSE 
    profile_prime_value = 0.0_wp
  END IF
END FUNCTION splProfile_eval_prime_at_rho

!===================================================================================================================================
!> evaluate the n-th derivative of the splProfile at position s
!!
!===================================================================================================================================
FUNCTION splProfile_eval_derivative( sf, rho, n ) RESULT(profile_prime_value)
  ! MODULES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_splProfile), INTENT(IN)  :: sf !! self
    REAL(wp)           , INTENT(IN)  :: rho !! evaluation point in the toroidal flux coordinate (s=sqrt(phi/phi_edge)=rho)
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
      CALL sf%bspl%eval_basis_and_n_derivs(rho,n,deriv_values,first_non_zero_bspl)
      profile_prime_value =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%deg)*deriv_values(n+1,:sf%deg+1))
    ELSE 
      profile_prime_value = 0.0_wp
    END IF
    SDEALLOCATE(deriv_values)
END FUNCTION splProfile_eval_derivative

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

