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
  
  CONTAINS
  
  PROCEDURE :: init                    => polyProfile_init
  PROCEDURE :: eval_at_phi_norm        => polyProfile_eval_at_phi_norm
  PROCEDURE :: eval_prime_at_phi_norm  => polyProfile_eval_prime_at_phi_norm
  PROCEDURE :: eval_derivative         => polyProfile_eval_derivative

  PROCEDURE :: free => polyProfile_free
  
END TYPE t_rProfile_poly

CONTAINS

SUBROUTINE polyProfile_init(sf, coefs, n_coefs)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL,    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients
    INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(OUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
    sf%profile_type = 0
    sf%deg   = n_coefs-1
    sf%n_coefs = n_coefs
    ALLOCATE(sf%coefs(1:n_coefs))
    sf%coefs = coefs
END SUBROUTINE polyProfile_init

!===================================================================================================================================
!> evaluate the polyProfile at position phi_norm
!!
!===================================================================================================================================
FUNCTION polyProfile_eval_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
    USE MODgvec_Globals    ,ONLY: EVAL1DPOLY
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(IN)  :: sf !! self
    REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

    profile_value = EVAL1DPOLY(sf%n_coefs,sf%coefs,phi_norm)

END FUNCTION polyProfile_eval_at_phi_norm

!===================================================================================================================================
!> evaluate the first derivative of the polyProfile at position phi_norm
!!
!===================================================================================================================================
FUNCTION polyProfile_eval_prime_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
! MODULES
    USE MODgvec_Globals    ,ONLY: EVAL1DPOLY_DERIV
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(IN)  :: sf !! self
    REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp)                         :: profile_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
    profile_value = EVAL1DPOLY_DERIV(sf%n_coefs,sf%coefs,phi_norm)

END FUNCTION polyProfile_eval_prime_at_phi_norm

!===================================================================================================================================
!> evaluate the n-th derivative of a power polynomial
!!
!===================================================================================================================================
FUNCTION polyProfile_eval_derivative(sf, phi_norm, n) RESULT(profile_prime_value)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(IN)  :: sf !! self
    REAL(wp)           , INTENT(IN)  :: phi_norm !! evaluation point in the toroidal flux coordinate (phi_norm=phi/phi_edge= spos^2)
    INTEGER            , OPTIONAL  :: n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
    REAL(wp)                         :: prefactor
    INTEGER                          :: d
!===================================================================================================================================
    IF (.NOT.PRESENT(n)) THEN
        n = 0
    END IF
    
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
!> finalize the type rProfile
!!
!===================================================================================================================================
SUBROUTINE polyProfile_free(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    CLASS(t_rProfile_poly), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
!=================================================================================================================================== 
    SDEALLOCATE(sf%coefs)
    sf%deg          =-1
    sf%n_coefs      =-1
    sf%profile_type =-1
END SUBROUTINE polyProfile_free

END MODULE MODgvec_rProfile_poly