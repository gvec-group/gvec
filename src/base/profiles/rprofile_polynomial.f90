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
  
  PROCEDURE :: eval_at_rho2        => polyProfile_eval_at_rho2
  !FINAL     :: polyProfile_free !! no finalizer necessary
  
END TYPE t_rProfile_poly

INTERFACE t_rProfile_poly
    MODULE PROCEDURE polyProfile_new
END INTERFACE t_rProfile_poly

CONTAINS

FUNCTION polyProfile_new(coefs, n_coefs) RESULT(sf)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL,    INTENT(IN) :: coefs(n_coefs)  !! B-Spline coefficients
    INTEGER, INTENT(IN) :: n_coefs !! number of coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    TYPE(t_rProfile_poly) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
    sf%deg   = n_coefs-1
    sf%n_coefs = n_coefs
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
    INTEGER               , OPTIONAL    :: deriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp)                         :: profile_prime_value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
    REAL(wp)                         :: prefactor
    INTEGER                          :: d
!===================================================================================================================================
    IF (.NOT.PRESENT(deriv)) THEN
        deriv = 0
    END IF
    
    IF (deriv>sf%deg) THEN
        profile_prime_value = 0.0_wp
    ELSE IF (deriv==0) THEN
        profile_prime_value = EVAL1DPOLY(sf%n_coefs, sf%coefs, rho2)
    ELSE IF (deriv==1) THEN
        profile_prime_value = EVAL1DPOLY_deriv(sf%n_coefs, sf%coefs, rho2)
    ELSE
        profile_prime_value = 0.0_wp
        DO d=deriv,sf%deg
            prefactor=poly_derivative_prefactor(d,deriv)
            profile_prime_value = profile_prime_value +prefactor*sf%coefs(d+1)*(rho2**(d-deriv))
        END DO
    END IF
END FUNCTION polyProfile_eval_at_rho2

END MODULE MODgvec_rProfile_poly
