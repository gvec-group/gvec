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
    INTEGER               :: profile_type
    INTEGER               :: deg = 0

    contains

    PROCEDURE(i_fun_eval_at_phi_norm      ), DEFERRED :: eval_at_phi_norm
    PROCEDURE(i_fun_eval_prime_at_phi_norm), DEFERRED :: eval_prime_at_phi_norm
    PROCEDURE(i_fun_eval_derivative       ), DEFERRED :: eval_derivative
    PROCEDURE(i_sub_free                  ), DEFERRED :: free

    PROCEDURE :: eval_rho_derivative => rProfile_eval_rho_derivative
    ! hard coded derivatives with respect to rho=sqrt(phi/phi_edge)
    PROCEDURE, PRIVATE :: rProfile_drho2
    PROCEDURE, PRIVATE :: rProfile_drho3
    PROCEDURE, PRIVATE :: rProfile_drho4

end type c_rProfile

ABSTRACT INTERFACE

    FUNCTION i_fun_eval_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
        IMPORT c_rProfile
        IMPORT wp
        CLASS(c_rProfile), INTENT(IN)  :: sf
        REAL(wp)           , INTENT(IN)  :: phi_norm 
        REAL(wp)                         :: profile_value
    END FUNCTION i_fun_eval_at_phi_norm

    FUNCTION i_fun_eval_prime_at_phi_norm( sf, phi_norm ) RESULT(profile_value)
        IMPORT c_rProfile
        IMPORT wp
        CLASS(c_rProfile), INTENT(IN) :: sf
        REAL(wp)         , INTENT(IN)   :: phi_norm 
        REAL(wp)                        :: profile_value
    END FUNCTION i_fun_eval_prime_at_phi_norm

    FUNCTION i_fun_eval_derivative(sf, phi_norm, n ) RESULT(profile_prime_value)
        IMPORT c_rProfile
        IMPORT wp
        CLASS(c_rProfile), INTENT(IN) :: sf
        REAL(wp)           , INTENT(IN) :: phi_norm
        INTEGER, OPTIONAL    :: n
        REAL(wp)                        :: profile_prime_value
    END FUNCTION i_fun_eval_derivative

    SUBROUTINE i_sub_free(sf)
        IMPORT c_rProfile
        CLASS(c_rProfile), INTENT(INOUT) :: sf
    END SUBROUTINE i_sub_free

END INTERFACE

CONTAINS
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
!> evaluate the 2nd derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! 
!===================================================================================================================================
FUNCTION rProfile_drho2(sf, spos) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CLASS(c_rProfile)  :: sf !! self
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
    CLASS(c_rProfile)  :: sf !! self
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
    CLASS(c_rProfile) :: sf !! self
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
!> evaluate the n-th derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! NOTE: n has to be in [0,4] due to a lazy implementation of the product rule.
!===================================================================================================================================
FUNCTION rProfile_eval_rho_derivative(sf, spos, n) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CLASS(c_rProfile)  :: sf !! self
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

END MODULE MODgvec_rProfile_base