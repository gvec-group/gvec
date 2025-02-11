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
    INTEGER               :: deg = 0

    contains

    PROCEDURE(i_fun_eval_at_rho2      ), DEFERRED :: eval_at_rho2

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
        INTEGER          , OPTIONAL    :: deriv
        REAL(wp)                       :: profile_value
    END FUNCTION i_fun_eval_at_rho2

END INTERFACE

CONTAINS
!===================================================================================================================================
!> calculate the prefactor for the d-th coefficient of the n-th derivative of a polynomial
!!
!===================================================================================================================================
FUNCTION poly_derivative_prefactor(D,deriv) RESULT(prefactor)
    INTEGER, INTENT(IN) :: D,deriv
    INTEGER :: i
    REAL(wp) :: prefactor
    prefactor = 1.0_wp
    DO i=D-deriv+1,D
        prefactor = prefactor*i
    END DO
END FUNCTION poly_derivative_prefactor

!===================================================================================================================================
!> evaluate the n-th derivative of rho2 with respect to rho=sqrt(phi/phi_edge)
!!
!===================================================================================================================================
FUNCTION rho2_derivative(spos,deriv) RESULT(rho2_prime)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL(wp), INTENT(IN) :: spos
    INTEGER, INTENT(IN)  :: deriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp) :: rho2_prime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
    IF (deriv>2) THEN
        rho2_prime = 0.0_wp
    ELSE
        rho2_prime = poly_derivative_prefactor(2,deriv)*spos**(2-deriv)
    END IF
END FUNCTION rho2_derivative

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
    REAL(wp) :: rho2
!===================================================================================================================================
    rho2 = rho2_derivative(spos,deriv=0)
    ! d^2/dx^2 f(g(x)) = [g'(x)]^2 * f''(g(x))+g''(x) * f'(g(x))
    derivative = rho2_derivative(spos,deriv=1)**2*sf%eval_at_rho2(rho2, deriv=2) &
                + rho2_derivative(spos,deriv=2)*sf%eval_at_rho2(rho2, deriv=1)
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
    REAL(wp) :: rho2
!===================================================================================================================================
    rho2 = rho2_derivative(spos,deriv=0)
    ! d^3/dx^3 f(g(x)) = 3g'(x) * f''(g(x)) * g''(x) + [g'(x)]^3*f'''(g(x)) + f'(x) * g'''(x)
    derivative = 3*rho2_derivative(spos,deriv=1)*sf%eval_at_rho2(rho2, deriv=2)*rho2_derivative(spos,deriv=2) &
                + rho2_derivative(spos,deriv=1)**3*sf%eval_at_rho2(rho2, deriv=3) &
                + rho2_derivative(spos,deriv=3)*sf%eval_at_rho2(rho2, deriv=1)
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
    REAL(wp) :: rho2
!===================================================================================================================================
    rho2 = rho2_derivative(spos,deriv=0)
    ! d^4/dx^4 f(g(x)) = f''''(g(x))g'(x)**4 
    !                   + 6f'''(g(x))g''(x)g'(x)^2
    !                   + 3f''(g(x))g''(x)^2 
    !                   + 4f''(g(x))g'''(x)g'(x)
    !                   + f'(g(x))g''''(x)
    derivative = sf%eval_at_rho2(rho2, deriv=4)*rho2_derivative(spos,deriv=1)**4 &
                + 6*sf%eval_at_rho2(rho2, deriv=3)*rho2_derivative(spos,deriv=2)*rho2_derivative(spos,deriv=1)**2 &
                + 3*sf%eval_at_rho2(rho2, deriv=2)*rho2_derivative(spos,deriv=2)**2 &
                + 4*sf%eval_at_rho2(rho2, deriv=2)*rho2_derivative(spos,deriv=3)*rho2_derivative(spos,deriv=1) &
                + sf%eval_at_rho2(rho2, deriv=1)*rho2_derivative(spos,deriv=4)
END FUNCTION rProfile_drho4

!===================================================================================================================================
!> evaluate the n-th derivative of a radial profile with respect to rho= sqrt(phi/phi-edge).
!! NOTE: n has to be in [0,4] due to a lazy implementation of the product rule.
!===================================================================================================================================
FUNCTION rProfile_eval_at_rho(sf, spos, deriv) RESULT(derivative)
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CLASS(c_rProfile)  :: sf !! self
    REAL(wp), INTENT(IN) :: spos
    INTEGER,  OPTIONAL   :: deriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp) :: derivative
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
    REAL(wp) :: rho2
!===================================================================================================================================
    IF (.NOT.PRESENT(deriv)) THEN
        deriv = 0
    END IF 

    rho2 = rho2_derivative(spos,deriv=0)
    SELECT CASE(deriv)
    CASE(0)
        derivative = sf%eval_at_rho2(rho2, deriv=0)
    CASE(1)
        derivative = sf%eval_at_rho2(rho2,deriv=1)*rho2_derivative(spos,deriv=1)
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
END FUNCTION rProfile_eval_at_rho

END MODULE MODgvec_rProfile_base