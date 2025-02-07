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
MODULE MODgvec_rProfile
! MODULES
USE MODgvec_Globals,       ONLY: wp,UNIT_stdOut
USE MODgvec_rProfile_base, ONLY: c_rProfile
USE MODgvec_rProfile_poly, ONLY: t_rProfile_poly
USE MODgvec_rProfile_bspl, ONLY: t_rProfile_bspl

IMPLICIT NONE

PUBLIC 

CONTAINS

SUBROUTINE rProfile_init(sf, profile_type, coefs, n_coefs, knots, n_knots)
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
    CLASS(c_rProfile), ALLOCATABLE, INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
    IF ((.NOT. PRESENT(knots)) .AND. (profile_type==1)) THEN
        WRITE(UNIT_stdOut,'(A)')'WARNING: Profile type "b_spline" specified without providing B-spline knots! Interpreting profile coefficients as for profile type "power_poly".'
        profile_type = 0
    END IF

    SELECT CASE(profile_type)
    CASE(0)
        ALLOCATE(t_rProfile_poly :: sf )
    CASE(1)
        ALLOCATE(t_rProfile_bspl :: sf )
    END SELECT

    SELECT TYPE( sf )
    TYPE IS( t_rProfile_poly )
        CALL sf%init(coefs, n_coefs)
    TYPE IS( t_rProfile_bspl )
        CALL sf% init(knots, n_knots, coefs, n_coefs)
    END SELECT
END SUBROUTINE rProfile_init


END MODULE MODgvec_rProfile
