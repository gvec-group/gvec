!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **Power polynomial to B-Spline representation**
!!
!! CONTAINS  all functions to transform the polynomial coefficients into B-Spline coefficients via a special case of the Marsden's
!! identity. (see e.g. https://www.uio.no/studier/emner/matnat/math/MAT4170/v22/undervisningsmateriale/marsden_lin_indep.pdf)
!!
!===================================================================================================================================
MODULE MODgvec_PolyToBspl
! MODULES
USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep
IMPLICIT NONE
PUBLIC

!===================================================================================================================================

CONTAINS

FUNCTION factorial(n) RESULT(n_fac)
    INTEGER, INTENT(IN)  :: n
    INTEGER :: n_fac, i
    n_fac = 1.0_wp
    DO i = 2, n
        n_fac = n_fac * i
    END DO
END FUNCTION factorial

!===================================================================================================================================
!> Calculate the binomial coefficient n over r => n!/(r!(n-r)!) 
!===================================================================================================================================
FUNCTION binomial(n,r) RESULT(coef)
    INTEGER, INTENT(IN) :: n, r 
    REAL :: coef

    coef = factorial(n)/(factorial(r)*factorial(n-r))
END FUNCTION binomial

!===================================================================================================================================
!> Calculate the additional prefactor (c_jr) for the polynomial coefficients originating from the derivative of Marsden's identity:
!! x^r = SUM_j^n c_jr*B_jd(x)
!! NOTE: Since we only consider 1-element splines with x in [0,1], we only have edge knots that are either 0 or 1. Therefore the
!! sum over the knots can be reduced to a combinatory problem that counts the number non-zero permutations in the sum
!! => binomial(n,r). This implementation does NOT work for internal knots or polynomials outside of [0,1].
!===================================================================================================================================
FUNCTION dual_factor(r, d, j, t)
    INTEGER, INTENT(IN) :: r, d, j
    REAL, INTENT(IN) :: t(:)

    INTEGER :: n
    REAL :: dual_factor

    n = INT(SUM(t(j:j+d)))
    IF (r>n) THEN
        dual_factor = 0.0_wp
    ELSE
        dual_factor = binomial(n,r)/binomial(d,r)
    END IF
END FUNCTION dual_factor

!===================================================================================================================================
!> Calculate the j-th B-Spline coefficient
!===================================================================================================================================
FUNCTION poly2bspl(c, j, t) RESULT(c_bspl)
    REAL, INTENT(IN)    :: c(:)
    INTEGER, INTENT(IN) :: j
    REAL, INTENT(IN)    :: t(:)

    REAL :: c_bspl, c_jr
    INTEGER :: d, r

    c_bspl = 0.0_wp

    d = SIZE(c)-1

    c_bspl = 0
    DO r=1,d+1
        c_jr = dual_factor(r-1,d,j,t)
        c_bspl = c_bspl + c(r)*c_jr
    END DO

END FUNCTIOn poly2bspl

END MODULE MODgvec_PolyToBspl