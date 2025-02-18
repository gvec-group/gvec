!===================================================================================================================================
! Copyright (C) 2025 Robert Koeberl <robert.koeberl@ipp.mpg.de>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================

#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** bspline_interpolation **
!!
!! Routines for interpolating 1D functions via B-splines.
!===================================================================================================================================
MODULE MODgvec_bspline_interpolation
! MODULES
USE MODgvec_Globals ,ONLY: wp, abort
USE MODgvec_LinAlg, ONLY: SOLVE
USE sll_m_bsplines  ,ONLY: sll_s_bsplines_new, sll_c_bsplines
IMPLICIT NONE

PUBLIC
 CONTAINS

! This tries to replicate the scipy implementation 
 FUNCTION not_a_knot(x, k) RESULT(knots)
    REAL(wp), INTENT(IN)  :: x(:) ! interpolation points x position
    INTEGER,  INTENT(IN)  :: k ! B-splines degree

    INTEGER               :: k2
    INTEGER               :: n, n_aux, n_internal, n_knots
    REAL(wp), ALLOCATABLE :: knots(:)
    REAL(wp), ALLOCATABLE :: knots_aux(:)

    n = SIZE(x)

    IF (MOD(k,2).EQ.1) THEN
        k2 = (k + 1) / 2
        ALLOCATE(knots_aux(n))
        knots_aux = x
    ELSE
        k2 = k/2
        ALLOCATE(knots_aux(n-1))
        knots_aux = (x(2:)-x(:n-1)) / 2
    END IF
 
    ALLOCATE(knots(SIZE(knots_aux)-2*k2 + 2*(k+1)))

    n_knots = size(knots)
    n_aux = size(knots_aux)
    n_internal = n_aux-2*k2
    knots(k+2:n_knots-k-1) =  knots_aux(k2+1:n_aux-k2)
    knots(:k+1) = x(1)
    knots(k+2+n_internal:) = x(n)
    DEALLOCATE(knots_aux)
END FUNCTION not_a_knot

SUBROUTINE SOLVE_LSTSQ(A, b, c)
    EXTERNAL DGELS 
    REAL(wp), INTENT(IN) :: A(:,:)
    REAL(wp), INTENT(IN) :: b(:)
    REAL(wp), INTENT(INOUT) :: c(:)
    INTEGER :: n_rows, n_cols, info, MN, LWORK
    REAL(wp), ALLOCATABLE :: WORK(:)

    n_rows = SIZE(A,1)
    n_cols = SIZE(A,2)
    MN = MIN(n_rows,n_cols)
    LWORK = MAX(1,MN+MAX(MN,1))
    ALLOCATE(WORK(LWORK))
    CALL DGELS("N", n_rows, n_cols, 1, A, n_rows, b,SIZE(b,1),WORK,LWORK,info)
    IF(info.NE.0)THEN
        CALL abort(__STAMP__,&
                   'Matrix solve for LSTSQ does not work!')
    END IF
    SDEALLOCATE(WORK)
    c=b

END SUBROUTINE SOLVE_LSTSQ

SUBROUTINE interpolate_not_a_knot(x, y, c, knots, deg)
    REAL(wp), INTENT(IN)                 :: x(:), y(:)
    REAL(wp), INTENT(INOUT), ALLOCATABLE :: c(:), knots(:)

    CLASS(sll_c_bsplines), ALLOCATABLE :: bspl
    INTEGER, OPTIONAL, INTENT(IN) :: deg
    INTEGER                  :: n_basis, n_x, n_knots, i, jmin, k
    REAL(wp), ALLOCATABLE    :: N(:,:)
    REAL(wp), ALLOCATABLE    :: basis_values(:)

    IF (PRESENT(deg)) THEN
        k = deg
    ELSE
        k = 3
    END IF
  
    ALLOCATE(basis_values(k+1))

    knots = not_a_knot(x,k)
    n_x = SIZE(x)
    n_basis = SIZE(knots)-k-1
    n_knots = SIZE(knots)
    IF (n_x .NE. n_basis) THEN
        CALL abort(__STAMP__,&
        'Profile interpolation sanity check failed!')
    END IF

    ALLOCATE(c(n_x))
    c = 0.0_wp
    ALLOCATE(N(n_basis,n_basis))
    N = 0.0_wp

    CALL sll_s_bsplines_new(bspl, k, .FALSE., xmin=x(1), xmax=x(n_x), &
    ncells=SIZE(knots(k+1:n_knots-k))-1, breaks=knots(k+1:n_knots-k))
    DO i=1,n_basis
        CALL bspl%eval_basis(x(i),basis_values, jmin)
        N(i,jmin:jmin+k) = basis_values
    END DO
    c = SOLVE(N,y)
    ! CALL SOLVE_LSTSQ(N,y,c)
    CALL bspl%free()
    SDEALLOCATE(bspl)
    SDEALLOCATE(basis_values)
    SDEALLOCATE(N)

END SUBROUTINE interpolate_not_a_knot

SUBROUTINE interpolate_complete_bspl(x, y, c, knots, y_BC)
    REAL(wp), INTENT(IN)                 :: x(:), y(:)
    REAL(wp), INTENT(INOUT), ALLOCATABLE :: c(:)
    REAL(wp), INTENT(INOUT), ALLOCATABLE :: knots(:)

    CLASS(sll_c_bsplines), ALLOCATABLE :: bspl
    REAL(wp)                           :: y_BC(2)
    INTEGER                            :: k, n_x, n_basis, n_knots, i, jmin
    REAL(wp), ALLOCATABLE              :: N(:,:)
    REAL(wp), ALLOCATABLE              :: basis_values(:), RHS(:)

    k = 3
    n_x = SIZE(x)
    n_knots = n_x+2*k

    ALLOCATE(c(n_x+2))
    c = 0.0_wp
    ALLOCATE(knots(n_knots))
    knots = 0.0_wp
    ALLOCATE(basis_values(k+1))
    basis_values = 0.0_wp
    ALLOCATE(RHS(n_x+2))
    RHS = 0.0_wp
    
    knots(k+1:n_x+k) = x(:)
    knots(:k) = x(1)
    knots(n_x+k+1:) = x(n_x)

    CALL sll_s_bsplines_new(bspl, k, .FALSE., xmin=x(1), xmax=x(n_x), &
    ncells=SIZE(knots(k+1:n_knots-k))-1, breaks=knots(k+1:n_knots-k))

    n_basis = SIZE(knots)-k-1
    ALLOCATE(N(n_basis,n_basis))
    N = 0.0_wp
    DO i=1,n_x
        CALL bspl%eval_basis(x(i),basis_values, jmin)
        N(i+1,jmin:jmin+k) = basis_values
    END DO

    CALL bspl%eval_deriv(x(1),basis_values, jmin)
    N(1,jmin:jmin+k) = basis_values

    CALL bspl%eval_deriv(x(n_x),basis_values, jmin)
    N(n_basis,jmin:jmin+k) = basis_values
    
    RHS(2:n_x+1) = y
    RHS(1) = y_BC(1)
    RHS(n_x+2) = y_BC(2)

    c = SOLVE(N, RHS)
    ! CALL SOLVE_LSTSQ(N,RHS,c)

    SDEALLOCATE(basis_values)
    SDEALLOCATE(N)
    SDEALLOCATE(RHS)
    CALL bspl%free()
    SDEALLOCATE(bspl)

END SUBROUTINE interpolate_complete_bspl

END MODULE MODgvec_bspline_interpolation