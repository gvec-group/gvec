!===================================================================================================================================
! Copyright (C) 2025 -  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module ** cubic spline **
!!
!! Define a cubic spline by points and function values, with boundary conditions and represent it as a B-Spline for evaluation
!!
!===================================================================================================================================
MODULE MODgvec_cubic_spline
  ! MODULES
  USE MODgvec_globals, ONLY: wp,abort
  use sll_m_bsplines, only: sll_c_bsplines
  IMPLICIT NONE

  PUBLIC

  TYPE :: t_cubspl
    !--------------------------------------------------------------------------------------------------------------------------------
    REAL(wp),ALLOCATABLE  :: coefs(:)   !! B-Spline coefficients
    REAL(wp),ALLOCATABLE  :: knots(:)   !! knots (=break points with repeated knots at the end)
    CLASS(sll_c_bsplines),ALLOCATABLE :: bspl !! b-spline class
    CONTAINS

    PROCEDURE :: eval  => cubspl_eval
    PROCEDURE :: eval_at_rho  => cubspl_eval  !! for testing
    FINAL :: cubspl_free

  END TYPE t_cubspl

  INTERFACE t_cubspl
    MODULE PROCEDURE cubspl_new
  END INTERFACE t_cubspl

  INTERFACE interpolate_cubic_spline
    MODULE PROCEDURE interpolate_cubic_spline
  END INTERFACE

  CONTAINS

  !===================================================================================================================================
  !> Interpolation of function values f(x_i)=f_i, i=1,n with a cubic spline, given left and right boundary condition
  !! types of boundary conditions:
  !! 0: not-a-knot
  !! 1: f'(x_boundary)=0
  !! 2: f''(x_boundary)=0
  !!
  !===============================================================================================================================
  FUNCTION cubspl_new(x,f,BC,BC_val) RESULT(sf)
    ! MODULES
    USE sll_m_bsplines, ONLY: sll_s_bsplines_new
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    REAL(wp) , INTENT(IN) :: x(:) !! x positions
    REAL(wp) , INTENT(IN) :: f(:) !! function values at x positions
    INTEGER  , INTENT(IN) :: BC(1:2) !! Boundary condition at x(1)/x(n): =0: not-a-knot, =1: first der. =BC_val(1)/BC_val(2), =2: second der. =BC_val(1)/BC_val(2)
    REAL(wp) , INTENT(IN),OPTIONAL :: BC_val(1:2) !! Boundary value for BC(1:2) >0,
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    TYPE(t_cubspl) :: sf !! self
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER, PARAMETER:: deg=3 !! degree of the spline
    INTEGER :: nbreaks
    !===================================================================================================================================

    CALL interpolate_cubic_spline(x,f,sf%coefs,sf%knots,BC,BC_val)
    nbreaks=SIZE(sf%knots)-2*deg
    CALL sll_s_bsplines_new(sf%bspl, deg, .FALSE., sf%knots(1+deg),sf%knots(nbreaks+deg),nbreaks-1,sf%knots(1+deg:nbreaks+deg))
  END FUNCTION cubspl_new


  SUBROUTINE interpolate_cubic_spline(x,f,coefs,knots,BC,BC_val)
    ! MODULES
    USE sll_m_bsplines, ONLY: sll_s_bsplines_new
    USE sll_m_spline_matrix,ONLY: sll_c_spline_matrix,sll_s_spline_matrix_new
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES

    REAL(wp) , INTENT(IN) :: x(:) !! x positions
    REAL(wp) , INTENT(IN) :: f(:) !! function values at x positions
    INTEGER  , INTENT(IN) :: BC(1:2) !! Boundary condition at x(1)/x(n): =0: not-a-knot, =1: first der. =BC_val(1)/BC_val(2), =2: second der. =BC_val(1)/BC_val(2)
    REAL(wp) , INTENT(IN),OPTIONAL :: BC_val(1:2) !! Boundary value for BC(1:2) >0,
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    REAL(wp),ALLOCATABLE,INTENT(INOUT) :: coefs(:)  !! B-Spline coefficients of interpolated cubic spline
    REAL(wp),ALLOCATABLE,INTENT(INOUT) :: knots(:)  !! B-Spline knots of interpolated cubic spline
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER :: innerpts(1:2),nbc_left,nbc_right,i,jmin,s,n,nbreaks,ncoefs
    INTEGER,PARAMETER :: deg=3
    REAL(wp):: base(0:deg),deriv_values(0:deg,0:deg)
    CLASS(sll_c_bsplines),ALLOCATABLE :: bspl
    CLASS(sll_c_spline_matrix),ALLOCATABLE :: solvemat
    !===================================================================================================================================
    n=SIZE(x)
    IF(SIZE(f).NE.n) THEN
      CALL abort(__STAMP__, &
                 "x and f must have the same size for cubic spline interpolation")
    END IF
    innerpts=(/2,n-1/)
    nbreaks=n
    ncoefs=n
    SELECT CASE(BC(1))
    CASE(0) !not-a-knot: leave out second knot
      nbreaks=nbreaks-1
      innerpts(1)=3
      nbc_left=0
    CASE(1,2)
      ncoefs=ncoefs+1
      nbc_left=1
    CASE DEFAULT
      CALL abort(__STAMP__, &
                 "BC_left not implemented")
    END SELECT
    SELECT CASE(BC(2))
    CASE(0) !not-a-knot: leave out second to last knot
      nbreaks=nbreaks-1
      innerpts(2)=n-2
      nbc_right=0
    CASE(1,2)
      ncoefs=ncoefs+1
      nbc_right=1
    CASE DEFAULT
      CALL abort(__STAMP__, &
                 "BC_right not implemented")
    END SELECT
    IF(n+nbc_left+nbc_right.LT.4) THEN
      CALL abort(__STAMP__, &
                 "minimum number of conditions of a cubic spline is 4")
    END IF
    ALLOCATE(coefs(ncoefs),knots(1:nbreaks+2*deg))
    knots = -42e5
    knots(1:1+deg)=x(1) !repreat first knot deg+1 times
    knots(nbreaks+deg:nbreaks+2*deg)=x(n) !repeat last knot deg+1 times
    IF(innerpts(2)-innerpts(1)+1.GT.0)THEN
      knots(2+deg:nbreaks+deg-1)=x(innerpts(1):innerpts(2))
    END IF
    CALL sll_s_bsplines_new(bspl, deg, .FALSE., x(1),x(n),nbreaks-1,knots(1+deg:nbreaks+deg))

    IF(bspl%nBasis.NE.ncoefs) CALL abort(__STAMP__, &
                            'problem with bspl basis in cubic spline')
    CALL sll_s_spline_matrix_new(solvemat , "banded",ncoefs,deg,deg)
    ! Interpolation points
    do i = 1, n
      call bspl % eval_basis( x(i), base, jmin )
      coefs(i+nbc_left)= f(i)  !! used as rhs, then overwritten in solve
      do s = 0, deg
        call solvemat % set_element( i+nbc_left, jmin+s, base(s) )
      end do
    end do

    ! Boundary conditions
    IF(BC(1).GT.0)THEN !deriv=BC(1)
      CALL bspl%eval_basis_and_n_derivs(x(1),BC(1),deriv_values(0:BC(1),0:deg),jmin)
      do s = 0, deg
        call solvemat % set_element( 1, jmin+s, deriv_values(BC(1),s) )
      end do
      IF(PRESENT(BC_val))THEN
        coefs(1)=BC_val(1)! used as rhs, then overwritten in solve
      ELSE
        coefs(1)=0.0_wp ! used as rhs, then overwritten in solve
      END IF
    END IF
    IF(BC(2).GT.0)THEN !deriv=BC(2)
      CALL bspl%eval_basis_and_n_derivs(x(n),BC(2),deriv_values(0:BC(2),0:deg),jmin)
      do s = 0, deg
        call solvemat % set_element( ncoefs, jmin+s, deriv_values(BC(2),s) )
      end do
      IF(PRESENT(BC_val))THEN
        coefs(ncoefs)=BC_val(2)! used as rhs, then overwritten in solve
      ELSE
        coefs(ncoefs)=0.0_wp ! used as rhs, then overwritten in solve
      END IF
    END IF
    CALL solvemat % factorize()
    CALL solvemat%solve_inplace(1,coefs)

  END SUBROUTINE interpolate_cubic_spline

!===================================================================================================================================
!> evaluate the n-th derivative of the bsplProfile at position s
!!
!===================================================================================================================================
  FUNCTION cubspl_eval( sf, xpos, deriv ) RESULT(y)
    ! MODULES
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
      CLASS(t_cubspl), INTENT(IN)  :: sf !! self
      REAL(wp)       , INTENT(IN)  :: xpos(:) !! position
      INTEGER        , INTENT(IN)  :: deriv !! derivative (=0: no derivative)

    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
      REAL(wp)                     :: y(size(xpos))
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
      REAL(wp) :: deriv_values(0:deriv,0:sf%bspl%degree) !! values of the (deg+1) B-splines that contribute at s_pos
      INTEGER  :: i,first_non_zero_bspl !! index offset for the coefficients
    !===================================================================================================================================
      IF(deriv.EQ.0)THEN
        DO i=1,size(xpos)
          CALL sf%bspl%eval_basis(xpos(i),deriv_values(0,0:sf%bspl%degree),first_non_zero_bspl)
          y(i) =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%bspl%degree)*deriv_values(0,0:sf%bspl%degree))
        END DO
      ELSEIF(deriv.LE.sf%bspl%degree) THEN
        DO i=1,size(xpos)
          CALL sf%bspl%eval_basis_and_n_derivs(xpos(i),deriv,deriv_values,first_non_zero_bspl)
          y(i) =  SUM(sf%coefs(first_non_zero_bspl:first_non_zero_bspl+sf%bspl%degree)*deriv_values(deriv,0:sf%bspl%degree))
        END DO
      ELSE
        y = 0.0_wp
      END IF
    END FUNCTION cubspl_eval

    !===================================================================================================================================
    !> finalize the type rProfile
    !!
    !===================================================================================================================================
    SUBROUTINE cubspl_free(sf)
    ! MODULES
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
      TYPE(t_cubspl), INTENT(INOUT) :: sf !! self
    !-----------------------------------------------------------------------------------------------------------------------------------
    !===================================================================================================================================
      IF (ALLOCATED(sf%bspl)) CALL sf%bspl%free()
    END SUBROUTINE cubspl_free


END MODULE MODgvec_cubic_spline
