!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
!
! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
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

!===================================================================================================================================
!>
!!# Module **NEWTON**
!!
!! Some simple Newton solvers
!!
!===================================================================================================================================
MODULE MOD_Newton
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PUBLIC

INTERFACE NewtonMin1D 
  MODULE PROCEDURE NewtonMin1D 
END INTERFACE

INTERFACE NewtonRoot1D 
  MODULE PROCEDURE NewtonRoot1D 
END INTERFACE

INTERFACE NewtonRoot1D_FdF
  MODULE PROCEDURE NewtonRoot1D_FdF
END INTERFACE

INTERFACE NewtonMin2D 
  MODULE PROCEDURE NewtonMin2D 
END INTERFACE

ABSTRACT INTERFACE 
  FUNCTION i_f1x1(x) RESULT (y1x1)
    IMPORT wp
    IMPLICIT NONE
    REAL(wp) :: x
    REAL(wp) :: y1x1
  END FUNCTION i_f1x1

  FUNCTION i_f2x1(x) RESULT (y2x1)
    IMPORT wp
    IMPLICIT NONE
    REAL(wp) :: x
    REAL(wp) :: y2x1(2)
  END FUNCTION i_f2x1

  FUNCTION i_f1x2(x) RESULT (y1x2)
    IMPORT wp
    IMPLICIT NONE
    REAL(wp) :: x(2)
    REAL(wp) :: y1x2
  END FUNCTION i_f1x2

  FUNCTION i_f2x2(x) RESULT (y2x2)
    IMPORT wp
    IMPLICIT NONE
    REAL(wp) :: x(2)
    REAL(wp) :: y2x2(2)
  END FUNCTION i_f2x2

  FUNCTION i_f22x2(x) RESULT (y22x2)
    IMPORT wp
    IMPLICIT NONE
    REAL(wp) :: x(2)
    REAL(wp) :: y22x2(2,2)
  END FUNCTION i_f22x2
END INTERFACE

!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Newton's iterative algorithm to find the minimimum of function f(x) in the interval [a,b], using df(x)=0 and the derivative 
!!
!===================================================================================================================================
FUNCTION NewtonMin1D(tol,a,b,x,FF,dFF,ddFF) RESULT (fmin)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)    :: tol  !! abort tolerance
REAL(wp),INTENT(IN)    :: a,b  !! search interval
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(INOUT) :: x    !! initial guess on input, result on output
PROCEDURE(i_f1x1)      :: FF   !! functional f(x) to minimize
PROCEDURE(i_f1x1)      :: dFF  !! d/dx f(x)
PROCEDURE(i_f1x1)      :: ddFF !! d^2/dx^2 f(x)
REAL(wp)               :: fmin !! on output =f(x)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)                :: x0
!===================================================================================================================================
x0=x
x=NewtonRoot1D(tol,a,b,x0,0.,dFF,ddFF)
fmin=FF(x)

END FUNCTION NewtonMin1D


!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function FR(x(:)) in the interval [a(:),b(:)], using d/dx(:)F(x)=0 and the derivative 
!!
!===================================================================================================================================
FUNCTION NewtonRoot1D(tol,a,b,xin,F0,FR,dFR) RESULT (xout)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: tol    !! abort tolerance
REAL(wp),INTENT(IN) :: a,b    !! search interval
REAL(wp),INTENT(IN) :: F0     !! function to find root is FR(x)-F0 
REAL(wp),INTENT(IN) :: xin    !! initial guess 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
PROCEDURE(i_f1x1)   :: FR      !! function to find root
PROCEDURE(i_f1x1)   :: dFR     !! multidimensional derivative d/dx f(x), size dim1
REAL(wp)            :: xout    !! on output =f(x)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL(wp)            :: x,dx
LOGICAL             :: converged
!===================================================================================================================================

converged=.FALSE.
x=xin
maxiter=50
DO iter=1,maxiter
  dx=-(FR(x)-F0)/dFR(x)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  converged=(ABS(dx).LT.tol).AND.(x.GT.a).AND.(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) STOP 'NewtonRoot1D not converged'
xout=x

END FUNCTION NewtonRoot1D


!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function FR(x(:)) in the interval [a(:),b(:)], using d/dx(:)F(x)=0 and the derivative 
!!
!===================================================================================================================================
FUNCTION NewtonRoot1D_FdF(tol,a,b,xin,F0,FRdFR) RESULT (xout)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: tol    !! abort tolerance
REAL(wp),INTENT(IN) :: a,b    !! search interval
REAL(wp),INTENT(IN) :: F0     !! function to find root is FR(x)-F0 
REAL(wp),INTENT(IN) :: xin    !! initial guess on input
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
PROCEDURE(i_f2x1)   :: FRdFR   !! function to find root f(x) & derivative d/dx f(x)
REAL(wp)            :: xout    !! output x for f(x)=0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL(wp)            :: x,dx
REAL(wp)            :: FRdFRx(2) !1: FR(x), 2: dFR(x)
LOGICAL             :: converged
!===================================================================================================================================

converged=.FALSE.
x=xin
maxiter=50
DO iter=1,maxiter
  FRdFRx=FRdFR(x)
  dx=-(FRdFRx(1)-F0)/FRdFRx(2)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  converged=(ABS(dx).LT.tol).AND.(x.GT.a).AND.(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) STOP 'NewtonRoot1D not converged'
xout=x

END FUNCTION NewtonRoot1D_FdF


!===================================================================================================================================
!> Newton's iterative algorithm to find the minimimum of function f(x,y) in the interval x(i)[a(i),b(i)],
!! using grad(f(x)=0 and the derivative 
!!
!===================================================================================================================================
FUNCTION NewtonMin2D(tol,a,b,x,FF,dFF,ddFF) RESULT (fmin)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)    :: tol        !! abort tolerance
REAL(wp),INTENT(IN)    :: a(2),b(2)  !! search interval (2D)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(INOUT) :: x(2) !! initial guess on input, result on output
PROCEDURE(i_f1x2)      :: FF   !! functional f(x,y) to minimize
PROCEDURE(i_f2x2)      :: dFF  !! d/dx f(x,y),d/dyf(x,y)
PROCEDURE(i_f22x2)     :: ddFF !! d^2/dx^2 f(x)
REAL(wp)               :: fmin !! on output =f(x,y)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL(wp)            :: dx(2)
REAL(wp)            :: det_Hess
REAL(wp)            :: gradF(2),Hess(2,2),HessInv(2,2)
LOGICAL             :: converged
!===================================================================================================================================
converged=.FALSE.
maxiter=50
DO iter=1,maxiter
  Hess=ddFF(x)
  det_Hess = Hess(1,1)*Hess(2,2)-Hess(1,2)*Hess(2,1)
  IF(det_Hess.LT.1.0E-12) STOP 'det Hessian=0 in NewtonMin'
  HessInv(1,1)= Hess(2,2)
  HessInv(1,2)=-Hess(1,2)
  HessInv(2,1)=-Hess(2,1)
  HessInv(2,2)= Hess(1,1)
  HessInv=HessInv/det_Hess
  gradF=dFF(x) 
  dx=-MATMUL(HessInv,gradF)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  converged=(SQRT(SUM(dx*dx)).LT.tol).AND.ALL(x.GT.a).AND.ALL(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) STOP 'NewtonMin not converged'
fmin=FF(x)

END FUNCTION NewtonMin2D

END MODULE MOD_Newton
