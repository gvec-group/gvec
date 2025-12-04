!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **NEWTON**
!!
!! Some simple Newton solvers
!!
!===================================================================================================================================
MODULE MODgvec_Newton
! MODULES
USE MODgvec_Globals, ONLY:wp,UNIT_stdout,abort
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

INTERFACE NewtonRoot2D
  MODULE PROCEDURE NewtonRoot2D
END INTERFACE

TYPE, ABSTRACT :: c_newton_Min1D
  CONTAINS
  PROCEDURE(i_newton_Min1D), DEFERRED :: FR
  PROCEDURE(i_newton_Min1D), DEFERRED :: dFR
  PROCEDURE(i_newton_Min1D), DEFERRED :: ddFR
END TYPE c_newton_Min1D

TYPE, ABSTRACT :: c_newton_Root1D
  CONTAINS
  PROCEDURE(i_newton_Root1D), DEFERRED :: FR
  PROCEDURE(i_newton_Root1D), DEFERRED :: dFR
END TYPE c_newton_Root1D

!===================================================================================================================================
!> This functional wraps a functional passed to Min1D for use with Root1D (returns the derivative: FR -> dFR, dFR -> ddFR)
!!
!===================================================================================================================================
TYPE, EXTENDS(c_newton_Root1D) :: t_newton_Root1D_wrap_Min1D
  CLASS(c_newton_Min1D), ALLOCATABLE :: parent
  CONTAINS
  PROCEDURE :: FR => newton_Root1D_wrap_Min1D_FR
  PROCEDURE :: dFR => newton_Root1D_wrap_Min1D_dFR
END TYPE t_newton_Root1D_wrap_Min1D

INTERFACE t_newton_Root1D_wrap_Min1D
    MODULE PROCEDURE newton_Root1D_wrap_Min1D_new
END INTERFACE t_newton_Root1D_wrap_Min1D

TYPE, ABSTRACT :: c_newton_Root1D_FdF
  CONTAINS
  PROCEDURE(i_newton_Root1D_FdF), DEFERRED :: FRdFR
END TYPE c_newton_Root1D_FdF

TYPE, ABSTRACT :: c_newton_Min2D
  CONTAINS
  PROCEDURE(i_newton_Min2D_FR), DEFERRED :: FR
  PROCEDURE(i_newton_Min2D_dFR), DEFERRED :: dFR
  PROCEDURE(i_newton_Min2D_ddFR), DEFERRED :: ddFR
END TYPE c_newton_Min2D

TYPE, ABSTRACT :: c_newton_Root2D
  CONTAINS
  PROCEDURE(i_newton_Root2D_FR), DEFERRED :: FR
  PROCEDURE(i_newton_Root2D_dFR), DEFERRED :: dFR
END TYPE c_newton_Root2D

ABSTRACT INTERFACE
  FUNCTION i_newton_Min1D(sf, x) RESULT (y1x1)
    IMPORT wp, c_newton_Min1D
    IMPLICIT NONE
    CLASS(c_newton_Min1D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: y1x1
  END FUNCTION i_newton_Min1D

  FUNCTION i_newton_Root1D(sf, x) RESULT (y1x1)
    IMPORT wp, c_newton_Root1D
    IMPLICIT NONE
    CLASS(c_newton_Root1D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: y1x1
  END FUNCTION i_newton_Root1D

  FUNCTION i_newton_Root1D_FdF(sf, x) RESULT (y2x1)
    IMPORT wp, c_newton_Root1D_FdF
    IMPLICIT NONE
    CLASS(c_newton_Root1D_FdF), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: y2x1(2)
  END FUNCTION i_newton_Root1D_FdF

  FUNCTION i_newton_Min2D_FR(sf, x) RESULT (y1x2)
    IMPORT wp, c_newton_Min2D
    IMPLICIT NONE
    CLASS(c_newton_Min2D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x(2)
    REAL(wp) :: y1x2
  END FUNCTION i_newton_Min2D_FR

  FUNCTION i_newton_Min2D_dFR(sf, x) RESULT (y2x2)
    IMPORT wp, c_newton_Min2D
    IMPLICIT NONE
    CLASS(c_newton_Min2D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x(2)
    REAL(wp) :: y2x2(2)
  END FUNCTION i_newton_Min2D_dFR

  FUNCTION i_newton_Min2D_ddFR(sf, x) RESULT (y22x2)
    IMPORT wp, c_newton_Min2D
    IMPLICIT NONE
    CLASS(c_newton_Min2D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x(2)
    REAL(wp) :: y22x2(2,2)
  END FUNCTION i_newton_Min2D_ddFR

  FUNCTION i_newton_Root2D_FR(sf, x) RESULT (y2x2)
    IMPORT wp, c_newton_Root2D
    IMPLICIT NONE
    CLASS(c_newton_Root2D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x(2)
    REAL(wp) :: y2x2(2)
  END FUNCTION i_newton_Root2D_FR

  FUNCTION i_newton_Root2D_dFR(sf, x) RESULT (y22x2)
    IMPORT wp, c_newton_Root2D
    IMPLICIT NONE
    CLASS(c_newton_Root2D), INTENT(IN) :: sf
    REAL(wp), INTENT(IN) :: x(2)
    REAL(wp) :: y22x2(2, 2)
  END FUNCTION i_newton_Root2D_dFR
END INTERFACE

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Newton's iterative algorithm to find the minimimum of function f(x) in the interval [a,b], using df(x)=0 and the derivative
!!
!===================================================================================================================================
FUNCTION NewtonMin1D(tol,a,b,maxstep,x,fobj) RESULT (fmin)
  ! MODULES
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  REAL(wp),INTENT(IN)               :: tol  !! abort tolerance
  REAL(wp),INTENT(IN)               :: a,b  !! search interval
  REAL(wp),INTENT(IN)               :: maxstep  !! max|dx| allowed
  CLASS(c_newton_Min1D),INTENT(IN)  :: fobj !! functional to minimize with f(x), d/dx f(x), d^2/dx^2 f(x)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),INTENT(INOUT) :: x    !! initial guess on input, result on output
  REAL(wp)               :: fmin !! on output =f(x)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  REAL(wp)                                      :: x0
  TYPE(t_newton_Root1D_wrap_Min1D), ALLOCATABLE :: fobj_wrap
  !=================================================================================================================================
  fobj_wrap = t_newton_Root1D_wrap_Min1D(fobj)
  x0=x
  x=NewtonRoot1D(tol,a,b,maxstep,x0,0.0_wp,fobj_wrap)
  fmin=fobj%FR(x)
  DEALLOCATE(fobj_wrap)
END FUNCTION NewtonMin1D


!===================================================================================================================================
!> constructor for the Min1D type wrapped for Root1D
!!
!===================================================================================================================================
FUNCTION newton_Root1D_wrap_Min1D_new(parent) RESULT(sf)
  ! MODULES
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(c_newton_Min1D), INTENT(IN) :: parent
  ! OUTPUT VARIABLES
  TYPE(t_newton_Root1D_wrap_Min1D) :: sf
  !=================================================================================================================================
  sf%parent = parent
END FUNCTION


!===================================================================================================================================
!> f(x) function of the Min1D type wrapped for Root1D
!!
!===================================================================================================================================
FUNCTION newton_Root1D_wrap_Min1D_FR(sf, x) RESULT(y1x1)
  ! MODULES
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_newton_Root1D_wrap_Min1D), INTENT(IN) :: sf
  REAL(wp), INTENT(IN) :: x
  ! OUTPUT VARIABLES
  REAL(wp) :: y1x1
  !=================================================================================================================================
  y1x1 = sf%parent%dFR(x)
END FUNCTION


!===================================================================================================================================
!> d/dx f(x) function of the Min1D type wrapped for Root1D
!!
!===================================================================================================================================
FUNCTION newton_Root1D_wrap_Min1D_dFR(sf, x) RESULT(y1x1)
  ! MODULES
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_newton_Root1D_wrap_Min1D), INTENT(IN) :: sf
  REAL(wp), INTENT(IN) :: x
  ! OUTPUT VARIABLES
  REAL(wp) :: y1x1
  !=================================================================================================================================
  y1x1 = sf%parent%ddFR(x)
END FUNCTION


!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function FR(x(:)) in the interval [a(:),b(:)], using d/dx(:)F(x)=0 and the derivative
!!
!===================================================================================================================================
FUNCTION NewtonRoot1D(tol,a,b,maxstep,xin,F0,fobj) RESULT (xout)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: tol    !! abort tolerance
REAL(wp),INTENT(IN) :: a,b    !! search interval
REAL(wp),INTENT(IN) :: maxstep !! max|dx| allowed
REAL(wp),INTENT(IN) :: xin    !! initial guess
REAL(wp),INTENT(IN) :: F0     !! function to find root is FR(x)-F0
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(c_newton_Root1D) :: fobj
REAL(wp)            :: xout    !! on output =f(x)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL(wp)            :: x,dx
LOGICAL             :: converged
LOGICAL             :: converged2
!===================================================================================================================================

converged=.FALSE.
x=xin
maxiter=20
DO iter=1,maxiter
  dx=-(fobj%FR(x)-F0)/fobj%dFR(x)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  IF(ABS(dx).GT.maxstep) dx=dx/ABS(dx)*maxstep
  converged=(ABS(dx).LT.tol).AND.(x.GT.a).AND.(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) THEN
  !repeat with maxstep /10 and a little change in the initial condition
  x=MIN(b,MAX(a,xin+0.01_wp*(b-a)))
  maxiter=200
  DO iter=1,maxiter
    dx=-(fobj%FR(x)-F0)/fobj%dFR(x)
    dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
    IF(ABS(dx).GT.maxstep) dx=dx/ABS(dx)*0.1_wp*maxstep
    x = x+dx
    converged2=(ABS(dx).LT.tol).AND.(x.GT.a).AND.(x.LT.b)
    IF(converged2) EXIT
  END DO !iter
  IF(converged2) THEN
    xout=x
    RETURN
  END IF
  WRITE(UNIT_stdout,*)'Newton abs(dx)<tol',ABS(dx),tol,ABS(dx).LT.tol
  WRITE(UNIT_stdout,*)'Newton x>a',x,a,(x.GT.a)
  WRITE(UNIT_stdout,*)'Newton x<b',x,b,(x.LT.b)
  WRITE(UNIT_stdout,*)'after iter',iter-1
  CALL abort(__STAMP__, &
             'NewtonRoot1D not converged')
END IF
xout=x

END FUNCTION NewtonRoot1D


!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function FR(x(:)) in the interval [a(:),b(:)], using d/dx(:)F(x)=0 and the derivative
!!
!===================================================================================================================================
FUNCTION NewtonRoot1D_FdF(tol,a,b,maxstep,xin,F0,fobj) RESULT (xout)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: tol     !! abort tolerance
REAL(wp),INTENT(IN) :: a,b     !! search interval
REAL(wp),INTENT(IN) :: maxstep !! max|dx| allowed
REAL(wp),INTENT(IN) :: xin     !! initial guess on input
REAL(wp),INTENT(IN) :: F0      !! function to find root is FR(x)-F0
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(c_newton_Root1D_FdF) :: fobj !! function to find root f(x) & derivative d/dx f(x) as FRdFR method
REAL(wp)            :: xout    !! output x for f(x)=0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL(wp)            :: x,dx
REAL(wp)            :: FRdFRx(2) !1: FR(x), 2: dFR(x)
LOGICAL             :: converged
LOGICAL             :: converged2
!===================================================================================================================================
converged=.FALSE.
x=xin
maxiter=20
DO iter=1,maxiter
  FRdFRx=fobj%FRdFR(x)
  dx=-(FRdFRx(1)-F0)/FRdFRx(2)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  IF(ABS(dx).GT.maxstep) dx=dx/ABS(dx)*maxstep
  x = x+dx
  converged=(ABS(dx).LT.tol).AND.(x.GE.a).AND.(x.LE.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) THEN
  !repeat with maxstep /10 and a little change in the initial condition
  converged2=.FALSE.
  x=MIN(b,MAX(a,xin+0.01_wp*(b-a)))
  maxiter=200
  DO iter=1,maxiter
    FRdFRx=fobj%FRdFR(x)
    dx=-(FRdFRx(1)-F0)/FRdFRx(2)
    dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
    IF(ABS(dx).GT.maxstep) dx=dx/ABS(dx)*0.1_wp*maxstep
    x = x+dx
    converged2=(ABS(dx).LT.tol).AND.(x.GE.a).AND.(x.LE.b)
    IF(converged2) EXIT
  END DO !iter
  IF(converged2) THEN
    xout=x
    RETURN
  END IF
  WRITE(UNIT_stdout,*)'Newton abs(dx)<tol',ABS(dx),tol,ABS(dx).LT.tol
  WRITE(UNIT_stdout,*)'Newton x>a',x,a,(x.GT.a)
  WRITE(UNIT_stdout,*)'Newton x<b',x,b,(x.LT.b)
  WRITE(UNIT_stdout,*)'after iter',iter-1
  CALL abort(__STAMP__,&
             'NewtonRoot1D_FdF not converged')
END IF
xout=x

END FUNCTION NewtonRoot1D_FdF


!===================================================================================================================================
!> Newton's iterative algorithm to find the minimimum of function f(x,y) in the interval x(i)[a(i),b(i)],
!! using grad(f(x)=0 and the derivative
!!
!===================================================================================================================================
FUNCTION NewtonMin2D(tol,a,b,x,fobj) RESULT (fmin)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)    :: tol        !! abort tolerance
REAL(wp),INTENT(IN)    :: a(2),b(2)  !! search interval (2D)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(INOUT) :: x(2) !! initial guess on input, result on output
CLASS(c_newton_Min2D)  :: fobj !! functional f(x,y) to minimize with f, f', f''
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
  Hess=fobj%ddFR(x)
  det_Hess = Hess(1,1)*Hess(2,2)-Hess(1,2)*Hess(2,1)
  IF(det_Hess.LT.1.0E-12) CALL abort(__STAMP__,&
                                     'det Hessian=0 in NewtonMin')
  HessInv(1,1)= Hess(2,2)
  HessInv(1,2)=-Hess(1,2)
  HessInv(2,1)=-Hess(2,1)
  HessInv(2,2)= Hess(1,1)
  HessInv=HessInv/det_Hess
  gradF=fobj%dFR(x)
  dx=-MATMUL(HessInv,gradF)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  converged=(SQRT(SUM(dx*dx)).LT.tol).AND.ALL(x.GT.a).AND.ALL(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) CALL abort(__STAMP__,&
                              'NewtonMin2D not converged')
fmin=fobj%FR(x)

END FUNCTION NewtonMin2D

!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function [f1(x1,x2),f2(x1,x2)]=[0,0] in the interval a(i)<=x(i)<=b(i),
!! using the Jacobian  dfi/dxj, i=1,2, j=1,2, such that fi(x1,x2)=fi(x1_0,x2_0)+  [dfi/dx1,dfi/dx2].[dx1,dx2]
!! in each step, we find dx1,dx2 st -[[dfi/dxj]] dxj =fi(x1_0,x2_0)
!!
!===================================================================================================================================
FUNCTION NewtonRoot2D(tol,a,b,maxstep,xin,fobj) RESULT (xout)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)    :: tol        !! abort tolerance
REAL(wp),INTENT(IN)    :: a(2),b(2)  !! search interval (2D)
REAL(wp),INTENT(IN) :: maxstep(2) !! max|dx| allowed
REAL(wp),INTENT(IN)    :: xin(2) !! initial guess
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(c_newton_Root2D) :: fobj !! function to find root f1(x1,x2)=0,f2(x1,x2)=0 and derivatives d fi(x1,x2)/dxj
REAL(wp)               :: xout(2) !! x1,x2 that have f1(x1,x2)=0 and f2(x1,x2)=0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL(wp)            :: x(2),dx(2)
REAL(wp)            :: det_Jac
REAL(wp)            :: F(2),Jac(2,2),JacInv(2,2)
LOGICAL             :: converged
!===================================================================================================================================
converged=.FALSE.
x=xin
maxiter=50
DO iter=1,maxiter
  Jac=fobj%dFR(x)
  det_Jac = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)
  IF(det_Jac.LT.1.0E-12) CALL abort(__STAMP__,&
                                    'det Jacobian<=0 in NewtonRoot2d')
  JacInv(1,1)= Jac(2,2)
  JacInv(1,2)=-Jac(1,2)
  JacInv(2,1)=-Jac(2,1)
  JacInv(2,2)= Jac(1,1)
  JacInv=JacInv/det_Jac
  F=fobj%FR(x)
  dx=-MATMUL(JacInv,F)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  IF(ABS(dx(1)).GT.maxstep(1)) dx(1)=dx(1)/ABS(dx(1))*maxstep(1)
  IF(ABS(dx(2)).GT.maxstep(2)) dx(2)=dx(2)/ABS(dx(2))*maxstep(2)

  x = x+dx
  converged=(SQRT(SUM(dx*dx)).LT.tol).AND.ALL(x.GT.a).AND.ALL(x.LT.b)
  IF(converged) EXIT
END DO !iter
xout=x
IF(.NOT.converged) THEN
  WRITE(UNIT_stdout,*)'Newton abs(dx)<tol',ABS(dx),tol,'F(x)',fobj%FR(xout)
  CALL abort(__STAMP__,&
             'NewtonRoot2D not converged')
END IF

END FUNCTION NewtonRoot2D

END MODULE MODgvec_Newton
