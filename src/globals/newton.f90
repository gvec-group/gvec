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
FUNCTION NewtonMin1D(tol,a,b,maxstep,x,FF,dFF,ddFF) RESULT (fmin)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)    :: tol  !! abort tolerance
REAL(wp),INTENT(IN)    :: a,b  !! search interval
REAL(wp),INTENT(IN)    :: maxstep  !! max|dx| allowed
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
x=NewtonRoot1D(tol,a,b,maxstep,x0,0.0_wp,dFF,ddFF)
fmin=FF(x)

END FUNCTION NewtonMin1D


!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function FR(x(:)) in the interval [a(:),b(:)], using d/dx(:)F(x)=0 and the derivative
!!
!===================================================================================================================================
FUNCTION NewtonRoot1D(tol,a,b,maxstep,xin,F0,FR,dFR) RESULT (xout)
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
PROCEDURE(i_f1x1)   :: FR      !! function to find root
PROCEDURE(i_f1x1)   :: dFR     !! multidimensional derivative d/dx f(x), size dim1
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
  dx=-(FR(x)-F0)/dFR(x)
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
    dx=-(FR(x)-F0)/dFR(x)
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
FUNCTION NewtonRoot1D_FdF(tol,a,b,maxstep,xin,F0,FRdFR) RESULT (xout)
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
PROCEDURE(i_f2x1)   :: FRdFR   !! function to find root f(x) & derivative d/dx f(x)
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
  FRdFRx=FRdFR(x)
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
    FRdFRx=FRdFR(x)
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
  IF(det_Hess.LT.1.0E-12) CALL abort(__STAMP__,&
                                     'det Hessian=0 in NewtonMin')
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
IF(.NOT.converged) CALL abort(__STAMP__,&
                              'NewtonMin2D not converged')
fmin=FF(x)

END FUNCTION NewtonMin2D

!===================================================================================================================================
!> Newton's iterative algorithm to find the root of function [f1(x1,x2),f2(x1,x2)]=[0,0] in the interval a(i)<=x(i)<=b(i),
!! using the Jacobian  dfi/dxj, i=1,2, j=1,2, such that fi(x1,x2)=fi(x1_0,x2_0)+  [dfi/dx1,dfi/dx2].[dx1,dx2]
!! in each step, we find dx1,dx2 st -[[dfi/dxj]] dxj =fi(x1_0,x2_0)
!!
!===================================================================================================================================
FUNCTION NewtonRoot2D(tol,a,b,maxstep,xin,FF,dFF) RESULT (xout)
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
PROCEDURE(i_f2x2)      :: FF  !! f1(x1,x2),f2(x1,x2) to be zero
PROCEDURE(i_f22x2)     :: dFF !! d fi(x1,x2) /dxj
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
  Jac=dFF(x)
  det_Jac = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)
  IF(det_Jac.LT.1.0E-12) CALL abort(__STAMP__,&
                                    'det Jacobian<=0 in NewtonRoot2d')
  JacInv(1,1)= Jac(2,2)
  JacInv(1,2)=-Jac(1,2)
  JacInv(2,1)=-Jac(2,1)
  JacInv(2,2)= Jac(1,1)
  JacInv=JacInv/det_Jac
  F=FF(x)
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
  WRITE(UNIT_stdout,*)'Newton abs(dx)<tol',ABS(dx),tol,'F(x)',FF(xout)
  CALL abort(__STAMP__,&
             'NewtonRoot2D not converged')
END IF

END FUNCTION NewtonRoot2D

END MODULE MODgvec_Newton
