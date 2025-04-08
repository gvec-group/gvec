!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **Linear Algebra**
!!
!! Provides the linear algebra wrapper routines using LAPACK.
!!
!!- matrix inverse
!!- solve linear system
!!
!===================================================================================================================================
MODULE MODgvec_LinAlg

USE MODgvec_Globals, ONLY:wp,abort
IMPLICIT NONE

PUBLIC

CONTAINS

!===================================================================================================================================
!> Computes matrix inverse using LAPACK
!! Input matrix should be a square matrix
!!
!===================================================================================================================================
FUNCTION INV(A) RESULT(Ainv)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(wp),INTENT(IN)  :: A(:,:)                      !! input matrix
REAL(wp)             :: Ainv(SIZE(A,1),SIZE(A,2))   !! result: inverse of A
!-----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI
! LOCAL VARIABLES
REAL(wp):: work(SIZE(A,1))  ! work array for LAPACK
INTEGER :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER :: n,info
!===================================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
CALL DGETRF(n, n, Ainv, n, ipiv, info)

IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
              'Matrix is numerically singular!')
END IF

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
              'Matrix inversion failed!')
END IF
END FUNCTION INV


!===================================================================================================================================
!> Solve  linear system of dimension dims and multiple RHS
!!
!===================================================================================================================================
FUNCTION SOLVE(A,RHS) RESULT(X)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: A(:,:) !! matrix
REAL(wp),INTENT(IN) :: RHS(:) !! RHS, sorting: (dimA,nRHS), two dimensions can be used in input
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)           :: X(SIZE(RHS,1))    !! result: solution of A X=RHS
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRS
! LOCAL VARIABLES
REAL(wp)    :: Atmp(SIZE(A,1), SIZE(A,1))
INTEGER     :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER     :: nRHS,n,info
!===================================================================================================================================
Atmp=A
X = RHS
n = SIZE(A,1)
nRHS=SIZE(RHS,1)/SIZE(A,1)

CALL DGETRF(n, n, Atmp, n, ipiv, info)

IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
              'Matrix is numerically singular!')
END IF

CALL DGETRS('N',n, nRHS,Atmp, n, ipiv,X,n, info)
IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
              'Matrix solve does not work!')
END IF
END FUNCTION SOLVE

!===================================================================================================================================
!> Solve  linear system of dimension dims and multiple RHS
!!
!===================================================================================================================================
FUNCTION SOLVEMAT(A,RHS) RESULT(X)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: A(:,:) !! matrix
REAL(wp),INTENT(IN) :: RHS(:,:) !! RHS, sorting: (dimA,nRHS), two dimensions can be used in input
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)           :: X(SIZE(A,1),SIZE(RHS,2))    !! result: solution of A X=RHS
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRS
! LOCAL VARIABLES
REAL(wp)    :: Atmp(SIZE(A,1), SIZE(A,1))
INTEGER     :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER     :: nRHS,n,info
!===================================================================================================================================
Atmp=A
X = RHS
n = SIZE(A,1)
nRHS=SIZE(RHS,2)

CALL DGETRF(n, n, Atmp, n, ipiv, info)

IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
              'Matrix is numerically singular!')
END IF

CALL DGETRS('N',n, nRHS,Atmp, n, ipiv,X,n, info)
IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
             'Matrix solve does not work!')
END IF
END FUNCTION SOLVEMAT

!===================================================================================================================================
!> Return P L U matrices of the LU decomposition, cmoputed from LAPACK Routine (if P is not passed, L=P*L)
!!
!===================================================================================================================================
SUBROUTINE getLU(dimA,A,L,U,P)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: dimA
REAL(wp),INTENT(IN) :: A(1:dimA,1:dimA) !! matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: L(1:dimA,1:dimA)   !! L or P*L (if P omitted)
REAL(wp),INTENT(OUT) :: U(1:dimA,1:dimA)
REAL(wp),INTENT(OUT),OPTIONAL :: P(1:dimA,1:dimA)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! External procedures defined in LAPACK
EXTERNAL DGETRF
! LOCAL VARIABLES
REAL(wp)    :: tmp,Atmp(dimA,dimA)
INTEGER     :: ipiv(dimA)  ! pivot indices
INTEGER     :: i,j,info
!===================================================================================================================================
Atmp=A

CALL DGETRF(dimA, dimA, Atmp, dimA, ipiv, info)

IF(info.NE.0)THEN
   CALL abort(__STAMP__,&
              'Matrix is numerically singular!')
END IF
!upper diagonal part
U=0.
DO i=1,dimA
  DO j=i,dimA
    U(i,j)=Atmp(i,j)
  END DO
END DO
!lower diagonal part, known: 1 on the diagonal
L=0.
DO i=1,dimA
  L(i,i)=1.
  DO j=1,i-1
    L(i,j)=Atmp(i,j)
  END DO
END DO

IF(PRESENT(P))THEN
  P=0.
  DO i=1,dimA
    P(i,i)=1.
  END DO

  !pivoting of rows in L is pivoting of columns in P
  DO i=1,dimA
    DO j=1,dimA
      tmp=P(j,i)
      P(j,i)=P(j,ipiv(i))
      P(j,ipiv(i))=tmp
    END DO
  END DO

  !CHECK
  IF(SUM(ABS(MATMUL(P,MATMUL(L,U))-A)).GT.1e-12*dimA*dimA) THEN
    CALL abort(__STAMP__,&
               'A=P*L*U decomposition not correct')
  END IF
ELSE
  !pivoting of rows in L, backwards
  DO i=dimA,1,-1
    DO j=1,dimA
      tmp=L(i,j)
      L(i,j)=L(ipiv(i),j)
      L(ipiv(i),j)=tmp
    END DO
  END DO

  !CHECK
  IF(SUM(ABS(MATMUL(L,U)-A)).GT.1e-12*dimA*dimA) THEN
    CALL abort(__STAMP__,&
               'A=L*U decomposition not correct')
  END IF
END IF

END SUBROUTINE getLU


END MODULE MODgvec_LinAlg
