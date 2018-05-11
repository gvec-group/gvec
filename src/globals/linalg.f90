!===================================================================================================================================
! Copyright (c) 2017 - 2018 Florian Hindenlang <hindenlang@gmail.com>
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

USE MODgvec_Globals, ONLY:wp
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
   STOP 'Matrix is numerically singular!'
END IF

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

IF(info.NE.0)THEN
   STOP 'Matrix inversion failed!'
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
   STOP 'Matrix is numerically singular!'
END IF

CALL DGETRS('N',n, nRHS,Atmp, n, ipiv,X,n, info)
IF(info.NE.0)THEN
   STOP 'Matrix solve does not work!'
END IF
END FUNCTION SOLVE


END MODULE MODgvec_LinAlg
