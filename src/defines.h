!===================================================================================================================================
! Copyright (c) 2017-2018 Florian Hindenlang <hindenlang@gmail.com>
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

! Abbrevations
#ifndef __FILENAME__ 
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#if MPI
#  define SWRITE IF(MPIRoot) WRITE
#else
#  define SWRITE WRITE
#endif
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)

!boundary condition for zero,odd and even m modes

!boundary types
#define NBC_TYPES    6

#define BC_TYPE_OPEN      1 
#define BC_TYPE_NEUMANN   2
#define BC_TYPE_DIRICHLET 3
#define BC_TYPE_SYMM      4
#define BC_TYPE_SYMMZERO  5
#define BC_TYPE_ANTISYMM  6

! index in BC arrays
#define BC_AXIS 1
#define BC_EDGE 2

!grid types
#define GRID_TYPE_UNIFORM 0 
#define GRID_TYPE_SQRT_S  1 
#define GRID_TYPE_S2      2 
#define GRID_TYPE_BUMP    3 
#define GRID_TYPE_BUMP_EDGE 4 

!fbase definitions
#define _SIN_    1
#define _COS_    2
#define _SINCOS_ 3

#define MN_ZERO  0
#define M_ODD    1
#define M_EVEN   2
#define M_ZERO   3
#define M_ODD_FIRST 4

#define DERIV_S    1
#define DERIV_THET 2
#define DERIV_ZETA 3


!!!!matvec with matmul
!!!#define __MATVEC_N(y,Mat,Vec)           y=MATMUL(Mat,Vec) 
!!!#define __MATVEC_T(y,Mat,Vec)           y=MATMUL(Vec,Mat)
!!!#define __PMATVEC_N(fy,y,Mat,Vec)       y=fy*y+MATMUL(Mat,Vec)
!!!#define __PMATVEC_T(fy,y,Mat,Vec)       y=fy*y+MATMUL(Vec,Mat)
!!!#define __AMATVEC_N(y,fMat,Mat,Vec)     y=fMat*MATMUL(Mat,Vec)
!!!#define __AMATVEC_T(y,fMat,Mat,Vec)     y=fMat*MATMUL(Vec,Mat)
!!!#define __PAMATVEC_N(fy,y,fMat,Mat,Vec) y=fy*y+fMat*MATMUL(Mat,Vec)
!!!#define __PAMATVEC_T(fy,y,fMat,Mat,Vec) y=fy*y+fMat*MATMUL(Vec,Mat)

! matvec with blas
#define __MATVEC_N(y,Mat,Vec)          __PAMATVEC('N',0.0_wp,y,1.0_wp,Mat,Vec) 
#define __MATVEC_T(y,Mat,Vec)          __PAMATVEC('T',0.0_wp,y,1.0_wp,Mat,Vec) 

#define __PMATVEC_N(fy,y,Mat,Vec)      __PAMATVEC('N',    fy,y,1.0_wp,Mat,Vec) 
#define __PMATVEC_T(fy,y,Mat,Vec)      __PAMATVEC('T',    fy,y,1.0_wp,Mat,Vec) 

#define __AMATVEC_N(y,fMat,Mat,Vec)    __PAMATVEC('N',0.0_wp,y,  fMat,Mat,Vec) 
#define __AMATVEC_T(y,fMat,Mat,Vec)    __PAMATVEC('T',0.0_wp,y,  fMat,Mat,Vec) 
!!!!
#define __PAMATVEC(NT,fy,y,fMat,Mat,Vec) CALL DGEMV(NT,SIZE(Mat,1),SIZE(Mat,2),fMat  ,Mat,SIZE(Mat,1),Vec,1,fy    ,y,1)


!!!!matmat with matmul
!!!#define __MATMAT_NN(Y,A,B)         Y=MATMUL(A,B)
!!!#define __MATMAT_TN(Y,A,B)         Y=MATMUL(TRANSPOSE(A),B)
!!!#define __MATMAT_NT(Y,A,B)         Y=MATMUL(A,TRANSPOSE(B))
!!!#define __MATMAT_TT(Y,A,B)         Y=TRANSPOSE(MATMUL(B,A))

!!!#define __PMATMAT_NN(fy,Y,A,B)     Y=fy*Y+MATMUL(A,B)
!!!#define __PMATMAT_TN(fy,Y,A,B)     Y=fy*Y+MATMUL(TRANSPOSE(A),B)
!!!#define __PMATMAT_NT(fy,Y,A,B)     Y=fy*Y+MATMUL(A,TRANSPOSE(B))
!!!#define __PMATMAT_TT(fy,Y,A,B)     Y=fy*Y+TRANSPOSE(MATMUL(B,A))

!!!#define __AMATMAT_NN(Y,fa,A,B)     Y=fa*MATMUL(A,B)
!!!#define __AMATMAT_TN(Y,fa,A,B)     Y=fa*MATMUL(TRANSPOSE(A),B)
!!!#define __AMATMAT_NT(Y,fa,A,B)     Y=fa*MATMUL(A,TRANSPOSE(B))
!!!#define __AMATMAT_TT(Y,fa,A,B)     Y=fa*TRANSPOSE(MATMUL(B,A))

!!!#define __PAMATMAT_NN(fy,Y,fa,A,B) Y=fy*Y+fa*MATMUL(A,B)
!!!#define __PAMATMAT_TN(fy,Y,fa,A,B) Y=fy*Y+fa*MATMUL(TRANSPOSE(A),B)
!!!#define __PAMATMAT_NT(fy,Y,fa,A,B) Y=fy*Y+fa*MATMUL(A,TRANSPOSE(B))
!!!#define __PAMATMAT_TT(fy,Y,fa,A,B) Y=fy*Y+fa*TRANSPOSE(MATMUL(B,A))


! matmat with blas
#define __MATMAT_NN(Y,A,B)     __PAMATMAT_NN(0.0_wp,Y,1.0_wp,A,B)
#define __MATMAT_TN(Y,A,B)     __PAMATMAT_TN(0.0_wp,Y,1.0_wp,A,B)
#define __MATMAT_NT(Y,A,B)     __PAMATMAT_NT(0.0_wp,Y,1.0_wp,A,B)
#define __MATMAT_TT(Y,A,B)     __PAMATMAT_TT(0.0_wp,Y,1.0_wp,A,B)

#define __PMATMAT_NN(fy,Y,A,B) __PAMATMAT_NN(    fy,Y,1.0_wp,A,B)
#define __PMATMAT_TN(fy,Y,A,B) __PAMATMAT_TN(    fy,Y,1.0_wp,A,B)
#define __PMATMAT_NT(fy,Y,A,B) __PAMATMAT_NT(    fy,Y,1.0_wp,A,B)
#define __PMATMAT_TT(fy,Y,A,B) __PAMATMAT_TT(    fy,Y,1.0_wp,A,B)

#define __AMATMAT_NN(Y,fa,A,B) __PAMATMAT_NN(0.0_wp,Y,    fa,A,B)
#define __AMATMAT_TN(Y,fa,A,B) __PAMATMAT_TN(0.0_wp,Y,    fa,A,B)
#define __AMATMAT_NT(Y,fa,A,B) __PAMATMAT_NT(0.0_wp,Y,    fa,A,B)
#define __AMATMAT_TT(Y,fa,A,B) __PAMATMAT_TT(0.0_wp,Y,    fa,A,B)

!!! GEMM does in general Y = fa A^?*B^? + fy Y
!!! with structure: (m x n) = (m x k) (k x n)  
!!! Y=A  *B   : DGEMM('N','N',m,n,k,fa,Amat ,m, Bmat,k, fy,Y,m)
!!! Y=A^T*B   : DGEMM('T','N',m,n,k,fa,Amat ,k, Bmat,k, fy,Y,m)
!!! Y=A  *B^T : DGEMM('N','T',m,n,k,fa,Amat ,m, Bmat,n, fy,Y,m)
!!! Y=A^T*B^T : DGEMM('T','T',m,n,k,fa,Amat ,k, Bmat,n, fy,Y,m)

#define __PAMATMAT_NN(fy,Y,fa,A,B) CALL DGEMM('N','N',SIZE(A,1),SIZE(B,2),SIZE(B,1),fa,A,SIZE(A,1),B,SIZE(B,1),fy,Y,SIZE(A,1))
#define __PAMATMAT_TN(fy,Y,fa,A,B) CALL DGEMM('T','N',SIZE(A,2),SIZE(B,2),SIZE(B,1),fa,A,SIZE(B,1),B,SIZE(B,1),fy,Y,SIZE(A,1))
#define __PAMATMAT_NT(fy,Y,fa,A,B) CALL DGEMM('N','T',SIZE(A,1),SIZE(B,1),SIZE(B,2),fa,A,SIZE(A,1),B,SIZE(B,1),fy,Y,SIZE(A,1))
#define __PAMATMAT_TT(fy,Y,fa,A,B) CALL DGEMM('T','T',SIZE(A,2),SIZE(B,1),SIZE(B,2),fa,A,SIZE(B,2),B,SIZE(B,1),fy,Y,SIZE(A,1))

