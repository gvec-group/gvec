!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
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
#  ifdef NO_MPI_F08
#    define USE_MPI USE mpi
#    define MPI_comm_TYPE INTEGER
#    define MPI_datatype_TYPE INTEGER
#    define MPI_request_TYPE INTEGER
#    define MPI_status_TYPE INTEGER, DIMENSION(MPI_STATUS_SIZE)
#    define MPI_op_TYPE INTEGER
#  else
#    define USE_MPI USE mpi_f08
#    define MPI_comm_TYPE TYPE(MPI_COMM)
#    define MPI_datatype_TYPE TYPE(MPI_DATATYPE)
#    define MPI_request_TYPE TYPE(MPI_REQUEST)
#    define MPI_status_TYPE TYPE(MPI_STATUS)
#    define MPI_op_TYPE TYPE(MPI_OP)
#  endif
#else
#  define USE_MPI
#  define SWRITE WRITE
#  define MPI_comm_TYPE INTEGER
#  define MPI_datatype_TYPE INTEGER
#  define MPI_request_TYPE INTEGER
#endif

#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)

#ifdef PP_FTIMINGS
#  define __PERFINIT   call perfinit()
#  define __PERFON(a)  call perfon(a)
#  define __PERFOFF(a) call perfoff(a)
#  define __PERFOUT(a) call perfout(a)
#else
#  define __PERFINIT
#  define __PERFON(a)
#  define __PERFOFF(a)
#  define __PERFOUT(a)
#endif

! if cmake option GVEC_FIX_HMAP is not off, sets PP_WHICH_HMAP
#if defined(PP_WHICH_HMAP)
#  if PP_WHICH_HMAP == 1
#    define  PP_MOD_HMAP MODgvec_hmap_RZ
#    define  PP_T_HMAP t_hmap_RZ
#    define  PP_T_HMAP_AUXVAR t_hmap_RZ_auxvar
#  elif PP_WHICH_HMAP == 3
#    define  PP_MOD_HMAP MODgvec_hmap_cyl
#    define  PP_T_HMAP t_hmap_cyl
#    define  PP_T_HMAP_AUXVAR t_hmap_cyl_auxvar
#  elif PP_WHICH_HMAP == 10
#    define  PP_MOD_HMAP MODgvec_hmap_knot
#    define  PP_T_HMAP t_hmap_knot
#    define  PP_T_HMAP_AUXVAR t_hmap_knot_auxvar
#  elif PP_WHICH_HMAP == 20
#    define  PP_MOD_HMAP MODgvec_hmap_frenet
#    define  PP_T_HMAP t_hmap_frenet
#    define  PP_T_HMAP_AUXVAR t_hmap_frenet_auxvar
#    define  PP_WHICH_HMAP 20
#  elif PP_WHICH_HMAP == 21
#    define  PP_MOD_HMAP MODgvec_hmap_axisNB
#    define  PP_T_HMAP t_hmap_axisNB
#    define  PP_T_HMAP_AUXVAR t_hmap_axisNB_auxvar
#  else
#    define  PP_MOD_HMAP MODgvec_hmap_?
#    define  PP_T_HMAP t_hmap_?
#    define  PP_T_HMAP_AUXVAR t_hmap_?_auxvar
#    define  PP_WHICH_HMAP ?
#  endif
#else
#  define  PP_T_HMAP c_hmap
#  define  PP_T_HMAP_AUXVAR c_hmap_auxvar
#endif /*PP_WHICH_HMAP defined*/


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
#define GRID_TYPE_CUSTOM -1
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
#define DERIV_S_S  2
#define DERIV_THET 2
#define DERIV_ZETA 3
#define DERIV_THET_THET 22
#define DERIV_THET_ZETA 23
#define DERIV_ZETA_ZETA 33


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
#define __MATVEC_N(y,Mat,Vec)           __GENERICMATVEC('N',0.0_wp,y,1.0_wp,Mat,Vec)
#define __MATVEC_T(y,Mat,Vec)           __GENERICMATVEC('T',0.0_wp,y,1.0_wp,Mat,Vec)

#define __PMATVEC_N(fy,y,Mat,Vec)       __GENERICMATVEC('N',fy,y,1.0_wp,Mat,Vec)
#define __PMATVEC_T(fy,y,Mat,Vec)       __GENERICMATVEC('T',fy,y,1.0_wp,Mat,Vec)

#define __AMATVEC_N(y,fMat,Mat,Vec)     __GENERICMATVEC('N',0.0_wp,y,fMat,Mat,Vec)
#define __AMATVEC_T(y,fMat,Mat,Vec)     __GENERICMATVEC('T',0.0_wp,y,fMat,Mat,Vec)

#define __PAMATVEC_N(fy,y,fMat,Mat,Vec) __GENERICMATVEC('N',fy,y,fMat,Mat,Vec)
#define __PAMATVEC_T(fy,y,fMat,Mat,Vec) __GENERICMATVEC('T',fy,y,fMat,Mat,Vec)

!!!!#define __GENERICMATVEC(NT,fy,y,fMat,Mat,Vec) CALL DGEMV(NT,SIZE(Mat,1),SIZE(Mat,2),fMat,Mat,SIZE(Mat,1),Vec,1,fy,y,1)

#define __GENERICMATVEC(NT,fy,y,fMat,Mat,Vec) __PADGEMV(NT,fy,y,fMat,SIZE(Mat,1),SIZE(Mat,2),Mat,Vec)

#define __PADGEMV(NT,fy,y,fMat,sz1,sz2,Mat,Vec) CALL DGEMV(NT,sz1,sz2,fMat,Mat,sz1,Vec,1,fy,y,1)

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
#define __MATMAT_NN(Y,A,B)     __GENERICMATMAT_NN(0.0_wp,Y,1.0_wp,A,B)
#define __MATMAT_TN(Y,A,B)     __GENERICMATMAT_TN(0.0_wp,Y,1.0_wp,A,B)
#define __MATMAT_NT(Y,A,B)     __GENERICMATMAT_NT(0.0_wp,Y,1.0_wp,A,B)
#define __MATMAT_TT(Y,A,B)     __GENERICMATMAT_TT(0.0_wp,Y,1.0_wp,A,B)

#define __PMATMAT_NN(fy,Y,A,B) __GENERICMATMAT_NN(fy,Y,1.0_wp,A,B)
#define __PMATMAT_TN(fy,Y,A,B) __GENERICMATMAT_TN(fy,Y,1.0_wp,A,B)
#define __PMATMAT_NT(fy,Y,A,B) __GENERICMATMAT_NT(fy,Y,1.0_wp,A,B)
#define __PMATMAT_TT(fy,Y,A,B) __GENERICMATMAT_TT(fy,Y,1.0_wp,A,B)

#define __AMATMAT_NN(Y,fa,A,B) __GENERICMATMAT_NN(0.0_wp,Y,fa,A,B)
#define __AMATMAT_TN(Y,fa,A,B) __GENERICMATMAT_TN(0.0_wp,Y,fa,A,B)
#define __AMATMAT_NT(Y,fa,A,B) __GENERICMATMAT_NT(0.0_wp,Y,fa,A,B)
#define __AMATMAT_TT(Y,fa,A,B) __GENERICMATMAT_TT(0.0_wp,Y,fa,A,B)

!!! GEMM does in general Y = fa A^?*B^? + fy Y
!!! with structure: (m x n) = (m x k) (k x n)
!!! Y=A  *B   : DGEMM('N','N',m,n,k,fa,Amat ,m, Bmat,k, fy,Y,m)
!!! Y=A^T*B   : DGEMM('T','N',m,n,k,fa,Amat ,k, Bmat,k, fy,Y,m)
!!! Y=A  *B^T : DGEMM('N','T',m,n,k,fa,Amat ,m, Bmat,n, fy,Y,m)
!!! Y=A^T*B^T : DGEMM('T','T',m,n,k,fa,Amat ,k, Bmat,n, fy,Y,m)

!!!#define __GENERICMATMAT_NN(fy,Y,fa,A,B) CALL DGEMM('N','N',SIZE(A,1),SIZE(B,2),SIZE(B,1),fa,A,SIZE(A,1),B,SIZE(B,1),fy,Y,SIZE(A,1))
!!!#define __GENERICMATMAT_TN(fy,Y,fa,A,B) CALL DGEMM('T','N',SIZE(A,2),SIZE(B,2),SIZE(B,1),fa,A,SIZE(B,1),B,SIZE(B,1),fy,Y,SIZE(A,2))
!!!#define __GENERICMATMAT_NT(fy,Y,fa,A,B) CALL DGEMM('N','T',SIZE(A,1),SIZE(B,1),SIZE(B,2),fa,A,SIZE(A,1),B,SIZE(B,1),fy,Y,SIZE(A,1))
!!!#define __GENERICMATMAT_TT(fy,Y,fa,A,B) CALL DGEMM('T','T',SIZE(A,2),SIZE(B,1),SIZE(B,2),fa,A,SIZE(B,2),B,SIZE(B,1),fy,Y,SIZE(A,2))

#define __GENERICMATMAT_NN(fy,Y,fa,A,B) __PADGEMM_NN(fy,Y,fa,SIZE(A,1),SIZE(A,2),A,SIZE(B,1),SIZE(B,2),B)
#define __GENERICMATMAT_TN(fy,Y,fa,A,B) __PADGEMM_TN(fy,Y,fa,SIZE(A,1),SIZE(A,2),A,SIZE(B,1),SIZE(B,2),B)
#define __GENERICMATMAT_NT(fy,Y,fa,A,B) __PADGEMM_NT(fy,Y,fa,SIZE(A,1),SIZE(A,2),A,SIZE(B,1),SIZE(B,2),B)
#define __GENERICMATMAT_TT(fy,Y,fa,A,B) __PADGEMM_TT(fy,Y,fa,SIZE(A,1),SIZE(A,2),A,SIZE(B,1),SIZE(B,2),B)

!!! SIMPLE INTERFACE FOR DGEMM, specifying nrows/ncols of mat A and nrows/ncols of mat B (for any transpose!)

#define __DGEMM_NN(Y,sza1,sza2,A,szb1,szb2,B)     __PADGEMM_NN(0.0_wp,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)
#define __DGEMM_TN(Y,sza1,sza2,A,szb1,szb2,B)     __PADGEMM_TN(0.0_wp,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)
#define __DGEMM_NT(Y,sza1,sza2,A,szb1,szb2,B)     __PADGEMM_NT(0.0_wp,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)
#define __DGEMM_TT(Y,sza1,sza2,A,szb1,szb2,B)     __PADGEMM_TT(0.0_wp,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)

#define __PDGEMM_NN(fy,Y,sza1,sza2,A,szb1,szb2,B) __PADGEMM_NN(fy,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)
#define __PDGEMM_TN(fy,Y,sza1,sza2,A,szb1,szb2,B) __PADGEMM_TN(fy,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)
#define __PDGEMM_NT(fy,Y,sza1,sza2,A,szb1,szb2,B) __PADGEMM_NT(fy,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)
#define __PDGEMM_TT(fy,Y,sza1,sza2,A,szb1,szb2,B) __PADGEMM_TT(fy,Y,1.0_wp,sza1,sza2,A,szb1,szb2,B)

#define __ADGEMM_NN(Y,fa,sza1,sza2,A,szb1,szb2,B) __PADGEMM_NN(0.0_wp,Y,fa,sza1,sza2,A,szb1,szb2,B)
#define __ADGEMM_TN(Y,fa,sza1,sza2,A,szb1,szb2,B) __PADGEMM_TN(0.0_wp,Y,fa,sza1,sza2,A,szb1,szb2,B)
#define __ADGEMM_NT(Y,fa,sza1,sza2,A,szb1,szb2,B) __PADGEMM_NT(0.0_wp,Y,fa,sza1,sza2,A,szb1,szb2,B)
#define __ADGEMM_TT(Y,fa,sza1,sza2,A,szb1,szb2,B) __PADGEMM_TT(0.0_wp,Y,fa,sza1,sza2,A,szb1,szb2,B)

!!! GEMM does in general Y = fa A^?*B^? + fy Y
!!! with structure: (m x n) = (m x k) (k x n)
!!! Y=A  *B   : DGEMM('N','N',m,n,k,fa,Amat ,m, Bmat,k, fy,Y,m)
!!! Y=A^T*B   : DGEMM('T','N',m,n,k,fa,Amat ,k, Bmat,k, fy,Y,m)
!!! Y=A  *B^T : DGEMM('N','T',m,n,k,fa,Amat ,m, Bmat,n, fy,Y,m)
!!! Y=A^T*B^T : DGEMM('T','T',m,n,k,fa,Amat ,k, Bmat,n, fy,Y,m)

#define __PADGEMM_NN(fy,Y,fa,sza1,sza2,A,szb1,szb2,B) CALL DGEMM('N','N',sza1,szb2,szb1,fa,A,sza1,B,szb1,fy,Y,sza1)
#define __PADGEMM_TN(fy,Y,fa,sza1,sza2,A,szb1,szb2,B) CALL DGEMM('T','N',sza2,szb2,szb1,fa,A,szb1,B,szb1,fy,Y,sza2)
#define __PADGEMM_NT(fy,Y,fa,sza1,sza2,A,szb1,szb2,B) CALL DGEMM('N','T',sza1,szb1,szb2,fa,A,sza1,B,szb1,fy,Y,sza1)
#define __PADGEMM_TT(fy,Y,fa,sza1,sza2,A,szb1,szb2,B) CALL DGEMM('T','T',sza2,szb1,szb2,fa,A,szb2,B,szb1,fy,Y,sza2)
