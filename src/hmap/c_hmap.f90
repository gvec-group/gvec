!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** c_hmap **
!!
!! contains only the abstract type to point to a specific map h (maps  omega_p x S^1 --> omega)
!!
!===================================================================================================================================
MODULE MODgvec_c_hmap
! MODULES
USE MODgvec_Globals    ,ONLY:wp,Unit_stdOut,abort
IMPLICIT NONE

PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE, ABSTRACT :: c_hmap
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: which_hmap         !! points to hmap (1: MHD3D)
  INTEGER              :: nfp=-1             !! number of field periods used in hmap. If =-1, its not used
  INTEGER              :: n_max=0            !! maximum number of toroidal modes needed to describe hmap. Used for estimating the number of integration points.
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS
    PROCEDURE(i_sub_hmap_eval_all   ),DEFERRED         :: eval_all
    ! eval?? is a generic name and can be called in three different ways, depending on the arguments :
    ! eval??_pw: pointwise evaluation, without precomputed auxiliary variables (slow)
    ! eval??_aux: pointwise evaluation, with precomputed auxiliary variables (fast)
    ! eval??_aux_all: evaluation on 1d array of points, with array of precomputed auxiliary variables (fast, with omp parallel loop)
    !%eval
    PROCEDURE(i_fun_hmap_eval       ),DEFERRED :: eval
    PROCEDURE                                  :: eval_aux     => hmap_eval_aux
    PROCEDURE                                  :: eval_aux_all => hmap_eval_aux_all
    !%eval_dxdq
    PROCEDURE(i_fun_hmap_eval_dxdq  ),DEFERRED :: eval_dxdq
    PROCEDURE                                  :: eval_dxdq_aux     => hmap_eval_dxdq_aux
    PROCEDURE                                  :: eval_dxdq_aux_all => hmap_eval_dxdq_aux_all    !%eval_Jh
    PROCEDURE(i_fun_hmap_eval_Jh    ),DEFERRED :: eval_Jh
    PROCEDURE                                  :: eval_Jh_aux     => hmap_eval_Jh_aux
    PROCEDURE                                  :: eval_Jh_aux_all => hmap_eval_Jh_aux_all
    !eval_Jh_dq1
    PROCEDURE(i_fun_hmap_eval_Jh_dq ),DEFERRED :: eval_Jh_dq
    PROCEDURE                                  :: eval_Jh_dq_aux     => hmap_eval_Jh_dq_aux
    PROCEDURE                                  :: eval_Jh_dq_aux_all => hmap_eval_Jh_dq_aux_all
    !eval_gij
    PROCEDURE(i_fun_hmap_eval_gij   ),DEFERRED :: eval_gij
    PROCEDURE                                  :: eval_gij_aux     => hmap_eval_gij_aux
    PROCEDURE                                  :: eval_gij_aux_all => hmap_eval_gij_aux_all
    !eval_gij_dq1
    PROCEDURE(i_fun_hmap_eval_gij_dq),DEFERRED :: eval_gij_dq
    PROCEDURE                                  :: eval_gij_dq_aux     => hmap_eval_gij_dq_aux
    PROCEDURE                                  :: eval_gij_dq_aux_all => hmap_eval_gij_dq_aux_all
    PROCEDURE(i_sub_hmap_get_dx_dqi) ,DEFERRED :: get_dx_dqi
    PROCEDURE                                  :: get_dx_dqi_aux => hmap_get_dx_dqi_aux
    PROCEDURE(i_sub_hmap_get_ddx_dqij),DEFERRED:: get_ddx_dqij
    PROCEDURE                                  :: get_ddx_dqij_aux => hmap_get_ddx_dqij_aux
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE c_hmap

TYPE :: c_hmap_auxvar
  REAL(wp) :: zeta
  LOGICAL  :: do_2nd_der
END TYPE c_hmap_auxvar


ABSTRACT INTERFACE

  SUBROUTINE i_sub_hmap_eval_all(sf,ndims,dim_zeta,xv,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                 Jh,    g_tt,    g_tz,    g_zz,&
                                 Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                 Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                 g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
    IMPORT c_hmap,c_hmap_auxvar,wp
    CLASS(c_hmap), INTENT(IN) :: sf
    INTEGER ,INTENT(IN)   :: ndims(3)
    INTEGER ,INTENT(IN)   :: dim_zeta
    CLASS(c_hmap_auxvar),INTENT(IN)   :: xv(ndims(dim_zeta))
    REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(IN) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz
    REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(OUT):: Jh,g_tt    ,g_tz    ,g_zz    , &
                                                                 Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                                                 Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                                                 g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
  END SUBROUTINE i_sub_hmap_eval_all

  !===============================================================================================================================
  !> evaluate the mapping h q=(X^1,X^2,zeta) -> (x,y,z)
  !!
  !===============================================================================================================================
  FUNCTION i_fun_hmap_eval( sf ,q_in) RESULT(x_out)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN) :: sf
    REAL(wp)     , INTENT(IN) :: q_in(3)
    REAL(wp)                  :: x_out(3)
  END FUNCTION i_fun_hmap_eval

  !===============================================================================================================================
  !> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
  !! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
  !!
  !===============================================================================================================================
  FUNCTION i_fun_hmap_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN) :: sf
    REAL(wp)     , INTENT(IN) :: q_in(3)
    REAL(wp)     , INTENT(IN) :: q_vec(3)
    REAL(wp)                  :: dxdq_qvec(3)
  END FUNCTION i_fun_hmap_eval_dxdq

  !===============================================================================================================================
  !> evaluate all first derivatives dx(1:3)/dq^i, i=1,2,3 , at q_in=(X^1,X^2,zeta),
  !!
  !===============================================================================================================================
  SUBROUTINE i_sub_hmap_get_dx_dqi( sf ,q_in,dx_dq1,dx_dq2,dx_dq3)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN)  :: sf
    REAL(wp)     , INTENT(IN)  :: q_in(3)
    REAL(wp)     , INTENT(OUT) :: dx_dq1(3)
    REAL(wp)     , INTENT(OUT) :: dx_dq2(3)
    REAL(wp)     , INTENT(OUT) :: dx_dq3(3)
  END SUBROUTINE i_sub_hmap_get_dx_dqi

  !===============================================================================================================================
  !> evaluate all second derivatives d^2x(1:3)/(dq^i dq^j), i,j=1,2,3 is evaluated at q_in=(X^1,X^2,zeta),
  !!
  !===============================================================================================================================
  SUBROUTINE i_sub_hmap_get_ddx_dqij( sf ,q_in,ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN)  :: sf
    REAL(wp)     , INTENT(IN)  :: q_in(3)
    REAL(wp)     , INTENT(OUT) :: ddx_dq11(3)
    REAL(wp)     , INTENT(OUT) :: ddx_dq12(3)
    REAL(wp)     , INTENT(OUT) :: ddx_dq13(3)
    REAL(wp)     , INTENT(OUT) :: ddx_dq22(3)
    REAL(wp)     , INTENT(OUT) :: ddx_dq23(3)
    REAL(wp)     , INTENT(OUT) :: ddx_dq33(3)
  END SUBROUTINE i_sub_hmap_get_ddx_dqij

  !===============================================================================================================================
  !> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(X^1,X^2,zeta)
  !!
  !===============================================================================================================================
  FUNCTION i_fun_hmap_eval_Jh( sf ,q_in) RESULT(Jh)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN) :: sf
    REAL(wp)     , INTENT(IN) :: q_in(3)
    REAL(wp)                  :: Jh
  END FUNCTION i_fun_hmap_eval_Jh

  !===============================================================================================================================
  !> evaluate derivative of Jacobian of mapping h: sum_k dJ_h(q)/dq^k q_vec^k, k=1,2 at q=(X^1,X^2,zeta)
  !!
  !===============================================================================================================================
  FUNCTION i_fun_hmap_eval_Jh_dq( sf ,q_in,q_vec) RESULT(Jh_dq)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN) :: sf
    REAL(wp)     , INTENT(IN) :: q_in(3)
    REAL(wp)     , INTENT(IN) :: q_vec(3)
    REAL(wp)                  :: Jh_dq
  END FUNCTION i_fun_hmap_eval_Jh_dq

  !===============================================================================================================================
  !>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,
  !! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
  !! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
  !!
  !===============================================================================================================================
  FUNCTION i_fun_hmap_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN) :: sf
    REAL(wp)     , INTENT(IN) :: qL_in(3)
    REAL(wp)     , INTENT(IN) :: q_G(3)
    REAL(wp)     , INTENT(IN) :: qR_in(3)
    REAL(wp)                  :: g_ab
  END FUNCTION i_fun_hmap_eval_gij


  !===============================================================================================================================
  !>  evaluate sum_k sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) q_vec^k, k=1,2,3
  !! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
  !! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
  !!
  !===============================================================================================================================
  FUNCTION i_fun_hmap_eval_gij_dq( sf ,qL_in,q_G,qR_in,q_vec) RESULT(g_ab_dq)
    IMPORT wp,c_hmap
    CLASS(c_hmap), INTENT(IN) :: sf
    REAL(wp)     , INTENT(IN) :: qL_in(3)
    REAL(wp)     , INTENT(IN) :: q_G(3)
    REAL(wp)     , INTENT(IN) :: qR_in(3)
    REAL(wp)     , INTENT(IN) :: q_vec(3)
    REAL(wp)                  :: g_ab_dq
  END FUNCTION i_fun_hmap_eval_gij_dq

END INTERFACE

!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> evaluate the mapping h (X^1,X^2,zeta) -> (x,y,z) cartesian
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the pointwise routine eval
!!
!===================================================================================================================================
FUNCTION hmap_eval_aux( sf ,q1,q2,xv) RESULT(x_out)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  REAL(wp)            ,INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: x_out(3)
!===================================================================================================================================
  x_out=sf%eval((/q1,q2,xv%zeta/))
END FUNCTION hmap_eval_aux

!===================================================================================================================================
!> call %eval_aux on 1d array of points of size np, using auxiliary variable array of same size
!!
!===================================================================================================================================
FUNCTION hmap_eval_aux_all( sf ,np,q1_in,q2_in,xv) RESULT(xyz)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  INTEGER             ,INTENT(IN) :: np
  REAL(wp)            ,INTENT(IN) :: q1_in(1:np),q2_in(1:np)
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: xyz(1:3,1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,np
    xyz(:,i)=sf%eval_aux(q1_in(i),q2_in(i),xv(i))
  END DO
  !$OMP END PARALLEL DO
END FUNCTION hmap_eval_aux_all


!===============================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the generic routine eval_dxdq_pw
!!
!===============================================================================================================================
FUNCTION hmap_eval_dxdq_aux(sf,q1,q2,q1_vec,q2_vec,q3_vec,xv) RESULT(dxdq_qvec)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  REAL(wp)            ,INTENT(IN) :: q1,q2
  REAL(wp)     , INTENT(IN)       :: q1_vec,q2_vec,q3_vec
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES

  REAL(wp)                  :: dxdq_qvec(3)
  !===================================================================================================================================
  dxdq_qvec=sf%eval_dxdq((/q1,q2,xv%zeta/),(/q1_vec,q2_vec,q3_vec/))
END FUNCTION hmap_eval_dxdq_aux

!===============================================================================================================================
!> evaluate all first derivatives dx(1:3)/dq^i, i=1,2,3 , at q_in=(X^1,X^2,zeta),
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the generic routine get_dx_dqi
!!
!===============================================================================================================================
SUBROUTINE hmap_get_dx_dqi_aux(sf,q1,q2,xv,dx_dq1,dx_dq2,dx_dq3)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  REAL(wp)            ,INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT) :: dx_dq1,dx_dq2,dx_dq3
  !===================================================================================================================================
  CALL sf%get_dx_dqi((/q1,q2,xv%zeta/),dx_dq1,dx_dq2,dx_dq3)
END SUBROUTINE hmap_get_dx_dqi_aux

!===============================================================================================================================
!> evaluate all second derivatives d^2x(1:3)/(dq^i dq^j), i,j=1,2,3 , at q_in=(X^1,X^2,zeta),
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the generic routine get_ddx_dqij
!!
!===============================================================================================================================
SUBROUTINE hmap_get_ddx_dqij_aux(sf,q1,q2,xv,ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  REAL(wp)            ,INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT) :: ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33
  !===================================================================================================================================
  CALL sf%get_ddx_dqij((/q1,q2,xv%zeta/),ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
END SUBROUTINE hmap_get_ddx_dqij_aux


!===================================================================================================================================
!> call %eval_dxdq_aux on 1d array of points of size np, using auxiliary variable array of same size
!!
!===================================================================================================================================
FUNCTION hmap_eval_dxdq_aux_all( sf ,np,q1,q2,q1_vec,q2_vec,q3_vec,xv) RESULT(dxdq_qvec)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  INTEGER             ,INTENT(IN) :: np
  REAL(wp)            ,INTENT(IN) :: q1(1:np),q2(1:np)
  REAL(wp)     , INTENT(IN)       :: q1_vec(1:np),q2_vec(1:np),q3_vec(1:np)
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(1:3,1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,np
    dxdq_qvec(:,i)=sf%eval_dxdq_aux(q1(i),q2(i),q1_vec(i),q2_vec(i),q3_vec(i),xv(i))
  END DO
  !$OMP END PARALLEL DO
END FUNCTION hmap_eval_dxdq_aux_all

!===============================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(X^1,X^2,zeta)
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the pointwise routine eval_Jh
!!
!===================================================================================================================================
FUNCTION hmap_eval_Jh_aux( sf ,q1,q2,xv) RESULT(Jh)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  REAL(wp)            ,INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh
!===================================================================================================================================
  Jh=sf%eval_Jh((/q1,q2,xv%zeta/))
END FUNCTION hmap_eval_Jh_aux

!===================================================================================================================================
!> call %eval_Jh_aux on 1d array of points of size np, using auxiliary variable array of same size
!!
!===================================================================================================================================
FUNCTION hmap_eval_Jh_aux_all( sf ,np,q1,q2,xv) RESULT(Jh)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  INTEGER             ,INTENT(IN) :: np
  REAL(wp)            ,INTENT(IN) :: q1(1:np),q2(1:np)
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,np
    Jh(i)=sf%eval_Jh_aux(q1(i),q2(i),xv(i))
  END DO
  !$OMP END PARALLEL DO
END FUNCTION hmap_eval_Jh_aux_all

!===============================================================================================================================
!> evaluate derivative of Jacobian of mapping h: sum_k dJ_h(q)/dq^k q_vec^k, k=1,2 at q=(X^1,X^2,zeta)
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the pointwise routine eval_Jh_dq
!!
!===================================================================================================================================
FUNCTION hmap_eval_Jh_dq_aux( sf ,q1,q2,q1_vec,q2_vec,q3_vec,xv) RESULT(Jh_dq)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  REAL(wp)            ,INTENT(IN) :: q1,q2
  REAL(wp)           , INTENT(IN) :: q1_vec,q2_vec,q3_vec
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh_dq
!===================================================================================================================================
  Jh_dq=sf%eval_Jh_dq((/q1,q2,xv%zeta/),(/q1_vec,q2_vec,q3_vec/))
END FUNCTION hmap_eval_Jh_dq_aux

!===================================================================================================================================
!> call %eval_Jh_dq1_aux on 1d array of points of size np, using auxiliary variable array of same size
!!
!===================================================================================================================================
FUNCTION hmap_eval_Jh_dq_aux_all( sf ,np,q1,q2,q1_vec,q2_vec,q3_vec,xv) RESULT(Jh_dq)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap)       ,INTENT(IN) :: sf
  INTEGER             ,INTENT(IN) :: np
  REAL(wp)            ,INTENT(IN) :: q1(1:np),q2(1:np)
  REAL(wp)     , INTENT(IN)       :: q1_vec(1:np),q2_vec(1:np),q3_vec(1:np)
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh_dq(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,np
    Jh_dq(i)=sf%eval_Jh_dq_aux(q1(i),q2(i),q1_vec(i),q2_vec(i),q3_vec(i),xv(i))
  END DO
  !$OMP END PARALLEL DO
END FUNCTION hmap_eval_Jh_dq_aux_all


!===============================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0.0 for ds and dtheta and 1.0 for dzeta
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the pointwise routine eval_gij
!!
!===============================================================================================================================
FUNCTION hmap_eval_gij_aux( sf ,qL1,qL2,qL3,q1,q2,qR1,qR2,qR3,xv) RESULT(g_ab)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap), INTENT(IN) :: sf
  REAL(wp)     , INTENT(IN) :: qL1,qL2,qL3
  REAL(wp)     , INTENT(IN) :: q1,q2
  REAL(wp)     , INTENT(IN) :: qR1,qR2,qR3
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                  :: g_ab
!===================================================================================================================================
  g_ab=sf%eval_gij((/qL1,qL2,qL3/),(/q1,q2,xv%zeta/),(/qR1,qR2,qR3/))
END FUNCTION hmap_eval_gij_aux

FUNCTION hmap_eval_gij_aux_all( sf ,np,qL1,qL2,qL3,q1,q2,qR1,qR2,qR3,xv) RESULT(g_ab)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap), INTENT(IN) :: sf
  INTEGER     , INTENT(IN) :: np
  REAL(wp)     , INTENT(IN) :: qL1(1:np),qL2(1:np),qL3(1:np)
  REAL(wp)     , INTENT(IN) :: q1(1:np),q2(1:np)
  REAL(wp)     , INTENT(IN) :: qR1(1:np),qR2(1:np),qR3(1:np)
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                  :: g_ab(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,np
    g_ab(i)=sf%eval_gij_aux(qL1(i),qL2(i),qL3(i),q1(i),q2(i),qR1(i),qR2(i),qR3(i),xv(i))
  END DO
  !$OMP END PARALLEL DO
END FUNCTION hmap_eval_gij_aux_all

!===============================================================================================================================
!>  evaluate sum_k sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) q_vec^k, k=1,2,3
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0.0 for ds and dtheta and 1.0 for dzeta
!! INFO: default routine that can be overwritten by specific hmap class,
!!       not using  additional hmap-dependent auxiliary variables,
!!       but calling the pointwise routine eval_gij_dq
!!
!===============================================================================================================================
FUNCTION hmap_eval_gij_dq_aux( sf ,qL1,qL2,qL3,q1,q2,qR1,qR2,qR3,q1_vec,q2_vec,q3_vec,xv) RESULT(g_ab_dq)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap), INTENT(IN) :: sf
  REAL(wp)     , INTENT(IN) :: qL1,qL2,qL3
  REAL(wp)     , INTENT(IN) :: q1,q2
  REAL(wp)     , INTENT(IN) :: qR1,qR2,qR3
  REAL(wp)     , INTENT(IN) :: q1_vec,q2_vec,q3_vec
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                  :: g_ab_dq
!===================================================================================================================================
  g_ab_dq=sf%eval_gij_dq((/qL1,qL2,qL3/),(/q1,q2,xv%zeta/),(/qR1,qR2,qR3/),(/q1_vec,q2_vec,q3_vec/))
END FUNCTION hmap_eval_gij_dq_aux

FUNCTION hmap_eval_gij_dq_aux_all( sf ,np,qL1,qL2,qL3,q1,q2,qR1,qR2,qR3,q1_vec,q2_vec,q3_vec,xv) RESULT(g_ab_dq)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap), INTENT(IN) :: sf
  INTEGER     , INTENT(IN) :: np
  REAL(wp)     , INTENT(IN) :: qL1(1:np),qL2(1:np),qL3(1:np)
  REAL(wp)     , INTENT(IN) :: q1(1:np),q2(1:np)
  REAL(wp)     , INTENT(IN) :: qR1(1:np),qR2(1:np),qR3(1:np)
  REAL(wp)     , INTENT(IN) :: q1_vec(1:np),q2_vec(1:np),q3_vec(1:np)
  CLASS(c_hmap_auxvar),INTENT(IN) :: xv(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                  :: g_ab_dq(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,np
    g_ab_dq(i)=sf%eval_gij_dq_aux(qL1(i),qL2(i),qL3(i),q1(i),q2(i),qR1(i),qR2(i),qR3(i),q1_vec(i),q2_vec(i),q3_vec(i),xv(i))
  END DO
  !$OMP END PARALLEL DO
END FUNCTION hmap_eval_gij_dq_aux_all

END MODULE MODgvec_c_hmap
