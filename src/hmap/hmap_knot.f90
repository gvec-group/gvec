!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** hmap_knot **
!!
!! contains the type that points to the routines of one chosen hmap_knot
!!
!===================================================================================================================================
MODULE MODgvec_hmap_knot
! MODULES
USE MODgvec_Globals, ONLY:PI,wp,Unit_stdOut,abort,MPIRoot
USE MODgvec_c_hmap,    ONLY:c_hmap, c_hmap_auxvar
IMPLICIT NONE

PUBLIC

TYPE,EXTENDS(c_hmap_auxvar) :: t_hmap_knot_auxvar
  !nothing more to add for knot hmap
END TYPE t_hmap_knot_auxvar

TYPE,EXTENDS(c_hmap) :: t_hmap_knot
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_knot:
  REAL(wp) :: k,  l    ! this map is based on the (k,l)-torus
  REAL(wp) :: R0       ! major radius
  REAL(wp) :: delta    ! shift of the axis
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  FINAL     :: hmap_knot_free
  PROCEDURE :: eval_all         => hmap_knot_eval_all
  PROCEDURE :: eval             => hmap_knot_eval
  PROCEDURE :: eval_dxdq        => hmap_knot_eval_dxdq
  PROCEDURE :: eval_Jh          => hmap_knot_eval_Jh
  PROCEDURE :: eval_Jh_dq       => hmap_knot_eval_Jh_dq
  PROCEDURE :: eval_gij         => hmap_knot_eval_gij
  PROCEDURE :: eval_gij_dq      => hmap_knot_eval_gij_dq
  PROCEDURE :: get_dx_dqi       => hmap_knot_get_dx_dqi
  PROCEDURE :: get_ddx_dqij     => hmap_knot_get_ddx_dqij
  !---------------------------------------------------------------------------------------------------------------------------------
  ! procedures for hmap_knot:
  PROCEDURE :: Rl            => hmap_knot_eval_Rl
  PROCEDURE :: Zl            => hmap_knot_eval_Zl
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap_knot

!INITIALIZATION FUNCTION:
INTERFACE t_hmap_knot
  MODULE PROCEDURE hmap_knot_init,hmap_knot_init_params
END INTERFACE t_hmap_knot

INTERFACE t_hmap_knot_auxvar
  MODULE PROCEDURE hmap_knot_init_aux
END INTERFACE t_hmap_knot_auxvar

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize the type hmap_knot, reading from parameter file and then call init_params
!!
!===================================================================================================================================
FUNCTION hmap_knot_init() RESULT(sf)
  ! MODULES
    USE MODgvec_ReadInTools, ONLY: GETINTARRAY, GETREAL
    IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    TYPE(t_hmap_knot) :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER                           :: knot_kl(1:2)         !parameters of the (k,l)-torus
    REAL(wp)                          :: knot_R0, knot_delta    !major radius and shift
  !===================================================================================================================================
    SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: KNOT ON A (k,l)-TORUS, GET PARAMETERS:'

    knot_kl=GETINTARRAY("hmap_knot_kl",2,proposal=(/2,3/))

    knot_R0=GETREAL("hmap_knot_major_radius",1.0_wp)


    knot_delta=GETREAL("hmap_knot_delta_shift",0.4_wp)

    sf=hmap_knot_init_params(knot_kl,knot_R0,knot_delta)
  END FUNCTION hmap_knot_init


!===================================================================================================================================
!> initialize the type hmap_knot, from given parameters as arguments
!!
!===================================================================================================================================
FUNCTION hmap_knot_init_params(knot_kl,knot_R0,knot_delta) RESULT(sf)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,INTENT(IN)  :: knot_kl(1:2)         !parameters of the (k,l)-torus
  REAL(wp),INTENT(IN) :: knot_R0, knot_delta    !major radius and shift
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_knot) :: sf !! self
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: KNOT ON A (k,l)-TORUS ...'

  sf%k=REAL(knot_kl(1), wp)
  sf%l=REAL(knot_kl(2), wp)
  sf%R0=knot_R0
  sf%delta=knot_delta
  sf%n_max=MAX(sf%k,sf%l)

  IF (.NOT.((sf%R0 - ABS(sf%delta)) > 0.0_wp)) THEN
     CALL abort(__STAMP__, &
          "hmap_knot init: condition R0 - |delta| > 0 not fulfilled!", &
          TypeInfo="InvalidParameterError")
  END IF

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_knot_test(sf)

END FUNCTION hmap_knot_init_params


!===================================================================================================================================
!> finalize the type hmap_knot
!!
!===================================================================================================================================
SUBROUTINE hmap_knot_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_knot), INTENT(INOUT) :: sf !! self
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN

  sf%initialized=.FALSE.
  sf%R0 = 0.0_wp
  sf%delta = -1.0_wp
  sf%k = 0.0_wp
  sf%l = 0.0_wp

END SUBROUTINE hmap_knot_free


!===================================================================================================================================
!> Allocate and initialize auxiliary variable at zeta position.
!!
!===================================================================================================================================
FUNCTION hmap_knot_init_aux( sf,zeta,do_2nd_der) RESULT(xv)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot),INTENT(IN) :: sf
  REAL(wp)          ,INTENT(IN) :: zeta
  LOGICAL           ,INTENT(IN) :: do_2nd_der !! compute second derivative and store second derivative terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_knot_auxvar)::xv
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  xv%do_2nd_der=do_2nd_der
  xv%zeta=zeta
END FUNCTION hmap_knot_init_aux

!===================================================================================================================================
!> evaluate all metrics necesseray for optimizer
!!
!===================================================================================================================================
SUBROUTINE hmap_knot_eval_all(sf,ndims,dim_zeta,xv,&
                                q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                Jh,    g_tt,    g_tz,    g_zz,&
                                Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot)  , INTENT(IN)   :: sf
  INTEGER             , INTENT(IN)   :: ndims(3)    !! 3D dimensions of input arrays
  INTEGER             , INTENT(IN)   :: dim_zeta    !! which dimension is zeta dependent
  CLASS(c_hmap_auxvar), INTENT(IN)   :: xv(ndims(dim_zeta))  !! zeta point positions
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(IN) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(OUT):: Jh    ,g_tt    ,g_tz    ,g_zz    , &
                                                               Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                                               Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                                               g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,j,k
  !===================================================================================================================================
  SELECT TYPE(xv)
  TYPE IS(t_hmap_knot_auxvar)
    SELECT CASE(dim_zeta)
    CASE(1)
      !$OMP PARALLEL DO &
      !$OMP COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        CALL hmap_knot_eval_all_e(sf%k,sf%l,sf%delta,sf%R0,xv(i)%zeta, &
                 q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
                 Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
                 Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
                 Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
                 g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
      END DO; END DO; END DO
      !$OMP END PARALLEL DO
    CASE(2)
!      !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        CALL hmap_knot_eval_all_e(sf%k,sf%l,sf%delta,sf%R0,xv(j)%zeta, &
                 q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
                 Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
                 Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
                 Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
                 g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
      END DO; END DO; END DO
!      !$OMP END PARALLEL DO
    CASE(3)
!      !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        CALL hmap_knot_eval_all_e(sf%k,sf%l,sf%delta,sf%R0,xv(k)%zeta, &
                 q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
                 Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
                 Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
                 Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
                 g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
      END DO; END DO; END DO
!      !$OMP END PARALLEL DO
    END SELECT !dim_zeta
  END SELECT !TYPE(xv)

END SUBROUTINE hmap_knot_eval_all

!===================================================================================================================================
!> evaluate all quantities at one given point (elemental)
!! NOTE: using calls to sf, not implemented/optimized for performance yet!
!!
!===================================================================================================================================
PURE SUBROUTINE hmap_knot_eval_all_e(k,l,delta,R0,zeta,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                       Jh,    g_tt,    g_tz,    g_zz,     &
                                       Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                       Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                       g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp),INTENT(IN)  :: k,l,delta,R0 !! hmap parameters
  REAL(wp),INTENT(IN)  :: zeta        !! zeta position
  REAL(wp),INTENT(IN)  :: q1,q2       !! solution variables q1,q2
  REAL(wp),INTENT(IN)  :: dX1_dt,dX2_dt  !! theta derivative of solution variables q1,q2
  REAL(wp),INTENT(IN)  :: dX1_dz,dX2_dz  !!  zeta derivative of solution variables q1,q2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),INTENT(OUT) :: Jh,g_tt,g_tz,g_zz              !! Jac,1/Jac,g_{ab} with a=theta/zeta b=theta/zeta
  REAL(wp),INTENT(OUT) :: Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1  !! and their variation vs q1
  REAL(wp),INTENT(OUT) :: Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2  !! and their variation vs q2
  REAL(wp),INTENT(OUT) :: g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  !! dq^{i}/dtheta*G^{i1}, dq^{i}/dtheta*G^{i2}, and dq^{i}/dzeta*G^{i1}, dq^{i}/dzeta*G^{i2} and G^{11},G^{22}
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: Rl,Gh31,Gh32,Gh33
!===================================================================================================================================
  Rl = R0 + delta*COS(l*zeta) + q1
  Jh = k*Rl
  Jh_dq1 = k
  Jh_dq2 = 0.0_wp
  Gh11 = 1.0_wp

  Gh22 = 1.0_wp
  Gh31 =-l*delta*SIN(l*zeta)
  Gh32 = l*delta*COS(l*zeta)
  Gh33 = (k * Rl)**2 + (l* delta)**2

  g_t1 = dX1_dt
  g_t2 = dX2_dt
  g_z1 = dX1_dz + Gh31
  g_z2 = dX2_dz + Gh32

  g_tt =   dX1_dt *  g_t1         +  dX2_dt *  g_t2
  g_tz =   dX1_dt *  g_z1         +  dX2_dt *  g_z2
  g_zz =   dX1_dz * (g_z1 + Gh31) +  dX2_dz * (g_z2 + Gh32)  + Gh33

  !Gh11/dq1 =0             Gh13/dq1 = 0
  !            Gh22/dq1 =0 Gh23/dq1 = 0
  !                        Gh33/dq1 = 2k**2 *Rl
  !Gh11/dq2 =0 Gh12/dq2 =0 Gh13/dq2 = 0
  !            Gh22/dq2 =0 Gh23/dq2 = 0
  !                        Gh33/dq2 = 0
  ! => g_t1 /dq1 =0, g_t1/dq2 =0, g_t2/dq1 =0, g_t2/dq2 =0
  ! => g_z1 /dq1 = 0, g_z1/dq2 =0, g_z2/dq1 =0, g_z2/dq2 =0

  g_tt_dq1 = 0.0_wp
  g_tz_dq1 = 0.0_wp
  g_zz_dq1 = 2.0_wp * k**2 * Rl

  g_tt_dq2 = 0.0_wp
  g_tz_dq2 = 0.0_wp
  g_zz_dq2 = 0.0_wp

END SUBROUTINE hmap_knot_eval_all_e

!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z)
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval( sf ,q_in) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN) :: q_in(3)
  CLASS(t_hmap_knot), INTENT(IN) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: x_out(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
 ! (k,l) are the indides of the (k,l)-torus
 ! q(:) = (q1,q2,zeta) are the variables in the domain of the map
 ! X(:) = (x,y,z) are the variables in the range of the map
 !
 !   Rl = R0 + delta * cos(l*zeta) + q1
 !   Zl = delta * sin(l*zeta) + q2
 !  |x |  | Rl*sin(k*zeta) |
 !  |y |= |-Rl*cos(k*zeta) |
 !  |z |  | Zl             |

 ASSOCIATE(zeta=>q_in(3))
 x_out(1:3)=(/ sf%Rl(q_in)*COS(sf%k*zeta), &
              -sf%Rl(q_in)*SIN(sf%k*zeta), &
               sf%Zl(q_in)                 /)
 END ASSOCIATE
END FUNCTION hmap_knot_eval

!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: q_in(3)
  REAL(wp)          , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: dxdq_qvec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: coskzeta,sinkzeta
!===================================================================================================================================
 ASSOCIATE(zeta=>q_in(3))
 coskzeta=COS(sf%k*zeta)
 sinkzeta=SIN(sf%k*zeta)
 dxdq_qvec(1:3)=   (/ coskzeta*q_vec(1),-sinkzeta*q_vec(1),q_vec(2)/) &
                  +(/ -sf%k*sf%Rl(q_in)*sinkzeta-sf%l*sf%delta*SIN(sf%l*zeta)*coskzeta, &
                      -sf%k*sf%Rl(q_in)*coskzeta+sf%l*sf%delta*SIN(sf%l*zeta)*sinkzeta, &
                                                 sf%l*sf%delta*COS(sf%l*zeta)      /)*q_vec(3)
 END ASSOCIATE


END FUNCTION hmap_knot_eval_dxdq

!===============================================================================================================================
!> evaluate all first derivatives dx(1:3)/dq^i, i=1,2,3 , at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_knot_get_dx_dqi( sf ,q_in,dx_dq1,dx_dq2,dx_dq3)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: q_in(3)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT):: dx_dq1,dx_dq2,dx_dq3
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  REAL(wp) :: coskzeta,sinkzeta
  !===================================================================================================================================
  ASSOCIATE(zeta=>q_in(3))
  coskzeta=COS(sf%k*zeta)
  sinkzeta=SIN(sf%k*zeta)
  dx_dq1(1:3)=   (/ coskzeta, -sinkzeta, 0.0_wp /)
  dx_dq2(1:3)=   (/ 0.0_wp,0.0_wp,1.0_wp /)
  dx_dq3(1:3)=   (/ -sf%k*sf%Rl(q_in)*sinkzeta-sf%l*sf%delta*SIN(sf%l*zeta)*coskzeta, &
                    -sf%k*sf%Rl(q_in)*coskzeta+sf%l*sf%delta*SIN(sf%l*zeta)*sinkzeta, &
                                               sf%l*sf%delta*COS(sf%l*zeta)      /)
 END ASSOCIATE
END SUBROUTINE hmap_knot_get_dx_dqi

!=================================================================================================================================
!> evaluate all second derivatives d^2x(1:3)/(dq^i dq^j), i,j=1,2,3 is evaluated at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_knot_get_ddx_dqij( sf ,q_in,ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: q_in(3)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT):: ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  REAL(wp) :: coskzeta,sinkzeta,dsinlzeta,dcoslzeta
  !===================================================================================================================================
  ASSOCIATE(zeta=>q_in(3))
  coskzeta=COS(sf%k*zeta)
  sinkzeta=SIN(sf%k*zeta)
  dsinlzeta=sf%delta*SIN(sf%l*zeta)
  dcoslzeta=sf%delta*COS(sf%l*zeta)
  !   Rl = R0 + delta * cos(l*zeta) + q1
  !   Zl = delta * sin(l*zeta) + q2
  ddx_dq11=0.0_wp
  ddx_dq12=0.0_wp
  ddx_dq13=(/ -sf%k*sinkzeta,-sf%k*coskzeta, 0.0_wp /)
  ddx_dq22=0.0_wp
  ddx_dq23=0.0_wp
  ddx_dq33(1)= -sf%k*sf%k*sf%Rl(q_in)*coskzeta+sf%l*sf%k*dsinlzeta*sinkzeta &
               +sf%k*sf%l*dsinlzeta  *sinkzeta-sf%l*sf%l*dcoslzeta*coskzeta
  ddx_dq33(2)=  sf%k*sf%k*sf%Rl(q_in)*sinkzeta+sf%l*sf%k*dsinlzeta*coskzeta  &
               +sf%k*sf%l*dsinlzeta  *coskzeta+sf%l*sf%l*dcoslzeta*sinkzeta
  ddx_dq33(3)= -sf%l*sf%l*dsinlzeta
  END ASSOCIATE
  END SUBROUTINE hmap_knot_get_ddx_dqij

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(q^1,q^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: Jh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Jh = sf%k*sf%Rl(q_in)   ! Jh = k * Rl
END FUNCTION hmap_knot_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: sum_k q_vec^k * dJ_h/dq^k, k=1,2,3 at q=(q^1,q^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh_dq( sf ,q_in,q_vec) RESULT(Jh_dq)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: q_in(3)
  REAL(wp)          , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: Jh_dq
!===================================================================================================================================
  Jh_dq = sf%k*(q_vec(1) -sf%delta*sf%l*SIN(sf%l*q_in(3))*q_vec(3)) ! dJh/dq^1 = d(kRl)/dq^1  dJh/dq^3 = dkRl/dzeta
END FUNCTION hmap_knot_eval_Jh_dq


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0.0 for ds and dtheta and 1.0 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: qL_in(3)
  REAL(wp)          , INTENT(IN) :: q_G(3)
  REAL(wp)          , INTENT(IN) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: g_ab
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                       :: A, B, C
!===================================================================================================================================
  ! A = - l * delta * sin(l*zeta),
  ! B = l * delta * cos(l*zeta)
  ! C = k**2 * Rl**2 + l**2 * delta**2
  !                       |q1  |   |1  0  A|        |q1  |
  !q_i G_ij q_j = (dalpha |q2  | ) |0  1  B| (dbeta |q2  | )
  !                       |q3  |   |A  B  C|        |q3  |
 ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
   A = - sf%l*sf%delta*SIN(sf%l*zeta)
   B = sf%l*sf%delta*COS(sf%l*zeta)
   C = sf%k**2 * sf%Rl(q_G)**2 + sf%l**2 * sf%delta**2
   g_ab=SUM(qL_in(:)*(/qR_in(1) + A*qR_in(3), qR_in(2) + B*qR_in(3), A*qR_in(1) + B*qR_in(2) + C*qR_in(3)/))
 END ASSOCIATE
END FUNCTION hmap_knot_eval_gij


!===================================================================================================================================
!>  evaluate sum_k sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) *q_vec^k, k=1,2,3
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0.0 for ds and dtheta and 1.0 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_gij_dq( sf ,qL_in,q_G,qR_in,q_vec) RESULT(g_ab_dq)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(IN) :: sf
  REAL(wp)          , INTENT(IN) :: qL_in(3)
  REAL(wp)          , INTENT(IN) :: q_G(3)
  REAL(wp)          , INTENT(IN) :: qR_in(3)
  REAL(wp)          , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: g_ab_dq
!===================================================================================================================================
  ! Rl = R0 + delta * cos(l*zeta) + q1
  ! A = - l * delta * sin(l*zeta),
  ! B = l * delta * cos(l*zeta)
  ! C = k**2 * Rl**2 + l**2 * delta**2
  !                       |q1  |   |1  0  A|        |q1  |
  !q_i G_ij q_j = (dalpha |q2  | ) |0  1  B| (dbeta |q2  | )
  !                       |q3  |   |A  B  C|        |q3  |
  ! dA_dq1 = 0 , dA_dq2 = 0 , dA_dq3 = - l*l * delta * cos(l*zeta)
  ! dB_dq1 = 0 , dB_dq2 = 0 , dB_dq3 = - l*l * delta * sin(l*zeta)
  ! dC_dq1 = 2*k**2 *Rl , dC_dq2 = 0 , dC_dq3 = -2*k**2 *Rl *l *delta * sin(l*zeta)
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  g_ab_dq =-sf%l*sf%l*sf%delta*( COS(sf%l*zeta)*(qL_in(1)*qR_in(3) + qL_in(3)*qR_in(1)) &
                                +SIN(sf%l*zeta)*(qL_in(2)*qR_in(3) + qL_in(3)*qR_in(2)) ) *q_vec(3) &
                  +2*sf%k**2 * sf%Rl(q_G) *qL_in(3)*qR_in(3)*(q_vec(1) - sf%l*sf%delta*SIN(sf%l*zeta)*q_vec(3))
  END ASSOCIATE
END FUNCTION hmap_knot_eval_gij_dq

!===================================================================================================================================
!> evaluate the effective major radius coordinate Rl(q)
!!
!===================================================================================================================================
PURE FUNCTION hmap_knot_eval_Rl( sf ,q_in) RESULT(Rl_out)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN) :: q_in(3)
  CLASS(t_hmap_knot), INTENT(IN) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: Rl_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
 !   Rl = R0 + delta * cos(l*zeta) + q1

 ASSOCIATE(q1=>q_in(1),zeta=>q_in(3))
   Rl_out = sf%R0 + sf%delta*COS(sf%l*zeta) + q1
 END ASSOCIATE
END FUNCTION hmap_knot_eval_Rl


!===================================================================================================================================
!> evaluate the effective vertical coordinate Zl(q)
!!
!===================================================================================================================================
PURE FUNCTION hmap_knot_eval_Zl( sf ,q_in) RESULT(Zl_out)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN) :: q_in(3)
  CLASS(t_hmap_knot), INTENT(IN) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: Zl_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  !   Zl = delta * sin(l*zeta) + q2

 ASSOCIATE(q2=>q_in(2),zeta=>q_in(3))
   Zl_out = sf%delta*SIN(sf%l*zeta) + q2
 END ASSOCIATE
END FUNCTION hmap_knot_eval_Zl


!===================================================================================================================================
!> test hmap_knot - evaluation of the map
!!
!===================================================================================================================================
SUBROUTINE hmap_knot_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,idir,jdir,kdir,qdir,izeta,i,j,k,ndims(1:3),ijk(3)
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3),q_test(3,3),x_eps(3),dxdq(3),gij,gij_eps,Jh_0,Jh_eps
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  REAL(wp),PARAMETER :: epsFD=1.0e-8
  REAL(wp),PARAMETER :: realtolFD=1.0e-4
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: a, Rl, Zl
  INTEGER,PARAMETER  :: nzeta=5
  INTEGER,PARAMETER  :: ns=2
  INTEGER,PARAMETER  :: nthet=3
  REAL(wp),ALLOCATABLE :: zeta(:)
  REAL(wp)           :: qloc(3),q_thet(3),q_zeta(3)
  REAL(wp)           :: dxdq_eps(3),dx_dqi(3,3),ddx_dqij(3,3,3)
  REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                     Jh,g_tt,    g_tz,    g_zz,     &
                                     Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                     Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                     g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
  TYPE(t_hmap_knot_auxvar),ALLOCATABLE :: xv(:)
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_knot TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    a = sf%R0 - ABS(sf%delta)
    q_in=(/0.5_wp*a, -0.2_wp*a, 0.5_wp*Pi/)
    Rl = sf%R0 + sf%delta*COS(sf%l*q_in(3)) + q_in(1)
    Zl = sf%delta*SIN(sf%l*q_in(3)) + q_in(2)
    x = sf%eval(q_in )
    checkreal=SUM((x-(/Rl*COS(sf%k*q_in(3)),-Rl*SIN(sf%k*q_in(3)),Zl/))**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    a = sf%R0 - ABS(sf%delta)
    q_in=(/0.3_wp*a,  0.1_wp*a, 0.4_wp*Pi/)
    Rl = sf%R0 + sf%delta*COS(sf%l*q_in(3)) + q_in(1)
    x = sf%eval_dxdq(q_in, (/1.1_wp,1.2_wp,1.3_wp/) )
    checkreal=SUM((x-( (/1.1_wp*COS(sf%k*q_in(3)),-1.1_wp*SIN(sf%k*q_in(3)),1.2_wp/) &
                      +1.3_wp*(/-(sf%k*Rl*SIN(sf%k*q_in(3))+sf%l*sf%delta*SIN(sf%l*q_in(3))*COS(sf%k*q_in(3))),&
                                -(sf%k*Rl*COS(sf%k*q_in(3))-sf%l*sf%delta*SIN(sf%l*q_in(3))*SIN(sf%k*q_in(3))),&
                                                            sf%l*sf%delta*COS(sf%l*q_in(3)) /)) )**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    q_test(1,:)=(/1.0_wp, 0.0_wp, 0.0_wp/)
    q_test(2,:)=(/0.0_wp, 1.0_wp, 0.0_wp/)
    q_test(3,:)=(/0.0_wp, 0.0_wp, 1.0_wp/)
    DO qdir=1,3
      !check dx/dq^i with FD
      iTest=102+qdir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      q_in=(/0.0_wp, 0.0_wp, 0.335_wp*PI/)
      x = sf%eval(q_in )
      x_eps = sf%eval(q_in+epsFD*q_test(qdir,:))
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      checkreal=SQRT(SUM((dxdq - (x_eps-x)/epsFD)**2)/SUM(x*x))
      refreal = 0.0_wp

      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be <',realtolFD,' : |dxdqFD-eval_dxdq|= ', checkreal,", dq=",qdir
      END IF !TEST
    END DO

    !! TEST G_ij
    DO idir=1,3; DO jdir=idir,3
      iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal= SUM(sf%eval_dxdq(q_in,q_test(idir,:))*sf%eval_dxdq(q_in,q_test(jdir,:)))
      refreal  =sf%eval_gij(q_test(idir,:),q_in,q_test(jdir,:))
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),2(A,I3))') &
       '\n =>  should be ', refreal,' : sum|Gij-eval_gij|= ', checkreal,', i=',idir,', j=',jdir
      END IF !TEST
    END DO; END DO
    !! TEST dJh_dqk with FD
    DO qdir=1,3
      iTest=iTest+1; IF(testdbg)WRITE(*,*)'iTest=',iTest
      Jh_0    = sf%eval_Jh(   q_in                     )
      Jh_eps  = sf%eval_Jh(   q_in+epsFD*q_test(qdir,:))
      refreal = sf%eval_Jh_dq(q_in                     ,q_test(qdir,:))
      checkreal=(Jh_eps-Jh_0)/epsFD-refreal
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be <', realtolFD,' : |dJh_dqFD-eval_Jh_dq|= ', checkreal,', dq=',qdir
      END IF !TEST
    END DO !qdir
    !! TEST dG_ij_dqk with FD
    DO qdir=1,3
    DO idir=1,3; DO jdir=1,3
      iTest=iTest+1; IF(testdbg)WRITE(*,*)'iTest=',iTest
      gij     = sf%eval_gij(   q_test(idir,:),q_in                     ,q_test(jdir,:))
      gij_eps = sf%eval_gij(   q_test(idir,:),q_in+epsFD*q_test(qdir,:),q_test(jdir,:))
      refreal = sf%eval_gij_dq(q_test(idir,:),q_in                     ,q_test(jdir,:),q_test(qdir,:))
      checkreal=(gij_eps-gij)/epsFD-refreal
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),3(A,I3))') &
       '\n =>  should be <', realtolFD,' : |dGij_dqFD-eval_gij_dq|= ', checkreal,', i=',idir,', j=',jdir,', dq=',qdir
      END IF !TEST
    END DO; END DO
    END DO !qdir


    CALL sf%get_dx_dqi(q_in,dx_dqi(1,:),dx_dqi(2,:),dx_dqi(3,:))
    DO qdir=1,3
      !check dx/dq^i with FD
      iTest=10+qdir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      checkreal=SUM(ABS(dxdq-dx_dqi(qdir,:)))
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be ', refreal,' : |dxdq-eval_dxdq|= ', checkreal,", dq=",qdir
      END IF !TEST
    END DO !qdir

    CALL sf%get_ddx_dqij(q_in,ddx_dqij(1,1,:),ddx_dqij(1,2,:),ddx_dqij(1,3,:), &
                              ddx_dqij(2,2,:),ddx_dqij(2,3,:),ddx_dqij(3,3,:))
    ddx_dqij(2,1,:)=ddx_dqij(1,2,:)
    ddx_dqij(3,1,:)=ddx_dqij(1,3,:)
    ddx_dqij(3,2,:)=ddx_dqij(2,3,:)
    DO qdir=1,3
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      DO idir=1,3
        iTest=20+idir+3*(qdir-1) ; IF(testdbg)WRITE(*,*)'iTest=',iTest
        dxdq_eps = sf%eval_dxdq(q_in+epsFD*q_test(idir,:),q_test(qdir,:))
        checkreal=SUM(ABS( (dxdq_eps-dxdq)/epsFD-ddx_dqij(qdir,idir,:)))
        refreal=0.0_wp
        IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
           nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
                '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
           nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),2(A,I3))') &
         '\n =>  should be <', realtolFD,' : |ddx_dqijFD-eval_ddx_dqij|= ', checkreal,", dqi=",qdir,", dqj=",idir
        END IF !TEST
      END DO !idir
    END DO !qdir
 END IF !testlevel >=1
 IF (testlevel .GE. 2) THEN
    DO idir=1,3
      SELECT CASE(idir)
      CASE(1)
        jdir=2; kdir=3
      CASE(2)
        jdir=1; kdir=3
      CASE(3)
        jdir=1; kdir=2
      END SELECT
      ndims(idir)=nzeta+idir
      ndims(jdir)=ns
      ndims(kdir)=nthet
      ALLOCATE(zeta(ndims(idir)),xv(ndims(idir)))
      DO izeta=1,ndims(idir)
        zeta(izeta)=0.333_wp+REAL(izeta-1,wp)/REAL(ndims(idir)-1,wp)*0.221_wp
        xv(izeta)=hmap_knot_init_aux(sf,zeta(izeta),.TRUE.)
      END DO
      ALLOCATE(q1(ndims(1),ndims(2),ndims(3)))
      ALLOCATE(q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz,Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2,g_t1,g_t2,g_z1,g_z2,Gh11,Gh22, &
               mold=q1)
      !assign somewhat randomly
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        q1(i,j,k) = 0.11_wp -0.22_wp *REAL((i+j)*k,wp)/REAL((ndims(idir)+ndims(jdir))*ndims(kdir),wp)
        q2(i,j,k) = 0.15_wp -0.231_wp*REAL((i+k)*j,wp)/REAL((ndims(idir)+ndims(kdir))*ndims(jdir),wp)
        dX1_dt(i,j,k)=-0.1_wp  +0.211_wp*REAL((i+2*j)*k,wp)/REAL((ndims(idir)+2*ndims(jdir))*ndims(kdir),wp)
        dX2_dt(i,j,k)= 0.231_wp-0.116_wp*REAL((2*i+k)*j,wp)/REAL((2*ndims(idir)+ndims(kdir))*ndims(jdir),wp)
        dX1_dz(i,j,k)=-0.024_wp+0.013_wp*REAL((3*i+2*j)*k,wp)/REAL((3*ndims(idir)+2*ndims(jdir))*ndims(kdir),wp)
        dX2_dz(i,j,k)=-0.06_wp +0.031_wp*REAL((2*k+3*k)*i,wp)/REAL((2*ndims(kdir)+3*ndims(kdir))*ndims(idir),wp)
      END DO; END DO; END DO
      CALL sf%eval_all(ndims,idir,xv,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
           Jh,g_tt,g_tz,g_zz,&
           Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,&
           Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2,&
           g_t1,g_t2,g_z1,g_z2,Gh11,Gh22)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        ijk=(/i,j,k/)
        izeta=ijk(idir)
        qloc=(/q1(i,j,k),q2(i,j,k),zeta(izeta)/)
        q_thet=(/dX1_dt(i,j,k),dX2_dt(i,j,k),0.0_wp/)
        q_zeta=(/dX1_dz(i,j,k),dX2_dz(i,j,k),1.0_wp/)
        Jh(i,j,k)       =Jh(i,j,k)       - sf%eval_Jh(qloc)
        g_tt(i,j,k)     =g_tt(i,j,k)     - sf%eval_gij(q_thet,qloc,q_thet)
        g_tz(i,j,k)     =g_tz(i,j,k)     - sf%eval_gij(q_thet,qloc,q_zeta)
        g_zz(i,j,k)     =g_zz(i,j,k)     - sf%eval_gij(q_zeta,qloc,q_zeta)
        Jh_dq1(i,j,k)   =Jh_dq1(i,j,k)   - sf%eval_Jh_dq(qloc,(/1.0_wp,0.0_wp,0.0_wp/))
        Jh_dq2(i,j,k)   =Jh_dq2(i,j,k)   - sf%eval_Jh_dq(qloc,(/0.0_wp,1.0_wp,0.0_wp/))
        g_tt_dq1(i,j,k) =g_tt_dq1(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_thet,(/1.0_wp,0.0_wp,0.0_wp/))
        g_tt_dq2(i,j,k) =g_tt_dq2(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_thet,(/0.0_wp,1.0_wp,0.0_wp/))
        g_tz_dq1(i,j,k) =g_tz_dq1(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_zeta,(/1.0_wp,0.0_wp,0.0_wp/))
        g_tz_dq2(i,j,k) =g_tz_dq2(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_zeta,(/0.0_wp,1.0_wp,0.0_wp/))
        g_zz_dq1(i,j,k) =g_zz_dq1(i,j,k) - sf%eval_gij_dq(q_zeta,qloc,q_zeta,(/1.0_wp,0.0_wp,0.0_wp/))
        g_zz_dq2(i,j,k) =g_zz_dq2(i,j,k) - sf%eval_gij_dq(q_zeta,qloc,q_zeta,(/0.0_wp,1.0_wp,0.0_wp/))
        g_t1(i,j,k)     =g_t1(i,j,k)     - sf%eval_gij(q_thet,qloc,(/1.0_wp,0.0_wp,0.0_wp/))
        g_t2(i,j,k)     =g_t2(i,j,k)     - sf%eval_gij(q_thet,qloc,(/0.0_wp,1.0_wp,0.0_wp/))
        g_z1(i,j,k)     =g_z1(i,j,k)     - sf%eval_gij(q_zeta,qloc,(/1.0_wp,0.0_wp,0.0_wp/))
        g_z2(i,j,k)     =g_z2(i,j,k)     - sf%eval_gij(q_zeta,qloc,(/0.0_wp,1.0_wp,0.0_wp/))
        Gh11(i,j,k)     =Gh11(i,j,k)     - sf%eval_gij((/1.0_wp,0.0_wp,0.0_wp/),qloc,(/1.0_wp,0.0_wp,0.0_wp/))
        Gh22(i,j,k)     =Gh22(i,j,k)     - sf%eval_gij((/0.0_wp,1.0_wp,0.0_wp/),qloc,(/0.0_wp,1.0_wp,0.0_wp/))
      END DO; END DO; END DO

      iTest=201+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Jh))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_all-eval_Jh(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=202+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tt))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tt_all-eval_g_tt(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=203+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tz))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tz_all-eval_g_tz(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=203+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_zz))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_zz_all-eval_g_zz(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=204+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Jh_dq1))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_dq1_all-eval_Jh_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=205+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Jh_dq2))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_dq2_all-eval_Jh_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=206+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tt_dq1))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tt_dq1_all-eval_g_tt_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=207+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tz_dq1))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tz_dq1_all-eval_g_tz_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=208+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_zz_dq1))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_zz_dq1_all-eval_g_zz_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=209+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tt_dq2))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tt_dq2_all-eval_g_tt_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=210+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tz_dq2))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tz_dq2_all-eval_g_tz_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=211+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_zz_dq2))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_zz_dq2_all-eval_g_zz_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=212+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_t1))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_t1_all-eval_g_t1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=213+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_t2))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_t2_all-eval_g_t2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=214+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_z1))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_z1_all-eval_g_z1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=215+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_z2))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_z2_all-eval_g_z2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=216+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Gh11))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Gh11_all-eval_Gh11(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=217+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Gh22))/REAL(PRODUCT(ndims),wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Gh22_all-eval_Gh22(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      DEALLOCATE(zeta,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                 Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,&
                 Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2,g_t1,g_t2,g_z1,g_z2,Gh11,Gh22)
      DEALLOCATE(xv)
    END DO !idir
 END IF

 test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_knot_test

END MODULE MODgvec_hmap_knot
