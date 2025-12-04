!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** hmap_cyl **
!!
!! contains the type that points to the routines of one chosen hmap_cyl
!!
!===================================================================================================================================
MODULE MODgvec_hmap_cyl
! MODULES
USE MODgvec_Globals, ONLY:PI,wp,Unit_stdOut,abort,MPIRoot
USE MODgvec_c_hmap,    ONLY:c_hmap, c_hmap_auxvar
IMPLICIT NONE

PUBLIC

TYPE,EXTENDS(c_hmap_auxvar) :: t_hmap_cyl_auxvar
  !nothing more to add for cyl hmap
END TYPE t_hmap_cyl_auxvar

TYPE,EXTENDS(c_hmap) :: t_hmap_cyl
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_cyl:
  REAL(wp) :: cyl_len
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  FINAL     :: hmap_cyl_free
  PROCEDURE :: eval_all         => hmap_cyl_eval_all
  PROCEDURE :: eval             => hmap_cyl_eval
  PROCEDURE :: eval_dxdq        => hmap_cyl_eval_dxdq
  PROCEDURE :: eval_Jh          => hmap_cyl_eval_Jh
  PROCEDURE :: eval_Jh_dq       => hmap_cyl_eval_Jh_dq
  PROCEDURE :: eval_gij         => hmap_cyl_eval_gij
  PROCEDURE :: eval_gij_dq      => hmap_cyl_eval_gij_dq
  PROCEDURE :: get_dx_dqi       => hmap_cyl_get_dx_dqi
  PROCEDURE :: get_ddx_dqij     => hmap_cyl_get_ddx_dqij

  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap_cyl

!INITIALIZATION FUNCTION:
INTERFACE t_hmap_cyl
  MODULE PROCEDURE hmap_cyl_init,hmap_cyl_init_params
END INTERFACE t_hmap_cyl

INTERFACE t_hmap_cyl_auxvar
  MODULE PROCEDURE hmap_cyl_init_aux
END INTERFACE t_hmap_cyl_auxvar

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type hmap_cyl, reading from parameterfile and then call init_params
!!
!===================================================================================================================================
FUNCTION hmap_cyl_init() RESULT(sf)
  ! MODULES
  USE MODgvec_Globals, ONLY:TWOPI
  USE MODgvec_ReadInTools, ONLY: GETREAL
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_cyl) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: cyl_len
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: CYLINDER WITH X1:=x, X2:=z, zeta := -2*pi*(y/cyl_len) , GET PARAMETERS:'

  cyl_len=GETREAL("hmap_cyl_len",TWOPI)
  sf= hmap_cyl_init_params(cyl_len)

END FUNCTION hmap_cyl_init

!===================================================================================================================================
!> initialize the type hmap_cyl, given the parameters as arguments
!!
!===================================================================================================================================
FUNCTION hmap_cyl_init_params(cyl_len) RESULT(sf)
  ! MODULES
  USE MODgvec_Globals, ONLY:TWOPI
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN) :: cyl_len !! total length of cylinder
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_cyl) :: sf !! self
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: CYLINDER WITH X1:=x, X2:=z, zeta := -2*pi*(y/cyl_len)  ...'

  sf%cyl_len=cyl_len/TWOPI

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_cyl_test(sf)

END FUNCTION hmap_cyl_init_params


!===================================================================================================================================
!> finalize the type hmap_cyl
!!
!===================================================================================================================================
SUBROUTINE hmap_cyl_free( sf )
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_cyl), INTENT(INOUT) :: sf !! self
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN
  sf%cyl_len=-1.
  sf%initialized=.FALSE.

END SUBROUTINE hmap_cyl_free

!===================================================================================================================================
!> Allocate and initialize auxiliary variable at zeta position.
!!
!===================================================================================================================================
FUNCTION hmap_cyl_init_aux( sf,zeta,do_2nd_der) RESULT(xv)
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl),INTENT(IN) :: sf
  REAL(wp)         ,INTENT(IN) :: zeta
  LOGICAL         , INTENT(IN) :: do_2nd_der !! compute second derivative and store second derivative terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_cyl_auxvar)::xv
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  xv%do_2nd_der=do_2nd_der
  xv%zeta=zeta
END FUNCTION hmap_cyl_init_aux

!===================================================================================================================================
!> evaluate all metrics necesseray for optimizer
!!
!===================================================================================================================================
SUBROUTINE hmap_cyl_eval_all(sf,ndims,dim_zeta,xv,&
                             q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                             Jh,    g_tt,    g_tz,    g_zz,&
                             Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                             Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                             g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl)   , INTENT(IN)    :: sf
  INTEGER             , INTENT(IN)    :: ndims(3)    !! 3D dimensions of input arrays
  INTEGER             , INTENT(IN)    :: dim_zeta    !! which dimension is zeta dependent
  CLASS(c_hmap_auxvar), INTENT(IN)    :: xv(ndims(dim_zeta))  !! zeta point positions
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(IN) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(OUT):: Jh,g_tt    ,g_tz    ,g_zz    , &
                                                               Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                                               Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                                               g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,j,k
  !===================================================================================================================================
  !xv not used here
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
  DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
    CALL hmap_cyl_eval_all_e(sf%cyl_len,&
             q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
             Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
             Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
             Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
             g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
  END DO; END DO; END DO
  !$OMP END PARALLEL DO

END SUBROUTINE hmap_cyl_eval_all

!===================================================================================================================================
!> evaluate all quantities at one given point (elemental)
!!
!===================================================================================================================================
PURE SUBROUTINE hmap_cyl_eval_all_e(cyl_len,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                    Jh,    g_tt,    g_tz,    g_zz,     &
                                    Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                    Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                    g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp),INTENT(IN)  :: cyl_len  !! input parameter
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
  !REAL(wp) :: Gh21,Gh31,Gh32,Gh33
!===================================================================================================================================

  Gh11=1.0_wp
  !Gh21=0.0_wp
  Gh22=1.0_wp
  !Gh31=0.0_wp
  !Gh32=0.0_wp
  !Gh33=(sf%cyl_len**2)

  Jh=cyl_len
  Jh_dq1=0.0_wp
  Jh_dq2=0.0_wp

  g_t1 = dX1_dt
  g_t2 = dX2_dt
  g_z1 = dX1_dz
  g_z2 = dX2_dz

  g_tt =   dX1_dt *  g_t1  +  dX2_dt *  g_t2
  g_tz =   dX1_dt *  g_z1  +  dX2_dt *  g_z2
  g_zz =   dX1_dz *  g_z1  +  dX2_dz *  g_z2  + (cyl_len**2)

  g_tt_dq1 = 0.0_wp
  g_tt_dq2 = 0.0_wp

  g_tz_dq1 = 0.0_wp
  g_tz_dq2 = 0.0_wp

  g_zz_dq1 = 0.0_wp
  g_zz_dq2 = 0.0_wp
END SUBROUTINE hmap_cyl_eval_all_e

!===================================================================================================================================
!> evaluate the mapping h (X^1,X^2,zeta) -> (x,y,z)
!!
!===================================================================================================================================
FUNCTION hmap_cyl_eval( sf ,q_in) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)         , INTENT(IN) :: q_in(3)
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: x_out(3)
!===================================================================================================================================
  !  q= (X1,X2,zeta), X1-> x, X2->z, zeta-> y [0,cyl_len]
  ! |x |  | X1           |
  ! |y |= |-zeta*cyl_len |
  ! |z |  | X2           |

  ASSOCIATE(xx=>q_in(1),zz=>q_in(2),zeta=>q_in(3))
  x_out(1:3)=(/ xx, &
               -zeta*sf%cyl_len , &
                zz               /)
  END ASSOCIATE
END FUNCTION hmap_cyl_eval

!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_cyl_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: q_in(3)
  REAL(wp)         , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: dxdq_qvec(3)
!===================================================================================================================================
  !  dxdq_qvec=
  ! | 1  0  0        | |q_vec(1) |
  ! | 0  0  -cyl_len | |q_vec(2) |
  ! | 0  1  0        | |q_vec(3) |

dxdq_qvec = (/q_vec(1),-sf%cyl_len*q_vec(3),q_vec(2)/)

END FUNCTION hmap_cyl_eval_dxdq


!===============================================================================================================================
!> evaluate all first derivatives dx(1:3)/dq^i, i=1,2,3 , at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_cyl_get_dx_dqi( sf ,q_in,dx_dq1,dx_dq2,dx_dq3)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: q_in(3)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp)     , INTENT(OUT) :: dx_dq1(3)
  REAL(wp)     , INTENT(OUT) :: dx_dq2(3)
  REAL(wp)     , INTENT(OUT) :: dx_dq3(3)
  !===================================================================================================================================
  dx_dq1(1:3) = (/ 1.0_wp,      0.0_wp, 0.0_wp /)
  dx_dq2(1:3) = (/ 0.0_wp,      0.0_wp, 1.0_wp /)
  dx_dq3(1:3) = (/ 0.0_wp, -sf%cyl_len, 0.0_wp /)
END SUBROUTINE hmap_cyl_get_dx_dqi

!=================================================================================================================================
!> evaluate all second derivatives d^2x(1:3)/(dq^i dq^j), i,j=1,2,3 is evaluated at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_cyl_get_ddx_dqij( sf ,q_in,ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: q_in(3)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp)     , INTENT(OUT) :: ddx_dq11(3)
  REAL(wp)     , INTENT(OUT) :: ddx_dq12(3)
  REAL(wp)     , INTENT(OUT) :: ddx_dq13(3)
  REAL(wp)     , INTENT(OUT) :: ddx_dq22(3)
  REAL(wp)     , INTENT(OUT) :: ddx_dq23(3)
  REAL(wp)     , INTENT(OUT) :: ddx_dq33(3)
  !===================================================================================================================================
  ddx_dq11(1:3) = 0.0_wp
  ddx_dq12(1:3) = 0.0_wp
  ddx_dq13(1:3) = 0.0_wp
  ddx_dq22(1:3) = 0.0_wp
  ddx_dq23(1:3) = 0.0_wp
  ddx_dq33(1:3) = 0.0_wp
END SUBROUTINE hmap_cyl_get_ddx_dqij

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(X^1,X^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_cyl_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: Jh
!===================================================================================================================================
  Jh=sf%cyl_len
END FUNCTION hmap_cyl_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(X^1,X^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_cyl_eval_Jh_dq( sf ,q_in, q_vec) RESULT(Jh_dq)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: q_in(3)
  REAL(wp)         , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: Jh_dq
!===================================================================================================================================
  Jh_dq = 0.0_wp !dJ_h / dR
END FUNCTION hmap_cyl_eval_Jh_dq


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0.0 for ds and dtheta and 1.0 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_cyl_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: qL_in(3)
  REAL(wp)         , INTENT(IN) :: q_G(3)
  REAL(wp)         , INTENT(IN) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: g_ab
!===================================================================================================================================
  g_ab=SUM(qL_in(:)*(/qR_in(1),qR_in(2),(sf%cyl_len**2)*qR_in(3)/))
END FUNCTION hmap_cyl_eval_gij


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0.0 for ds and dtheta and 1.0 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_cyl_eval_gij_dq( sf ,qL_in,q_G,qR_in,q_vec) RESULT(g_ab_dq)
  CLASS(t_hmap_cyl), INTENT(IN) :: sf
  REAL(wp)         , INTENT(IN) :: qL_in(3)
  REAL(wp)         , INTENT(IN) :: q_G(3)
  REAL(wp)         , INTENT(IN) :: qR_in(3)
  REAL(wp)         , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: g_ab_dq
!===================================================================================================================================
  g_ab_dq=0.0_wp
END FUNCTION hmap_cyl_eval_gij_dq


!===================================================================================================================================
!> test hmap_cyl
!!
!===================================================================================================================================
SUBROUTINE hmap_cyl_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cyl), INTENT(INOUT) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,idir,jdir,kdir,qdir,izeta,i,j,k,ndims(1:3),ijk(3)
  INTEGER,PARAMETER  :: nzeta=5
  INTEGER,PARAMETER  :: ns=2
  INTEGER,PARAMETER  :: nthet=3
  REAL(wp),ALLOCATABLE  :: zeta(:)
  REAL(wp)           :: qloc(3),q_thet(3),q_zeta(3)
  REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                     Jh,g_tt,    g_tz,    g_zz,     &
                                     Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                     Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                     g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3)
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
  TYPE(t_hmap_cyl_auxvar),ALLOCATABLE :: xv(:)
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_cyl TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.1_wp,-0.2_wp,0.5_wp*Pi/)
    x = sf%eval(q_in )
    checkreal=SUM((x-(/q_in(1),-q_in(3)*sf%cyl_len,q_in(2)/))**2)
    refreal  =0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
      '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.3_wp, 0.1_wp,0.4_wp*Pi/)
    x = sf%eval_dxdq(q_in, (/1.1_wp,1.2_wp,1.3_wp/) )
    checkreal=SUM((x-(/1.1_wp,-sf%cyl_len*1.3_wp,1.2_wp/))**2)
    refreal  =0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
      '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST
  END IF !testlevel>=1
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
      ndims(idir)=nzeta
      ndims(jdir)=ns
      ndims(kdir)=nthet
      ALLOCATE(zeta(ndims(idir)),xv(ndims(idir)))
      DO izeta=1,ndims(idir)
        zeta(izeta)=0.333_wp+REAL(izeta-1,wp)/REAL(ndims(idir)-1,wp)*0.221_wp
        xv(izeta)=hmap_cyl_init_aux(sf,zeta(izeta),.TRUE.)
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
      checkreal=SUM(ABS(Jh))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_all-eval_Jh(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=202+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tt))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tt_all-eval_g_tt(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=203+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tz))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tz_all-eval_g_tz(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=203+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_zz))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_zz_all-eval_g_zz(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=204+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Jh_dq1))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_dq1_all-eval_Jh_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=205+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Jh_dq2))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_dq2_all-eval_Jh_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=206+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tt_dq1))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tt_dq1_all-eval_g_tt_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=207+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tz_dq1))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tz_dq1_all-eval_g_tz_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=208+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_zz_dq1))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_zz_dq1_all-eval_g_zz_dq1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=209+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tt_dq2))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tt_dq2_all-eval_g_tt_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=210+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_tz_dq2))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_tz_dq2_all-eval_g_tz_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=211+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_zz_dq2))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_zz_dq2_all-eval_g_zz_dq2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=212+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_t1))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_t1_all-eval_g_t1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=213+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_t2))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_t2_all-eval_g_t2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=214+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_z1))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_z1_all-eval_g_z1(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=215+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(g_z2))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|g_z2_all-eval_g_z2(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=216+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Gh11))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Gh11_all-eval_Gh11(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      iTest=217+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Gh22))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_cyl TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Gh22_all-eval_Gh22(xall)|)|= ', checkreal, " ,idir=",idir
      END IF

      DEALLOCATE(zeta,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz,Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                 Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2,g_t1,g_t2,g_z1,g_z2,Gh11,Gh22)
      DEALLOCATE(xv)
    END DO !idir
 END IF

  test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_cyl_test

END MODULE MODgvec_hmap_cyl
