!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** hmap_RphiZ **
!!
!! contains the type that points to the routines of one chosen hmap_RphiZ
!!
!===================================================================================================================================
MODULE MODgvec_hmap_RphiZ
! MODULES
USE MODgvec_Globals, ONLY:PI,wp,Unit_stdOut,abort,MPIRoot
USE MODgvec_c_hmap,    ONLY:c_hmap
IMPLICIT NONE

PUBLIC


TYPE,EXTENDS(c_hmap) :: t_hmap_RphiZ
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_RphiZ:

  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  PROCEDURE :: init          => hmap_RphiZ_init
  PROCEDURE :: free          => hmap_RphiZ_free
  PROCEDURE :: eval          => hmap_RphiZ_eval
  PROCEDURE :: eval_dxdq     => hmap_RphiZ_eval_dxdq
  PROCEDURE :: eval_Jh       => hmap_RphiZ_eval_Jh
  PROCEDURE :: eval_Jh_dq1   => hmap_RphiZ_eval_Jh_dq1
  PROCEDURE :: eval_Jh_dq2   => hmap_RphiZ_eval_Jh_dq2
  PROCEDURE :: eval_gij      => hmap_RphiZ_eval_gij
  PROCEDURE :: eval_gij_dq1  => hmap_RphiZ_eval_gij_dq1
  PROCEDURE :: eval_gij_dq2  => hmap_RphiZ_eval_gij_dq2
  ! --- Not used
  PROCEDURE :: init_aux      => dummy_sub_hmap_init_aux
  PROCEDURE :: free_aux      => dummy_sub_hmap
  PROCEDURE :: eval_aux      => dummy_sub_hmap
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap_RphiZ

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS
!===============================================================================================================================
!> dummy routine that does noting
!!
SUBROUTINE dummy_sub_hmap( sf )
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
END SUBROUTINE dummy_sub_hmap

!===============================================================================================================================
!> dummy routine that does noting
!!
SUBROUTINE dummy_sub_hmap_init_aux( sf ,nzeta_aux,zeta_aux)
  INTEGER,INTENT(IN)   :: nzeta_aux
  REAL(wp),INTENT(IN)  :: zeta_aux(1:nzeta_aux)
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
END SUBROUTINE dummy_sub_hmap_init_aux

!===================================================================================================================================
!> initialize the type hmap_RphiZ with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_RphiZ_init( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: (R,phi,Z) TORUS WITH X1:=R, X2:=Z, -zeta := toroidal angle  ...'

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_RphiZ_test(sf)

END SUBROUTINE hmap_RphiZ_init


!===================================================================================================================================
!> finalize the type hmap_RphiZ
!!
!===================================================================================================================================
SUBROUTINE hmap_RphiZ_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN

  sf%initialized=.FALSE.

END SUBROUTINE hmap_RphiZ_free

!===================================================================================================================================
!> evaluate the mapping h (X^1,X^2,zeta) -> (x,y,z) cartesian
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval( sf ,q_in) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: x_out(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  !  q= (R,Z,zeta)
  ! |x |  | R*cos(zeta) |
  ! |y |= | R*sin(zeta) |
  ! |z |  | Z           |

  ASSOCIATE(R=>q_in(1),Z=>q_in(2),zeta=>q_in(3))
  x_out(1:3)=(/ R*COS(zeta), &
                R*SIN(zeta), &
                Z           /)
  END ASSOCIATE
END FUNCTION hmap_RphiZ_eval

!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
  REAL(wp)        , INTENT(IN   ) :: q_vec(3)
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: coszeta,sinzeta
!===================================================================================================================================
  !  dxdq_qvec=
  ! |  cos(zeta)  0  -q^1 sin(zeta) | |q_vec(1) |
  ! |  sin(zeta)  0  +q^1 cos(zeta) | |q_vec(2) |
  ! |     0       1        0        | |q_vec(3) |

sinzeta=SIN(q_in(3))
coszeta=COS(q_in(3))
dxdq_qvec(1:3) = (/ q_vec(1)*coszeta-q_vec(3)*q_in(1)*sinzeta, &
                    q_vec(1)*sinzeta+q_vec(3)*q_in(1)*coszeta, &
                    q_vec(2) /)


END FUNCTION hmap_RphiZ_eval_dxdq

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(X^1,X^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh
!===================================================================================================================================
  !  q= (R,Z,zeta)
  Jh=-q_in(1)
END FUNCTION hmap_RphiZ_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(X^1,X^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_Jh_dq1( sf ,q_in) RESULT(Jh_dq1)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh_dq1
!===================================================================================================================================
  !  q= (R,Z,zeta)
  Jh_dq1 =-1.0_wp !dJ_h / dR
END FUNCTION hmap_RphiZ_eval_Jh_dq1

!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(X^1,X^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_Jh_dq2( sf ,q_in) RESULT(Jh_dq2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh_dq2
!===================================================================================================================================
  !  q= (R,Z,zeta)
  Jh_dq2 = 0.0_wp !dJ_h / dZ
END FUNCTION hmap_RphiZ_eval_Jh_dq2


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: qL_in(3)
  REAL(wp)        , INTENT(IN   ) :: q_G(3)
  REAL(wp)        , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab
!===================================================================================================================================
  !                       |R   |   |1  0  0   |        |R   |
  !q_i G_ij q_j = (dalpha |Z   | ) |0  1  0   | (dbeta |Z   | )
  !                       |zeta|   |0  0  R^2 |        |zeta|
  g_ab=SUM(qL_in(:)*(/qR_in(1),qR_in(2),q_G(1)**2*(-qR_in(3))/))
END FUNCTION hmap_RphiZ_eval_gij


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_gij_dq1( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq1)
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: qL_in(3)
  REAL(wp)        , INTENT(IN   ) :: q_G(3)
  REAL(wp)        , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab_dq1
!===================================================================================================================================
  !                            |R   |   |0  0  0   |        |R   |
  !q_i dG_ij/dq1 q_j = (dalpha |Z   | ) |0  0  0   | (dbeta |Z   | )
  !                            |zeta|   |0  0  2*R |        |zeta|
  g_ab_dq1=(-qL_in(3))*2.0_wp*q_G(1)*(-qR_in(3))
END FUNCTION hmap_RphiZ_eval_gij_dq1


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_RphiZ_eval_gij_dq2( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq2)
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: qL_in(3)
  REAL(wp)        , INTENT(IN   ) :: q_G(3)
  REAL(wp)        , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab_dq2
!===================================================================================================================================
  !                            |R   |   |0  0  0  |        |R   |
  !q_i dG_ij/dq1 q_j = (dalpha |Z   | ) |0  0  0  | (dbeta |Z   | ) =0
  !                            |zeta|   |0  0  0  |        |zeta|
  g_ab_dq2=0.0_wp
END FUNCTION hmap_RphiZ_eval_gij_dq2


!===================================================================================================================================
!> test hmap_RphiZ
!!
!===================================================================================================================================
SUBROUTINE hmap_RphiZ_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_RphiZ), INTENT(INOUT) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3)
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_RphiZ TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.1_wp,-0.2_wp,0.5_wp*Pi/)
    x = sf%eval(q_in )
    checkreal=SUM((x-(/q_in(1)*COS(q_in(3)),-q_in(1)*SIN(q_in(3)),q_in(2)/))**2)
    refreal  =0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! hmap_RphiZ TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
      '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.3_wp, 0.1_wp,0.4_wp*Pi/)
    x = sf%eval_dxdq(q_in, (/1.1_wp,1.2_wp,1.3_wp/) )
    checkreal=SUM((x-(/ 1.1_wp*COS(q_in(3))-1.3_wp*q_in(1)*SIN(q_in(3)), &
                       -1.1_wp*SIN(q_in(3))-1.3_wp*q_in(1)*COS(q_in(3)),1.2_wp/))**2)
    refreal  =0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! hmap_RphiZ TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
      '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST
  END IF !testlevel >=1

  test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_RphiZ_test

END MODULE MODgvec_hmap_RphiZ
