!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** hmap_knot **
!!
!! contains the type that points to the routines of one chosen hmap_knot
!!
!===================================================================================================================================
MODULE MOD_hmap_knot
! MODULES
USE MOD_Globals, ONLY:PI,wp,Unit_stdOut,abort
USE MOD_c_hmap,    ONLY:c_hmap
IMPLICIT NONE

PUBLIC
 

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

  PROCEDURE :: init          => hmap_knot_init
  PROCEDURE :: free          => hmap_knot_free
  PROCEDURE :: eval          => hmap_knot_eval          
  PROCEDURE :: eval_dxdq     => hmap_knot_eval_dxdq
  PROCEDURE :: eval_Jh       => hmap_knot_eval_Jh       
  PROCEDURE :: eval_Jh_dq1   => hmap_knot_eval_Jh_dq1    
  PROCEDURE :: eval_Jh_dq2   => hmap_knot_eval_Jh_dq2    
  PROCEDURE :: eval_gij      => hmap_knot_eval_gij      
  PROCEDURE :: eval_gij_dq1  => hmap_knot_eval_gij_dq1  
  PROCEDURE :: eval_gij_dq2  => hmap_knot_eval_gij_dq2  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! procedures for hmap_knot:
  PROCEDURE :: Rl            => hmap_knot_eval_Rl
  PROCEDURE :: Zl            => hmap_knot_eval_Zl
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap_knot

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type hmap_knot with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_knot_init( sf )
! MODULES
USE MOD_ReadInTools, ONLY: GETINT, GETREAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                           :: knot_k, knot_l         !parameters of the (k,l)-torus
  REAL(wp)                          :: knot_R0, knot_delta    !major radius and shift
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: KNOT ON A (k,l)-TORUS ...'

  knot_k=GETINT("hmap_knot_k","2")
  sf%k=REAL(knot_k, wp)

  knot_l=GETINT("hmap_knot_l","3")
  sf%l=REAL(knot_l, wp)

  knot_R0=GETREAL("hmap_knot_major_radius","1.0")
  sf%R0=knot_R0

  knot_delta=GETREAL("hmap_knot_delta_shift","0.4")
  sf%delta=knot_delta

  IF (.NOT.((sf%R0 - ABS(sf%delta)) > 0.0_wp)) THEN
     CALL abort(__STAMP__, &
          "hmap_knot init: condition R0 - |delta| > 0 not fulfilled!")
  END IF

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_knot_test(sf)

END SUBROUTINE hmap_knot_init


!===================================================================================================================================
!> finalize the type hmap_knot
!!
!===================================================================================================================================
SUBROUTINE hmap_knot_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN

  sf%initialized=.FALSE.
  sf%R0 = 0.0_wp
  sf%delta = -1.0_wp
  sf%k = 0.0_wp
  sf%l = 0.0_wp

END SUBROUTINE hmap_knot_free


!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval( sf ,q_in) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: x_out(3)
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
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_vec(3)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(3)
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

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Jh = sf%k*sf%Rl(q_in)   ! Jh = k * Rl
END FUNCTION hmap_knot_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh_dq1( sf ,q_in) RESULT(Jh_dq1)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh_dq1
!===================================================================================================================================
  Jh_dq1 = sf%k ! dJh/dq^1 = d(kRl)/dq^1 = p, since dRl/dq^1 = 1.
END FUNCTION hmap_knot_eval_Jh_dq1


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh_dq2( sf ,q_in) RESULT(Jh_dq2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh_dq2
!===================================================================================================================================
  Jh_dq2 = 0.0_wp ! dJh/dq^2 = d(kRl)/dq^2 = 0, Rl is independent of q^2
END FUNCTION hmap_knot_eval_Jh_dq2


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: qL_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_G(3)
  REAL(wp)          , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: g_ab
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                          :: A, B, C
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
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and 
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_gij_dq1( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq1)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)           , INTENT(IN   ) :: qL_in(3)
  REAL(wp)           , INTENT(IN   ) :: q_G(3)
  REAL(wp)           , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                           :: g_ab_dq1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  !                       |q1  |   |0  0  0        |        |q1  |  
  !q_i G_ij q_j = (dalpha |q2  | ) |0  0  0        | (dbeta |q2  | )
  !                       |q3  |   |0  0  2k**2 *Rl|        |q3  |  
  g_ab_dq1 = qL_in(3) * 2.0_wp * sf%k**2 * sf%Rl(q_G) * qR_in(3)  
END FUNCTION hmap_knot_eval_gij_dq1


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and 
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_gij_dq2( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq2)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: qL_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_G(3)
  REAL(wp)          , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: g_ab_dq2
!===================================================================================================================================
  !                            |q1  |   |0  0  0  |        |q1   |  
  !q_i dG_ij/dq1 q_j = (dalpha |q2  | ) |0  0  0  | (dbeta |q1   | ) =0
  !                            |q3  |   |0  0  0  |        |q3   |  
  g_ab_dq2=0.0_wp 
END FUNCTION hmap_knot_eval_gij_dq2


!===================================================================================================================================
!> evaluate the effective major radius coordinate Rl(q) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Rl( sf ,q_in) RESULT(Rl_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Rl_out
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
FUNCTION hmap_knot_eval_Zl( sf ,q_in) RESULT(Zl_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Zl_out
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
USE MOD_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3)
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: a, Rl, Zl
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
  IF(testlevel.LE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    a = sf%R0 - ABS(sf%delta)
    q_in=(/0.5_wp*a, -0.2_wp*a, 0.5_wp*Pi/)
    Rl = sf%R0 + sf%delta*COS(sf%l*q_in(3)) + q_in(1)
    Zl = sf%delta*SIN(sf%l*q_in(3)) + q_in(2)
    x = sf%eval(q_in )
    checkreal=SUM((x-(/Rl*COS(sf%k*q_in(3)),-Rl*SIN(sf%k*q_in(3)),Zl/))**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
            '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(2(A,E11.3))') &
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
       nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
            '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST
    
 END IF !testlevel <1
 
 test_called=.FALSE. ! to prevent infinite loop in this routine
 

END SUBROUTINE hmap_knot_test

END MODULE MOD_hmap_knot

