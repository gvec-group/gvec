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
  INTEGER  :: p,  q    ! this map is based on the (p,q)-torus
  REAL(wp) :: R0       ! major radius
  REAL(wp) :: delta    ! shift of the axis
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_knot:

  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  PROCEDURE :: init          => hmap_knot_init
  PROCEDURE :: free          => hmap_knot_free
  PROCEDURE :: eval          => hmap_knot_eval          
  PROCEDURE :: eval_Jh       => hmap_knot_eval_Jh       
  PROCEDURE :: eval_Jh_dq1   => hmap_knot_eval_Jh_dq1    
  PROCEDURE :: eval_Jh_dq2   => hmap_knot_eval_Jh_dq2    
  PROCEDURE :: eval_gij      => hmap_knot_eval_gij      
  PROCEDURE :: eval_gij_dq1  => hmap_knot_eval_gij_dq1  
  PROCEDURE :: eval_gij_dq2  => hmap_knot_eval_gij_dq2  
  PROCEDURE :: Rq            => hmap_knot_eval_Rq
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap_knot

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type hmap_knot with number of elements
!!
!===================================================================================================================================
!!!!SUBROUTINE hmap_knot_init( sf, R0, delta, p, q)
SUBROUTINE hmap_knot_init( sf )
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
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: KNOT ...'

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

END SUBROUTINE hmap_knot_free

!===================================================================================================================================
!> evaluate the mapping h (X^1,X^2,zeta) -> (x,y,z) 
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
  REAL(wp)                          :: Rq, Zq
!===================================================================================================================================
 ! (p,q) are the indides of the (p,q)-torus
 ! q(:) = (X1,X2,zeta) are the variables in the domain of the map
 ! X(:) = (x,y,z) are the variables in the range of the map
 !
 !   Rq = R0 + delta * cos(q*zeta) + X1
 !   Zq = delta * sin(q*zeta) + X2
 !  |x |  | Rq*sin(p*zeta) |
 !  |y |= |-Rq*cos(p*zeta) |
 !  |z |  | Zq             |

 ASSOCIATE(X1=>q_in(1),X2=>q_in(2),zeta=>q_in(3))
 Rq = sf%R0 + sf%delta*COS(sf%q*zeta) + X1
 Zq = sf%delta*SIN(sf%q*zeta) + X2
 x_out(1:3)=(/ Rq*COS(sf%p*zeta), &
              -Rq*SIN(sf%p*zeta), &
               Zq                 /)
 END ASSOCIATE
END FUNCTION hmap_knot_eval

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(X^1,X^2,zeta) 
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
  REAL(wp)                          :: Rq
!===================================================================================================================================
 ASSOCIATE(X1=>q_in(1),X2=>q_in(2),zeta=>q_in(3))
   Rq = sf%R0 + sf%delta*COS(sf%q*zeta) + X1
   Jh=sf%p*Rq
 END ASSOCIATE
END FUNCTION hmap_knot_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(X^1,X^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh_dq1( sf ,q_in) RESULT(Jh_dq1)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh_dq1
!===================================================================================================================================
  Jh_dq1 = sf%p ! dJh/dX2 = d(pRq)/dX^1 = p, since dRq/dX^1 = 1.
END FUNCTION hmap_knot_eval_Jh_dq1

!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(X^1,X^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Jh_dq2( sf ,q_in) RESULT(Jh_dq2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: Jh_dq2
!===================================================================================================================================
  Jh_dq2 = 0.0_wp ! dJh/dX2 = d(pRq)/dX^2 = 0, Rq is independent of X^2
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
  REAL(wp)        , INTENT(IN   ) :: qL_in(3)
  REAL(wp)        , INTENT(IN   ) :: q_G(3)
  REAL(wp)        , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                        :: Rq, A, B, C
!===================================================================================================================================
  ! A = - q*delta*sin(q*zeta), 
  ! B = q*delta*cos(q*zeta)
  ! C = p**2 * Rq**2 + q**2 * delta**2
  !                       |X1  |   |1  0  A|        |X1  |  
  !q_i G_ij q_j = (dalpha |X2  | ) |0  1  B| (dbeta |X2  | )
  !                       |zeta|   |A  B  C|        |zeta|  
 ASSOCIATE(X1=>q_G(1),X2=>q_G(2),zeta=>q_G(3))
   Rq = sf%R0 + sf%delta*COS(sf%q*zeta) + X1
   A = - sf%q*sf%delta*SIN(sf%q*zeta)
   B = sf%q*sf%delta*COS(sf%q*zeta)
   C = sf%p**2 * Rq**2 + sf%q**2 * sf%delta**2
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
  REAL(wp)        , INTENT(IN   ) :: qL_in(3)
  REAL(wp)        , INTENT(IN   ) :: q_G(3)
  REAL(wp)        , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab_dq1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                        :: Rq
!===================================================================================================================================
  !                       |X1  |   |0  0  0        |        |X1  |  
  !q_i G_ij q_j = (dalpha |X2  | ) |0  0  0        | (dbeta |X2  | )
  !                       |zeta|   |0  B  2p**2 *Rq|        |zeta|  
 ASSOCIATE(X1=>q_G(1),X2=>q_G(2),zeta=>q_G(3))
   Rq = sf%R0 + sf%delta*COS(sf%q*zeta) + X1
   g_ab_dq1=SUM(qL_in(:)*(/0.0_wp, 0.0_wp, 2.0_wp * sf%p**2 * Rq*qR_in(3)/))
 END ASSOCIATE
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
END FUNCTION hmap_knot_eval_gij_dq2


!===================================================================================================================================
!> evaluate the mapping h (X^1,X^2,zeta) -> (x,y,z) 
!!
!===================================================================================================================================
FUNCTION hmap_knot_eval_Rq( sf ,q_in) RESULT(Rq_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
  CLASS(t_hmap_knot), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Rq_out
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
 !   Rq = R0 + delta * cos(q*zeta) + X1

 ASSOCIATE(X1=>q_in(1),X2=>q_in(2),zeta=>q_in(3))
   Rq_out = sf%R0 + sf%delta*COS(sf%q*zeta) + X1
 END ASSOCIATE
END FUNCTION hmap_knot_eval_Rq


!===================================================================================================================================
!> test hmap_knot 
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
    q_in=(/0.1_wp,-0.2_wp,0.5_wp*Pi/)
    x = sf%eval(q_in )
!   checkreal=SUM((x-(/q_in(1)*COS(q_in(3)),-q_in(1)*SIN(q_in(3)),q_in(2)/))**2)
!    refreal  =0.0_wp
!
!    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
!      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
!      '\n!! hmap_knot TEST ID',nTestCalled ,': TEST ',iTest,Fail
!      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(2(A,E11.3))') &
!      '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
!    END IF !TEST

  END IF !testlevel <1

  test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_knot_test

END MODULE MOD_hmap_knot

