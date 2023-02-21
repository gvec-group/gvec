!===================================================================================================================================
! Copyright (C) 2022 - 2023  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module ** hmap_frenet **
!!
!! This map uses the Frenet frame of a given periodic input curve X0(zeta) along the curve parameter zeta in [0,2pi]. 
!! It uses the signed orthonormal Frenet-Serret frame (TNB frame) that can be computed from derivatives of X0  in zeta. 
!! h:  X_0(\zeta) + q_1 \sigma N(\zeta) + q_2 \sigma B(\zeta)
!! with a sign switching function \sigma(\zeta), that accounts for points of zero curvature.
!! the tangent is T=X_0' / |X_0'|, the bi-normal is B= (X_0' x X_0'') / |X_0' x X_0''|, and the normal N= B X T
!! Derivatives use the Frenet-Serret formulas:
!!
!! dT/dl = k N
!! dN/dl = -kappa T + tau B
!! dB/dl = -tau N
!!
!! With  l(\zeta) being the arc-length, and l' = |X_0'|. 
!! the curvature is kappa=  |X_0' x  X_0''| / (l')^3, 
!! and the torsion tau= (X_0' x X_0'').X_0''' /  |X_0' x X_0''|^2
!!
!! As a first representation of the curve X0(\zeta), we choose zeta to be the geometric toroidal angle zeta=phi, such that
!!             R0(zeta)*cos(zeta)
!!  X0(zeta)=( R0(zeta)*sin(zeta)  )
!!             Z0(zeta)
!! and both R0,Z0 are represented as a real Fourier series with modes 0... n_max and number of Field periods Nfp
!! R0(zeta) = sum_{n=0}^{n_{max}} R0c(n)*cos(n*Nfp*zeta) + R0s(n)*sin(n*Nfp*zeta)
!===================================================================================================================================
MODULE MODgvec_hmap_frenet
! MODULES
USE MODgvec_Globals, ONLY:PI,wp,Unit_stdOut,abort
USE MODgvec_c_hmap,    ONLY:c_hmap
IMPLICIT NONE

PUBLIC
 

TYPE,EXTENDS(c_hmap) :: t_hmap_frenet
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_frenet:
  INTEGER              :: nfp     !! input number of field periods
  !curve description
  INTEGER              :: n_max   !! input maximum mode number (without nfp), 0...n_max, 
  REAL(wp),ALLOCATABLE :: R0c(:)  !! input cosine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  REAL(wp),ALLOCATABLE :: R0s(:)  !! input   sine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp  
  REAL(wp),ALLOCATABLE :: Z0c(:)  !! input cosine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  REAL(wp),ALLOCATABLE :: Z0s(:)  !! input   sine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  INTEGER,ALLOCATABLE  :: Xn(:)   !! array of mode numbers,  local variable =(0,1,...,n_max)*nfp 
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  PROCEDURE :: init          => hmap_frenet_init
  PROCEDURE :: free          => hmap_frenet_free
  PROCEDURE :: eval          => hmap_frenet_eval          
  PROCEDURE :: eval_dxdq     => hmap_frenet_eval_dxdq
  PROCEDURE :: eval_Jh       => hmap_frenet_eval_Jh       
  PROCEDURE :: eval_Jh_dq1   => hmap_frenet_eval_Jh_dq1    
  PROCEDURE :: eval_Jh_dq2   => hmap_frenet_eval_Jh_dq2    
  PROCEDURE :: eval_gij      => hmap_frenet_eval_gij      
  PROCEDURE :: eval_gij_dq1  => hmap_frenet_eval_gij_dq1  
  PROCEDURE :: eval_gij_dq2  => hmap_frenet_eval_gij_dq2  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! procedures for hmap_frenet:
  PROCEDURE :: eval_X0       => hmap_frenet_eval_X0_fromRZ
END TYPE t_hmap_frenet

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type hmap_frenet with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_init( sf )
! MODULES
USE MODgvec_ReadInTools, ONLY: GETINT, GETREALARRAY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: n
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: FRENET FRAME OF A CLOSED CURVE ...'

  sf%nfp=GETINT("hmap_nfp")
  sf%n_max=GETINT("hmap_n_max")
  ALLOCATE(sf%Xn(0:sf%n_max))
  DO n=0,sf%n_max
    sf%Xn(n)=n*sf%nfp
  END DO
  ALLOCATE(sf%R0c(0:sf%n_max));sf%R0c=0.
  ALLOCATE(sf%R0s(0:sf%n_max));sf%R0s=0.
  ALLOCATE(sf%Z0c(0:sf%n_max));sf%Z0c=0.
  ALLOCATE(sf%Z0s(0:sf%n_max));sf%Z0s=0.


  sf%R0c=GETREALARRAY("hmap_R0c",sf%n_max+1,sf%R0c)
  sf%R0s=GETREALARRAY("hmap_R0s",sf%n_max+1,sf%R0s)
  sf%Z0c=GETREALARRAY("hmap_Z0c",sf%n_max+1,sf%Z0c)
  sf%Z0s=GETREALARRAY("hmap_Z0s",sf%n_max+1,sf%Z0s)


  IF (.NOT.(sf%R0c(0) > 0.0_wp)) THEN
     CALL abort(__STAMP__, &
          "hmap_frenet init: condition R0c(n=0) > 0 not fulfilled!")
  END IF

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_frenet_test(sf)

END SUBROUTINE hmap_frenet_init


!===================================================================================================================================
!> finalize the type hmap_frenet
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN
  DEALLOCATE(sf%R0c)
  DEALLOCATE(sf%R0s)
  DEALLOCATE(sf%Z0c)
  DEALLOCATE(sf%Z0s)

  sf%initialized=.FALSE.

END SUBROUTINE hmap_frenet_free


!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z) 
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval( sf ,q_in) RESULT(x_out)
! MODULES
USE MODgvec_Globals, ONLY:CROSS
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: x_out(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),PARAMETER:: sigma=1.0_wp  !! placeholder, no sign change for now
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
!===================================================================================================================================
 ! q(:) = (q1,q2,zeta) are the variables in the domain of the map
 ! X(:) = (x,y,z) are the variables in the range of the map
 !
 !  |x |  
 !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
 !  |z |  

 ASSOCIATE(zeta=>q_in(3))
 CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
 T=X0p
 T=T/SQRT(SUM(T*T))
 B=CROSS(X0p,X0pp)
 B=B/SQRT(SUM(B*B))
 N=CROSS(B,T)
 x_out=X0 +sigma*(q_in(1)*N + q_in(2)*B)
 END ASSOCIATE
END FUNCTION hmap_frenet_eval

!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_vec(3)
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp
!===================================================================================================================================
 ASSOCIATE(zeta=>q_in(3))
 CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
 dxdq_qvec(1:3)= X0p !placeholder
                                                 
 END ASSOCIATE
 

END FUNCTION hmap_frenet_eval_dxdq

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Jh = 0.
END FUNCTION hmap_frenet_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh_dq1( sf ,q_in) RESULT(Jh_dq1)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh_dq1
!===================================================================================================================================
  Jh_dq1 = 0. 
END FUNCTION hmap_frenet_eval_Jh_dq1


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh_dq2( sf ,q_in) RESULT(Jh_dq2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh_dq2
!===================================================================================================================================
  Jh_dq2 = 0.0_wp 
END FUNCTION hmap_frenet_eval_Jh_dq2


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
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
  ! A = 
  ! B = 
  ! C = 
  !                       |q1  |   |1  0  A|        |q1  |  
  !q_i G_ij q_j = (dalpha |q2  | ) |0  1  B| (dbeta |q2  | )
  !                       |q3  |   |A  B  C|        |q3  |  
 ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
   A = 0. 
   B = 0.
   C = 1.
   g_ab=SUM(qL_in(:)*(/qR_in(1) + A*qR_in(3), qR_in(2) + B*qR_in(3), A*qR_in(1) + B*qR_in(2) + C*qR_in(3)/))
 END ASSOCIATE
END FUNCTION hmap_frenet_eval_gij


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and 
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij_dq1( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq1)
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
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
  !                       |q3  |   |0  0           |        |q3  |  
  g_ab_dq1 = 2.0_wp
END FUNCTION hmap_frenet_eval_gij_dq1


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and 
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij_dq2( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq2)
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
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
END FUNCTION hmap_frenet_eval_gij_dq2


!===================================================================================================================================
!> evaluate curve X0(zeta), position and first three derivatives, from given R0,Z0 Fourier 
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_eval_X0_fromRZ( sf,zeta,X0,X0p,X0pp,X0ppp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN ) :: sf
  REAL(wp)            , INTENT(IN ) :: zeta       !! position along closed curve parametrized in [0,2pi]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)            , INTENT(OUT) :: X0(1:3)      !! curve position in cartesian coordinates
  REAL(wp)            , INTENT(OUT) :: X0p(1:3)     !! 1st derivative in zeta
  REAL(wp)            , INTENT(OUT) :: X0pp(1:3)    !! 2nd derivative in zeta
  REAL(wp)            , INTENT(OUT) :: X0ppp(1:3)   !! 3rd derivative in zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: R0,R0p,R0pp,R0ppp
  REAL(wp) :: coszeta,sinzeta
!===================================================================================================================================
  CALL eval_fourier1d(sf%n_max,sf%Xn,sf%R0c,sf%R0s,zeta,R0,R0p,R0pp,R0ppp)
  CALL eval_fourier1d(sf%n_max,sf%Xn,sf%Z0c,sf%Z0s,zeta,X0(3),X0p(3),X0pp(3),X0ppp(3)) !=Z0,Z0p,Z0pp,Z0ppp
  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  ASSOCIATE(x   =>X0(1)   ,y   =>X0(2)   , &
            xp  =>X0p(1)  ,yp  =>X0p(2)  , &
            xpp =>X0pp(1) ,ypp =>X0pp(2) , &
            xppp=>X0ppp(1),yppp=>X0ppp(2))
    !! angle zeta=geometric toroidal angle phi=atan(y/x)
    x=R0*coszeta
    x=R0*sinzeta
    
    !xp = R0p*coszeta  - R0*sinzeta
    !yp = R0p*sinzeta  + R0*coszeta
    xp  = R0p*coszeta  -y
    yp  = R0p*sinzeta  +x
    
    !xpp = R0pp*coszeta - 2*R0p*sinzeta - R0*coszeta  =  ... -2*yp + x
    !ypp = R0pp*sinzeta + 2*R0p*coszeta - R0*sinzeta  = ...  +2*xp + y
    xpp  = R0pp*coszeta -2.0_wp*yp + x
    ypp  = R0pp*sinzeta +2.0_wp*xp + y
    
    !xppp = R0ppp*coszeta - 3*R0pp*sinzeta - 3*R0p*coszeta + R0*sinzeta  = ... -3*ypp +3*xp +y 
    !yppp = R0ppp*sinzeta + 3*R0pp*coszeta - 3*R0p*sinzeta - R0*coszeta  = ... +3*xpp +3*yp -x 
    xppp  = R0ppp*coszeta +3.0_wp*(xp-ypp) + y
    yppp  = R0ppp*sinzeta +3.0_wp*(yp+xpp) + x 

  END ASSOCIATE !x,y,xp,yp,...

END SUBROUTINE hmap_frenet_eval_X0_fromRZ


!===================================================================================================================================
!> evaluate 1d fourier series from given cos/sin coefficients and mode numbers xn
!! SUM(xc(0:n_max)*COS(xn(0:n_max)*zeta)+xs(0:n_max)*SIN(xn(0:n_max)*zeta)
!! evaluate all derivatives 1,2,3 alongside
!!
!===================================================================================================================================
SUBROUTINE eval_fourier1d(n_max,xn,xc,xs,zeta,x,xp,xpp,xppp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER  , INTENT(IN ) :: n_max        !! number of modes is n_max+1  (0...n_max)
  INTEGER  , INTENT(IN ) :: xn(0:n_max)  !! array of mode numbers  
  REAL(wp) , INTENT(IN ) :: xc(0:n_max)  !! cosine coefficients
  REAL(wp) , INTENT(IN ) :: xs(0:n_max)  !!   sine coefficients
  REAL(wp) , INTENT(IN ) :: zeta         !! angular position [0,2pi] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) , INTENT(OUT) :: x      !! value at zeta 
  REAL(wp) , INTENT(OUT) :: xp     !! 1st derivative in zeta
  REAL(wp) , INTENT(OUT) :: xpp    !! 2nd derivative in zeta
  REAL(wp) , INTENT(OUT) :: xppp   !! 3rd derivative in zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(0:n_max) :: cos_nzeta,sin_nzeta,xtmp,xptmp
!===================================================================================================================================
  cos_nzeta=COS(REAL(xn,wp)*zeta)
  sin_nzeta=SIN(REAL(xn,wp)*zeta)
  xtmp = xc*cos_nzeta+xs*sin_nzeta
  xptmp= REAL(xn,wp)*(-xc*sin_nzeta+xs*cos_nzeta)
  x    = SUM(xtmp)
  xp   = SUM(xptmp)
  xpp  = SUM(REAL(-xn*xn,wp)*xtmp)
  xppp = SUM(REAL(-xn*xn,wp)*xptmp)

END SUBROUTINE eval_fourier1d


!===================================================================================================================================
!> test hmap_frenet - evaluation of the map
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3)
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: a, R0, Z0
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_frenet TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    !evaluate on the axis q1=q2=0
    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.0_wp, 0.0_wp, 0.335_wp*PI/)
    R0 = SUM(sf%R0c(:)*COS(sf%Xn(:)*q_in(3)) + sf%R0s(:)*SIN(sf%Xn(:)*q_in(3)))
    Z0 = SUM(sf%Z0c(:)*COS(sf%Xn(:)*q_in(3)) + sf%Z0s(:)*SIN(sf%Xn(:)*q_in(3)))
    x = sf%eval(q_in )
    checkreal=SUM((x-(/R0*COS(q_in(3)),R0*SIN(q_in(3)),Z0/))**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    
 END IF !testlevel >=1
 
 test_called=.FALSE. ! to prevent infinite loop in this routine
 

END SUBROUTINE hmap_frenet_test

END MODULE MODgvec_hmap_frenet

