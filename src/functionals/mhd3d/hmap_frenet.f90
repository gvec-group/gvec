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
!! R0(zeta) = sum_{n=0}^{n_{max}} rc(n)*cos(n*Nfp*zeta) + rs(n)*sin(n*Nfp*zeta)
!===================================================================================================================================
MODULE MODgvec_hmap_frenet
! MODULES
USE MODgvec_Globals, ONLY:PI,TWOPI,CROSS,wp,Unit_stdOut,abort,MPIroot
USE MODgvec_c_hmap,    ONLY:c_hmap
IMPLICIT NONE

PUBLIC
 

TYPE,EXTENDS(c_hmap) :: t_hmap_frenet
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_frenet:
  !curve description
  !INTEGER             :: nfp  !already part of c_hmap. Is overwritten in init
  INTEGER              :: n_max=0  !! input maximum mode number (without nfp), 0...n_max, 
  REAL(wp),ALLOCATABLE :: rc(:)  !! input cosine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  REAL(wp),ALLOCATABLE :: rs(:)  !! input   sine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp  
  REAL(wp),ALLOCATABLE :: zc(:)  !! input cosine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  REAL(wp),ALLOCATABLE :: zs(:)  !! input   sine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  INTEGER,ALLOCATABLE  :: Xn(:)   !! array of mode numbers,  local variable =(0,1,...,n_max)*nfp 
  LOGICAL              :: omnig=.FALSE.   !! omnigenity. True: sign change of frame at pi/nfp , False: no sign change
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
  PROCEDURE :: sigma         => hmap_frenet_sigma
  ! --- Not used
  PROCEDURE :: init_aux      => dummy_sub_hmap_init_aux
  PROCEDURE :: free_aux      => dummy_sub_hmap
  PROCEDURE :: eval_aux      => dummy_sub_hmap   
END TYPE t_hmap_frenet

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS
!===============================================================================================================================
!> dummy routine that does noting
!!
SUBROUTINE dummy_sub_hmap( sf )
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
END SUBROUTINE dummy_sub_hmap

!===============================================================================================================================
!> dummy routine that does noting
!!
SUBROUTINE dummy_sub_hmap_init_aux( sf ,nzeta_aux,zeta_aux)
  INTEGER,INTENT(IN)   :: nzeta_aux
  REAL(wp),INTENT(IN)  :: zeta_aux(1:nzeta_aux)
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
END SUBROUTINE dummy_sub_hmap_init_aux

!===================================================================================================================================
!> initialize the type hmap_frenet with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_init( sf )
! MODULES
USE MODgvec_ReadInTools, ONLY: GETLOGICAL,GETINT, GETREALARRAY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: n
  INTEGER :: nvisu
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: FRENET FRAME OF A CLOSED CURVE ...'

  sf%nfp=GETINT("hmap_nfp")
  IF(sf%nfp.LE.0) &
     CALL abort(__STAMP__, &
          "hmap_frenet init: nfp > 0 not fulfilled!")

  sf%n_max=GETINT("hmap_n_max")
  ALLOCATE(sf%Xn(0:sf%n_max))
  DO n=0,sf%n_max
    sf%Xn(n)=n*sf%nfp
  END DO
  ALLOCATE(sf%rc(0:sf%n_max));sf%rc=0.
  ALLOCATE(sf%rs(0:sf%n_max));sf%rs=0.
  ALLOCATE(sf%zc(0:sf%n_max));sf%zc=0.
  ALLOCATE(sf%zs(0:sf%n_max));sf%zs=0.


  sf%rc=GETREALARRAY("hmap_rc",sf%n_max+1,sf%rc)
  sf%rs=GETREALARRAY("hmap_rs",sf%n_max+1,sf%rs)
  sf%zc=GETREALARRAY("hmap_zc",sf%n_max+1,sf%zc)
  sf%zs=GETREALARRAY("hmap_zs",sf%n_max+1,sf%zs)
  sf%omnig=GETLOGICAL("hmap_omnig",.FALSE.) !omnigenity 


  IF (.NOT.(sf%rc(0) > 0.0_wp)) THEN
     CALL abort(__STAMP__, &
          "hmap_frenet init: condition rc(n=0) > 0 not fulfilled!")
  END IF

  nvisu=GETINT("hmap_nvisu",0) 

  IF(MPIroot)THEN
    IF(nvisu.GT.0) CALL VisuFrenet(sf,nvisu)
  
    CALL CheckZeroCurvature(sf)
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
  DEALLOCATE(sf%rc)
  DEALLOCATE(sf%rs)
  DEALLOCATE(sf%zc)
  DEALLOCATE(sf%zs)

  sf%initialized=.FALSE.

END SUBROUTINE hmap_frenet_free

!===================================================================================================================================
!> Sample axis and check for zero (<1.e-12) curvature 
!!
!===================================================================================================================================
SUBROUTINE checkZeroCurvature( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iz,nz
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)              :: lp,absB
  REAL(wp),DIMENSION((sf%n_max+1)*8) :: zeta,kappa
  LOGICAL ,DIMENSION((sf%n_max+1)*8) :: checkzero
!===================================================================================================================================
  nz=(sf%n_max+1)*8
  DO iz=1,nz
    zeta(iz)=REAL(iz-1,wp)/REAL(nz,wp)*TWOPI/sf%nfp  !0...2pi/nfp without endpoint
    CALL sf%eval_X0(zeta(iz),X0,X0p,X0pp,X0ppp) 
    lp=SQRT(SUM(X0p*X0p))
    B=CROSS(X0p,X0pp)
    absB=SQRT(SUM(B*B))
    kappa(iz)=absB/(lp**3)
  END DO !iz
  checkzero=(kappa.LT.1.0e-8)
  IF(ANY(checkzero))THEN
    IF(sf%omnig)THEN
      !omnig=True: kappa can only be zero once, at 0,pi/nfp,[2pi/nfp...]
      IF(.NOT.(checkzero(1).AND.checkzero(nz/2+1).AND.(COUNT(checkzero).EQ.2)))THEN
        DO iz=1,nz  
          IF(checkzero(iz)) WRITE(UNIT_StdOut,'(A,E15.5)')'         ...curvature <1e-8 at zeta/(2pi/nfp)=',zeta(iz)*sf%nfp/TWOPI
        END DO
        CALL abort(__STAMP__, &
             "hmap_frenet checkZeroCurvature with omnig=True: found additional points with zero curvature")
      END IF
    ELSE
      DO iz=1,nz  
        IF(checkzero(iz)) WRITE(UNIT_StdOut,'(A,E15.5)')'         ...curvature <1e-8 at zeta/(2pi/nfp)=',zeta(iz)*sf%nfp/TWOPI
      END DO
      CALL abort(__STAMP__, &
           "hmap_frenet checkZeroCurvature with omnig=False: found points with zero curvature")
    END IF
  END IF
END SUBROUTINE CheckZeroCurvature

!===================================================================================================================================
!> Write evaluation of the axis and signed frenet frame to file
!!
!===================================================================================================================================
SUBROUTINE VisuFrenet( sf ,nvisu)
! MODULES
USE MODgvec_Output_CSV,     ONLY: WriteDataToCSV
USE MODgvec_Output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_Output_netcdf,     ONLY: WriteDataToNETCDF
USE MODgvec_Analyze_vars,     ONLY: outfileType
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
  INTEGER             , INTENT(IN   ) :: nvisu     !!
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)              :: zeta,sigma,lp,absB,kappa,tau,eps
  INTEGER               :: iVar,ivisu,itest
  INTEGER,PARAMETER     :: nVars=26
  CHARACTER(LEN=20)     :: VarNames(1:nVars)
  REAL(wp)              :: values(1:nVars,1:nvisu*sf%nfp+1) 
!===================================================================================================================================
  IF(nvisu.LE.0) RETURN
  iVar=0
  VarNames(ivar+1:iVar+3)=(/ "x", "y", "z"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"TX","TY","TZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"NX","NY","NZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"BX","BY","BZ"/);iVar=iVar+3
  VarNames(iVar+1       )="zeta_norm"       ;iVar=iVar+1
  VarNames(iVar+1       )="sigma_sign"      ;iVar=iVar+1
  VarNames(iVar+1       )="lprime"          ;iVar=iVar+1
  VarNames(iVar+1       )="kappa"           ;iVar=iVar+1
  VarNames(iVar+1       )="tau"             ;iVar=iVar+1
  VarNames(ivar+1:iVar+3)=(/ "X0pX", "X0pY", "X0pZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/ "X0ppX", "X0ppY", "X0ppZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/ "X0pppX", "X0pppY", "X0pppZ"/);iVar=iVar+3
  
!  values=0.
  DO ivisu=1,nvisu*sf%nfp+1
    eps=0.
    kappa=0.
    itest=0
    DO WHILE(kappa.LT.1.0e-6)! for tau being meaningful
      zeta=(REAL(ivisu-1,wp)+eps)/REAL(nvisu*sf%nfp,wp)*TWOPI
      CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
      lp=SQRT(SUM(X0p*X0p))
      T=X0p/lp
      B=CROSS(X0p,X0pp)
      absB=SQRT(SUM(B*B))
      kappa=absB/(lp**3)
      itest=itest+1
      eps=10**REAL(-16+itest,wp)
      IF(itest.EQ.15)THEN !meaningful kappa not found
        B=0.
        kappa=1.0e-6 !-12
        absB=1.
      END IF   
    END DO
    tau=SUM(X0ppp*B)/(absB**2)
    B=B/absB
    N=CROSS(B,T)
    sigma=sf%sigma(zeta)
    iVar=0
    values(ivar+1:iVar+3,ivisu)=X0                ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=T                 ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=N*sigma           ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=B*sigma           ;iVar=iVar+3
    values(iVar+1       ,ivisu)=zeta*sf%nfp/TWOPI ;iVar=iVar+1
    values(iVar+1       ,ivisu)=sigma             ;iVar=iVar+1
    values(iVar+1       ,ivisu)=lp                ;iVar=iVar+1
    values(iVar+1       ,ivisu)=kappa             ;iVar=iVar+1
    values(iVar+1       ,ivisu)=tau               ;iVar=iVar+1
    values(ivar+1:iVar+3,ivisu)=X0p               ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=X0pp              ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=X0ppp             ;iVar=iVar+3
  END DO !ivisu
  IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
    CALL WriteDataToVTK(1,3,nVars-3,(/nvisu*sf%nfp/),1,VarNames(4:nVars),values(1:3,:),values(4:nVars,:),"visu_hmap_frenet.vtu")
  END IF
  IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
#if NETCDF
    CALL WriteDataToNETCDF(1,3,nVars-3,(/nvisu*sf%nfp/),(/"dim_zeta"/),VarNames(4:nVars),values(1:3,:),values(4:nVars,:), &
         "visu_hmap_frenet")
#else
    CALL WriteDataToCSV(VarNames(:) ,values, TRIM("out_visu_hmap_frenet.csv") ,append_in=.FALSE.)
#endif
  END IF
END SUBROUTINE VisuFrenet

!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z) 
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval( sf ,q_in) RESULT(x_out)
! MODULES
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
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)          :: lp,absB,kappa,sigma,Jh
!===================================================================================================================================
  ! q(:) = (q1,q2,zeta) are the variables in the domain of the map
  ! X(:) = (x,y,z) are the variables in the range of the map
  !
  !  |x |  
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |  
 
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)
  sigma=sf%sigma(zeta)
  Jh=lp*(1.0_wp-sigma*kappa*q1)
  IF(Jh.LT.1.0e-12) &
      CALL abort(__STAMP__, &
           "hmap_frenet, evaluation outside curvature radius (sigma*q1 >= 1./(kappa))",RealInfo=zeta*sf%nfp/TWOPI)
  B=B/absB
  N=CROSS(B,T)
  x_out=X0 +sigma*(q1*N + q2*B)
  END ASSOCIATE
END FUNCTION hmap_frenet_eval

!===================================================================================================================================
!> sign function depending on zeta, 
!! if omnig=False, sigma=1
!! if omnig=True, sigma=+1 for 0<=zeta<=pi/nfp, and -1 for pi/nfp<zeta<2pi
!!
!===================================================================================================================================
FUNCTION hmap_frenet_sigma(sf,zeta) RESULT(sigma)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   )   :: zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: sigma
!===================================================================================================================================
  sigma=MERGE(SIGN(1.0_wp,SIN(sf%nfp*zeta)),1.0_wp,sf%omnig)
END FUNCTION hmap_frenet_sigma

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
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)          :: lp,absB,kappa,tau,sigma,Jh
!===================================================================================================================================
  !  |x |  
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |  
  !  dh/dq1 =sigma*N , dh/dq2=sigma*B 
  !  dh/dq3 = l' [(1-sigma*kappa*q1)T + sigma*tau*(B*q1-N*q2) ]
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)

  sigma=sf%sigma(zeta)
  Jh=lp*(1.0_wp-sigma*kappa*q1)
  IF(Jh.LT.1.0e-12) &
      CALL abort(__STAMP__, &
           "hmap_frenet, evaluation outside curvature radius (sigma*q1 >= 1/(kappa))",RealInfo=zeta*sf%nfp/TWOPI)

  tau=SUM(X0ppp*B)/(absB**2)
  B=B/absB
  N=CROSS(B,T)
  dxdq_qvec(1:3)= sigma*(N*q_vec(1)+B*q_vec(2))+(Jh*T +sigma*lp*tau*(B*q1-N*q2))*q_vec(3)
                                                  
  END ASSOCIATE !zeta
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
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)          :: lp,absB,kappa,sigma
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)
  sigma=sf%sigma(zeta)

  Jh=lp*(1.0_wp-sigma*kappa*q1)
  IF(Jh .LT. 1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet, evaluation outside curvature radius, Jh<0",RealInfo=zeta*sf%nfp/TWOPI)

  END ASSOCIATE !zeta
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
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)          :: lp,absB,kappa,sigma
!===================================================================================================================================
  !  |x |  
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |  
  ASSOCIATE(zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)

  Jh_dq1=-lp*sigma*kappa

  END ASSOCIATE !zeta
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
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)              :: lp,absB,kappa,tau,sigma
  REAL(wp)              :: Ga, Gb, Gc
!===================================================================================================================================
  ! A = -q2*l' * tau 
  ! B =  q1*l' * tau
  ! C = Jh^2 + (l'*tau)^2(q1^2+q2^2) 
  !                       |q1  |   |1   0   Ga|        |q1  |  
  !q_i G_ij q_j = (dalpha |q2  | ) |0   1   Gb| (dbeta |q2  | )
  !                       |q3  |   |Ga  Gb  Gc|        |q3  |  
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  tau=SUM(X0ppp*B)/(absB**2)
 
  Ga = -lp*tau*q2 
  Gb =  lp*tau*q1
  Gc = (lp**2)*((1.0_wp-sigma*kappa*q1)**2+tau**2*(q1**2+q2**2))
  g_ab=      qL_in(1)*qR_in(1) &
            +qL_in(2)*qR_in(2) &
       + Gc* qL_in(3)*qR_in(3) &
       + Ga*(qL_in(1)*qR_in(3)+qL_in(3)*qR_in(1)) &
       + Gb*(qL_in(2)*qR_in(3)+qL_in(3)*qR_in(2))  
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
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)              :: lp,absB,kappa,tau,sigma
!===================================================================================================================================
  !                       |q1  |   |0  0        0           |        |q1  |  
  !q_i G_ij q_j = (dalpha |q2  | ) |0  0      l'*tau        | (dbeta |q2  | )
  !                       |q3  |   |0  l'*tau  dG33/dq1     |        |q3  |  
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  tau=SUM(X0ppp*B)/(absB**2)

  g_ab_dq1 = lp*tau*(qL_in(2)*qR_in(3)+ qL_in(3)*qR_in(2)) &
            +2.0_wp*(lp**2)*((tau**2+kappa**2)*q1-sigma*kappa)*(qL_in(3)*qR_in(3))

  END ASSOCIATE
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
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)              :: lp,absB,tau
!===================================================================================================================================
  !                            |q1  |   |0       0  -l'*tau  |        |q1   |  
  !q_i dG_ij/dq1 q_j = (dalpha |q2  | ) |0       0        0  | (dbeta |q1   | ) =0
  !                            |q3  |   |-l'*tau 0   dG33/dq2|        |q3   |  
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  tau=SUM(X0ppp*B)/(absB**2)

  g_ab_dq2=-lp*tau*(qL_in(1)*qR_in(3)+qL_in(3)*qR_in(1)) + 2.0_wp*(lp*tau)**2*q2*(qL_in(3)*qR_in(3))

  END ASSOCIATE
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
  CALL eval_fourier1d(sf%n_max,sf%Xn,sf%rc,sf%rs,zeta,R0,R0p,R0pp,R0ppp)
  CALL eval_fourier1d(sf%n_max,sf%Xn,sf%zc,sf%zs,zeta,X0(3),X0p(3),X0pp(3),X0ppp(3)) !=Z0,Z0p,Z0pp,Z0ppp
  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  ASSOCIATE(x   =>X0(1)   ,y   =>X0(2)   , &
            xp  =>X0p(1)  ,yp  =>X0p(2)  , &
            xpp =>X0pp(1) ,ypp =>X0pp(2) , &
            xppp=>X0ppp(1),yppp=>X0ppp(2))
    !! angle zeta=geometric toroidal angle phi=atan(y/x)
    x=R0*coszeta
    y=R0*sinzeta
    
    xp = R0p*coszeta  - R0*sinzeta
    yp = R0p*sinzeta  + R0*coszeta
    !xp  = R0p*coszeta  -y
    !yp  = R0p*sinzeta  +x
    
    xpp = R0pp*coszeta - 2*R0p*sinzeta - R0*coszeta 
    ypp = R0pp*sinzeta + 2*R0p*coszeta - R0*sinzeta
    !xpp  = R0pp*coszeta -2.0_wp*yp + x
    !ypp  = R0pp*sinzeta +2.0_wp*xp + y
    
    xppp = R0ppp*coszeta - 3*R0pp*sinzeta - 3*R0p*coszeta + R0*sinzeta
    yppp = R0ppp*sinzeta + 3*R0pp*coszeta - 3*R0p*sinzeta - R0*coszeta
    !xppp  = R0ppp*coszeta +3.0_wp*(xp-ypp) + y
    !yppp  = R0ppp*sinzeta +3.0_wp*(yp+xpp) + x 

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
  INTEGER            :: iTest,idir,jdir,qdir
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3),q_test(3,3),x_eps(3),dxdq(3),gij,gij_eps
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  REAL(wp),PARAMETER :: epsFD=1.0e-8
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: R0, Z0
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
    R0 = SUM(sf%rc(:)*COS(sf%Xn(:)*q_in(3)) + sf%rs(:)*SIN(sf%Xn(:)*q_in(3)))
    Z0 = SUM(sf%zc(:)*COS(sf%Xn(:)*q_in(3)) + sf%zs(:)*SIN(sf%Xn(:)*q_in(3)))
    x = sf%eval(q_in )
    checkreal=SUM((x-(/R0*COS(q_in(3)),R0*SIN(q_in(3)),Z0/))**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    q_test(1,:)=(/1.0_wp, 0.0_wp, 0.0_wp/)
    q_test(2,:)=(/0.0_wp, 1.0_wp, 0.0_wp/)
    q_test(3,:)=(/0.0_wp, 0.0_wp, 1.0_wp/)
    DO qdir=1,3
      !check dx/dq^i with FD
      iTest=101+qdir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      q_in=(/0.0_wp, 0.0_wp, 0.335_wp*PI/)
      x = sf%eval(q_in )
      x_eps = sf%eval(q_in+epsFD*q_test(qdir,:))
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      checkreal=SQRT(SUM((dxdq - (x_eps-x)/epsFD)**2)/SUM(x*x))
      refreal = 0.0_wp
      
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. 100*epsFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be <',100*epsFD,' : |dxdqFD-eval_dxdq|= ', checkreal,", dq=",qdir
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
    !! TEST dG_ij_dq1 with FD 
    DO qdir=1,2
    DO idir=1,3; DO jdir=idir,3
      iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      gij  =sf%eval_gij(q_test(idir,:),q_in,q_test(jdir,:))
      gij_eps = sf%eval_gij(q_test(idir,:),q_in+epsFD*q_test(qdir,:),q_test(jdir,:))
      IF(qdir.EQ.1) refreal = sf%eval_gij_dq1(q_test(idir,:),q_in,q_test(jdir,:))
      IF(qdir.EQ.2) refreal = sf%eval_gij_dq2(q_test(idir,:),q_in,q_test(jdir,:))
      checkreal=(gij_eps-gij)/epsFD-refreal
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. 100*epsFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),3(A,I3))') &
       '\n =>  should be < ', 100*epsFD,' : |dGij_dqFD-eval_gij_dq|= ', checkreal,', i=',idir,', j=',jdir,', dq=',qdir
      END IF !TEST
    END DO; END DO
    END DO

    
      
    
 END IF !testlevel >=1
 
 test_called=.FALSE. ! to prevent infinite loop in this routine
 

END SUBROUTINE hmap_frenet_test

END MODULE MODgvec_hmap_frenet

