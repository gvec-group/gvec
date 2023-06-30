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
!!# Module ** hmap_axisNB **
!!
!! This map uses a given periodic input curve X0(zeta) along the curve parameter zeta in [0,2pi]. 
!! It is very similar to a Frenet-Serret frame (TNB frame), but with normal N and binormal B vectors as input functions as well.
!! h:  X_0(\zeta) + q_1 N(\zeta) + q_2 B(\zeta)
!! the tangent is T=X_0', and N and B are assumed to be continous along zeta, so no switches.
!! Note that since N, B are input functions, they are not assumed to be unit length nor
!! orthogonal, but together with the tangent of the curve T = X_0'  , (T, N, B) should form a
!! linearly independent set of basis vectors, with T.(B x N) > 0.
!!
!! As a representation of the curve X0(\zeta), zeta is the curve parametrization in [0,2pi]
!! and the 3 cartesian coordinates of X0,N,B are given at a set of zeta points over one full turn, with nzeta*Nfp number of points.
!! these will then be fourier transformed to compute derivatives
!===================================================================================================================================
MODULE MODgvec_hmap_axisNB
! MODULES
USE MODgvec_Globals, ONLY:PI,TWOPI,CROSS,wp,Unit_stdOut,abort
USE MODgvec_c_hmap,    ONLY:c_hmap
IMPLICIT NONE

PUBLIC
 

TYPE,EXTENDS(c_hmap) :: t_hmap_axisNB
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_axisNB:
  INTEGER              :: nfp   !! input number of field periods
  !curve description
  !INTEGER              :: n_max=0  !! input maximum mode number (without nfp), 0...n_max, 
  REAL(wp),ALLOCATABLE :: rc(:)  !! input cosine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  REAL(wp),ALLOCATABLE :: rs(:)  !! input   sine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp  
  REAL(wp),ALLOCATABLE :: zc(:)  !! input cosine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  REAL(wp),ALLOCATABLE :: zs(:)  !! input   sine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp 
  INTEGER,ALLOCATABLE  :: Xn(:)   !! array of mode numbers,  local variable =(0,1,...,n_max)*nfp 
  LOGICAL              :: omnig=.FALSE.   !! omnigenity. True: sign change of frame at pi/nfp , False: no sign change
  !curve description
  INTEGER              :: nzeta=0       !! number of points in zeta direction of the input axis 
  INTEGER              :: n_max=0       !! maximum number of fourier coefficients
  REAL(wp),ALLOCATABLE :: xyz(:,:)      !! cartesian coordinates of the axis for a full turn, (1:3,1:NFP*nzeta), zeta is on 'half' grid: zeta(i)=(i-0.5)*(2pi)/(NFP*nzeta)
  REAL(wp),ALLOCATABLE :: Nxyz(:,:)     !! "normal" vector of axis frame in cartesian coordinates for a full turn (1:3,1:NFP*nzeta). NOT ASSUMED TO BE ORTHOGONAL to tangent of curve
  REAL(wp),ALLOCATABLE :: Bxyz(:,:)      !! "Bi-normal" vector of axis frame in cartesian coordinates for a full turn (1:3,1:NFP*nzeta). NOT ASSUMED TO BE ORTHOGONAL to tangent of curve or Nxyz
  CHARACTER(LEN=100)   :: axis_ncfile=" " !! name of netcdf file with axis information
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp)             :: rot_origin(1:3)=(/0.,0.,0./)   !! origin of rotation, needed for checking field periodicity
  REAL(wp)             :: rot_axis(1:3)=(/0.,0.,1./)   !! rotation axis (unit length), needed for checking field periodicity
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  PROCEDURE :: init          => hmap_axisNB_init
  PROCEDURE :: free          => hmap_axisNB_free
  PROCEDURE :: eval          => hmap_axisNB_eval          
  PROCEDURE :: eval_dxdq     => hmap_axisNB_eval_dxdq
  PROCEDURE :: eval_Jh       => hmap_axisNB_eval_Jh       
  PROCEDURE :: eval_Jh_dq1   => hmap_axisNB_eval_Jh_dq1    
  PROCEDURE :: eval_Jh_dq2   => hmap_axisNB_eval_Jh_dq2    
  PROCEDURE :: eval_gij      => hmap_axisNB_eval_gij      
  PROCEDURE :: eval_gij_dq1  => hmap_axisNB_eval_gij_dq1  
  PROCEDURE :: eval_gij_dq2  => hmap_axisNB_eval_gij_dq2  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! procedures for hmap_axisNB:
  PROCEDURE :: eval_X0      => hmap_axisNB_eval_X0_fromRZ
  PROCEDURE :: eval_TNB     => hmap_axisNB_eval_TNB_from_frenet
END TYPE t_hmap_axisNB

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

SUBROUTINE init_dummy( sf )
IMPLICIT NONE
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
  CALL abort(__STAMP__, &
             "dummy init in hmap_axisNB should not be used")
END SUBROUTINE init_dummy

!===================================================================================================================================
!> initialize the type hmap_axisNB with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_init( sf )
! MODULES
USE MODgvec_ReadInTools, ONLY: GETLOGICAL,GETINT, GETREALARRAY,GETSTR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: n
  INTEGER :: nvisu,error_nfp
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: axisNB FRAME OF A CLOSED CURVE ...'

  !sf%nfp=GETINT("hmap_nfp") !<= already set in hmap_new, before init!
  WRITE(UNIT_stdOut,*)'nfp in hmap is ', sf%nfp

  sf%n_max=GETINT("hmap_n_max")
  ALLOCATE(sf%Xn(0:sf%n_max))
  DO n=0,sf%n_max
    sf%Xn(n)=n*sf%nfp
  END DO
  ALLOCATE(sf%rc(0:sf%n_max));sf%rc=0.
  ALLOCATE(sf%rs(0:sf%n_max));sf%rs=0.
  ALLOCATE(sf%zc(0:sf%n_max));sf%zc=0.
  ALLOCATE(sf%zs(0:sf%n_max));sf%zs=0.
  sf%omnig=GETLOGICAL("hmap_omnig",.FALSE.) !omnigenity 


  sf%rc=GETREALARRAY("hmap_rc",sf%n_max+1,sf%rc)
  sf%rs=GETREALARRAY("hmap_rs",sf%n_max+1,sf%rs)
  sf%zc=GETREALARRAY("hmap_zc",sf%n_max+1,sf%zc)
  sf%zs=GETREALARRAY("hmap_zs",sf%n_max+1,sf%zs)

#if NETCDF
    ! read axis from netcdf
    sf%axis_ncfile=GETSTR("hmap_ncfile")
    CALL ReadAxis_NETCDF(sf,sf%axis_ncfile)
#else
    CALL abort(__STAMP__,&
        "cannot read axis netcdf file, since code is compiled with BUILD_NETCDF=OFF")
#endif
  


  IF (.NOT.(sf%rc(0) > 0.0_wp)) THEN
     CALL abort(__STAMP__, &
          "hmap_axisNB init: condition rc(n=0) > 0 not fulfilled!")
  END IF

  nvisu=GETINT("hmap_nvisu",0) 

  IF(nvisu.GT.0) CALL Visu_axisNB(sf,nvisu)
  
  CALL CheckFieldPeriodicity(sf,2*(sf%n_max)+1,error_nfp)
  IF(error_nfp.LT.0) &
     CALL abort(__STAMP__, &
          "hmap_axisNB check Field Periodicity failed!")
  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_axisNB_test(sf)

END SUBROUTINE hmap_axisNB_init

!===================================================================================================================================
!> finalize the type hmap_axisNB
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN
  DEALLOCATE(sf%rc)
  DEALLOCATE(sf%rs)
  DEALLOCATE(sf%zc)
  DEALLOCATE(sf%zs)

  sf%initialized=.FALSE.

END SUBROUTINE hmap_axisNB_free

#if NETCDF
!===================================================================================================================================
!> READ axis from netcdf file, needs netcdf library!
!! 
!===================================================================================================================================
SUBROUTINE ReadAxis_NETCDF(sf,fileName)
  IMPLICIT NONE
  INCLUDE "netcdf.inc"
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
  CHARACTER(LEN = *), INTENT(IN) :: fileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: ioError, ncid, id
!===================================================================================================================================
!======= HEADER ===================================================================================
!--------- GENERAL -------------
!* NFP: number of field periods
!* m_max: maximum mode number in theta
!* n_max: maximum mode number in zeta (in one field period)
!* lasym: asymmetry, logical. 
!          if lasym=0, boundary surface position X,Y in the N-B plane of the axis frame can be represented only with
!            X(theta,zeta)=sum X_mn*cos(m*theta-n*NFP*zeta), with {m=0,n=0...n_max},{m=1...m_max,n=-n_max...n_max}
!            Y(theta,zeta)=sum Y_mn*sin(m*theta-n*NFP*zeta), with {m=0,n=1...n_max},{m=1...m_max,n=-n_max...n_max}
!          if lasym=1, full fourier series is taken for X,Y

! ntheta: number of points in theta (>=2*m_max+1)
! nzeta: number of points along the axis, in one field period (>=2*n_max+1)
!
!---------- GEOMETRY ---------------
!* theta(:): theta positions, 1D array of size ntheta, must exclude the end point (last point <  (first point + 2pi), preferred on "half grid" points
!* zeta(:) :  zeta positions, 1D array of size nzeta,  must exclude the end point (last point < (first point + 2pi/NFP)), preferred on "half grid" points
!
!definition of the axis-following frame in cartesian coordinates ( boundary surface at rho=1):
!
!    {x,y,z}(rho,theta,zeta)=axis_{x,y,z}(zeta) + X(rho,theta,zeta)*N_{x,y,z}(zeta)+Y(rho,theta,zeta)*B_{x,y,z}(zeta)  
!
!* axis_xyz(::) : cartesian positions along the axis for ONE FULL TURN, 2D array of size (3,NFP*nzeta), sampled at zeta positions, must exclude the endpoint
!                 axis_xyz[i,j]=axis_i(zeta[j]),        
!* axis_Nxyz(::): cartesian components of the normal vector of the axis frame, 2D array of size (3, NFP*nzeta), evaluated analogously to the axis
!* axis_Bxyz(::): cartesian components of the bi-normal vector of the axis frame, 2D array of size (3, NFP*nzeta), evaluated analogously to the axis

!* X(::),Y(::): boundary position X,Y in the N-B plane of the axis frame, in one field period, 2D array of size(ntheta, nzeta),  with
!                  X[i, j]=X(theta[i],zeta[j])
!                  Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1
!---- PLASMA PARAMETERS:
! ....
!======= END HEADER,START DATA ===================================================================================


  SWRITE(UNIT_stdOut,'(4X,A)')'READ AXIS FILE "'//TRIM(sf%axis_ncfile)//'" in NETCDF format ...'



  !! open NetCDF input file
  ioError = NF_OPEN(TRIM(fileName), NF_NOWRITE, ncid)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
        " Cannot open "//TRIM(fileName)//" in ReadAxis_NETCDF!")

  !ioError = NF_INQ_VARID(ncid, "nfp", sf%nfp)  !OVERWRITE!
  ioError = NF_INQ_VARID(ncid, "nzeta", sf%nzeta)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                          'NETCDF READIN: problem reading nzeta' )

WRITE(*,*)'DEBUG, nzeta=',sf%nzeta
  ioError = NF_CLOSE(ncid)
  SWRITE(*,'(4X,A)')'...DONE.'
END SUBROUTINE ReadAxis_NETCDF
#endif /*NETCDF*/

!===================================================================================================================================
!> Check that the TNB frame  really has the field periodicity of NFP: 
!> assumption for now is that the origin is fixed at rot_origin=(/0.,0.,0./)
!> and the rotation axis is fixed at rot_axis=(/0.,0.,1./)
!> the sign of rotation is being checked as well
!!
!===================================================================================================================================
SUBROUTINE CheckFieldPeriodicity( sf ,nzeta,error_nfp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  INTEGER, INTENT(IN)  :: nzeta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(out)  :: error_nfp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)              :: zeta(nzeta),dzeta_fp
  REAL(wp),DIMENSION(3) :: X0_a,T,N_a,B_a,Np,Bp
  REAL(wp),DIMENSION(3) :: X0_b,N_b,B_b,P_a,P_b
  INTEGER               :: ifp,iz,error_x,error_n,error_b
!===================================================================================================================================
  error_nfp=0
  DO iz=1,nzeta
    zeta(iz)=(REAL(iz-0.5,wp))/REAL(nzeta*sf%nfp,wp)*TWOPI
  END DO

  dzeta_fp=TWOPI/REAL(sf%nfp,wp)
  CALL sf%eval_TNB(zeta(1),X0_a,T,N_a,B_a,Np,Bp) 
  P_a=X0_a-sf%rot_origin
  CALL sf%eval_TNB(zeta(1)+dzeta_fp,X0_b,T,N_b,B_b,Np,Bp) 
  P_b=X0_b-sf%rot_origin
  !account for possible sign change
  IF(SUM((P_b -  rodrigues(P_a,sf%rot_axis,dzeta_fp))**2).LT.1.0e-12)THEN
    dzeta_fp=dzeta_fp
  ELSEIF(SUM((P_b -  rodrigues(P_a,sf%rot_axis,-dzeta_fp))**2).LT.1.0e-12)THEN
    dzeta_fp=-dzeta_fp
  ELSE
    WRITE(UNIT_stdOut,*)'problem in CheckFieldPeriodicity: axis point at zeta=0 does not rotate to point at zeta= +/-2pi/nfp.'
    error_nfp=-1
    RETURN
  END IF
  DO ifp=0,sf%nfp
    error_x=0
    error_n=0
    error_b=0
    DO iz=1,nzeta
      CALL sf%eval_TNB(zeta(iz)+ ifp*dzeta_fp   ,X0_a,T,N_a,B_a,Np,Bp) 
      CALL sf%eval_TNB(zeta(iz)+(ifp+1)*dzeta_fp,X0_b,T,N_b,B_b,Np,Bp) 
      P_a=X0_a-sf%rot_origin
      P_b=X0_b-sf%rot_origin
      IF(.NOT.(SUM((P_b -  rodrigues(P_a,sf%rot_axis,dzeta_fp))**2).LT.1.0e-12))THEN
        error_x=error_x+1
      END IF
      P_a=X0_a-sf%rot_origin+N_a
      P_b=X0_b-sf%rot_origin+N_b
      IF(.NOT.(SUM((P_b -  rodrigues(P_a,sf%rot_axis,dzeta_fp))**2).LT.1.0e-12))THEN
        error_n=error_n+1
      END IF
      P_a=X0_a-sf%rot_origin+B_a
      P_b=X0_b-sf%rot_origin+B_b
      IF(.NOT.(SUM((P_b -  rodrigues(P_a,sf%rot_axis,dzeta_fp))**2).LT.1.0e-12))THEN
        error_b=error_b+1
      END IF
    END DO !iz
    IF(error_x.NE.0)THEN
      WRITE(UNIT_stdOut,*)'problem in CheckFieldPeriodicity: at least one pair of two axis points do not match at NFP' 
      error_nfp=-10
      IF(error_n.NE.0)THEN
        WRITE(UNIT_stdOut,*)'problem in CheckFieldPeriodicity: at least one pair of two normal vectors do not match at NFP' 
        error_nfp=-20
      ELSEIF(error_b.NE.0)THEN
        WRITE(UNIT_stdOut,*)'problem in CheckFieldPeriodicity: at least one pair of two bi-normal vectors do not match at NFP' 
        error_nfp=-30
      END IF
    END IF
  END DO !ifp
        
  CONTAINS 
  ! Rodrigues' rotation formula
  FUNCTION rodrigues(vec,axis,angle) RESULT(vec_rot)
     REAL(wp) :: vec(3),axis(3),angle
     REAL(wp) :: vec_rot(3)
     vec_rot = vec*COS(angle)+CROSS(axis,vec)*SIN(angle)+axis*SUM(axis*vec)*(1.0_wp-COS(angle)) 
  END FUNCTION rodrigues
END SUBROUTINE CheckFieldPeriodicity
  

!===================================================================================================================================
!> Write evaluation of the axis and signed axisNB frame to file
!!
!===================================================================================================================================
SUBROUTINE Visu_axisNB( sf ,nvisu)
! MODULES
USE MODgvec_Output_CSV,     ONLY: WriteDataToCSV
USE MODgvec_Output_vtk,     ONLY: WriteDataToVTK
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  INTEGER             , INTENT(IN   ) :: nvisu     !!
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp
  REAL(wp)              :: zeta,lp
  INTEGER               :: iVar,ivisu
  INTEGER,PARAMETER     :: nVars=14
  CHARACTER(LEN=20)     :: VarNames(1:nVars)
  REAL(wp)              :: values(1:nVars,1:nvisu*sf%nfp+1) 
!===================================================================================================================================
  IF(nvisu.LE.0) RETURN
  iVar=0
  VarNames(ivar+1:iVar+3)=(/ "x", "y", "z"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"TX","TY","TZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"NX","NY","NZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"BX","BY","BZ"/);iVar=iVar+3
  VarNames(iVar+1       )="zeta/(2pi/nfp)"  ;iVar=iVar+1
  VarNames(iVar+1       )="lprime"          ;iVar=iVar+1
  
!  values=0.
  DO ivisu=1,nvisu*sf%nfp+1
    zeta=(REAL(ivisu-0.5,wp))/REAL(nvisu*sf%nfp,wp)*TWOPI
    CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
    lp=SQRT(SUM(T*T))
    iVar=0
    values(ivar+1:iVar+3,ivisu)=X0                ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=T                 ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=N                 ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=B                 ;iVar=iVar+3
    values(iVar+1       ,ivisu)=zeta*sf%nfp/TWOPI ;iVar=iVar+1
    values(iVar+1       ,ivisu)=lp                ;iVar=iVar+1
  END DO !ivisu
  CALL WriteDataToCSV(VarNames(:) ,values, TRIM("out_visu_hmap_axisNB.csv") ,append_in=.FALSE.)
  CALL WriteDataToVTK(1,3,nVars-3,(/nvisu*sf%nfp/),1,VarNames(4:nVars),values(1:3,:),values(4:nVars,:),"visu_hmap_axisNB.vtu")
END SUBROUTINE Visu_axisNB

!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z) 
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval( sf ,q_in) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: x_out(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp
!===================================================================================================================================
  ! q(:) = (q1,q2,zeta) are the variables in the domain of the map
  ! X(:) = (x,y,z) are the variables in the range of the map
  !
  !  |x |  
  !  |y |=  X0(zeta) + (N(zeta)*q1 + B(zeta)*q2)
  !  |z |  
 
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  x_out=X0 +(q1*N + q2*B)
  END ASSOCIATE
END FUNCTION hmap_axisNB_eval


!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_vec(3)
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp
!===================================================================================================================================
  !  |x |  
  !  |y |=  X0(zeta) + (N(zeta)*q1 + B(zeta)*q2)
  !  |z |  
  !  dh/dq1 =N , dh/dq2=B 
  !  dh/dq3 = t + q1 N' + q2 * B'
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  dxdq_qvec(1:3)= N(:)*q_vec(1)+B(:)*q_vec(2)+(T(:)+q1*Np(:)+q2*Bp(:))*q_vec(3)
                                                  
  END ASSOCIATE !zeta
END FUNCTION hmap_axisNB_eval_dxdq

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  REAL(wp)        , INTENT(IN   )   :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp,Tq
  REAL(wp)              :: TqTq,NN,BB,NTq,BTq,NB,Jh2 
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  Tq=(T+q1*Np+q2*Bp)
  TqTq=SUM(Tq(:)*Tq(:))
  NTq =SUM(N(:)*Tq(:))
  BTq =SUM(B(:)*Tq(:))
  NN  =SUM(N(:)*N(:))
  BB  =SUM(B(:)*B(:))
  NB  =SUM(N(:)*B(:))
  Jh2=TqTq*(NN*BB-NB*NB)  +2.0_wp*NB*Btq*Ntq - NTq*NTq*BB - BTq*BTq*NN   
  IF(Jh2 .LT. 1.0e-4) &
      CALL abort(__STAMP__, &
           "hmap_axisNB,  Jh<0",RealInfo=zeta*sf%nfp/TWOPI)
  Jh=SQRT(Jh2)
  END ASSOCIATE !zeta
END FUNCTION hmap_axisNB_eval_Jh


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_Jh_dq1( sf ,q_in) RESULT(Jh_dq1)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh_dq1
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp,Tq
  REAL(wp)              :: TqTq,NN,BB,NTq,BTq,NB,Jh2 
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  Tq=(T+q1*Np+q2*Bp)
  TqTq=SUM(Tq(:)*Tq(:))
  NTq =SUM(N(:)*Tq(:))
  BTq =SUM(B(:)*Tq(:))
  NN  =SUM(N(:)*N(:))
  BB  =SUM(B(:)*B(:))
  NB  =SUM(N(:)*B(:))
  Jh2=TqTq*(NN*BB-NB*NB)  +2.0_wp*NB*Btq*Ntq - NTq*NTq*BB - BTq*BTq*NN
  IF(Jh2 .LT. 1.0e-4) &
      CALL abort(__STAMP__, &
           "hmap_axisNB,  Jh<0",RealInfo=zeta*sf%nfp/TWOPI)

  Jh_dq1=(SUM(Tq*Np)*(NN*BB-NB*NB) + SUM(B*Np)*(NB*NTq-NN*BTq) + SUM(N*Np)*(NB*BTq-BB*Ntq)  )/SQRT(Jh2) 

  END ASSOCIATE !zeta
END FUNCTION hmap_axisNB_eval_Jh_dq1


!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: dJ_h/dq^k, k=1,2 at q=(q^1,q^2,zeta) 
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_Jh_dq2( sf ,q_in) RESULT(Jh_dq2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: Jh_dq2
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp,Tq
  REAL(wp)              :: TqTq,NN,BB,NTq,BTq,NB,Jh2 
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  Tq=(T+q1*Np+q2*Bp)
  TqTq=SUM(Tq(:)*Tq(:))
  NTq =SUM(N(:)*Tq(:))
  BTq =SUM(B(:)*Tq(:))
  NN  =SUM(N(:)*N(:))
  BB  =SUM(B(:)*B(:))
  NB  =SUM(N(:)*B(:))
  Jh2=TqTq*(NN*BB-NB*NB)  +2.0_wp*NB*Btq*Ntq - NTq*NTq*BB - BTq*BTq*NN
  IF(Jh2 .LT. 1.0e-4) &
      CALL abort(__STAMP__, &
           "hmap_axisNB,  Jh<0",RealInfo=zeta*sf%nfp/TWOPI)

  Jh_dq2=(SUM(Tq*Bp)*(NN*BB-NB*NB) + SUM(B*Bp)*(NB*NTq-NN*BTq) + SUM(N*Bp)*(NB*BTq-BB*Ntq)  )/SQRT(Jh2) 
  END ASSOCIATE !zeta
END FUNCTION hmap_axisNB_eval_Jh_dq2


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: qL_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_G(3)
  REAL(wp)          , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: g_ab
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp,Tq
!===================================================================================================================================
  !                       |q1  |   |N.N   B.N   Tq.N |        |q1  |  
  !q_i G_ij q_j = (dalpha |q2  | ) |N.B   B.B   Tq.B | (dbeta |q2  | )
  !                       |q3  |   |N.Tq  B.Tq  Tq.Tq|        |q3  |  
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  Tq=(T+q1*Np+q2*Bp)
  g_ab=    SUM( N* N)*qL_in(1)*qR_in(1) &
         + SUM( B* B)*qL_in(2)*qR_in(2) &
         + SUM(Tq*Tq)*qL_in(3)*qR_in(3) &
         + SUM( N* B)*(qL_in(1)*qR_in(2)+qL_in(2)*qR_in(1)) &
         + SUM( N*Tq)*(qL_in(1)*qR_in(3)+qL_in(3)*qR_in(1)) &
         + SUM( B*Tq)*(qL_in(2)*qR_in(3)+qL_in(3)*qR_in(2))  

  END ASSOCIATE
END FUNCTION hmap_axisNB_eval_gij


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and 
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_gij_dq1( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq1)
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  REAL(wp)           , INTENT(IN   ) :: qL_in(3)
  REAL(wp)           , INTENT(IN   ) :: q_G(3)
  REAL(wp)           , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                           :: g_ab_dq1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp,Tq
!===================================================================================================================================
  !                            |q1  |   |0     0      N'.N  |        |q1  |  
  !q_i dG_ij/dq1 q_j = (dalpha |q2  | ) |0     0      N'.B  | (dbeta |q2  | )
  !                            |q3  |   |N.N'  B.N'  2N'.Tq |        |q3  |  
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  Tq=(T+q1*Np+q2*Bp)

  g_ab_dq1 =         SUM(N *Np)*(qL_in(1)*qR_in(3)+ qL_in(3)*qR_in(1)) &
             +       SUM(B *Np)*(qL_in(2)*qR_in(3)+ qL_in(3)*qR_in(2)) & 
             +2.0_wp*SUM(Tq*Np)*(qL_in(3)*qR_in(3))

  END ASSOCIATE
END FUNCTION hmap_axisNB_eval_gij_dq1


!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and 
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and 
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_axisNB_eval_gij_dq2( sf ,qL_in,q_G,qR_in) RESULT(g_ab_dq2)
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  REAL(wp)          , INTENT(IN   ) :: qL_in(3)
  REAL(wp)          , INTENT(IN   ) :: q_G(3)
  REAL(wp)          , INTENT(IN   ) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                          :: g_ab_dq2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp,Tq
!===================================================================================================================================
  !                            |q1  |   |0     0      B'.N  |        |q1  |  
  !q_i dG_ij/dq2 q_j = (dalpha |q2  | ) |0     0      B'.B  | (dbeta |q2  | )
  !                            |q3  |   |N.B'  B.B'  2B'.Tq |        |q3  |  
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp) 
  Tq=(T+q1*Np+q2*Bp)

  g_ab_dq2 =         SUM(N *Bp)*(qL_in(1)*qR_in(3)+ qL_in(3)*qR_in(1)) &
             +       SUM(B *Bp)*(qL_in(2)*qR_in(3)+ qL_in(3)*qR_in(2)) & 
             +2.0_wp*SUM(Tq*Bp)*(qL_in(3)*qR_in(3))

  END ASSOCIATE
END FUNCTION hmap_axisNB_eval_gij_dq2


!===================================================================================================================================
!> evaluate curve X0(zeta), and T,N,B,N',B'. Here still using the frenet formulas 
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_eval_TNB_from_frenet( sf,zeta,X0,T,N,B,Np,Bp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(IN ) :: sf
  REAL(wp)            , INTENT(IN ) :: zeta       !! position along closed curve parametrized in [0,2pi]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)            , INTENT(OUT) :: X0(1:3)      !! curve position in cartesian coordinates
  REAL(wp)            , INTENT(OUT) :: T(1:3)       !! tangent X0'
  REAL(wp)            , INTENT(OUT) :: N(1:3)       !! Normal
  REAL(wp)            , INTENT(OUT) :: B(1:3)       !! bi-Normal
  REAL(wp)            , INTENT(OUT) :: Np(1:3)      !! derivative of Normal in zeta (N')
  REAL(wp)            , INTENT(OUT) :: Bp(1:3)      !! derivative of bi-Normal in zeta  (B')
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0p,X0pp,X0ppp,Tp,N_add,Np_add,B_add,Bp_add
  REAL(wp)          :: lp,absB,kappa,tau,sigma
!===================================================================================================================================
  !USING FRENET FRAME HERE, FOR DEBUGGING THE REST...
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp) 
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_axisNB cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)

  tau=SUM(X0ppp*B)/(absB**2)
  B=B/absB
  N=CROSS(B,T)
  ! END FRENET FRAME

  !scaling with lprime, so changing to the derivative in zeta
  Np=lp*(-kappa*T+tau*B)
  Bp=lp*(-tau*N)

  sigma=MERGE(SIGN(1.0_wp,SIN(sf%nfp*zeta)),1.0_wp,sf%omnig)
  !add NFP periodic perturbation to make T,N,B non-orthogonal (T'=lp*kappa*N), for DEBUGGING
  Tp    = lp*kappa*N
  N_add =  0.05*sigma*(SIN(zeta*sf%nfp-0.5)*T)+0.0333*(COS(zeta*sf%nfp+0.2)*B)
  Np_add=  0.05  *sigma*( sf%nfp*COS(zeta*sf%nfp-0.5)*T+SIN(zeta*sf%nfp-0.5)*Tp) &
          +0.0333*(-sf%nfp*SIN(zeta*sf%nfp+0.2)*T+COS(zeta*sf%nfp+0.2)*Bp)
  B_add =  0.044*sigma*(COS(zeta*sf%nfp+0.1)*T)+0.023*(SIN(zeta*sf%nfp-0.1)*N)
  Bp_add=  0.044*sigma*(-sf%nfp*SIN(zeta*sf%nfp+0.1)*T+COS(zeta*sf%nfp+0.1)*Tp) &
          +0.023*( sf%nfp*SIN(zeta*sf%nfp-0.1)*N+SIN(zeta*sf%nfp-0.1)*Np)

  N=sigma*(N+N_add)
  B=sigma*(B+B_add)
  Np=sigma*(Np+Np_add)
  Bp=sigma*(Bp+Bp_add)
  !scaling with lprime, so changing to the derivative in zeta
  T =T*lp

END SUBROUTINE hmap_axisNB_eval_TNB_from_frenet

!===================================================================================================================================
!> evaluate curve X0(zeta), position and first three derivatives, from given R0,Z0 Fourier 
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_eval_X0_fromRZ( sf,zeta,X0,X0p,X0pp,X0ppp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(IN ) :: sf
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

END SUBROUTINE hmap_axisNB_eval_X0_fromRZ


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
!> test hmap_axisNB - evaluation of the map
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf  !!self
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
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_axisNB TEST ID',nTestCalled,'    >>>>>>>>>'
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
            '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
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
              '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be <',100*epsFD,' : |dxdqFD-eval_dxdq|= ', checkreal,", dq=",qdir
      END IF !TEST
    END DO

    !! TEST G_ij
    DO idir=1,3; DO jdir=1,3
      iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal= SUM(sf%eval_dxdq(q_in,q_test(idir,:))*sf%eval_dxdq(q_in,q_test(jdir,:))) 
      refreal  =sf%eval_gij(q_test(idir,:),q_in,q_test(jdir,:))
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),2(A,I3))') &
       '\n =>  should be ', refreal,' : sum|Gij-eval_gij|= ', checkreal,', i=',idir,', j=',jdir
      END IF !TEST
    END DO; END DO
    !! TEST dG_ij_dq1 with FD 
    DO qdir=1,2
    DO idir=1,3; DO jdir=1,3
      iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      gij  =sf%eval_gij(q_test(idir,:),q_in,q_test(jdir,:))
      gij_eps = sf%eval_gij(q_test(idir,:),q_in+epsFD*q_test(qdir,:),q_test(jdir,:))
      IF(qdir.EQ.1) refreal = sf%eval_gij_dq1(q_test(idir,:),q_in,q_test(jdir,:))
      IF(qdir.EQ.2) refreal = sf%eval_gij_dq2(q_test(idir,:),q_in,q_test(jdir,:))
      checkreal=(gij_eps-gij)/epsFD-refreal
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. 100*epsFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),3(A,I3))') &
       '\n =>  should be <', 100*epsFD,' : |dGij_dqFD-eval_gij_dq|= ', checkreal,', i=',idir,', j=',jdir,', dq=',qdir
      END IF !TEST
    END DO; END DO
    END DO

    
      
    
 END IF !testlevel >=1
 
 test_called=.FALSE. ! to prevent infinite loop in this routine
 

END SUBROUTINE hmap_axisNB_test

END MODULE MODgvec_hmap_axisNB

