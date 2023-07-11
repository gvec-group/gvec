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
USE MODgvec_fBase   ,ONLY: t_fbase
IMPLICIT NONE

PUBLIC
 

TYPE,EXTENDS(c_hmap) :: t_hmap_axisNB
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_axisNB:
  !INTEGER              :: nfp   !! already part of c_hmap. Is overwritten in init!
  !curve description
  INTEGER              :: nzeta=0       !! number of points in zeta direction of the input axis 
  REAL(wp),ALLOCATABLE :: zeta(:)       !! zeta positions in one field period (1:nzeta),  on 'half' grid: zeta(i)=(i-0.5)/nzeta*(2pi/nfp)
  INTEGER              :: n_max=0       !! maximum number of fourier coefficients
  REAL(wp),ALLOCATABLE :: xyz(:,:)      !! cartesian coordinates of the axis for a full turn, (1:NFP*nzeta,1:3), zeta is on 'half' grid: zeta(i)=(i-0.5)/(NFP*nzeta)*(2pi)
  REAL(wp),ALLOCATABLE :: Nxyz(:,:)     !! "normal" vector of axis frame in cartesian coordinates for a full turn (1:NFP*nzeta,1:3). NOT ASSUMED TO BE ORTHOGONAL to tangent of curve
  REAL(wp),ALLOCATABLE :: Bxyz(:,:)      !! "Bi-normal" vector of axis frame in cartesian coordinates for a full turn (1:NFP*nzeta,1:3). NOT ASSUMED TO BE ORTHOGONAL to tangent of curve or Nxyz
  REAL(wp),ALLOCATABLE :: xyz_modes(:,:)   !! fourier modes of xyz
  REAL(wp),ALLOCATABLE :: Nxyz_modes(:,:)   !! 1d fourier modes of Nxyz
  REAL(wp),ALLOCATABLE :: Bxyz_modes(:,:)   !! 1d fourier modes of Bxyz
  CHARACTER(LEN=100)   :: axis_ncfile=" " !! name of netcdf file with axis information
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp)             :: rot_origin(1:3)=(/0.,0.,0./)   !! origin of rotation, needed for checking field periodicity
  REAL(wp)             :: rot_axis(1:3)=(/0.,0.,1./)   !! rotation axis (unit length), needed for checking field periodicity
  !---------------------------------------------------------------------------------------------------------------------------------
  CLASS(t_fbase),ALLOCATABLE  :: fb  !! container for 1d fourier base on full turn

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
  PROCEDURE :: eval_TNB     => hmap_axisNB_eval_TNB
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
USE MODgvec_fbase,ONLY:fbase_new

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: n,i
  INTEGER :: nvisu,error_nfp
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: axisNB FRAME OF A CLOSED CURVE ...'


#if NETCDF
  ! read axis from netcdf
  nvisu=GETINT("hmap_nvisu",2*(sf%n_max+1)) 
  sf%axis_ncfile=GETSTR("hmap_ncfile")
  CALL ReadAxis_NETCDF(sf,sf%axis_ncfile)
#else
  CALL abort(__STAMP__,&
      "cannot read axis netcdf file, since code is compiled with BUILD_NETCDF=OFF")
#endif
  
  !initialize DOFS by projection.  nzeta*nfp >= 2*n_max*nfp+1  
  sf%n_max=MIN(sf%n_max,(sf%nzeta*sf%nfp-1)/(2*sf%nfp))

  CALL fbase_new(sf%fb,(/0,sf%n_max*sf%nfp/),(/1,sf%nzeta*sf%nfp/),1,"_sincos_",.FALSE.)

  IF(MAXVAL(ABS(sf%fb%X_IP(2,1:sf%nzeta)-sf%zeta)).GT.1.0e-14*sf%nzeta) &
     CALL abort(__STAMP__,&
          "zeta positions from axis file do not coincide with zeta positions in fbase.")

  WRITE(*,*)'DEBUG,MAXVAL(|N.B|)',MAXVAL(ABS(sf%Nxyz(:,1)*sf%Bxyz(:,1)+sf%Nxyz(:,2)*sf%Bxyz(:,2)+sf%Nxyz(:,3)*sf%Bxyz(:,3)))

  ALLOCATE(sf%xyz_modes(sf%fb%modes,3))
  ALLOCATE(sf%Nxyz_modes(sf%fb%modes,3))
  ALLOCATE(sf%Bxyz_modes(sf%fb%modes,3))
  DO i=1,3
     sf%xyz_modes(:,i) =sf%fb%initDOF(sf%xyz(:,i))
     sf%Nxyz_modes(:,i)=sf%fb%initDOF(sf%Nxyz(:,i))
     sf%Bxyz_modes(:,i)=sf%fb%initDOF(sf%Bxyz(:,i))
  END DO



  IF(nvisu.GT.0) CALL Visu_axisNB(sf,nvisu)
  
  CALL CheckFieldPeriodicity(sf,error_nfp)
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
  DEALLOCATE(sf%zeta)
  DEALLOCATE(sf%xyz)
  DEALLOCATE(sf%Nxyz)
  DEALLOCATE(sf%Bxyz)
  DEALLOCATE(sf%xyz_modes)
  DEALLOCATE(sf%Nxyz_modes)
  DEALLOCATE(sf%Bxyz_modes)
  CALL sf%fb%free()
  DEALLOCATE(sf%fb)

  sf%initialized=.FALSE.

END SUBROUTINE hmap_axisNB_free

#if NETCDF
!===================================================================================================================================
!> READ axis from netcdf file, needs netcdf library!
!> ======= HEADER OF THE NETCDF FILE VERSION 3.0 ===================================================================================
!> === FILE DESCRIPTION:
!>   * axis, normal and binormal of the frame are given in cartesian coordinates along the curve parameter zeta [0,2pi].
!>   * The curve is allowed to have a field periodicity NFP, but the curve must be provided on a full turn.
!>   * The adata is given in real space, sampled along equidistant zeta point positions:
!>       zeta(i)=(i+0.5)/nzeta * (2pi/NFP), i=0,...,nzeta-1
!>     always shifted by (2pi/NFP) for the next field period.
!>     Thus the number of points along the axis for a full turn is NFP*nzeta
!>   * definition of the axis-following frame in cartesian coordinates ( boundary surface at rho=1):
!>
!>      {x,y,z}(rho,theta,zeta)=axis_{x,y,z}(zeta) + X(rho,theta,zeta)*N_{x,y,z}(zeta)+Y(rho,theta,zeta)*B_{x,y,z}(zeta)  
!>
!> === DATA DESCRIPTION
!> - general data
!>   * NFP: number of field periods
!>   * VERSION: version number as integer: V3.0 => 300
!> - axis data group:
!>   * 'axis_nzeta'   : number of points along the axis, in one field period (>=2*n_max+1)
!>   * 'axis_zeta(:)' : zeta positions, 1D array of size 'axis_nzeta',  must exclude the end point, on half grid! zeta(i)=(i+0.5)/nzeta*(2pi/nfp), i=0,...nzeta-1
!>   * 'axis_xyz(::)' : cartesian positions along the axis for ONE FULL TURN, 2D array of size (3,NFP* axis_nzeta ), sampled at zeta positions, must exclude the endpoint
!>                      xyz[i,j]=axis_i(zeta[j]),        
!>   * 'axis_Nxyz(::)': cartesian components of the normal vector of the axis frame, 2D array of size (3, NFP* axis_nzeta), evaluated analogously to the axis
!>   * 'axis_Bxyz(::)': cartesian components of the bi-normal vector of the axis frame, 2D array of size (3, NFP*axis_nzeta), evaluated analogously to the axis
!> - boundary data group:
!>   * 'boundary_m_max'    : maximum mode number in theta 
!>   * 'boundary_n_max'    : maximum mode number in zeta (in one field period)
!>   * 'boundary_lasym'    : asymmetry, logical. 
!>                            if lasym=0, boundary surface position X,Y in the N-B plane of the axis frame can be represented only with
!>                              X(theta,zeta)=sum X_mn*cos(m*theta-n*NFP*zeta), with {m=0,n=0...n_max},{m=1...m_max,n=-n_max...n_max}
!>                              Y(theta,zeta)=sum Y_mn*sin(m*theta-n*NFP*zeta), with {m=0,n=1...n_max},{m=1...m_max,n=-n_max...n_max}
!>                            if lasym=1, full fourier series is taken for X,Y
!>   * 'boundary_ntheta'    : number of points in theta (>=2*m_max+1)
!>   * 'boundary_nzeta'     : number of points in zeta  (>=2*n_max+1)
!>   * 'boundary_theta(:)'  : theta positions, 1D array of size 'boundary_ntheta', on half grid! theta(i)=(i+0.5)/ntheta*(2pi), i=0,...ntheta-1
!>   * 'boundary_zeta(:)'   : zeta positions, 1D array of size 'boundary_nzeta', on half grid for one field period! zeta(i)=(i+0.5)/nzeta*(2pi/nfp), i=0,...nzeta-1
!>   * 'boundary_X(::)',
!>     'boundary_Y(::)'     : boundary position X,Y in the N-B plane of the axis frame, in one field period, 2D array of size(ntheta, nzeta),  with
!>                               X[i, j]=X(theta[i],zeta[j])
!>                               Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1                                         
!> 
!> ---- PLASMA PARAMETERS:
!>  ....
!> ======= END HEADER,START DATA ===================================================================================
!! 
!> NOTE THAT ONLY THE AXIS DATA IS NEEDED FOR THE AXIS DEFINITION 
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
  INTEGER :: nc_id,ioError
  INTEGER :: tmp_int
!===================================================================================================================================


  SWRITE(UNIT_stdOut,'(4X,A)')'READ AXIS FILE "'//TRIM(fileName)//'" in NETCDF format ...'

  !! open NetCDF input file
  ioError = NF_OPEN(TRIM(fileName), NF_NOWRITE, nc_id)
  IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"opening file",filename)

  CALL GETSCALAR_NC(nc_id,"NFP",intout=sf%nfp)
  CALL GETSCALAR_NC(nc_id,"axis_nzeta",intout=sf%nzeta)
  sf%n_max= (sf%nzeta-1)/2

  ALLOCATE(sf%zeta(sf%nzeta))
  CALL GETARR1D_NC(nc_id,"axis_zeta(:)",sf%nzeta,realout=sf%zeta)

  ALLOCATE(sf%xyz(sf%nfp*sf%nzeta,3))
  CALL GETARR2D_NC(nc_id,"axis_xyz(::)",sf%nfp*sf%nzeta,3,realout=sf%xyz)
  WRITE(*,*)'DEBUG,xyz(1:3,1)',sf%xyz(1,1:3)

  ALLOCATE(sf%Nxyz(sf%nfp*sf%nzeta,3))
  CALL GETARR2D_NC(nc_id,"axis_Nxyz(::)",sf%nfp*sf%nzeta,3,realout=sf%Nxyz)
  WRITE(*,*)'DEBUG,Nxyz(1:3,1)',sf%Nxyz(1,1:3)

  ALLOCATE(sf%Bxyz(sf%nfp*sf%nzeta,3))
  CALL GETARR2D_NC(nc_id,"axis_Bxyz(::)",sf%nfp*sf%nzeta,3,realout=sf%Bxyz)
  WRITE(*,*)'DEBUG,Bxyz(1:3,1)',sf%Bxyz(1,1:3)

  ioError = NF_CLOSE(nc_id)
  SWRITE(*,'(4X,A)')'...DONE.'

  CONTAINS


  !using nc_id,iError
  SUBROUTINE enter_groups(grpid_in,varname_in,grpid,varname) 
    INTEGER,INTENT(IN)          :: grpid_in
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
    CHARACTER(LEN=255),INTENT(OUT) :: varname
    INTEGER,INTENT(OUT)          :: grpid
    CHARACTER(LEN=255) :: grpname
    INTEGER          :: grpid_old,id
    !split the varname at first occurence of "/" to get the group name. Then get the group id. 
    ! repeat until no "/" is found anymore.
    ! output the final groupid and the variable name without the group names.
    grpid=grpid_in 
    varname=varname_in
    id=INDEX(varname,"/")
    DO WHILE (id.NE.0)
      grpname=varname(1:id-1)
      varname=varname(id+1:)
      grpid_old=grpid
      ioError = NF_INQ_NCID(grpid_old, TRIM(grpname), grpid) 
      IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"finding group",grpname)
      id=INDEX(varname,"/")
    END DO
  END SUBROUTINE Enter_groups


  SUBROUTINE GETSCALAR_NC(grpid_in,varname,intout,realout) 
    INTEGER,INTENT(IN)          :: grpid_in
    CHARACTER(LEN=*),INTENT(IN) :: varname
    INTEGER,INTENT(OUT),OPTIONAL:: intout
    REAL(wp),INTENT(OUT),OPTIONAL:: realout
    CHARACTER(LEN=255) :: tmpname
    INTEGER          :: grpid,id
    CALL enter_groups(grpid_in,varname,grpid,tmpname)
    ioError = NF_INQ_VARID(grpid, TRIM(tmpname), id) 
    IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"finding",varname)
    IF(PRESENT(intout))THEN
      ioError = NF_GET_VAR(grpid, id, intout)
      IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"reading",varname)
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,I6)')'read integer ',TRIM(varname),' :: ',intout
    ELSEIF(PRESENT(realout))THEN
      ioError = NF_GET_VAR(grpid, id, realout)
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,E21.11)')'read double  ',TRIM(varname),' :: ',realout
    END IF
  END SUBROUTINE GETSCALAR_NC


  SUBROUTINE GETARR1D_NC(grpid_in,varname,dim_1,intout,realout) 
    INTEGER,INTENT(IN)          :: grpid_in
    CHARACTER(LEN=*),INTENT(IN) :: varname
    INTEGER          :: dim_1
    INTEGER,INTENT(OUT),OPTIONAL:: intout(1:dim_1)
    REAL(wp),INTENT(OUT),OPTIONAL:: realout(1:dim_1)
    CHARACTER(LEN=255) :: tmpname
    INTEGER          :: grpid,id
    CALL enter_groups(grpid_in,varname,grpid,tmpname)
    ioError = NF_INQ_VARID(grpid, TRIM(tmpname), id) 
    IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"finding",varname)
    IF(PRESENT(intout))THEN
      ioError = NF_GET_VAR(grpid, id, intout )
      IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"reading",varname)
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,I6)')'read 1d integer array ',TRIM(varname)
    ELSEIF(PRESENT(realout))THEN
      ioError = NF_GET_VAR(grpid, id, realout )
      IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"reading",varname)
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,I6)')'read 1d double array ',TRIM(varname)
    END IF
  END SUBROUTINE GETARR1D_NC

  SUBROUTINE GETARR2D_NC(grpid_in,varname,dim_1,dim_2,intout,realout) 
    INTEGER,INTENT(IN)          :: grpid_in
    CHARACTER(LEN=*),INTENT(IN) :: varname
    INTEGER          :: dim_1,dim_2
    INTEGER,INTENT(OUT),OPTIONAL:: intout(1:dim_1,1:dim_2)
    REAL(wp),INTENT(OUT),OPTIONAL:: realout(1:dim_1,1:dim_2)
    CHARACTER(LEN=255) :: tmpname
    INTEGER          :: grpid,id
    CALL enter_groups(grpid_in,varname,grpid,tmpname)
    ioError = NF_INQ_VARID(grpid, TRIM(tmpname), id) 
    IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"finding",varname)
    IF(PRESENT(intout))THEN
      ioError = NF_GET_VAR(grpid, id, intout )
      IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"reading",varname)
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,I6)')'read 2d integer array ',TRIM(varname)
    ELSEIF(PRESENT(realout))THEN
      ioError = NF_GET_VAR(grpid, id, realout )
      IF (ioError .NE. NF_NOERR) CALL handle_error(ioError,"reading",varname)
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,I6)')'read 2d double array ',TRIM(varname)
    END IF
  END SUBROUTINE GETARR2D_NC

  SUBROUTINE handle_error(ioerr,act,varname)
    INTEGER, intent(in) :: ioerr
    CHARACTER(LEN=*) :: act,varname
    IF (ioerr .ne. NF_NOERR) THEN
       WRITE(*,'(6X,A)')"A netCDF error has occurred when "//TRIM(act)//" variable '"//TRIM(varname)//"'"
       CALL abort(__STAMP__,&
                 NF_STRERROR(ioerr))
    END IF
  END SUBROUTINE handle_error

END SUBROUTINE ReadAxis_NETCDF
#endif /*NETCDF*/

!===================================================================================================================================
!> Check that the TNB frame  really has the field periodicity of NFP: 
!> assumption for now is that the origin is fixed at rot_origin=(/0.,0.,0./)
!> and the rotation axis is fixed at rot_axis=(/0.,0.,1./)
!> the sign of rotation is being checked as well
!!
!===================================================================================================================================
SUBROUTINE CheckFieldPeriodicity( sf ,error_nfp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(out)  :: error_nfp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)              :: dzeta_fp
  REAL(wp),DIMENSION(3) :: X0_a,T,N_a,B_a,Np,Bp
  REAL(wp),DIMENSION(3) :: X0_b,N_b,B_b,P_a,P_b
  INTEGER               :: ifp,iz,error_x,error_n,error_b
!===================================================================================================================================
  error_nfp=0

  dzeta_fp=TWOPI/REAL(sf%nfp,wp)
  CALL sf%eval_TNB(sf%zeta(1),X0_a,T,N_a,B_a,Np,Bp) 
  P_a=X0_a-sf%rot_origin
  CALL sf%eval_TNB(sf%zeta(1)+dzeta_fp,X0_b,T,N_b,B_b,Np,Bp) 
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
    DO iz=1,sf%nzeta
      CALL sf%eval_TNB(sf%zeta(iz)+ ifp*dzeta_fp   ,X0_a,T,N_a,B_a,Np,Bp) 
      CALL sf%eval_TNB(sf%zeta(iz)+(ifp+1)*dzeta_fp,X0_b,T,N_b,B_b,Np,Bp) 
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
    zeta=(REAL(ivisu-1,wp))/REAL(nvisu*sf%nfp,wp)*TWOPI
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
!> evaluate curve X0(zeta), and T=X0',N,B,N',B', using the fourier series of X0,N and B in cartesian coordinates
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_eval_TNB( sf,zeta,X0,T,N,B,Np,Bp)
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
  INTEGER :: i
  REAL(wp)                      :: base_x(sf%fb%modes) 
  REAL(wp)                      :: base_dxdz(sf%fb%modes) 
!===================================================================================================================================
  base_x =sf%fb%eval(         0,(/0.,zeta/))
  base_dxdz=sf%fb%eval(DERIV_ZETA,(/0.,zeta/))

  X0(:)=MATMUL(base_x   ,sf%xyz_modes( :,:))
  T( :)=MATMUL(base_dxdz,sf%xyz_modes( :,:))

  N( :)=MATMUL(base_x   ,sf%Nxyz_modes(:,:))
  Np(:)=MATMUL(base_dxdz,sf%Nxyz_modes(:,:))

  B( :)=MATMUL(base_x   ,sf%Bxyz_modes(:,:))
  Bp(:)=MATMUL(base_dxdz,sf%Bxyz_modes(:,:))

END SUBROUTINE hmap_axisNB_eval_TNB


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
    q_in=(/0.0_wp, 0.0_wp, sf%zeta(sf%nzeta/2+1)/)
    x = sf%eval(q_in )
    checkreal=SUM((x-sf%xyz(sf%nzeta/2+1,:))**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    !evaluate at q1=0.44,q2=-0.33 (= x+0.44*N-0.33*B)
    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.44_wp,-0.33_wp, sf%zeta(sf%nzeta/2+1)/)
    x = sf%eval(q_in )
    checkreal=SUM((x-(sf%xyz(sf%nzeta/2+1,:)+0.44_wp*sf%Nxyz(sf%nzeta/2+1,:)-0.33_wp*sf%Bxyz(sf%nzeta/2+1,:)))**2)
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
      iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
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

