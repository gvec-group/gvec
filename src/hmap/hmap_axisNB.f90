!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
USE MODgvec_Globals, ONLY:PI,TWOPI,CROSS,wp,Unit_stdOut,abort,MPIroot
USE MODgvec_c_hmap,    ONLY:c_hmap
USE MODgvec_fBase   ,ONLY: t_fbase
USE MODgvec_io_netcdf   ,ONLY: t_ncfile
IMPLICIT NONE

PUBLIC


!---------------------------------------------------------------------------------------------------------------------------------
!> Store data that can be precomputed on a set ot zeta points
!> depends on hmap_axisNB, but could be used for different point sets in zeta
!
TYPE :: t_hmap_axisNB_aux
  INTEGER :: nzeta    !! size of zeta point positions
  REAL(wp),ALLOCATABLE :: zeta(:)       !! zeta positions, size(1:nzeta_eval)
  REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: X0,T,N,B,Np,Bp    !! Position,Tangent,Normal,Bi-Normal and N',B' size (3,1:nzeta_eval)
  REAL(wp),ALLOCATABLE,DIMENSION(:) :: BB,NN,NB,BpN,BpB,NpN,NpB !!dot-products of above vectors, size(nzeta_eval)
END TYPE t_hmap_axisNB_aux

TYPE,EXTENDS(c_hmap) :: t_hmap_axisNB
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_axisNB:
  !INTEGER              :: nfp   !! already part of c_hmap. Is overwritten in init!
  !curve description
  INTEGER              :: nzeta=0       !! number of points in zeta direction of the input axis
  INTEGER              :: sgn_rot       !! sign of rotation around Z axis,  either +1 or -1: positive means that from one field period to the next,
                                        !! xyz rotate counterclockwise around the Z-axis (right hand rule), negative then clockwise.
  REAL(wp),ALLOCATABLE :: zeta(:)       !! zeta positions in one field period (1:nzeta),  on 'half' grid: zeta(i)=(i-0.5)/nzeta*(2pi/nfp)
  INTEGER              :: n_max=0       !! maximum number of fourier coefficients on a field period (<=2*nzeta-1)
  REAL(wp),ALLOCATABLE :: xyz(:,:)      !! cartesian coordinates of the axis for a full turn, (1:NFP*nzeta,1:3), zeta is on 'half' grid: zeta(i)=(i-0.5)/(NFP*nzeta)*(2pi)
  REAL(wp),ALLOCATABLE :: Nxyz(:,:)     !! "normal" vector of axis frame in cartesian coordinates for a full turn (1:NFP*nzeta,1:3). NOT ASSUMED TO BE ORTHOGONAL to tangent of curve
  REAL(wp),ALLOCATABLE :: Bxyz(:,:)      !! "Bi-normal" vector of axis frame in cartesian coordinates for a full turn (1:NFP*nzeta,1:3). NOT ASSUMED TO BE ORTHOGONAL to tangent of curve or Nxyz
  !fourier modes
  REAL(wp),ALLOCATABLE :: xyz_hat_modes(:,:)   !! fourier modes of xhat,yhat,zhat on one field period,
                                               !! x=cos(zeta)xhat-sgn_rot*sin(zeta)yhat,
                                               !! y=sin(zeta)xhat+sgn_rot*cos(zeta)yhat, z=zhat
  REAL(wp),ALLOCATABLE :: Nxyz_hat_modes(:,:)   !! 1d fourier modes of Nxyz, one field period
  REAL(wp),ALLOCATABLE :: Bxyz_hat_modes(:,:)   !! 1d fourier modes of Bxyz, one field period

  CHARACTER(LEN=100)   :: ncfile=" " !! name of netcdf file with axis information
  !---------------------------------------------------------------------------------------------------------------------------------
  TYPE(t_hmap_axisNB_aux)      :: aux !! container for preevaluations
  TYPE(t_fbase),ALLOCATABLE   :: fb_hat  !! container for 1d fourier base of xhat
  CLASS(t_ncfile),ALLOCATABLE  :: nc  !! container for netcdf-file


  CONTAINS

  PROCEDURE :: init          => hmap_axisNB_init
  PROCEDURE :: free          => hmap_axisNB_free
  PROCEDURE :: eval          => hmap_axisNB_eval
  PROCEDURE :: init_aux      => hmap_axisNB_init_aux
  PROCEDURE :: free_aux      => hmap_axisNB_free_aux
  PROCEDURE :: eval_aux      => hmap_axisNB_eval_aux
  PROCEDURE :: eval_dxdq     => hmap_axisNB_eval_dxdq
  PROCEDURE :: eval_Jh       => hmap_axisNB_eval_Jh
  PROCEDURE :: eval_Jh_dq1   => hmap_axisNB_eval_Jh_dq1
  PROCEDURE :: eval_Jh_dq2   => hmap_axisNB_eval_Jh_dq2
  PROCEDURE :: eval_gij      => hmap_axisNB_eval_gij
  PROCEDURE :: eval_gij_dq1  => hmap_axisNB_eval_gij_dq1
  PROCEDURE :: eval_gij_dq2  => hmap_axisNB_eval_gij_dq2
  !---------------------------------------------------------------------------------------------------------------------------------
  ! procedures for hmap_axisNB:
  PROCEDURE :: eval_TNB     => hmap_axisNB_eval_TNB_hat
END TYPE t_hmap_axisNB

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

! SUBROUTINE init_dummy( sf )
! IMPLICIT NONE
!   CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!   CALL abort(__STAMP__, &
!              "dummy init in hmap_axisNB should not be used")
! END SUBROUTINE init_dummy

!===================================================================================================================================
!> initialize the type hmap_axisNB with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_init( sf )
! MODULES
USE MODgvec_ReadInTools,ONLY: GETLOGICAL,GETINT, GETREALARRAY,GETSTR
USE MODgvec_fbase      ,ONLY: fbase_new
USE MODgvec_io_netcdf  ,ONLY: ncfile_init
USE MODgvec_MPI        ,ONLY: par_BCast,par_barrier
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
  INTEGER :: nvisu,error_nfp
  REAL(wp),ALLOCATABLE :: cosz(:),sinz(:)
!===================================================================================================================================
  CALL par_Barrier(beforeScreenOut='INIT HMAP :: axisNB FRAME OF A CLOSED CURVE ...')

  sf%ncfile=GETSTR("hmap_ncfile")
  IF(MPIroot)THEN
    ! read axis from netcdf
    CALL ncfile_init(sf%nc,sf%ncfile,"r")
    CALL ReadNETCDF(sf)


    sf%sgn_rot=1
    IF(sf%nfp.GT.2)THEN !no need to check for nfp=1 and nfp=2 is either a positive or negative rotation by pi
      IF(SUM((sf%xyz(:,1+sf%nzeta) -  rodrigues(sf%xyz(:,1),TWOPI/REAL(sf%nfp,wp)))**2).LT.1.0e-12)THEN
         sf%sgn_rot=1
      ELSEIF(SUM((sf%xyz(:,1+sf%nzeta) -  rodrigues(sf%xyz(:,1),-TWOPI/REAL(sf%nfp,wp)))**2).LT.1.0e-12)THEN
         sf%sgn_rot=-1
      ELSE
        CALL abort(__STAMP__, &
           'problem with check of field period rotation: point at zeta=0 does not rotate to point at zeta= +/-2pi/nfp.')
      END IF
    END IF
    WRITE(UNIT_stdOut,'(4X,A,I2)')'INFO: sign of the rotation from zeta=0 to zeta=2pi/nfp is: ',sf%sgn_rot

    sf%n_max=(sf%nzeta-1)/2 ! maximum mode number on a field period


  END IF !MPIroot
  CALL par_BCast(sf%nzeta,0)
  CALL par_BCast(sf%nfp,0)
  CALL par_BCast(sf%n_max,0)
  CALL par_Bcast(sf%sgn_rot,0)
  IF(.NOT.MPIroot) CALL allocate_readin_vars(sf)
  CALL par_Bcast(sf%zeta,0)
  CALL par_Bcast(sf%xyz,0)
  CALL par_Bcast(sf%Nxyz,0)
  CALL par_Bcast(sf%Bxyz,0)
  !Fourier 1D base on one field period for "hat" coordinates
  CALL fbase_new(sf%fb_hat,(/0,sf%n_max/),(/1,sf%nzeta/),sf%nfp,"_sincos_",.FALSE.)
  ALLOCATE(sf%xyz_hat_modes( 3,sf%fb_hat%modes),&
           sf%Nxyz_hat_modes(3,sf%fb_hat%modes),&
           sf%Bxyz_hat_modes(3,sf%fb_hat%modes))
  IF(MPIroot)THEN
    ALLOCATE(cosz(sf%nzeta*sf%nfp),sinz(sf%nzeta*sf%nfp))
    !build cos(zeta),sin(zeta) on full torus
    DO i=1,sf%nfp
      cosz(sf%nzeta*(i-1)+1:sf%nzeta*i)=COS(sf%zeta(1:sf%nzeta)+TWOPI*(REAL(i-1,wp)/REAL(sf%nfp,wp)))
      sinz(sf%nzeta*(i-1)+1:sf%nzeta*i)=SIN(sf%zeta(1:sf%nzeta)+TWOPI*(REAL(i-1,wp)/REAL(sf%nfp,wp)))
    END DO
    sf%xyz_hat_modes =transform_to_hat(sf%xyz, "xyz to xyzhat"  )
    sf%Nxyz_hat_modes=transform_to_hat(sf%Nxyz,"Nxyz to Nxyzhat")
    sf%Bxyz_hat_modes=transform_to_hat(sf%Bxyz,"Bxyz to Bxyzhat")
    DEALLOCATE(cosz,sinz)

    nvisu=GETINT("hmap_nvisu",2*(sf%n_max*sf%nfp+1))

    IF(nvisu.GT.0) CALL Visu_axisNB(sf,nvisu)

    CALL CheckFieldPeriodicity(sf,sf%sgn_rot,error_nfp)
    IF(error_nfp.LT.0) &
       CALL abort(__STAMP__, &
          "hmap_axisNB check Field Periodicity failed!")
  END IF !MPIroot

  CALL par_BCast(sf%xyz_hat_modes,0)
  CALL par_BCast(sf%Nxyz_hat_modes,0)
  CALL par_BCast(sf%Bxyz_hat_modes,0)

  sf%aux%nzeta=0
  sf%initialized=.TRUE.
  CALL par_barrier(afterScreenOut='...DONE')

  IF(.NOT.test_called) CALL hmap_axisNB_test(sf)

  !------------------
  CONTAINS
  !------------------
  !! transform from full period cartesian coordinates to one-field period "hat" cartesian coordinates,
  !! by rotating around the zaxis with the local angle zeta,
  !! with the sign given by which direction xyz(zeta=0) rotates to xyz(zeta=2pi/nfp)
  !! INVERSION OF: x=cos(zeta)*xhat - sgn_rot*sin(zeta)*yhat,
  !!               y=cos(zeta)*yhat + sgn_rot*sin(zeta)*xhat, z <=> zhat
  !! ==>           xhat=cos(zeta)x+sgn_rot*sin(zeta)y,
  !!               yhat=cos(zeta)y-sgn_rot*sin(zeta)x,
  !! check that all points on full period are the same in the xhat,yhat,zhat coordinates
  !! NOTE THIS FUNCTION IS USING sinz=sin(zeta),cosz=cos(zeta) and sgn_rot computed above!!!!
  FUNCTION transform_to_hat(xyz_in,msg)  RESULT(to_hat_modes)
    IMPLICIT NONE

    REAL(wp),INTENT(IN) :: xyz_in(3,sf%nzeta*sf%nfp)
    CHARACTER(*),INTENT(IN) :: msg
    REAL(wp) :: to_hat_modes(3,sf%fb_hat%modes)
    !------------------
    REAL(wp) :: to_hat(3,sf%nzeta*sf%nfp)
    REAL(wp) :: check_xhat(1:3,2:sf%nfp)
    !------------------
    to_hat(1,:)=xyz_in(1,:)*cosz(:)+sf%sgn_rot*xyz_in(2,:)*sinz(:)
    to_hat(2,:)=xyz_in(2,:)*cosz(:)-sf%sgn_rot*xyz_in(1,:)*sinz(:)
    to_hat(3,:)=xyz_in(3,:)
    !check periodicity for remaining nfp
    DO i=2,sf%nfp
      check_xhat(1,i)=MAXVAL(ABS(to_hat(1,1:sf%nzeta)-to_hat(1,sf%nzeta*(i-1)+1:sf%nzeta*i)))
      check_xhat(2,i)=MAXVAL(ABS(to_hat(2,1:sf%nzeta)-to_hat(2,sf%nzeta*(i-1)+1:sf%nzeta*i)))
      check_xhat(3,i)=MAXVAL(ABS(to_hat(3,1:sf%nzeta)-to_hat(3,sf%nzeta*(i-1)+1:sf%nzeta*i)))
    END DO
    IF(ANY(check_xhat.GT.1.0e-12))THEN
      DO i=2,sf%nfp
        WRITE(UNIT_stdout,'(A,I4,A,3E9.2)')'fp=1 vs fp=',i,', check_xyz=',check_xhat(1:3,i)
      END DO
      CALL abort(__STAMP__,&
            "transform from cartesian to hat coordinates"//TRIM(msg)//" yields non-field periodic data!")
    END IF
    DO i=1,3
      to_hat_modes(i,:)=sf%fb_hat%initDOF(to_hat(i,1:sf%nzeta),thet_zeta_start=(/0.,sf%zeta(1)/))
    END DO
  END FUNCTION transform_to_hat

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
  CALL sf%free_aux()

  SDEALLOCATE(sf%zeta)
  SDEALLOCATE(sf%xyz)
  SDEALLOCATE(sf%Nxyz)
  SDEALLOCATE(sf%Bxyz)
  SDEALLOCATE(sf%xyz_hat_modes)
  SDEALLOCATE(sf%Nxyz_hat_modes)
  SDEALLOCATE(sf%Bxyz_hat_modes)
  IF(ALLOCATED(sf%fb_hat))THEN
    CALL sf%fb_hat%free()
    DEALLOCATE(sf%fb_hat)
  END IF
  IF(ALLOCATED(sf%nc))THEN
    CALL sf%nc%free()
    DEALLOCATE(sf%nc)
  END IF

  sf%initialized=.FALSE.

END SUBROUTINE hmap_axisNB_free

!===================================================================================================================================
!> initialize the aux variable
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_init_aux( sf ,nzeta_aux,zeta_aux)
  ! MODULES
    IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER, INTENT(IN)  :: nzeta_aux
    REAL(wp),INTENT(IN)  :: zeta_aux(1:nzeta_aux)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
  !===================================================================================================================================
  IF(sf%aux%nzeta.NE.0) CALL sf%free_aux()
  sf%aux%nzeta=nzeta_aux
  ALLOCATE(sf%aux%zeta(nzeta_aux))
  sf%aux%zeta=zeta_aux
  ALLOCATE(sf%aux%X0(3,nzeta_aux))
  ALLOCATE(sf%aux%T  ,sf%aux%N  ,sf%aux%B  ,sf%aux%Np ,sf%aux%Bp ,mold=sf%aux%X0)
  ALLOCATE(sf%aux%NN ,sf%aux%BB ,sf%aux%NB ,sf%aux%NpN, &
           sf%aux%NpB,sf%aux%BpN,sf%aux%BpB,    mold=sf%aux%zeta)

END SUBROUTINE hmap_axisNB_init_aux

!===================================================================================================================================
!> free the aux variable
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_free_aux( sf)
  ! MODULES
    IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
  !===================================================================================================================================
    IF(sf%aux%nzeta.NE.0)THEN !ALLOCATED
      sf%aux%nzeta=0
      DEALLOCATE(sf%aux%X0,sf%aux%T,sf%aux%N,sf%aux%B,sf%aux%Np,sf%aux%Bp)
      DEALLOCATE(sf%aux%zeta,sf%aux%NN,sf%aux%BB,sf%aux%NB,sf%aux%NpN,sf%aux%NpB,sf%aux%BpN,sf%aux%BpB)
    END IF
  END SUBROUTINE hmap_axisNB_free_aux


!===================================================================================================================================
!> evaluate the auxiliar variables on the given zeta points
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_eval_aux(sf)
  ! MODULES
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! sf%aux variables
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: izeta
  !===================================================================================================================================
!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE) SHARED(sf) PRIVATE(izeta)
  DO izeta=1,sf%aux%nzeta
    CALL sf%eval_TNB(sf%aux%zeta(izeta),sf%aux%X0(:,izeta),sf%aux%T(:,izeta),sf%aux%N(:,izeta), &
                     sf%aux%B(:,izeta),sf%aux%Np(:,izeta),sf%aux%Bp(:,izeta))
    sf%aux%NN( izeta)=SUM(sf%aux%N( :,izeta)* sf%aux%N(:,izeta))
    sf%aux%BB( izeta)=SUM(sf%aux%B( :,izeta)* sf%aux%B(:,izeta))
    sf%aux%NB( izeta)=SUM(sf%aux%N( :,izeta)* sf%aux%B(:,izeta))
    sf%aux%NpN(izeta)=SUM(sf%aux%Np(:,izeta)* sf%aux%N(:,izeta))
    sf%aux%NpB(izeta)=SUM(sf%aux%Np(:,izeta)* sf%aux%B(:,izeta))
    sf%aux%BpN(izeta)=SUM(sf%aux%Bp(:,izeta)* sf%aux%N(:,izeta))
    sf%aux%BpB(izeta)=SUM(sf%aux%Bp(:,izeta)* sf%aux%B(:,izeta))
  END DO !izeta
!$OMP END PARALLEL DO
END SUBROUTINE hmap_axisNB_eval_aux

!===================================================================================================================================
!> READ axis from netcdf file, needs netcdf library!
!! ======= HEADER OF THE NETCDF FILE VERSION 3.1 ===================================================================================
!! === FILE DESCRIPTION:
!!   * axis, normal and binormal of the frame are given in cartesian coordinates along the curve parameter zeta [0,2pi].
!!   * The curve is allowed to have a field periodicity NFP, but the curve must be provided on a full turn.
!!   * The adata is given in real space, sampled along equidistant zeta point positions:
!!       zeta(i)=(i+0.5)/nzeta * (2pi/NFP), i=0,...,nzeta-1
!!     always shifted by (2pi/NFP) for the next field period.
!!     Thus the number of points along the axis for a full turn is NFP*nzeta
!!   * definition of the axis-following frame in cartesian coordinates ( boundary surface at rho=1):
!!
!!      {x,y,z}(rho,theta,zeta)={x,y,z}(zeta) + X(rho,theta,zeta)*N_{x,y,z}(zeta)+Y(rho,theta,zeta)*B_{x,y,z}(zeta)
!!
!! === DATA DESCRIPTION
!! - general data
!!   * NFP: number of field periods
!!   * VERSION: version number as integer: V3.0 => 300
!! - axis/ data group:
!!   * 'axis/n_max'   : maximum mode number in zeta (in one field period)
!!   * 'axis/nzeta'   : number of points along the axis, in one field period (>=2*n_max+1)
!!   * 'axis/zeta(:)' : zeta positions, 1D array of size 'axis/nzeta', for one field period. zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp), i=1,..nzeta, zeta[1] is arbitrary
!!   * 'axis/xyz(::)' : cartesian positions along the axis for ONE FULL TURN, 2D array of size (3,NFP* nzeta ), sampled at zeta positions, must exclude the endpoint
!!                      xyz[:,j+fp*nzeta]=axis(zeta[j]+fp*2pi/NFP), for j=0,..nzeta-1 and  fp=0,...,NFP-1
!!   * 'axis/Nxyz(::)': cartesian components of the normal vector of the axis frame, 2D array of size (3, NFP* nzeta), evaluated analogously to the axis
!!   * 'axis/Bxyz(::)': cartesian components of the bi-normal vector of the axis frame, 2D array of size (3, NFP*nzeta), evaluated analogously to the axis
!! - boundary data group:
!!   * 'boundary/m_max'    : maximum mode number in theta
!!   * 'boundary/n_max'    : maximum mode number in zeta (in one field period)
!!   * 'boundary/lasym'    : asymmetry, logical.
!!                            if lasym=0, boundary surface position X,Y in the N-B plane of the axis frame can be represented only with
!!                              X(theta,zeta)=sum X_mn*cos(m*theta-n*NFP*zeta), with {m=0,n=0...n_max},{m=1...m_max,n=-n_max...n_max}
!!                              Y(theta,zeta)=sum Y_mn*sin(m*theta-n*NFP*zeta), with {m=0,n=1...n_max},{m=1...m_max,n=-n_max...n_max}
!!                            if lasym=1, full fourier series is taken for X,Y
!!   * 'boundary/ntheta'    : number of points in theta (>=2*m_max+1)
!!   * 'boundary/nzeta'     : number of points in zeta  (>=2*n_max+1), can be different to 'axis/nzeta'
!!   * 'boundary/theta(:)'  : theta positions, 1D array of size 'boundary/ntheta',  theta[i]=theta[1] + (i-1)/ntheta*(2pi), starting value arbitrary
!!   * 'boundary/zeta(:)'   : zeta positions, 1D array of size 'boundary/nzeta', for one field period! zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp), i=1,..nzeta, zeta[1] is arbitrary
!!   * 'boundary/X(::)',
!!     'boundary/Y(::)'     : boundary position X,Y in the N-B plane of the axis frame, in one field period, 2D array of size(ntheta, nzeta),  with
!!                               X[i, j]=X(theta[i],zeta[j])
!!                               Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1
!!
!! ---- PLASMA PARAMETERS:
!!  ....
!! ======= END HEADER,START DATA ===================================================================================
!!
!! NOTE THAT ONLY THE AXIS DATA IS NEEDED FOR THE AXIS DEFINITION
!===================================================================================================================================
SUBROUTINE ReadNETCDF(sf)
  USE MODgvec_io_netcdf
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.MPIroot)RETURN
  WRITE(UNIT_stdOut,'(4X,A)')'READ AXIS FILE "'//TRIM(sf%nc%fileName)//'" in NETCDF format ...'

  CALL sf%nc%get_scalar("NFP",intout=sf%nfp)
  !sf%nzeta=sf%nc%get_dimension("axis/nzeta")
  CALL sf%nc%get_scalar("axis/nzeta",intout=sf%nzeta)

  CALL allocate_readin_vars(sf)
  CALL sf%nc%get_array("axis/zeta(:)",realout_1d=sf%zeta)

  CALL sf%nc%get_array("axis/xyz(::)",realout_2d=sf%xyz)

  CALL sf%nc%get_array("axis/Nxyz(::)",realout_2d=sf%Nxyz)

  CALL sf%nc%get_array("axis/Bxyz(::)",realout_2d=sf%Bxyz)
  !SWRITE(*,*)'DEBUG,zeta(1)',sf%zeta(1)
  !SWRITE(*,*)'DEBUG,xyz(1:3,1)',sf%xyz(1:3,1)
  !SWRITE(*,*)'DEBUG,Nxyz(1:3,1)',sf%Nxyz(1:3,1)
  !SWRITE(*,*)'DEBUG,Bxyz(1:3,1)',sf%Bxyz(1:3,1)
  CALL sf%nc%closefile()
  WRITE(UNIT_stdout,'(4X,A)')'...DONE.'

END SUBROUTINE ReadNETCDF

SUBROUTINE allocate_readin_vars(sf)
  IMPLICIT NONE
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  IF(sf%nzeta.EQ.0) CALL abort(__STAMP__, &
       'sf%nzeta must be set before allocation')
  ALLOCATE(sf%zeta(sf%nzeta))
  ALLOCATE(sf%xyz(3,sf%nfp*sf%nzeta))
  ALLOCATE(sf%Nxyz(3,sf%nfp*sf%nzeta))
  ALLOCATE(sf%Bxyz(3,sf%nfp*sf%nzeta))
END SUBROUTINE allocate_readin_vars

!===================================================================================================================================
!> Check that the TNB frame  really has the field periodicity of NFP:
!! assumption for now is that the origin is fixed at rot_origin=(/0.,0.,0./)
!! and the rotation axis is fixed at rot_axis=(/0.,0.,1./)
!! sign of the rotation 'sgn_rot' is now accounted for in the transformation to xhat, so it has to be passed here.
!!
!===================================================================================================================================
SUBROUTINE CheckFieldPeriodicity( sf ,sgn_rot,error_nfp)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT) :: sf
  INTEGER, INTENT(in)   :: sgn_rot !! sign of rotation
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(out)  :: error_nfp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)              :: dzeta_fp
  REAL(wp),DIMENSION(3) :: X0_a,T,N_a,B_a,Np,Bp
  REAL(wp),DIMENSION(3) :: X0_b,N_b,B_b
  INTEGER               :: ifp,iz,error_x,error_n,error_b
!===================================================================================================================================
  error_nfp=0

  dzeta_fp=TWOPI/REAL(sf%nfp,wp)
  CALL sf%eval_TNB(sf%zeta(1),X0_a,T,N_a,B_a,Np,Bp)
  CALL sf%eval_TNB(sf%zeta(1)+dzeta_fp,X0_b,T,N_b,B_b,Np,Bp)
  DO ifp=0,sf%nfp
    error_x=0
    error_n=0
    error_b=0
    DO iz=1,sf%nzeta
      CALL sf%eval_TNB(sf%zeta(iz)+ ifp*dzeta_fp   ,X0_a,T,N_a,B_a,Np,Bp)
      CALL sf%eval_TNB(sf%zeta(iz)+(ifp+1)*dzeta_fp,X0_b,T,N_b,B_b,Np,Bp)
      IF(.NOT.(SUM((X0_b -  rodrigues(X0_a,sgn_rot*dzeta_fp))**2).LT.1.0e-12))THEN
        error_x=error_x+1
      END IF
      IF(.NOT.(SUM((X0_b+N_b -  rodrigues(X0_a+N_a,sgn_rot*dzeta_fp))**2).LT.1.0e-12))THEN
        error_n=error_n+1
      END IF
      IF(.NOT.(SUM((X0_b+B_b -  rodrigues(X0_a+B_a,sgn_rot*dzeta_fp))**2).LT.1.0e-12))THEN
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

END SUBROUTINE CheckFieldPeriodicity


!===================================================================================================================================
!> Rodrigues' rotation formula
!> assumption for now is that the origin is fixed at rot_origin=(/0.,0.,0./)
!> and the rotation axis is fixed at rot_axis=(/0.,0.,1./)
!!
!===================================================================================================================================
FUNCTION rodrigues(pos,angle) RESULT(pos_rot)
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  REAL(wp) :: pos(3),angle
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp) :: pos_rot(3)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  REAL(wp) :: vec(3),vec_rot(3)
  REAL(wp),PARAMETER :: origin(1:3)=(/0.,0.,0./)   !! origin of rotation, needed for checking field periodicity
  REAL(wp),PARAMETER :: axis(1:3)=(/0.,0.,1./)   !! rotation axis (unit length), needed for checking field periodicity
  !=================================================================================================================================
  vec=pos-origin
  vec_rot = vec*COS(angle)+CROSS(axis,vec)*SIN(angle)+axis*SUM(axis*vec)*(1.0_wp-COS(angle))
  pos_rot=origin+vec_rot
END FUNCTION rodrigues

!===================================================================================================================================
!> Write evaluation of the axis and signed axisNB frame to file
!!
!===================================================================================================================================
SUBROUTINE Visu_axisNB( sf ,nvisu)
! MODULES
USE MODgvec_Output_CSV,     ONLY: WriteDataToCSV
USE MODgvec_Output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_Output_netcdf,     ONLY: WriteDataToNETCDF
USE MODgvec_Analyze_vars,     ONLY: outfileType
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
  INTEGER,PARAMETER     :: nVars=21
  CHARACTER(LEN=20)     :: VarNames(1:nVars)
  REAL(wp)              :: values(1:nVars,1:nvisu*sf%nfp+1)
!===================================================================================================================================
  IF(nvisu.LE.0) RETURN
  iVar=0
  VarNames(ivar+1:iVar+3)=(/ "x", "y", "z"/)   ;iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"TX","TY","TZ"/)   ;iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"NX","NY","NZ"/)   ;iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"BX","BY","BZ"/)   ;iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"NpX","NpY","NpZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"BpX","BpY","BpZ"/);iVar=iVar+3
  VarNames(iVar+1       )="zeta"               ;iVar=iVar+1
  VarNames(iVar+1       )="zeta_norm"          ;iVar=iVar+1
  VarNames(iVar+1       )="lprime"             ;iVar=iVar+1

!  values=0.
  DO ivisu=1,nvisu*sf%nfp+1
    zeta=(0.5_wp+REAL(ivisu-1,wp))/REAL(nvisu*sf%nfp,wp)*TWOPI
    CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp)
    lp=SQRT(SUM(T*T))
    iVar=0
    values(ivar+1:iVar+3,ivisu)=X0                 ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=T                  ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=N                  ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=B                  ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=Np                 ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=Bp                 ;iVar=iVar+3
    values(iVar+1       ,ivisu)=zeta               ;iVar=iVar+1
    values(iVar+1       ,ivisu)=zeta*sf%nfp/TWOPI  ;iVar=iVar+1
    values(iVar+1       ,ivisu)=lp                 ;iVar=iVar+1
  END DO !ivisu
  IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
    CALL WriteDataToVTK(1,3,nVars-3,(/nvisu*sf%nfp/),1,VarNames(4:nVars),values(1:3,:),values(4:nVars,:),"visu_hmap_axisNB.vtu")
  END IF
  IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
#if NETCDF
    CALL WriteDataToNETCDF(1,3,nVars-3,(/nvisu*sf%nfp/),(/"dim_zeta"/),VarNames(4:nVars),values(1:3,:),values(4:nVars,:), &
        "visu_hmap_axisNB")
#else
    CALL WriteDataToCSV(VarNames(:) ,values, TRIM("out_visu_hmap_axisNB.csv") ,append_in=.FALSE.)
#endif
  END IF
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
!> evaluate all metrics necesseray for optimizer
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_eval_all(sf,ndims,dim_zeta,q1_in,q2_in,q1_thet,q2_thet,q1_zeta,q2_zeta, &
                                Jh,g_tt    ,g_tz    ,g_zz    , &
                                Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2  )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_axisNB), INTENT(INOUT)                  :: sf
  INTEGER                           ,INTENT(IN)   :: ndims(3)    !! 3D dimensions
  INTEGER                           ,INTENT(IN)   :: dim_zeta    !! which dimension is zeta dependent
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(IN)   :: q1_in,q2_in,q1_thet,q2_thet,q1_zeta,q2_zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(OUT)   :: Jh,g_tt    ,g_tz    ,g_zz    , &
                                                                  Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                                                  Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,j,k
  !===================================================================================================================================
#define OMPDEF PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
#define TEMPLATE_LOOP(IZT) \
  DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1);\
    CALL eval_all_e(sf%aux,IZT, \
                    q1_in(i,j,k),q2_in(i,j,k),q1_thet(i,j,k),q2_thet(i,j,k),q1_zeta(i,j,k),q2_zeta(i,j,k), \
                    Jh(i,j,k),g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k)    , \
                    Jh_dq1(i,j,k)     ,g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), \
                    Jh_dq2(i,j,k)     ,g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k) ); \
  END DO; END DO; END DO \

  SELECT CASE(dim_zeta)
  CASE(1)
#define IZT i
!$OMP OMPDEF
    TEMPLATE_LOOP(IZT)
!$OMP END PARALLEL DO
#undef IZT
  CASE(2)
#define IZT j
!$OMP OMPDEF
    TEMPLATE_LOOP(IZT)
!$OMP END PARALLEL DO
#undef IZT
  CASE(3)
#define IZT k
!$OMP OMPDEF
    TEMPLATE_LOOP(IZT)
!$OMP END PARALLEL DO
#undef IZT
  END SELECT
#undef OMPDEF
#undef TEMPLATE_LOOP

END SUBROUTINE hmap_axisNB_eval_all

!===================================================================================================================================
!> evaluate all quantities at one given point (elemental)
!!
!===================================================================================================================================
PURE SUBROUTINE eval_all_e(aux,iz,q1,q2,q1t,q2t,q1z,q2z, &
                                     Jh,g_tt,    g_tz,    g_zz,     &
                                     Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                     Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  TYPE(t_hmap_axisNB_aux),INTENT(IN) :: aux
  INTEGER ,INTENT(IN)  :: iz
  REAL(wp),INTENT(IN)  :: q1,q2       !! solution variables q1,q2
  REAL(wp),INTENT(IN)  :: q1t,q2t     !! theta derivative of solution variables q1,q2
  REAL(wp),INTENT(IN)  :: q1z,q2z     !!  zeta derivative of solution variables q1,q2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),INTENT(OUT) :: Jh,g_tt,g_tz,g_zz              !! Jac,1/Jac,g_{ab} with a=theta/zeta b=theta/zeta
  REAL(wp),INTENT(OUT) :: Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1  !! and their variation vs q1
  REAL(wp),INTENT(OUT) :: Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2  !! and their variation vs q2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  !REAL(wp)            :: TqTq,NTq,BTq,NpTq,BpTq    !! solution dependent dot products
  !REAL(wp)            :: NN,BB,NB,NpN,NpB,BpN,BpB  !! dot products of frame dependent vectors
  REAL(wp) :: Tq(3),NxB(3),TqTq,NTq,BTq,NpTq,BpTq
!===================================================================================================================================
  ASSOCIATE(T=>aux%T(:,iz),N=>aux%N(:,iz),B=>aux%B(:,iz),Np=>aux%Np(:,iz),Bp=>aux%Bp(:,iz), &
            NN=>aux%NN(iz),BB=>aux%BB(iz),NB=>aux%NB(iz),NpN=>aux%NpN(iz),NpB=>aux%NpB(iz), &
            BpN=>aux%BpN(iz),BpB=>aux%BpB(iz))
  Tq=(T+q1*Np+q2*Bp)
  NxB=CROSS(N,B)
  TqTq =SUM(Tq(:)*Tq(:))
  NTq  =SUM(N( :)*Tq(:))
  BTq  =SUM(B( :)*Tq(:))
  NpTq =SUM(Np(:)*Tq(:))
  BpTq =SUM(Bp(:)*Tq(:))
  Jh=SUM(Tq*NxB)
  Jh_dq1=SUM(Np*NxB)
  Jh_dq2=SUM(Bp*NxB)
  !! template for g_ab
  !! q1q1_ab = q1a*q1b
  !! q2q2_ab = q2a*q2b
  !! q3q3_ab = q3a*q3b
  !! q1q2_ab  = (q1a*q2b+q2a*q1b)
  !! q1q3_ab  = (q1a*q3b+q3a*q1b)
  !! q2q3_ab  = (q2a*q3b+q3a*q2b)
  !! g_ab     = NN*q1q1_ab + BB*q2q2_ab + TqTq*q3q3_ab + NB *q1q2_ab + NTq*q1q3_ab + BTq*q2q3_ab
  !! g_ab_dq1 =                    2.0_wp*NpTq*q3q3_ab               + NpN*q1q3_ab + NpB*q2q3_ab
  !! g_ab_dq2 =                    2.0_wp*BpTq*q3q3_ab               + BpN*q1q3_ab + BpB*q2q3_ab
  !q3t=dzeta/dtheta=0,q3z=dzeta/dzeta=1
  !g_tt
  !q1q1_tt = q1t*q1t
  !q2q2_tt = q2t*q2t
  !q3q3_tt = q3t*q3t =0
  !q1q2_tt  = (q1t*q2t+q2t*q1t)
  !q1q3_tt  = (q1t*q3t+q3t*q1t) =0
  !q2q3_tt  = (q2t*q3t+q3t*q2t) =0
  g_tt     = NN*q1t*q1t + BB*q2t*q2t + NB *(q1t*q2t+q2t*q1t)
  g_tt_dq1 = 0.
  g_tt_dq2 = 0.
  !g_tz
  !q1q1_tz = q1t*q1z
  !q2q2_tz = q2t*q2z
  !q3q3_tz = q3t*q3z =0
  !q1q2_tz  = (q1t*q2z+q2t*q1z)
  !q1q3_tz  = q1t  !(q1t*q3z+q3t*q1z)
  !q2q3_tz  = q2t  !(q2t*q3z+q3t*q2z)
  g_tz     = NN*q1t*q1z + BB*q2t*q2z + NB *(q1t*q2z+q2t*q1z)+ NTq*q1t + BTq*q2t
  g_tz_dq1 =                                                  NpN*q1t + NpB*q2t
  g_tz_dq2 =                                                  BpN*q1t + BpB*q2t
  !g_zz
  !q1q1_zz = q1z*q1z
  !q2q2_zz = q2z*q2z
  !q3q3_zz = q3z*q3z =1
  !q1q2_zz  = (q1z*q2z+q2z*q1z)
  !q1q3_zz  = (q1z+q1z)  !(q1z*q3z+q3z*q1z)
  !q2q3_zz  = (q2z+q2z)  !(q2z*q3z+q3z*q2z)
  g_zz     = NN*q1z*q1z + BB*q2z*q2z + TqTq + NB *(q1z*q2z+q2z*q1z) + NTq*(q1z+q1z) + BTq*(q2z+q2z)
  g_zz_dq1 =                   2.0_wp*(NpTq                         + NpN* q1z      + NpB* q2z      )
  g_zz_dq2 =                   2.0_wp*(BpTq                         + BpN* q1z      + BpB* q2z      )
  END ASSOCIATE
END SUBROUTINE eval_all_e



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
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp)
  Jh=SUM((T+q1*Np+q2*Bp)*CROSS(N,B))  ! Tq. (N x B)
  IF(Jh .LT. 1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_axisNB,  Jh<0",RealInfo=zeta*sf%nfp/TWOPI)
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
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp)
  Jh_dq1=SUM(Np*CROSS(N,B))  ! Tq. (N x B)
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
  REAL(wp),DIMENSION(3) :: X0,T,N,B,Np,Bp
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_TNB(zeta,X0,T,N,B,Np,Bp)
  Jh_dq2=SUM(Bp*CROSS(N,B))
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
!> evaluate curve X0(zeta), and T=X0',N,B,N',B', using the fourier series of X0_hat,N_hat and B_hat and transform from "hat"
!! coordinates to cartesian coordinates:
!! x=xhat*cos(zeta)-sgn_rot*yhat*sin(zeta), y=yhat*cos(zeta)+sgn_rot*sin(zeta)*xhat, z=zhat
!!
!===================================================================================================================================
SUBROUTINE hmap_axisNB_eval_TNB_hat( sf,zeta,X0,T,N,B,Np,Bp)
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
  REAL(wp)                      :: base_x(sf%fb_hat%modes)
  REAL(wp)                      :: base_dxdz(sf%fb_hat%modes)
  REAL(wp)                      :: cosz,sinz
  REAL(wp),DIMENSION(3)         :: X0_hat,T_hat,N_hat,B_hat,Np_hat,Bp_hat
!===================================================================================================================================
  base_x =sf%fb_hat%eval(           0,(/0.,zeta/))
  base_dxdz=sf%fb_hat%eval(DERIV_ZETA,(/0.,zeta/))

  !hat coordinates (in one field period)
  __MATVEC_N(X0_hat,sf%xyz_hat_modes( :,:),base_x)
  __MATVEC_N(N_hat ,sf%Nxyz_hat_modes(:,:),base_x)
  __MATVEC_N(B_hat ,sf%Bxyz_hat_modes(:,:),base_x)
  __MATVEC_N(T_hat ,sf%xyz_hat_modes( :,:),base_dxdz)
  __MATVEC_N(Np_hat,sf%Nxyz_hat_modes(:,:),base_dxdz)
  __MATVEC_N(Bp_hat,sf%Bxyz_hat_modes(:,:),base_dxdz)

  ! apply transform to x,y,z: x=xhat*cos(zeta)-sgn*yhat*sin(zeta), y=yhat*cos(zeta)+xhat*sgn*sin(zeta), z=zhat
  ! sgn_rot was used in the transform from x to xhat, so it has to be used here as well!!
  cosz=COS(zeta); sinz=SIN(zeta)
  X0=(/X0_hat(1)*cosz - sf%sgn_rot*X0_hat(2)*sinz, &
       X0_hat(2)*cosz + sf%sgn_rot*X0_hat(1)*sinz , &
       X0_hat(3)/)
  !dX0/dzeta:
  T= (/(T_hat(1)-sf%sgn_rot*X0_hat(2))*cosz - (sf%sgn_rot*T_hat(2)+X0_hat(1))*sinz, &
       (sf%sgn_rot*T_hat(1)-X0_hat(2))*sinz + (T_hat(2)+sf%sgn_rot*X0_hat(1))*cosz, &
       T_hat(3)/)
  !transform N to x,y,z
  N=(/N_hat(1)*cosz - sf%sgn_rot*N_hat(2)*sinz, &
      N_hat(2)*cosz + sf%sgn_rot*N_hat(1)*sinz, &
      N_hat(3)/)
  !dN/dzeta:
  Np= (/(Np_hat(1)-sf%sgn_rot*N_hat(2))*cosz - (sf%sgn_rot*Np_hat(2)+N_hat(1))*sinz, &
        (sf%sgn_rot*Np_hat(1)-N_hat(2))*sinz + (Np_hat(2)+sf%sgn_rot*N_hat(1))*cosz, &
         Np_hat(3)/)
  !transform B to x,y,z
  B=(/B_hat(1)*cosz - sf%sgn_rot*B_hat(2)*sinz, &
      B_hat(2)*cosz + sf%sgn_rot*B_hat(1)*sinz, &
      B_hat(3)/)
  !dB/dzeta:
  Bp= (/(Bp_hat(1)-sf%sgn_rot*B_hat(2))*cosz - (sf%sgn_rot*Bp_hat(2)+B_hat(1))*sinz, &
        (sf%sgn_rot*Bp_hat(1)-B_hat(2))*sinz + (Bp_hat(2)+sf%sgn_rot*B_hat(1))*cosz, &
         Bp_hat(3)/)

END SUBROUTINE hmap_axisNB_eval_TNB_hat


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
  INTEGER            :: iTest,idir,jdir,kdir,qdir,izeta,i,j,k,ndims(1:3),ijk(3)
  INTEGER,PARAMETER  :: nzeta=5
  INTEGER,PARAMETER  :: ns=2
  INTEGER,PARAMETER  :: nthet=3
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3),q_test(3,3),x_eps(3),dxdq(3),gij,gij_eps
  REAL(wp)           :: zeta(nzeta)
  REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: q1,q2,q1t,q2t,q1z,q2z, &
                                     Jh,g_tt,    g_tz,    g_zz,     &
                                     Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                     Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  REAL(wp),PARAMETER :: epsFD=1.0e-8
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
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_axisNB TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN
    ijk=(/1,4*sf%nzeta/7+1,2*sf%nzeta/3+1/)
    DO i=1,3
      izeta=ijk(i)
    !evaluate on the axis q1=q2=0
    iTest=100+i ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.0_wp, 0.0_wp, sf%zeta(izeta)/)
    x = sf%eval(q_in )
    checkreal=SQRT(SUM((x-sf%xyz(:,izeta))**2))
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|= ', checkreal, " ,izeta=",izeta
    END IF !TEST


    !evaluate at q1=0.01,q2=0. (= x+0.01*N)
    iTest=110+i ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.01_wp,0.0_wp, sf%zeta(izeta)/)
    x = sf%eval(q_in )
    checkreal=SQRT(SUM((x-(sf%xyz(:,izeta)+0.01_wp*sf%Nxyz(:,izeta)))**2))
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|= ', checkreal, " ,izeta=",izeta
    END IF !TEST

    !evaluate at q1=0.,q2=0.2 (= x+0.*N+0.2*B)
    iTest=120+i ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.0_wp,0.2_wp, sf%zeta(izeta)/)
    x = sf%eval(q_in )
    checkreal=SQRT(SUM((x-(sf%xyz(:,izeta)+0.2_wp*sf%Bxyz(:,izeta)))**2))
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|= ', checkreal, " ,izeta=",izeta
    END IF !TEST

    !evaluate at q1=-0.44,q2=0.33 (= x-0.44*N+0.33*B)
    iTest=130+i ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/-0.44_wp,0.33_wp, sf%zeta(izeta)/)
    x = sf%eval(q_in )
    checkreal=SQRT(SUM((x-(sf%xyz(:,izeta)-0.44_wp*sf%Nxyz(:,izeta)+0.33_wp*sf%Bxyz(:,izeta)))**2))
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|= ', checkreal, " ,izeta=",izeta
    END IF !TEST

    END DO ! izeta

    q_test(1,:)=(/1.0_wp, 0.0_wp, 0.0_wp/)
    q_test(2,:)=(/0.0_wp, 1.0_wp, 0.0_wp/)
    q_test(3,:)=(/0.0_wp, 0.0_wp, 1.0_wp/)
    q_in=(/0.1_wp, -0.15_wp, 0.335_wp*PI/)
    DO qdir=1,3
      !check dx/dq^i with FD
      iTest=140+qdir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
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
      iTest=150+idir+3*(jdir-1) ; IF(testdbg)WRITE(*,*)'iTest=',iTest
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
      iTest=160+idir+3*(jdir-1); IF(testdbg)WRITE(*,*)'iTest=',iTest
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
 IF (testlevel .GE. 2) THEN
    DO izeta=1,nzeta
      zeta(izeta)=0.333_wp+REAL(izeta-1,wp)/REAL(nzeta-1,wp)*0.221_wp
    END DO
    CALL sf%init_aux(nzeta,zeta)
    CALL sf%eval_aux()
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
      ALLOCATE(q1(ndims(1),ndims(2),ndims(3)))
      ALLOCATE(q2,q1t,q2t,q1z,q2z,Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
               mold=q1)
      !assign somewhat randomly
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        q1(i,j,k) = 0.11_wp -0.22_wp *REAL((i+j)*k,wp)/REAL((ndims(idir)+ndims(jdir))*ndims(kdir),wp)
        q2(i,j,k) = 0.15_wp -0.231_wp*REAL((i+k)*j,wp)/REAL((ndims(idir)+ndims(kdir))*ndims(jdir),wp)
        q1t(i,j,k)=-0.1_wp  +0.211_wp*REAL((i+2*j)*k,wp)/REAL((ndims(idir)+2*ndims(jdir))*ndims(kdir),wp)
        q2t(i,j,k)= 0.231_wp-0.116_wp*REAL((2*i+k)*j,wp)/REAL((2*ndims(idir)+ndims(kdir))*ndims(jdir),wp)
        q1z(i,j,k)=-0.024_wp+0.013_wp*REAL((3*i+2*j)*k,wp)/REAL((3*ndims(idir)+2*ndims(jdir))*ndims(kdir),wp)
        q2z(i,j,k)=-0.06_wp +0.031_wp*REAL((2*k+3*k)*i,wp)/REAL((2*ndims(kdir)+3*ndims(kdir))*ndims(idir),wp)
      END DO; END DO; END DO
      CALL hmap_axisNB_eval_all(sf,ndims,idir,q1,q2,q1t,q2t,q1z,q2z, &
           Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        ijk=(/i,j,k/)
        izeta=ijk(idir)
        Jh(i,j,k)       =Jh(i,j,k)       - sf%eval_Jh((/q1(i,j,k),q2(i,j,k),zeta(izeta)/))
        !g_tt(i,j,k)     =g_tt(i,j,k)     - sf%eval_ ()
        !g_tz(i,j,k)     =g_tz(i,j,k)     - sf%eval_ ()
        !g_zz(i,j,k)     =g_zz(i,j,k)     - sf%eval_ ()
        !Jh_dq1(i,j,k)   =Jh_dq1(i,j,k)   - sf%eval_ ()
        !g_tt_dq1(i,j,k) =g_tt_dq1(i,j,k) - sf%eval_ ()
        !g_tz_dq1(i,j,k) =g_tz_dq1(i,j,k) - sf%eval_ ()
        !g_zz_dq1(i,j,k) =g_zz_dq1(i,j,k) - sf%eval_ ()
        !Jh_dq2(i,j,k)   =Jh_dq2(i,j,k)   - sf%eval_ ()
        !g_tt_dq2(i,j,k) =g_tt_dq2(i,j,k) - sf%eval_ ()
        !g_tz_dq2(i,j,k) =g_tz_dq2(i,j,k) - sf%eval_ ()
        !g_zz_dq2(i,j,k) =g_zz_dq2(i,j,k) - sf%eval_ ()
      END DO; END DO; END DO
      iTest=200+idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal=SUM(ABS(Jh))/REAL(ns*nthet*nzeta,wp)
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
             '\n!! hmap_axisNB TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
      '\n =>  should be ', refreal,' : |sum(|Jh_all-eval_Jh(xall)|)|= ', checkreal, " ,idir=",idir
      END IF
      DEALLOCATE(q1,q2,q1t,q2t,q1z,q2z,Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2)
    END DO !idir
    CALL sf%free_aux()
 END IF

 test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_axisNB_test

END MODULE MODgvec_hmap_axisNB
