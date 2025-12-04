!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **Read in the boundary from a specifically defined NETCDF file **
!!
!!
!===================================================================================================================================
MODULE MODgvec_boundaryFromFile
! MODULES
USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut
USE MODgvec_io_netcdf, ONLY:t_ncfile
IMPLICIT NONE
PRIVATE

TYPE                 :: t_boundaryFromFile
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: nfp           !! number of field periods
  INTEGER              :: m_max         !! maximum number of fourier modes of the boundary in theta direction
  INTEGER              :: n_max         !! maximum number of fourier modes of the boundary in zeta direction (one field period only!)
  INTEGER              :: ntheta,nzeta  !! number of interpolation points in theta and zeta (one field period, half grid!!)
  INTEGER              :: lasym         !! =0: symmetric, =1: asymmetric fourier series
  REAL(wp),ALLOCATABLE :: theta(:)      !! theta positions [0,2pi), should be half grid!
  REAL(wp),ALLOCATABLE :: zeta(:)       !! zeta positions [0,2pi/nfp) should be on half grid!
  REAL(wp),ALLOCATABLE :: X(:,:) ,Y(:,:)!! boundary data X/Y positions X[i, j]=X(theta[i],zeta[j]),
                                        !!    Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1
  CLASS(t_ncfile),ALLOCATABLE  :: nc  !! container for netcdf-file
  CHARACTER(LEN=255)   :: ncfile=" " !! name of netcdf file with axis information
  CONTAINS

  PROCEDURE :: init       => bff_init
  PROCEDURE :: convert_to_modes => bff_convert_to_modes
  PROCEDURE :: free        => bff_free
END TYPE t_boundaryFromFile


INTERFACE boundaryFromFile_new
  MODULE PROCEDURE boundaryFromFile_new
END INTERFACE


PUBLIC::t_boundaryFromFile,boundaryFromFile_new
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Allocate class and call init
!!
!===================================================================================================================================
SUBROUTINE boundaryFromFile_new( sf,fileString)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)    , INTENT(IN   ) :: fileString
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_boundaryFromFile), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!===================================================================================================================================
ALLOCATE(t_boundaryFromFile :: sf)
CALL sf%init(fileString)
END SUBROUTINE boundaryFromFile_new


!===================================================================================================================================
!> initialize class: read file and save data to class structure
!!
!===================================================================================================================================
SUBROUTINE bff_init(sf,fileString)
  ! MODULES
  USE MODgvec_io_netcdf, ONLY:ncfile_init
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CHARACTER(LEN=*)    , INTENT(IN   ) :: fileString
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(t_boundaryFromFile), INTENT(INOUT) :: sf !! self
  !===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)')'   READ BOUNDARY FROM NETCDF FILE "'//TRIM(FileString)//'" ...'
  sf%ncfile=TRIM(FileString)
  CALL ncfile_init(sf%nc,sf%ncfile,"r")
  CALL READNETCDF(sf)
  sf%initialized=.TRUE.
END SUBROUTINE bff_init

!===================================================================================================================================
!> READ FROM SPECIFIC NETCDF FILE: general data and "boundary" group
!! ======= HEADER OF THE NETCDF FILE VERSION 3.1 ===================================================================================
!! === FILE DESCRIPTION:
!!   * axis, normal and binormal of the frame are given in cartesian coordinates along the curve parameter zeta [0,2pi].
!!   * The curve is allowed to have a field periodicity NFP, but the curve must be provided on a full turn.
!!   * The data is given in REAL SPACE, sampled along equidistant zeta point positions:
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
!!   * 'boundary/theta(:)'  : theta positions, 1D array of size 'boundary_ntheta', on half grid! theta(i)=(i+0.5)/ntheta*(2pi), i=0,...ntheta-1
!!   * 'boundary/zeta(:)'   : zeta positions, 1D array of size 'boundary/nzeta', for one field period! zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp), i=1,..nzeta, zeta[1] is arbitrary
!!   * 'boundary/X(::)',
!!     'boundary/Y(::)'     : boundary position X,Y in the N-B plane of the axis frame, in one field period, 2D array of size(ntheta, nzeta),  with
!!                               X[i, j]=X(theta[i],zeta[j])
!!                               Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1
!===================================================================================================================================
  SUBROUTINE ReadNETCDF(sf)
    USE MODgvec_io_netcdf
    IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT/OUTPUT VARIABLES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    CLASS(t_boundaryFromFile), INTENT(INOUT) :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  !===================================================================================================================================

  CALL sf%nc%get_scalar("NFP",intout=sf%nfp)
  CALL sf%nc%get_scalar("boundary/lasym",intout=sf%lasym)
  CALL sf%nc%get_scalar("boundary/ntheta",intout=sf%ntheta)
  CALL sf%nc%get_scalar("boundary/nzeta",intout=sf%nzeta)
   IF(sf%nc%var_exists("boundary/m_max"))THEN
    CALL sf%nc%get_scalar("boundary/m_max",intout=sf%m_max)
    sf%m_max=MIN(sf%m_max,(sf%ntheta-1)/2)  !maximum mode number based on number of interpolation points ntheta>=2*m_max+1
  ELSE
    sf%m_max=(sf%ntheta-1)/2  !maximum mode number based on number of interpolation points ntheta>=2*m_max+1
    WRITE(UNIT_stdOut,'(6X,A,I8)')'"boundary/m_max" not found, set to: ',sf%m_max
  END IF
  IF(sf%nc%var_exists("boundary/n_max"))THEN
    CALL sf%nc%get_scalar("boundary/n_max",intout=sf%n_max)
    sf%n_max = MIN(sf%n_max,(sf%nzeta-1)/2)  !maximum mode number based on number of interpolation points nzeta>=2*n_max+1
  ELSE
    sf%n_max=(sf%nzeta-1)/2   !maximum mode number based on number of interpolation points nzeta>=2*n_max+1
  END IF
  WRITE(UNIT_stdOut,'(6X,A,2I4)')'" boundary (m_max,n_max)" is set to: ',sf%m_max,sf%n_max
  ALLOCATE(sf%theta(sf%ntheta))
  CALL sf%nc%get_array("boundary/theta(:)",realout_1d=sf%theta)
  ALLOCATE(sf%zeta(sf%nzeta))
  CALL sf%nc%get_array("boundary/zeta(:)",realout_1d=sf%zeta)
  ALLOCATE(sf%X(sf%ntheta,sf%nzeta))
  CALL sf%nc%get_array("boundary/X(::)",realout_2d=sf%X)
  ALLOCATE(sf%Y(sf%ntheta,sf%nzeta))
  CALL sf%nc%get_array("boundary/Y(::)",realout_2d=sf%Y)
  CALL sf%nc%closefile()
END SUBROUTINE READNETCDF


!===================================================================================================================================
!> convert from interpolation points X=> X1_b, Y=> X2_b to fourier modes, given from the input fbase
!! convert to maximum allowable number of modes (ntheta>=2*m_max+1, nzeta>=2*n_max+1)
!! the final m_max/n_max can be smaller or larger. If larger, a change of base is necessary
!!
!===================================================================================================================================

SUBROUTINE bff_convert_to_modes(sf,x1_fbase_in,x2_fbase_in,X1_b,X2_b,scale_minor_radius)
! MODULES
  USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(t_fbase), INTENT(IN):: x1_fbase_in,x2_fbase_in
  REAL(wp), INTENT(IN)  :: scale_minor_radius
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp), INTENT(INOUT)  :: X1_b(x1_fbase_in%modes),X2_b(x2_fbase_in%modes)
  CLASS(t_boundaryFromFile), INTENT(INOUT) :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  TYPE(t_fBase),ALLOCATABLE        :: X_fbase,Y_fbase
  INTEGER                           :: i,nIP,mIP,mn_max_pts(2)
  REAL(wp)                          :: xn(2,sf%ntheta*sf%nzeta)
  REAL(wp),ALLOCATABLE              :: xydofs(:,:),X12dofs(:,:)
  !===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)')'  CONVERT BOUNDARY FROM POINTS TO MODES:'

  IF(sf%nfp.NE.x1_fbase_in%nfp)CALL abort(__STAMP__,&
                 " error in convert boundary to modes, nfp from GVEC parameterfile does not match to nfp of boundary file",&
                 TypeInfo="InvalidParameterError")
  mn_max_pts(1:2)=(/(sf%ntheta-1)/2,(sf%nzeta-1)/2/)
  WRITE(UNIT_stdOut,'(6X,2(A,I4),2A)')' from X,Y(ntheta= ',sf%ntheta,', nzeta= ', &
                                               sf%nzeta, '), lasym=',MERGE("symmetric","asymetric",(sf%lasym.EQ.0))
  WRITE(UNIT_stdOut,'(6x,2(3A,2(2I4,A)))') ' => X1 ',sin_cos_map(x1_fbase_in%sin_cos), &
                                                   ', (m_max,n_max)= (',mn_max_pts,')=>(',x1_fbase_in%mn_max,')', &
                                            ' , X2 ',sin_cos_map(x2_fbase_in%sin_cos), &
                                                  ', (m_max,n_max)= (',mn_max_pts,')=>(',x2_fbase_in%mn_max,')'
  IF(ALL(x1_fbase_in%mn_max.LE.mn_max_pts))THEN  !X1_base is smaller/equal
    CALL fbase_new( X_fbase, x1_fbase_in%mn_max,  (/sf%ntheta,sf%nzeta/), &
                    sf%nfp, sin_cos_map(x1_fbase_in%sin_cos), x1_fbase_in%exclude_mn_zero)
    X1_b = X_fbase%initDOF(RESHAPE(sf%X*scale_minor_radius,(/sf%ntheta*sf%nzeta/)) ,thet_zeta_start=(/sf%theta(1),sf%zeta(1)/))
  ELSE
    CALL fbase_new( X_fbase, mn_max_pts,  (/sf%ntheta,sf%nzeta/), &
                    sf%nfp, sin_cos_map(x1_fbase_in%sin_cos), x1_fbase_in%exclude_mn_zero)
    ALLOCATE(xydofs(1,1:X_fbase%modes),X12dofs(1,1:x1_fbase_in%modes))
    xydofs(1,:) = X_fbase%initDOF(RESHAPE(sf%X*scale_minor_radius,(/sf%ntheta*sf%nzeta/)),thet_zeta_start=(/sf%theta(1),sf%zeta(1)/))
    CALL x1_fbase_in%change_base(X_fbase,1,xydofs,X12dofs)
    X1_b=X12dofs(1,:)
    DEALLOCATE(xydofs,X12dofs)
  END IF
  IF(ALL(x2_fbase_in%mn_max.LE.mn_max_pts))THEN  !X2_base is smaller/equal
    CALL fbase_new( Y_fbase, x2_fbase_in%mn_max,  (/sf%ntheta,sf%nzeta/), &
                    sf%nfp,  sin_cos_map(x2_fbase_in%sin_cos),  x2_fbase_in%exclude_mn_zero)
    X2_b = Y_fbase%initDOF(RESHAPE(sf%Y*scale_minor_radius,(/sf%ntheta*sf%nzeta/)) ,thet_zeta_start=(/sf%theta(1),sf%zeta(1)/))
  ELSE
    CALL fbase_new( Y_fbase, mn_max_pts,  (/sf%ntheta,sf%nzeta/), &
                    sf%nfp,  sin_cos_map(x2_fbase_in%sin_cos),  x2_fbase_in%exclude_mn_zero)
    ALLOCATE(xydofs(1,1:Y_fbase%modes),X12dofs(1,1:x2_fbase_in%modes))
    xydofs(1,:) = Y_fbase%initDOF(RESHAPE(sf%Y*scale_minor_radius,(/sf%ntheta*sf%nzeta/)),thet_zeta_start=(/sf%theta(1),sf%zeta(1)/))
    CALL x2_fbase_in%change_base(Y_fbase,1,xydofs,X12dofs)
    X2_b=X12dofs(1,:)
  END IF
  !evaluate at interpolation points and check the error
  i=0
  DO nIP=1,sf%nzeta
    DO mIP=1,sf%ntheta
      i=i+1
      xn(1,i)=sf%theta(mIP)
      xn(2,i)=sf%zeta(nIP)
    END DO !m
  END DO !n

  WRITE(UNIT_stdOut,'(6X,A)')      ' => APPROXIMATION ERROR COMPARED TO INPUT POINTS:'
  WRITE(UNIT_stdOut,'(6X,A,E11.3)')'     max(|X1_fourier-X_input|)=',&
                      MAXVAL(ABS(x1_fbase_in%evalDOF_xn(sf%ntheta*sf%nzeta,xn,0,X1_b)/scale_minor_radius &
                             -RESHAPE(sf%X,(/sf%ntheta*sf%nzeta/))))
  WRITE(UNIT_stdOut,'(6X,A,E11.3)')'     max(|X2_fourier-Y_input|)', &
                      MAXVAL(ABS(x2_fbase_in%evalDOF_xn(sf%ntheta*sf%nzeta,xn,0,X2_b)/scale_minor_radius &
                             -RESHAPE(sf%Y,(/sf%ntheta*sf%nzeta/))))


  CALL X_fbase%free()
  CALL Y_fbase%free()
  DEALLOCATE(X_fbase,Y_fbase)
  WRITE(UNIT_stdOut,'(A)')'  ... CONVERT BOUNDARY DONE.'
END SUBROUTINE bff_convert_to_modes

!===================================================================================================================================
!! deallocate everything
!===================================================================================================================================
SUBROUTINE bff_free(sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_boundaryFromFile), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.sf%initialized)RETURN
SDEALLOCATE(sf%theta)
SDEALLOCATE(sf%zeta)
SDEALLOCATE(sf%X)
SDEALLOCATE(sf%Y)
IF(ALLOCATED(sf%nc))THEN
  CALL sf%nc%free()
  DEALLOCATE(sf%nc)
END IF
sf%initialized=.FALSE.
END SUBROUTINE bff_free

END MODULE MODgvec_boundaryFromFile
