!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **Output to netcdf**
!!
!! Write structured visualization data to multidimensional arrays of a netcdf file
!!
!===================================================================================================================================
MODULE MODgvec_Output_netcdf
! MODULES
USE MODgvec_Globals, ONLY: wp
IMPLICIT NONE
PRIVATE

!INTERFACE WriteDataToNETCDF
!  MODULE PROCEDURE WriteDataToNETCDF
!END INTERFACE

PUBLIC::WriteDataToNETCDF
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write multidimensional data to netCDF format
!!
!===================================================================================================================================
SUBROUTINE WriteDataToNETCDF(dim1,vecdim,nVal,ndims,Dimnames,VarNames,Coord,Values,FileString,coord1,coord2,coord3, CoordNames, attr_values, attr_coords)
! MODULES
USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut
USE MODgvec_io_netcdf, ONLY:t_ncfile,ncfile_init
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    !! dimension of the data (either 1D,2D or 3D)
INTEGER,INTENT(IN)            :: vecdim                  !! dimension of coordinates
INTEGER,INTENT(IN)            :: nVal                    !! Number of nodal output variables
INTEGER,INTENT(IN)            :: ndims(1:dim1)           !! size of the data in each dimension
CHARACTER(LEN=*),INTENT(IN)   :: DimNames(1:dim1)        !! Names of dimensions of multi-dimensional array
REAL(wp),INTENT(IN)           :: Coord(vecdim,1:PRODUCT(ndims))      ! CoordinatesVector
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          !! Names of all variables that will be written out
REAL(wp),INTENT(IN)           :: Values(nVal,1:PRODUCT(ndims))   !! Statevector
CHARACTER(LEN=*),INTENT(IN)   :: FileString              !! Output file name (without .nc ending)
REAL(wp),INTENT(IN),OPTIONAL  :: coord1(:),coord2(:),coord3(:) !! Netcdf coordinate values e.g. rho, theta and zeta
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: CoordNames(1:dim1) !! Names of the dimensions
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: attr_values(nVal,2)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: attr_coords(dim1,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CLASS(t_ncfile),ALLOCATABLE  :: nc  !! container for netcdf-file
CHARACTER(LEN=255) :: tmpVarName,tmpVarNameY,tmpVarNameZ, coord_attr
CHARACTER(LEN=255) :: tmpCoordNames(1:dim1)
INTEGER :: def_put_mode,i,iVal,StrLen,dimids(1:dim1),vecdimid
LOGICAL            :: isVector,maybeVector
!===================================================================================================================================
#if NETCDF
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITE DATA TO NETCDF FILE "'//TRIM(FileString)//'.nc" ...'
#else
WRITE(UNIT_stdOut,'(A)')'   OMIT WRITING DATA TO NETCDF FILE "'//TRIM(FileString)//'.nc" (not compiled with netcdf)'
RETURN
#endif
CALL ncfile_init(nc,TRIM(fileString)//".nc","o")
CALL nc%def_dim("xyz",vecdim,vecdimid)
DO i=1,dim1
  CALL nc%def_dim(TRIM(DimNames(i)),ndims(i),dimids(i))
END DO
IF (PRESENT(CoordNames)) THEN
  tmpCoordNames = CoordNames
ELSE
  tmpCoordNames = DimNames
END IF
DO def_put_mode=1,2
  IF(def_put_mode.EQ.2) CALL nc%end_def_mode()
  CALL nc%put_array("pos",1+dim1,(/vecdim,ndims/),(/vecdimid,dimids/),def_put_mode,real_in=Coord)
  IF(def_put_mode.EQ.1) THEN
    coord_attr = "xyz"
    DO i=1,dim1
      coord_attr = TRIM(coord_attr)//" "//TRIM(tmpCoordNames(i))
    END DO
    CALL nc%put_attr_char("pos",1,(/"coordinates"/),(/coord_attr/))
    CALL nc%put_attr_char("pos",2,(/"long_name","symbol   "/),(/"position vector", "\mathbf{x}     "/))
  END IF
  !accout for vectors:
  ! if Variable Name ends with an X and the following have the same name with Y and Z
  ! then it forms a vector variable (X is omitted for the name)

  iVal=0 !scalars
  IF(def_put_mode.EQ.1) THEN
    coord_attr = ""
    DO i=1,dim1
      coord_attr = TRIM(coord_attr)//" "//TRIM(tmpCoordNames(i))
    END DO
  END IF
  DO WHILE(iVal.LT.nVal)
    iVal=iVal+1
    tmpVarName=TRIM(VarNames(iVal))
    StrLen=LEN(TRIM(tmpVarName))
    maybeVector=(iVal+vecdim-1.LE.nVal)
    isVector=.FALSE.
    IF(maybeVector)THEN
      SELECT CASE(vecdim)
      CASE(2)
        tmpVarNameY=TRIM(VarNames(iVal+1))
        isVector=((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0) &
                                  .AND.(INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0))
      CASE(3)
        tmpVarNameY=TRIM(VarNames(iVal+1))
        tmpVarNameZ=TRIM(VarNames(iVal+2))
        isVector=((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0) &
                                  .AND.(INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0) &
                                  .AND.(INDEX(tmpVarNameZ(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Z").NE.0))

      END SELECT
    END IF !maybevector

    IF(isvector)THEN !variable is a vector!
      tmpVarName=tmpVarName(:StrLen-1)
      CALL nc%put_array(TRIM(tmpVarName),1+dim1,(/vecdim,ndims/),(/vecdimid,dimids/),def_put_mode, &
                        real_in=Values(iVal:ival-1+vecdim,:))
      iVal=iVal+vecdim-1 !skip the Y (& Z) components
    ELSE
      CALL nc%put_array(TRIM(tmpVarName),dim1,ndims,dimids,def_put_mode,real_in=Values(iVal,:))
    END IF !isvector
    IF (def_put_mode.EQ.1) THEN
      CALL nc%put_attr_char(TRIM(tmpVarName),1,(/"coordinates"/),(/coord_attr/))
      IF (PRESENT(attr_values)) THEN
        CALL nc%put_attr_char(TRIM(tmpVarName),2,(/"long_name","symbol   "/),attr_values(iVal,:))
      END IF
    END IF
  END DO !iVal <=nVal
  IF (PRESENT(coord1)) THEN
    CALL nc%put_array(TRIM(tmpCoordNames(1)),1,(/ndims(1)/),(/dimids(1)/),def_put_mode,real_in=coord1)
  END IF
  IF (PRESENT(coord2)) THEN
    CALL nc%put_array(TRIM(tmpCoordNames(2)),1,(/ndims(2)/),(/dimids(2)/),def_put_mode,real_in=coord2)
  END IF
  IF (PRESENT(coord3)) THEN
    CALL nc%put_array(TRIM(tmpCoordNames(3)),1,(/ndims(3)/),(/dimids(3)/),def_put_mode,real_in=coord3)
  END IF
  IF (def_put_mode.EQ.1) THEN
    IF (PRESENT(attr_coords)) THEN
      DO i=1,dim1
        CALL nc%put_attr_char(TRIM(tmpCoordNames(i)),2,(/"long_name","symbol   "/),attr_coords(i,:))
      END DO
    END IF
  END IF
  CALL nc%put_char("xyz",vecdimid,def_put_mode,char_in="xyz")
  IF (def_put_mode.EQ.1) THEN
    CALL nc%put_attr_char("xyz",2,(/"long_name","symbol   "/),(/"cartesian vector components","\mathbf{x}                 "/))
  END IF
END DO !mode
CALL nc%free()
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToNETCDF

END MODULE MODgvec_Output_netcdf
