!===================================================================================================================================
! Copyright (c) 2017 - 2018 Florian Hindenlang <hindenlang@gmail.com>
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/hopr)
!
! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
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
SUBROUTINE WriteDataToNETCDF(dim1,vecdim,nVal,ndims,Dimnames,VarNames,Coord,Values,FileString,coord1,coord2,coord3)
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CLASS(t_ncfile),ALLOCATABLE  :: nc  !! container for netcdf-file
CHARACTER(LEN=255) :: tmpVarName,tmpVarNameY,tmpVarNameZ
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
CALL nc%def_dim("vecdim",vecdim,vecdimid) 
DO i=1,dim1
  CALL nc%def_dim(TRIM(DimNames(i)),ndims(i),dimids(i)) 
END DO
DO def_put_mode=1,2
  IF(def_put_mode.EQ.2) CALL nc%end_def_mode()
  CALL nc%put_array("Coord",1+dim1,(/vecdim,ndims/),(/vecdimid,dimids/),def_put_mode,real_in=Coord)
  !accout for vectors: 
  ! if Variable Name ends with an X and the following have the same name with Y and Z 
  ! then it forms a vector variable (X is omitted for the name) 
  
  iVal=0 !scalars
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
  END DO !iVal <=nVal
  IF (PRESENT(coord1)) THEN
    CALL nc%put_array(TRIM(DimNames(1)),1,(/ndims(1)/),(/dimids(1)/),def_put_mode,real_in=coord1)
  END IF
  IF (PRESENT(coord2)) THEN
    CALL nc%put_array(TRIM(DimNames(2)),1,(/ndims(2)/),(/dimids(2)/),def_put_mode,real_in=coord2)
  END IF
  IF (PRESENT(coord3)) THEN
    CALL nc%put_array(TRIM(DimNames(3)),1,(/ndims(3)/),(/dimids(3)/),def_put_mode,real_in=coord3)
  END IF
END DO !mode
CALL nc%free()
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToNETCDF 

END MODULE MODgvec_Output_netcdf
