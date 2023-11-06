!===================================================================================================================================
! Copyright (c) 2023 Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **IO_NETCDF: SIMPLE NETCDF INTERFACE**
!!
!! Provides simplified read routines for netcdf files, via a class "t_ncfile"
!! start with defining a variable as 
!!  CLASS(t_ncfile),ALLOCATABLE  :: nc
!! and to allocate and initialize  
!!   CALL ncfile_init(nc,Filename,rw_mode)
!!
!!
!===================================================================================================================================
MODULE MODgvec_IO_NETCDF
USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut
#if NETCDF
USE netcdf
#endif /*NETCDF*/
IMPLICIT NONE

PUBLIC 

TYPE :: t_ncfile
  INTEGER  :: nc_id
  INTEGER  :: ioError
  LOGICAL  :: isopen
  CHARACTER(LEN=1)  :: rw_mode
  CHARACTER(LEN=255) :: Filename 
  CONTAINS
  PROCEDURE :: openfile       => ncfile_openfile
  PROCEDURE :: closefile      => ncfile_closefile
  PROCEDURE :: var_exists     => ncfile_var_exists
  PROCEDURE :: get_var_ndims  => ncfile_get_var_ndims
  PROCEDURE :: get_var_dims   => ncfile_get_var_dims
  PROCEDURE :: get_scalar     => ncfile_get_scalar
  PROCEDURE :: get_array      => ncfile_get_array
  PROCEDURE :: enter_groups   => ncfile_enter_groups
  PROCEDURE :: handle_error   => ncfile_handle_error
  PROCEDURE :: free   => ncfile_free

END TYPE t_ncfile

CONTAINS

  !=================================================================================================================================
  !> allocate and initialize class and open/close the netcdf file and define  read ("r") or write ("w" includes read) mode
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_init(sf,FileName,rw_mode) 
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: Filename
    CHARACTER(LEN=1),INTENT(IN) :: rw_mode        !either read "r" or write "w"
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===============================================================================================================================
    ALLOCATE(t_ncfile :: sf)
    sf%isopen=.FALSE.
    sf%nc_id=0
    sf%filename=TRIM(FileName)
    sf%rw_mode=rw_mode
    CALL sf%openfile() 
    CALL sf%closefile()

  END SUBROUTINE ncfile_init

  !=================================================================================================================================
  !> open netcdf file 
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_openfile( sf)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===============================================================================================================================
    IF(sf%isopen) RETURN
#if NETCDF
    SELECT CASE(sf%rw_mode)
    CASE("r")
      sf%ioError = nf90_OPEN(TRIM(sf%fileName), nf90_NOWRITE, sf%nc_id)
    CASE("w")
      sf%ioError = nf90_OPEN(TRIM(sf%fileName), nf90_WRITE, sf%nc_id)
    END SELECT
      CALL sf%handle_error("opening file '"//TRIM(sf%filename)//"' in '"//TRIM(sf%rw_mode)//"' mode")
    sf%isopen=.TRUE.
#else
  CALL abort(__STAMP__,&
      "cannot open netcdf file, since code is compiled with BUILD_NETCDF=OFF")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_openfile

  !=================================================================================================================================
  !> open netcdf file 
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_closefile( sf)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===============================================================================================================================
    IF(.NOT.sf%isopen) RETURN
#if NETCDF
    sf%ioError = nf90_CLOSE(sf%nc_id)
    CALL sf%handle_error("closing file ")
    sf%isopen=.FALSE.
#endif /*NETCDF*/
  END SUBROUTINE ncfile_closefile


  !=================================================================================================================================
  !> if variable name contains "/", these are interpreted as groups/subgroups.
  !> split the varname at first occurence of "/" to get the first group name on the file level. Then get the group id. 
  !>   repeat until no "/" is found anymore.
  !>   output the final groupid and the variable name without the group names.
  !=================================================================================================================================
  SUBROUTINE ncfile_enter_groups(sf,varname_in,grpid,varname,exists) 
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the variable (can include "/" for groups)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)  :: sf !! self
    CHARACTER(LEN=255),INTENT(OUT) :: varname  !! name of the variable without groups
    INTEGER,INTENT(OUT)            :: grpid    !! id of the last group found
    LOGICAL,INTENT(OUT)            :: exists 
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    CHARACTER(LEN=255) :: grpname
    INTEGER          :: grpid_old,id
    !===============================================================================================================================
    IF(.NOT.sf%isopen) CALL sf%openfile()
    grpid=sf%nc_id 
    varname=varname_in
#if NETCDF
    id=INDEX(varname,"/")
    DO WHILE (id.NE.0)
      grpname=varname(1:id-1)
      varname=varname(id+1:)
      grpid_old=grpid
      sf%ioError = nf90_INQ_NCID(grpid_old, TRIM(grpname), grpid) 
      exists=(sf%ioError .NE. nf90_NOERR)
      IF(.NOT.exists) RETURN
      id=INDEX(varname,"/")
    END DO
#endif /*NETCDF*/
  END SUBROUTINE ncfile_enter_groups

  !=================================================================================================================================
  !> check if variable name exists (also including groups separated with "/")
  !!
  !=================================================================================================================================
  FUNCTION ncfile_var_exists(sf,varname_in) RESULT(exists)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the variable (can include "/" for groups)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    LOGICAL                              :: exists  
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    CHARACTER(LEN=255) :: varname
    INTEGER :: grpid,varid
    !===============================================================================================================================
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(exists)THEN
      sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid) 
    END IF 
    exists=(sf%ioError.NE.nf90_NOERR)
#endif /*NETCDF*/
  END FUNCTION ncfile_var_exists

  !=================================================================================================================================
  !> get the number of dimensions of a variable
  !!
  !=================================================================================================================================
  FUNCTION ncfile_get_var_ndims(sf,varname_in) RESULT(ndims_out)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in !! name of the variable (can include "/" for groups)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    INTEGER                              :: ndims_out !0: scalar, 1: vector, 2: matrix...
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    CHARACTER(LEN=255) :: varname
    INTEGER :: grpid,varid
    LOGICAL :: exists
    !===============================================================================================================================
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(.NOT.exists) CALL sf%handle_error("finding group in '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid) 
    CALL sf%handle_error("finding of variable '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_inquire_variable(grpid,  varid, ndims=ndims_out) 
    CALL sf%handle_error("finding ndims of variable '"//TRIM(varname_in)//"'")
#endif /*NETCDF*/
  END FUNCTION ncfile_get_var_ndims


  !=================================================================================================================================
  !> get the size of a ulti-dimensional  array for all dimensions ndims
  !!
  !=================================================================================================================================
  FUNCTION ncfile_get_var_dims(sf,varname_in,ndims_in,transpose_in) RESULT(dims_out)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the multi-dimensional array (can include "/" for groups)
    INTEGER,INTENT(IN)          :: ndims_in    !! number of dimensions in the array
    LOGICAL,INTENT(IN),OPTIONAL :: transpose_in !! transpose the data array, default is true, because of fortran ordering
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    INTEGER                              :: dims_out(ndims_in)  !! size of each dimension of the array
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    CHARACTER(LEN=255) :: varname,dimname
    INTEGER :: grpid,varid,ndims_var,dim_ids(1:ndims_in),i
    LOGICAL :: exists,transpose
    !===============================================================================================================================
    IF(PRESENT(transpose_in))THEN
      transpose=transpose_in
    ELSE
      transpose=.TRUE.
    END IF 
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(.NOT.exists) CALL sf%handle_error("finding group in '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid) 
    CALL sf%handle_error("finding of variable '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_inquire_variable(grpid,  varid, ndims=ndims_var) 
    CALL sf%handle_error("finding ndims & dimids of variable '"//TRIM(varname_in)//"'")
    IF(ndims_var.NE.ndims_in) & 
      CALL sf%handle_error("ndims_in not correct for variable '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_inquire_variable(grpid,  varid, dimids=dim_ids) 
    DO i=1,ndims_var
      sf%ioError = nf90_inquire_dimension(grpid, dim_ids(i),name=dimname, len=dims_out(i))
      CALL sf%handle_error("finding size of dimension  '"//TRIM(dimname)//"'")
    END DO
    IF(transpose) dims_out=dims_out(ndims_var:1:-1)
#endif /*NETCDF*/
  END FUNCTION ncfile_get_var_dims


  !=================================================================================================================================
  !> get integer or real scalar (depends on optional argument)
  !! abort if variable does not exist. USE var_exists for checking
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_get_scalar( sf,varname_in,intout,realout) 
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in !! name of the variable (can include "/" for groups)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    INTEGER ,INTENT(OUT),OPTIONAL        :: intout   !! choose for integer out
    REAL(wp),INTENT(OUT),OPTIONAL        :: realout  !! choose for real(wp) out (double)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    CHARACTER(LEN=255) :: varname
    INTEGER :: grpid,varid
    LOGICAL :: exists
    !===============================================================================================================================
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(.NOT.exists) CALL sf%handle_error("finding group in '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid) 
    CALL sf%handle_error("finding scalar variable '"//TRIM(varname_in)//"'")
    IF(PRESENT(intout))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout)
      CALL sf%handle_error("reading scalar variable '"//TRIM(varname_in)//"'")
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,I8)')"read scalar  ",TRIM(varname_in),' :: ',intout
    ELSEIF(PRESENT(realout))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout)
      CALL sf%handle_error("reading scalar variable '"//TRIM(varname_in)//"'")
      SWRITE(UNIT_stdOut,'(6X,A,A50,A,E21.11)')'read scalar  ',TRIM(varname_in),' :: ',realout
    END IF
#endif /*NETCDF*/
  END SUBROUTINE ncfile_get_scalar

  !=================================================================================================================================
  !> get integer or real array of dimension 1d,2d,3d,4d (depends on optional argument)
  !> netcdf call get_var knows type and dimensions directly from argument
  !! abort if variable does not exist. USE var_exists for checking
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_get_array( sf,varname_in,transpose_in, &
                                             intout_1d,realout_1d, &
                                             intout_2d,realout_2d, &
                                             intout_3d,realout_3d, &
                                             intout_4d,realout_4d) 
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the variable (can include "/" for groups)
    LOGICAL,INTENT(IN),OPTIONAL :: transpose_in !! transpose the data array, default is true, because of fortran ordering
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    INTEGER ,INTENT(OUT),OPTIONAL        :: intout_1d(:)         !! choose for integer out 1d array
    REAL(wp),INTENT(OUT),OPTIONAL        :: realout_1d(:)        !! choose for real(wp) out (double)  1d array
    INTEGER ,INTENT(OUT),OPTIONAL        :: intout_2d(:,:)       !! choose for integer out 2d array
    REAL(wp),INTENT(OUT),OPTIONAL        :: realout_2d(:,:)      !! choose for real(wp) out (double) 2d array 
    INTEGER ,INTENT(OUT),OPTIONAL        :: intout_3d(:,:,:)     !! choose for integer out 3d array
    REAL(wp),INTENT(OUT),OPTIONAL        :: realout_3d(:,:,:)    !! choose for real(wp) out (double)  3d array
    INTEGER ,INTENT(OUT),OPTIONAL        :: intout_4d(:,:,:,:)   !! choose for integer out 4d array
    REAL(wp),INTENT(OUT),OPTIONAL        :: realout_4d(:,:,:,:)  !! choose for real(wp) out (double) 4darray
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    CHARACTER(LEN=255) :: varname,dimname
    INTEGER :: grpid,varid,i,ndims_var,dim_ids(1:4),dims(1:4)
    INTEGER,ALLOCATABLE :: tmpint2d(:,:),tmpint3d(:,:,:),tmpint4d(:,:,:,:)
    REAL(wp),ALLOCATABLE :: tmpreal2d(:,:),tmpreal3d(:,:,:),tmpreal4d(:,:,:,:)
    LOGICAL :: exists,transpose
    !===============================================================================================================================
    IF(PRESENT(transpose_in))THEN
      transpose=transpose_in
    ELSE
      transpose=.TRUE.
    END IF 
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(.NOT.exists) CALL sf%handle_error("finding group in '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid) 
    CALL sf%handle_error("finding array '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_inquire_variable(grpid,  varid, ndims=ndims_var) 
    CALL sf%handle_error("finding ndims & dimids of variable '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_inquire_variable(grpid,  varid, dimids=dim_ids(1:ndims_var)) 
    DO i=1,ndims_var
      sf%ioError = nf90_inquire_dimension(grpid, dim_ids(i),name=dimname, len=dims(i))
      CALL sf%handle_error("finding size of dimension  '"//TRIM(dimname)//"'")
    END DO
    IF(transpose) dims(1:ndims_var)=dims(ndims_var:1:-1)
    IF(PRESENT(intout_1d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout_1d)
    ELSEIF(PRESENT(realout_1d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout_1d)
    ELSEIF(PRESENT(intout_2d))THEN
      IF(transpose)THEN
        ALLOCATE(tmpint2d(dims(2),dims(1)))
        sf%ioError = nf90_GET_VAR(grpid, varid, tmpint2d)
        intout_2d=RESHAPE(tmpint2d,shape(intout_2d),order=[2,1])
        DEALLOCATE(tmpint2d)
      ELSE
        sf%ioError = nf90_GET_VAR(grpid, varid, intout_2d)
      END IF
    ELSEIF(PRESENT(realout_2d))THEN
      IF(transpose)THEN
        ALLOCATE(tmpreal2d(dims(2),dims(1)))
        sf%ioError = nf90_GET_VAR(grpid, varid, tmpreal2d)
        realout_2d=RESHAPE(tmpreal2d,shape(realout_2d),order=[2,1])
        DEALLOCATE(tmpreal2d)
      ELSE
        sf%ioError = nf90_GET_VAR(grpid, varid, realout_2d)
      END IF
    ELSEIF(PRESENT(intout_3d))THEN
      IF(transpose)THEN
        ALLOCATE(tmpint3d(dims(3),dims(2),dims(1)))
        sf%ioError = nf90_GET_VAR(grpid, varid, tmpint3d)
        intout_3d=RESHAPE(tmpint3d,shape(intout_3d),order=[3,2,1])
        DEALLOCATE(tmpint3d)
      ELSE
        sf%ioError = nf90_GET_VAR(grpid, varid, intout_3d)
      END IF
    ELSEIF(PRESENT(realout_3d))THEN
      IF(transpose)THEN
        ALLOCATE(tmpreal3d(dims(3),dims(2),dims(1)))
        sf%ioError = nf90_GET_VAR(grpid, varid, tmpreal3d)
        realout_3d=RESHAPE(tmpreal3d,shape(realout_3d),order=[3,2,1])
        DEALLOCATE(tmpreal3d)
      ELSE
        sf%ioError = nf90_GET_VAR(grpid, varid, realout_3d)
      END IF
    ELSEIF(PRESENT(intout_4d))THEN
      IF(transpose)THEN
        ALLOCATE(tmpint4d(dims(4),dims(3),dims(2),dims(1)))
        sf%ioError = nf90_GET_VAR(grpid, varid, tmpint4d)
        intout_4d=RESHAPE(tmpint4d,shape(intout_4d),order=[4,3,2,1])
        DEALLOCATE(tmpint4d)
      ELSE
        sf%ioError = nf90_GET_VAR(grpid, varid, intout_4d)
      END IF
    ELSEIF(PRESENT(realout_4d))THEN
      IF(transpose)THEN
        ALLOCATE(tmpreal4d(dims(4),dims(3),dims(2),dims(1)))
        sf%ioError = nf90_GET_VAR(grpid, varid, tmpreal4d)
        realout_4d=RESHAPE(tmpreal4d,shape(realout_4d),order=[4,3,2,1])
        DEALLOCATE(tmpreal4d)
      ELSE
        sf%ioError = nf90_GET_VAR(grpid, varid, realout_4d)
      END IF
    END IF
    CALL sf%handle_error("reading array '"//TRIM(varname_in)//"'")
    SWRITE(UNIT_stdOut,'(6X,A,A50,A,*(I4,:,","))')'read array  ',TRIM(varname_in),dims(ndims_var)
#endif /*NETCDF*/
  END SUBROUTINE ncfile_get_array

  !=================================================================================================================================
  !> netcdf error handling via sf%ioError variable
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_handle_error(sf,errmsg)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: errmsg
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    !===============================================================================================================================
#if NETCDF
    IF (sf%ioError .ne. nf90_NOERR) THEN
       SWRITE(UNIT_stdOut,'(6X,A)')"A netCDF error has occurred:  "//TRIM(errmsg)
       CALL abort(__STAMP__,&
                 nf90_STRERROR(sf%ioError))
    END IF
#endif
  END SUBROUTINE ncfile_handle_error

  !=================================================================================================================================
  !> closes file and frees variable 
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_free(sf)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile), INTENT(INOUT)        :: sf !! self
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===============================================================================================================================
    IF(sf%isopen) CALL sf%closefile()
    sf%nc_id=0
    sf%filename=""
    sf%rw_mode=""
  END SUBROUTINE ncfile_free

END MODULE MODgvec_IO_NETCDF
