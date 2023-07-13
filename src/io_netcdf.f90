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
!! Provides simplified read routines for netcdf files, via a the class "t_ncfile"
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
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)  :: sf !! self
    CHARACTER(LEN=255),INTENT(OUT) :: varname
    INTEGER,INTENT(OUT)            :: grpid
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
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
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
  !> get integer or real scalar (depends on optional argument)
  !! abort if variable does not exist. USE var_exists for checking
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_get_scalar( sf,varname_in,intout,realout) 
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
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
  SUBROUTINE ncfile_get_array( sf,varname_in,intout_1d,realout_1d, &
                                             intout_2d,realout_2d, &
                                             intout_3d,realout_3d, &
                                             intout_4d,realout_4d) 
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
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
    CHARACTER(LEN=255) :: varname
    INTEGER :: grpid,varid
    LOGICAL :: exists
    !===============================================================================================================================
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(.NOT.exists) CALL sf%handle_error("finding group in '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid) 
    CALL sf%handle_error("finding array '"//TRIM(varname_in)//"'")
    IF(PRESENT(intout_1d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout_1d)
    ELSEIF(PRESENT(realout_1d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout_1d)
    ELSEIF(PRESENT(intout_2d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout_2d)
    ELSEIF(PRESENT(realout_2d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout_2d)
    ELSEIF(PRESENT(intout_3d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout_3d)
    ELSEIF(PRESENT(realout_3d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout_3d)
    ELSEIF(PRESENT(intout_4d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout_4d)
    ELSEIF(PRESENT(realout_4d))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout_4d)
    END IF
    CALL sf%handle_error("reading array '"//TRIM(varname_in)//"'")
    SWRITE(UNIT_stdOut,'(6X,A,A50,A)')'read array  ',TRIM(varname_in)
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
