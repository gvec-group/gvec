!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
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
!!   CALL ncfile_init(nc,Filename,rwo_mode)
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
  CHARACTER(LEN=1)  :: rwo_mode
  CHARACTER(LEN=255) :: Filename
  CONTAINS
  PROCEDURE :: openfile       => ncfile_openfile
  PROCEDURE :: closefile      => ncfile_closefile
  PROCEDURE :: var_exists     => ncfile_var_exists
  PROCEDURE :: get_var_ndims  => ncfile_get_var_ndims
  PROCEDURE :: get_var_dims   => ncfile_get_var_dims
  PROCEDURE :: get_scalar     => ncfile_get_scalar
  PROCEDURE :: get_array      => ncfile_get_array
  PROCEDURE :: def_dim        => ncfile_def_dim
  PROCEDURE :: end_def_mode   => ncfile_end_def_mode
  PROCEDURE :: put_scalar     => ncfile_put_scalar
  PROCEDURE :: put_attr_char  => ncfile_put_attributes_char
  PROCEDURE :: put_char     => ncfile_put_char
  PROCEDURE :: put_array      => ncfile_put_array
  PROCEDURE :: enter_groups   => ncfile_enter_groups
  PROCEDURE :: handle_error   => ncfile_handle_error
  PROCEDURE :: free   => ncfile_free

END TYPE t_ncfile

CONTAINS

  SUBROUTINE mpi_check_single_access()
    USE MODgvec_Globals,ONLY: abort,MPIroot
    IMPLICIT NONE
    !===============================================================================================================================
#if MPI
    IF(.NOT. MPIroot) &
       CALL abort(__STAMP__,&
                  "netcdf routines are supposed to be called by MPIroot only")
#endif /*MPI*/
  END SUBROUTINE mpi_check_single_access

  !=================================================================================================================================
  !> allocate and initialize class and open/close the netcdf file and define  read ("r") or write ("w" includes read) mode
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_init(sf,FileName,rwo_mode)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: Filename
    CHARACTER(LEN=1),INTENT(IN) :: rwo_mode        !either read "r" or write "w" (existing file) or "o" createnew or overwrite
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===============================================================================================================================
    CALL mpi_check_single_access()
    ALLOCATE(t_ncfile :: sf)
    sf%isopen=.FALSE.
    sf%nc_id=0
    sf%filename=TRIM(FileName)
    sf%rwo_mode=rwo_mode
    CALL sf%openfile()

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
    CALL mpi_check_single_access()
    IF(sf%isopen) RETURN
#if NETCDF
    SELECT CASE(sf%rwo_mode)
    CASE("r")
      sf%ioError = nf90_OPEN(TRIM(sf%fileName), nf90_NOWRITE, sf%nc_id)
      CALL sf%handle_error("opening file '"//TRIM(sf%filename)//"' in read  mode",&
                           TypeInfo="FileNotFoundError")
    CASE("w")
      sf%ioError = nf90_OPEN(TRIM(sf%fileName), nf90_WRITE, sf%nc_id)
      CALL sf%handle_error("opening file '"//TRIM(sf%filename)//"' in write mode",&
                           TypeInfo="FileNotFoundError")
    CASE("o")
      sf%ioError = nf90_CREATE(TRIM(sf%fileName), NF90_64BIT_OFFSET, sf%nc_id)
      CALL sf%handle_error("creating or overwriting file '"//TRIM(sf%filename))
    END SELECT
    sf%isopen=.TRUE.
#else
  CALL abort(__STAMP__,&
      "cannot open netcdf file, since code is compiled with BUILD_NETCDF=OFF")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_openfile

  !=================================================================================================================================
  !> close netcdf file
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
    CALL mpi_check_single_access()
    IF(.NOT.sf%isopen) RETURN
#if NETCDF
    sf%ioError = nf90_CLOSE(sf%nc_id)
    CALL sf%handle_error("closing file")
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
    CALL mpi_check_single_access()
    IF(.NOT.sf%isopen) CALL sf%openfile()
    grpid=sf%nc_id
    varname=varname_in
    exists=.TRUE.
#if NETCDF
    id=INDEX(varname,"/")
    DO WHILE (id.NE.0)
      grpname=varname(1:id-1)
      varname=varname(id+1:)
      grpid_old=grpid
      sf%ioError = nf90_INQ_NCID(grpid_old, TRIM(grpname), grpid)
      exists=(sf%ioError .EQ. nf90_NOERR)
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
    CALL mpi_check_single_access()
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(exists)THEN
      sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid)
    END IF
    exists=(sf%ioError.EQ.nf90_NOERR)
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
    CALL mpi_check_single_access()
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
    CALL mpi_check_single_access()
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
    CALL mpi_check_single_access()
    CALL sf%enter_groups(varname_in,grpid,varname,exists)
#if NETCDF
    IF(.NOT.exists) CALL sf%handle_error("finding group in '"//TRIM(varname_in)//"'")
    sf%ioError = nf90_INQ_VARID(grpid, TRIM(varname), varid)
    CALL sf%handle_error("finding scalar variable '"//TRIM(varname_in)//"'")
    IF(PRESENT(intout))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, intout)
      CALL sf%handle_error("reading scalar variable '"//TRIM(varname_in)//"'")
      WRITE(UNIT_stdOut,'(6X,A,A30,A,I8)')'read netCDF scalar ','"'//TRIM(varname_in)//'"',' :: ',intout
    ELSEIF(PRESENT(realout))THEN
      sf%ioError = nf90_GET_VAR(grpid, varid, realout)
      CALL sf%handle_error("reading scalar variable '"//TRIM(varname_in)//"'")
      WRITE(UNIT_stdOut,'(6X,A,A30,A,E21.11)')'read netCDF scalar ','"'//TRIM(varname_in)//'"',' :: ',realout
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
    CHARACTER(LEN=100) :: fmtstr
    INTEGER :: grpid,varid,i,ndims_var,dim_ids(1:4),dims(1:4)
    INTEGER,ALLOCATABLE :: tmpint2d(:,:),tmpint3d(:,:,:),tmpint4d(:,:,:,:)
    REAL(wp),ALLOCATABLE :: tmpreal2d(:,:),tmpreal3d(:,:,:),tmpreal4d(:,:,:,:)
    LOGICAL :: exists,transpose
    !===============================================================================================================================
    CALL mpi_check_single_access()
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
    WRITE(fmtstr,'(A,I4,A)')'(6X,A,A30,A,"["',ndims_var,'(I4),"]")'
    WRITE(UNIT_stdOut,fmtstr)'read netCDF array ','"'//TRIM(varname_in)//'"',', shape: ',dims(1:ndims_var)
#endif /*NETCDF*/
  END SUBROUTINE ncfile_get_array


  !=================================================================================================================================
  !> define a dimension to the netCDF file
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_def_dim(sf,dimname_in,dimlen,dimid)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: dimname_in  !! name of the dimension
    INTEGER,INTENT(IN)          :: dimlen      !! length of the dimension
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT):: sf    !! self
    INTEGER,INTENT(OUT)          :: dimid !! id of the dimension
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF
    sf%ioError = nf90_def_dim(sf%nc_id, dimname_in, dimlen, dimid)
    CALL sf%handle_error("define dimension '"//TRIM(dimname_in)//"'")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_def_dim


  !=================================================================================================================================
  !> after creating a new file and making all definitions, one has to call end_def_mode
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_end_def_mode(sf)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT):: sf    !! self
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF
    sf%ioError = nf90_enddef(sf%nc_id)
    CALL sf%handle_error("finalize definition mode")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_end_def_mode

  !=================================================================================================================================
  !> define and put a scalar value to the netCDF file
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_put_scalar(sf,varname_in,def_put_mode,int_in,real_in)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the variable
    INTEGER,INTENT(IN) :: def_put_mode  !! 1:"def" or 2:"put" mode
    INTEGER,INTENT(IN),OPTIONAL :: int_in      !! scalar integer input
    REAL(wp),INTENT(IN),OPTIONAL :: real_in      !! scalar double input
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT):: sf    !! self
    ! LOCAL VARIABLES
    INTEGER    :: varid
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF
    IF(PRESENT(int_in))THEN
      SELECT CASE(def_put_mode)
      CASE(1) !def
        sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_INT,varid)
        CALL sf%handle_error("define scalar integer '"//TRIM(varname_in)//"'")
      CASE(2) !put
        sf%ioError = nf90_INQ_VARID(sf%nc_id, TRIM(varname_in), varid)
        CALL sf%handle_error("find varid of real array'"//TRIM(varname_in)//"'")
        sf%ioError = nf90_put_var(sf%nc_id, varid,int_in)
        CALL sf%handle_error("write scalar integer '"//TRIM(varname_in)//"'")
      END SELECT
    END IF
    IF(PRESENT(real_in))THEN
      SELECT CASE(def_put_mode)
      CASE(1) !def
        sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_DOUBLE,varid)
        CALL sf%handle_error("define scalar real '"//TRIM(varname_in)//"'")
      CASE(2) !put
        sf%ioError = nf90_INQ_VARID(sf%nc_id, TRIM(varname_in), varid)
        CALL sf%handle_error("find varid of real array'"//TRIM(varname_in)//"'")
        sf%ioError = nf90_put_var(sf%nc_id, varid,real_in)
        CALL sf%handle_error("write scalar real '"//TRIM(varname_in)//"'")
      END SELECT
    END IF
#else
  CALL abort(__STAMP__,&
      "cannot write scalar, BUILD_NETCDF=OFF")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_put_scalar

  SUBROUTINE ncfile_put_attributes_char(sf,varname_in,n_attr,attrs_names,attr_values)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in !! name of the variable
    INTEGER, INTENT(IN) :: n_attr           !! number of attributes
    CHARACTER(LEN=*),INTENT(IN) :: attrs_names(:)   !! Array of attribute names
    CHARACTER(LEN=*),INTENT(IN) :: attr_values(:)   !!  double input
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT):: sf    !! self
    ! LOCAL VARIABLES
    INTEGER :: varid !! ID of the variable
    INTEGER :: i !! iterable
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF
    DO i=1,n_attr
      sf%ioError = nf90_INQ_VARID(sf%nc_id, TRIM(varname_in), varid)
      CALL sf%handle_error("find varid during attribute write for '"//TRIM(varname_in)//"'")
      sf%ioError = nf90_put_att(sf%nc_id, varid, TRIM(attrs_names(i)), TRIM(attr_values(i)))
      CALL sf%handle_error("Putting attribute '"//TRIM(attrs_names(i))//"'")
    END DO
#else
  CALL abort(__STAMP__,&
      "cannot write array, BUILD_NETCDF=OFF")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_put_attributes_char

  !=================================================================================================================================
  !> define and put a char to the netCDF file as a string
  !! NOTE: This is only used for naming coordinate directions with a single character
  !=================================================================================================================================
  SUBROUTINE ncfile_put_char(sf,varname_in,dimid,def_put_mode,char_in)
        ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the variable
    INTEGER,INTENT(IN)          :: dimid !! ids of dimensions, must be created before by put_dim
    INTEGER,INTENT(IN)          :: def_put_mode  !! 1:"def" or 2:"put" mode
    CHARACTER(LEN=*),INTENT(IN) :: char_in   !!  double input
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT):: sf    !! self
    ! LOCAL VARIABLES
    INTEGER    :: varid
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF

    SELECT CASE(def_put_mode)
    CASE(1) !def
      sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_CHAR,dimid,varid)
      CALL sf%handle_error("define String'"//TRIM(varname_in)//"'")
    CASE(2) !put
      sf%ioError = nf90_INQ_VARID(sf%nc_id, TRIM(varname_in), varid)
      CALL sf%handle_error("find varid of string'"//TRIM(varname_in)//"'")
      sf%ioError = nf90_put_var(sf%nc_id, varid,char_in)
      CALL sf%handle_error("write string '"//TRIM(varname_in)//"'")
    END SELECT !CASE(def_put_mode)
#else
  CALL abort(__STAMP__,&
      "cannot write string, BUILD_NETCDF=OFF")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_put_char
  !=================================================================================================================================
  !> define and put an array value to the netCDF file
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_put_array(sf,varname_in,ndims_var,dims,dimids,def_put_mode,transpose_in, int_in,real_in)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in  !! name of the variable
    INTEGER,INTENT(IN)          :: ndims_var       !! number of dimensions
    INTEGER,INTENT(IN)          :: dims(1:ndims_var)   !! number of dimensions
    INTEGER,INTENT(IN)          :: dimids(1:ndims_var) !! ids of dimensions, must be created before by put_dim
    INTEGER,INTENT(IN) :: def_put_mode  !! 1:"def" or 2:"put" mode
    LOGICAL,INTENT(IN),OPTIONAL :: transpose_in !! transpose the data array, default is true, because of fortran ordering
    INTEGER,INTENT(IN),OPTIONAL  :: int_in(PRODUCT(dims))     !! integer input
    REAL(wp),INTENT(IN),OPTIONAL :: real_in(PRODUCT(dims))   !!  double input
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT):: sf    !! self
    ! LOCAL VARIABLES
    INTEGER    :: varid
    LOGICAL ::transpose
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF
    IF(PRESENT(transpose_in))THEN
      transpose=transpose_in
    ELSE
      transpose=.TRUE.
    END IF
    IF(.NOT.sf%isopen) CALL sf%openfile()
    IF(PRESENT(int_in))THEN
      SELECT CASE(def_put_mode)
      CASE(1) !def
        IF(transpose)THEN
          sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_INT,dimids(ndims_var:1:-1),varid)
        ELSE
          sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_INT,dimids(1:ndims_var),varid)
        END IF
        CALL sf%handle_error("define integer array'"//TRIM(varname_in)//"'")
      CASE(2) !put
        SELECT CASE(ndims_var)
        CASE(1)
          sf%ioError = nf90_put_var(sf%nc_id, varid,int_in)
        CASE(2)
          IF(transpose)THEN
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(int_in,dims(2:1:-1),order=[2,1]))
          ELSE
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(int_in,dims(1:2)))
          END IF
        CASE(3)
          IF(transpose)THEN
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(int_in,dims(3:1:-1),order=[3,2,1]))
          ELSE
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(int_in,dims(1:3)))
          END IF
        CASE(4)
          IF(transpose)THEN
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(int_in,dims(4:1:-1),order=[4,3,2,1]))
          ELSE
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(int_in,dims(1:4)))
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__,&
                     "ndims_var>4 not implemented yet in io_netcdf")
        END SELECT
        CALL sf%handle_error("write integer array '"//TRIM(varname_in)//"'")
      END SELECT !CASE(def_put_mode)
    END IF
    IF(PRESENT(real_in))THEN
      SELECT CASE(def_put_mode)
      CASE(1) !def
        IF(transpose)THEN
          sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_DOUBLE,dimids(ndims_var:1:-1),varid)
        ELSE
          sf%ioError = nf90_def_var(sf%nc_id, varname_in,NF90_DOUBLE,dimids(1:ndims_var),varid)
        END IF
        CALL sf%handle_error("define real array'"//TRIM(varname_in)//"'")
      CASE(2) !put
        sf%ioError = nf90_INQ_VARID(sf%nc_id, TRIM(varname_in), varid)
        CALL sf%handle_error("find varid of real array'"//TRIM(varname_in)//"'")
        SELECT CASE(ndims_var)
        CASE(1)
          sf%ioError = nf90_put_var(sf%nc_id, varid,real_in)
        CASE(2)
          IF(transpose)THEN
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(real_in,dims(2:1:-1),order=[2,1]))
          ELSE
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(real_in,dims(1:2)))
          END IF
        CASE(3)
          IF(transpose)THEN
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(real_in,dims(3:1:-1),order=[3,2,1]))
          ELSE
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(real_in,dims(1:3)))
          END IF
        CASE(4)
          IF(transpose)THEN
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(real_in,dims(4:1:-1),order=[4,3,2,1]))
          ELSE
            sf%ioError = nf90_put_var(sf%nc_id, varid,RESHAPE(real_in,dims(1:4)))
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__,&
                     "ndims_var>4 not implemented yet in io_netcdf")
        END SELECT
        CALL sf%handle_error("write real array '"//TRIM(varname_in)//"'")
      END SELECT !CASE(def_put_mode)
    END IF
#else
  CALL abort(__STAMP__,&
      "cannot write array, BUILD_NETCDF=OFF")
#endif /*NETCDF*/
  END SUBROUTINE ncfile_put_array

  !=================================================================================================================================
  !> netcdf error handling via sf%ioError variable
  !!
  !=================================================================================================================================
  SUBROUTINE ncfile_handle_error(sf,errmsg,TypeInfo)
    ! MODULES
    IMPLICIT NONE
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: errmsg
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: TypeInfo
    !-------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    CLASS(t_ncfile),INTENT(INOUT)        :: sf !! self
    CHARACTER(LEN=50)                    :: errtype
    !===============================================================================================================================
    CALL mpi_check_single_access()
#if NETCDF
    IF (sf%ioError .NE. nf90_NOERR) THEN
      IF(PRESENT(TypeInfo))THEN
        errtype=TRIM(TypeInfo)
      ELSE
        errtype="netCDF_error"
      END IF
       CALL abort(__STAMP__,&
                 "netCDF error: '"//TRIM(nf90_STRERROR(sf%ioError))//"' when "//TRIM(errmsg),&
                 IntInfo=sf%ioError,&
                 TypeInfo=errtype)
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
    CALL mpi_check_single_access()
    IF(sf%isopen) CALL sf%closefile()
    sf%nc_id=0
    sf%filename=""
    sf%rwo_mode=""
  END SUBROUTINE ncfile_free

END MODULE MODgvec_IO_NETCDF
