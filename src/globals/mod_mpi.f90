!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!>
!!# Module **MOD_MPI**
!!
!!  MPI related stuff, including communication
!!
!===================================================================================================================================
MODULE MODgvec_MPI
  USE_MPI
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: par_Init, par_Finalize, par_Barrier
  PUBLIC :: par_Bcast, par_Reduce, par_AllReduce
  PUBLIC :: par_IReduce, par_IBcast, par_Wait
  PUBLIC :: dType, req, req1, req2, req3

  !interfaces here

  ! allows same call for both scalar and different array ranks
  INTERFACE par_AllReduce
    MODULE PROCEDURE par_AllReduce_scalar,par_AllReduce_scalar_int, par_AllReduce_array1D, par_AllReduce_array2D
  END INTERFACE par_AllReduce

  ! allows same call for both scalar and different array ranks
  INTERFACE par_Reduce
    MODULE PROCEDURE par_Reduce_scalar, par_Reduce_scalar_int, par_Reduce_array1D, par_Reduce_array2D
  END INTERFACE par_Reduce

  ! allows same call for both scalar and different array ranks (nonblocking)
  INTERFACE par_IReduce
    MODULE PROCEDURE par_IReduce_array1D, par_IReduce_array2D
  END INTERFACE par_IReduce

  ! allows same call for both scalar and different array ranks
  INTERFACE par_Bcast
    MODULE PROCEDURE par_BCast_scalar, par_BCast_scalar_int, par_BCast_scalar_str, &
                     par_BCast_array1D, par_BCast_array1D_int, par_BCast_array1D_str, par_BCast_array2D
  END INTERFACE par_Bcast

  ! allows same call for both scalar and different array ranks (nonblocking)
  INTERFACE par_IBcast
    MODULE PROCEDURE par_IBCast_array1D, par_IBCast_array2D
  END INTERFACE par_IBcast

  ! allows same call for different scalar and array MPI request handles
  INTERFACE par_Wait
    MODULE PROCEDURE par_Wait, par_WaitAll
  END INTERFACE par_Wait

  ! Variables: whole module scope
  MPI_comm_TYPE :: worldComm
  MPI_datatype_TYPE :: dType
  MPI_request_TYPE, ALLOCATABLE :: req(:), req1(:), req2(:), req3(:)

CONTAINS

  !================================================================================================================================
  !> Initialization of MPI.
  !================================================================================================================================
  SUBROUTINE par_Init(comm_in)
  ! MODULES
    USE MODgvec_Globals, ONLY : MPIRoot,myRank,nRanks
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLE
  INTEGER, INTENT(IN),OPTIONAL :: comm_in
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: ierr,provided
  !================================================================================================================================
  ! BODY
#   if MPI
    IF(PRESENT(comm_in)) THEN
      worldComm = TRANSFER(comm_in,worldcomm)
    ELSE
      !CALL MPI_INIT(ierr)
      CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,provided,ierr)
      !CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      worldComm = MPI_COMM_WORLD
    END IF
    CALL MPI_COMM_SIZE(worldComm, nRanks, ierr)
    CALL MPI_COMM_RANK(worldComm, myRank, ierr)
    dType=MPI_DOUBLE_PRECISION
#   endif
    MPIRoot=(myRank.EQ.0)
    ALLOCATE(req (0:nRanks-1))
    ALLOCATE(req1(0:nRanks-1))
    ALLOCATE(req2(0:nRanks-1))
    ALLOCATE(req3(0:nRanks-1))
  END SUBROUTINE par_Init

  !================================================================================================================================
  !> Deinitialization of MPI.
  !================================================================================================================================
  SUBROUTINE par_Finalize()
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: ierr
  !================================================================================================================================
  ! BODY
    SDEALLOCATE(req )
    SDEALLOCATE(req1)
    SDEALLOCATE(req2)
    SDEALLOCATE(req3)
#   if MPI
    CALL MPI_FINALIZE(ierr)
#   endif
    !STOP !'END OF PROGRAM'
  END SUBROUTINE par_Finalize

  !================================================================================================================================
  !> Barrier for specified communicator, or world-communicator otherwise.
  !================================================================================================================================
  SUBROUTINE par_Barrier(Comm,beforeScreenOut,afterScreenOut)
  ! MODULES
    USE MODgvec_Globals, ONLY : MPIRoot,UNIT_Stdout
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    MPI_comm_TYPE, INTENT(IN), OPTIONAL :: Comm
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: beforeScreenOut,afterScreenOut
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: ierr
    MPI_comm_TYPE :: Communicator
  !================================================================================================================================
  ! BODY
#   if MPI
    IF (PRESENT(Comm)) THEN
      Communicator=Comm
    ELSE
      Communicator=worldComm
    END IF
    IF(PRESENT(beforeScreenOut).AND.MPIroot) WRITE(UNIT_StdOut,'(A)') beforeScreenOut
    CALL MPI_BARRIER(Communicator, ierr)
    IF(PRESENT(afterScreenOut) .AND.MPIroot) WRITE(UNIT_StdOut,'(A)') afterScreenOut
#   else
    IF(PRESENT(beforeScreenOut).AND.MPIroot) WRITE(UNIT_StdOut,'(A)') beforeScreenOut
    IF(PRESENT(afterScreenOut) .AND.MPIroot) WRITE(UNIT_StdOut,'(A)') afterScreenOut
#   endif
  END SUBROUTINE par_Barrier

  !================================================================================================================================
  !> Find MAX/MIN/SUM scalar value across MPI ranks and bradcast result back to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_AllReduce_scalar(scalar,parOP)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: scalar
    CHARACTER(LEN=3), INTENT(IN) :: parOP
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER     :: ierr
#   if MPI
    MPI_op_TYPE :: mpiOP
  !================================================================================================================================
  ! BODY
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    CALL MPI_AllReduce(MPI_IN_PLACE, scalar, 1, dType, mpiOP, worldComm, ierr)
#   endif
  END SUBROUTINE par_AllReduce_scalar

  !================================================================================================================================
  !> Find MAX/MIN/SUM scalar value across MPI ranks and bradcast result back to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_AllReduce_scalar_int(scalar_int,parOP)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER, INTENT(INOUT)       :: scalar_int
    CHARACTER(LEN=3), INTENT(IN) :: parOP
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER     :: ierr
#   if MPI
    MPI_op_TYPE :: mpiOP
  !================================================================================================================================
  ! BODY
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    CALL MPI_AllReduce(MPI_IN_PLACE, scalar_int, 1, MPI_INTEGER, mpiOP, worldComm, ierr)
#   endif
  END SUBROUTINE par_AllReduce_scalar_int

  !================================================================================================================================
  !> Find MAX/MIN/SUM of 1D array (assumed-shape) across all MPI ranks and bradcast result back to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_AllReduce_array1D(arr,parOP)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:)
    CHARACTER(LEN=3), INTENT(IN) :: parOP
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    MPI_op_TYPE :: mpiOP
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    CALL MPI_AllReduce(MPI_IN_PLACE, arr, sz, dType, mpiOP, worldComm, ierr)
#   endif
  END SUBROUTINE par_AllReduce_array1D

  !================================================================================================================================
  !> Find MAX/MIN/SUM of 2D array (assumed-shape) across all MPI ranks and bradcast result back to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_AllReduce_array2D(arr,parOP)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:,:)
    CHARACTER(LEN=3), INTENT(IN) :: parOP
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    MPI_op_TYPE :: mpiOP
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    CALL MPI_AllReduce(MPI_IN_PLACE, arr, sz, dType, mpiOP, worldComm, ierr)
#   endif
  END SUBROUTINE par_AllReduce_array2D

  !================================================================================================================================
  !> Find on MPI rank 'toRank' MAX/MIN/SUM scalar value across MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Reduce_scalar_int(scalar_int,parOP,toRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp,myRank
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER, INTENT(INOUT)       :: scalar_int
    CHARACTER(LEN=3), INTENT(IN) :: parOP
    INTEGER, INTENT(IN)          :: toRank  ! =0 by default
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER     :: ierr
#   if MPI
    MPI_op_TYPE :: mpiOP
  !================================================================================================================================
  ! BODY
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    IF (myRank.EQ.toRank) THEN
      CALL MPI_Reduce(MPI_IN_PLACE, scalar_int, 1, dType, mpiOP, toRank, worldComm, ierr)
    ELSE
      CALL MPI_Reduce(scalar_int, scalar_int, 1, dType, mpiOP, toRank, worldComm, ierr)
    END IF
#   endif
  END SUBROUTINE par_Reduce_scalar_int

  !================================================================================================================================
  !> Find on MPI rank 'toRank' MAX/MIN/SUM scalar value across MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Reduce_scalar(scalar,parOP,toRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp,myRank
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: scalar
    CHARACTER(LEN=3), INTENT(IN) :: parOP
    INTEGER, INTENT(IN)          :: toRank  ! =0 by default
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER     :: ierr
#   if MPI
    MPI_op_TYPE :: mpiOP
  !================================================================================================================================
  ! BODY
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    IF (myRank.EQ.toRank) THEN
      CALL MPI_Reduce(MPI_IN_PLACE, scalar, 1, dType, mpiOP, toRank, worldComm, ierr)
    ELSE
      CALL MPI_Reduce(scalar, scalar, 1, dType, mpiOP, toRank, worldComm, ierr)
    END IF
#   endif
  END SUBROUTINE par_Reduce_scalar

  !================================================================================================================================
  !> Find on MPI rank 'toRank' MAX/MIN/SUM of 1D array (assumed-shape) across all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Reduce_array1D(arr,parOP,toRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp,myRank
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:)
    CHARACTER(LEN=3), INTENT(IN) :: parOP
    INTEGER, INTENT(IN)          :: toRank  ! =0 by default
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    MPI_op_TYPE :: mpiOP
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    IF (myRank.EQ.toRank) THEN
      CALL MPI_Reduce(MPI_IN_PLACE, arr, sz, dType, mpiOP, toRank, worldComm, ierr)
    ELSE
      CALL MPI_Reduce(arr, arr, sz, dType, mpiOP, toRank, worldComm, ierr)
    END IF
#   endif
  END SUBROUTINE par_Reduce_array1D

  !================================================================================================================================
  !> Find on MPI rank 'toRank' MAX/MIN/SUM of 2D array (assumed-shape) across all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Reduce_array2D(arr,parOP,toRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp,myRank
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:,:)
    CHARACTER(LEN=3), INTENT(IN) :: parOP
    INTEGER, INTENT(IN)          :: toRank  ! =0 by default
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    MPI_op_TYPE :: mpiOP
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    IF (myRank.EQ.toRank) THEN
      CALL MPI_Reduce(MPI_IN_PLACE, arr, sz, dType, mpiOP, toRank, worldComm, ierr)
    ELSE
      CALL MPI_Reduce(arr, arr, sz, dType, mpiOP, toRank, worldComm, ierr)
    END IF
#   endif
  END SUBROUTINE par_Reduce_array2D

  !================================================================================================================================
  !> Find on MPI rank 'toRank' MAX/MIN/SUM of 1D array (assumed-shape) across all MPI ranks (nonblocking).
  !================================================================================================================================
  SUBROUTINE par_IReduce_array1D(arr,parOP,toRank,req_out)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp,myRank
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)        :: arr(:)
    CHARACTER(LEN=3), INTENT(IN)   :: parOP
    INTEGER, INTENT(IN)            :: toRank  ! =0 by default
    MPI_request_TYPE, INTENT(OUT) :: req_out
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    MPI_op_TYPE  :: mpiOP
    INTEGER      :: ierr
    INTEGER      :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    IF (myRank.EQ.toRank) THEN
      CALL MPI_IReduce(MPI_IN_PLACE, arr, sz, dType, mpiOP, toRank, worldComm, req_out, ierr)
    ELSE
      CALL MPI_IReduce(arr, arr, sz, dType, mpiOP, toRank, worldComm, req_out, ierr)
    END IF
#   endif
  END SUBROUTINE par_IReduce_array1D

  !================================================================================================================================
  !> Find on MPI rank 'toRank' MAX/MIN/SUM of 2D array (assumed-shape) across all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_IReduce_array2D(arr,parOP,toRank,req_out)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp,myRank
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)       :: arr(:,:)
    CHARACTER(LEN=3), INTENT(IN)  :: parOP
    INTEGER, INTENT(IN)           :: toRank  ! =0 by default
    MPI_request_TYPE, INTENT(OUT) :: req_out
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    MPI_op_TYPE :: mpiOP
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    SELECT CASE(parOP)
      CASE('MAX')
        mpiOP=MPI_MAX
      CASE('MIN')
        mpiOP=MPI_MIN
      CASE('SUM')
        mpiOP=MPI_SUM
    END SELECT
    IF (myRank.EQ.toRank) THEN
      CALL MPI_IReduce(MPI_IN_PLACE, arr, sz, dType, mpiOP, toRank, worldComm, req_out, ierr)
    ELSE
      CALL MPI_IReduce(arr, arr, sz, dType, mpiOP, toRank, worldComm, req_out, ierr)
    END IF
#   endif
  END SUBROUTINE par_IReduce_array2D

  !================================================================================================================================
  !> Broadcast a scalar from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_scalar_int(scalar_int,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER, INTENT(INOUT)       :: scalar_int
    INTEGER                      :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
  !================================================================================================================================
  ! BODY
    CALL MPI_Bcast(scalar_int, 1, MPI_INTEGER, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_scalar_int

  !================================================================================================================================
  !> Broadcast a scalar from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_scalar_str(scalar_str,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CHARACTER(LEN=*), INTENT(INOUT) :: scalar_str
    INTEGER                         :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
  !================================================================================================================================
  ! BODY
    CALL MPI_Bcast(scalar_str, LEN(scalar_str), MPI_CHARACTER, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_scalar_str

  !================================================================================================================================
  !> Broadcast a scalar from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_scalar(scalar,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: scalar
    INTEGER                      :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
  !================================================================================================================================
  ! BODY
    CALL MPI_Bcast(scalar, 1, dType, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_scalar

  !================================================================================================================================
  !> Broadcast a 1D array (assumed-shape) from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_array1D(arr,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:)
    INTEGER                      :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    CALL MPI_Bcast(arr, sz, dType, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_array1D

  !================================================================================================================================
  !> Broadcast a 1D array (assumed-shape) from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_array1D_int(arr_int,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER, INTENT(INOUT)       :: arr_int(:)
    INTEGER                      :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr_int)
    CALL MPI_Bcast(arr_int, sz, MPI_INTEGER, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_array1D_int

  !================================================================================================================================
  !> Broadcast a 1D array (assumed-shape) from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_array1D_str(arr_str,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CHARACTER(LEN=*), INTENT(INOUT) :: arr_str(:)
    INTEGER                         :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr_str)*LEN(arr_str(1))
    CALL MPI_Bcast(arr_str, sz, MPI_CHARACTER, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_array1D_str

  !================================================================================================================================
  !> Broadcast a 2D array (assumed-shape) from MPI rank 'fromRank' to all MPI ranks.
  !================================================================================================================================
  SUBROUTINE par_Bcast_array2D(arr,fromRank)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:,:)
    INTEGER                      :: fromRank
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    CALL MPI_Bcast(arr, sz, dType, fromRank, worldComm, ierr)
#   endif
  END SUBROUTINE par_Bcast_array2D

  !================================================================================================================================
  !> Broadcast a 1D array (assumed-shape) from MPI rank 'fromRank' to all MPI ranks (nonblocking)
  !================================================================================================================================
  SUBROUTINE par_IBcast_array1D(arr,fromRank,req_out)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:)
    INTEGER                      :: fromRank
    MPI_request_TYPE, INTENT(OUT) :: req_out
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    CALL MPI_IBcast(arr, sz, dType, fromRank, worldComm, req_out, ierr)
#   endif
  END SUBROUTINE par_IBcast_array1D

  !================================================================================================================================
  !> Broadcast a 2D array (assumed-shape) from MPI rank 'fromRank' to all MPI ranks (nonblocking)
  !================================================================================================================================
  SUBROUTINE par_IBcast_array2D(arr,fromRank,req_out)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    REAL(wp), INTENT(INOUT)      :: arr(:,:)
    INTEGER                      :: fromRank
    MPI_request_TYPE, INTENT(OUT) :: req_out
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: ierr
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(arr)
    CALL MPI_IBcast(arr, sz, dType, fromRank, worldComm, req_out, ierr)
#   endif
  END SUBROUTINE par_IBcast_array2D

  !================================================================================================================================
  !> Wait for completion of a single nonblocking communication
  !================================================================================================================================
  SUBROUTINE par_Wait(req_in)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  ! INPUT VARIABLES
    MPI_request_TYPE, INTENT(INOUT) :: req_in
#   if MPI
    CALL MPI_Wait(req_in, MPI_STATUS_IGNORE)
#   endif
  END SUBROUTINE par_Wait

  !================================================================================================================================
  !> Wait for completion of all nonblocking communications for req(:)
  !================================================================================================================================
  SUBROUTINE par_WaitAll(req_in)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    MPI_request_TYPE, INTENT(INOUT) :: req_in(:)
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER     :: sz
  !================================================================================================================================
  ! BODY
    sz=SIZE(req_in)
    CALL MPI_WaitAll(sz, req_in(:), MPI_STATUSES_IGNORE)
#   endif
  END SUBROUTINE par_WaitAll

  !================================================================================================================================
  !> Sum an array across MPI ranks: explicit-shape with implicit reshaping Multi-D->1D.
  !================================================================================================================================
  SUBROUTINE parSumArrayES(arr,sz)
  ! MODULES
    USE MODgvec_Globals, ONLY : wp
    IMPLICIT NONE
  !--------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER, INTENT(IN) :: sz
    REAL(wp), DIMENSION(sz), INTENT(INOUT) :: arr  !implicit array reshaping
  !--------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
#   if MPI
    INTEGER :: ierr
  !================================================================================================================================
  ! BODY
    CALL MPI_AllReduce(MPI_IN_PLACE, arr, sz, dType, MPI_SUM, worldComm, ierr)
#   endif
  END SUBROUTINE parSumArrayES


END MODULE MODgvec_MPI
