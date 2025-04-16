!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** sGrid **
!!
!! 1D grid in radial coordinate "s": Contains sgrid type definition and associated routines
!!
!===================================================================================================================================
MODULE MODgvec_sGrid
! MODULES
USE MODgvec_Globals,    ONLY : wp,Unit_stdOut,abort,MPIRoot
IMPLICIT NONE

PRIVATE
PUBLIC c_sgrid,t_sgrid
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES
TYPE, ABSTRACT :: c_sgrid
  CONTAINS
    PROCEDURE(i_sub_sgrid_init     ),DEFERRED :: init
    PROCEDURE(i_sub_sgrid_free     ),DEFERRED :: free
    PROCEDURE(i_sub_sgrid_copy     ),DEFERRED :: copy
    PROCEDURE(i_sub_sgrid_compare  ),DEFERRED :: compare
    PROCEDURE(i_fun_sgrid_find_elem),DEFERRED :: find_elem

END TYPE c_sgrid

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sgrid_init( sf , nElems_in, grid_type_in,sp_in)
    IMPORT wp,c_sgrid
    INTEGER       , INTENT(IN   ) :: nElems_in
    INTEGER       , INTENT(IN   ) :: grid_type_in
    REAL(wp),INTENT(IN),OPTIONAL  :: sp_in(0:nElems_in)
    CLASS(c_sgrid), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sgrid_init

  SUBROUTINE i_sub_sgrid_free( sf )
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sgrid_free

  SUBROUTINE i_sub_sgrid_copy( sf, tocopy )
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(INOUT) :: sf
    CLASS(c_sgrid), INTENT(IN   ) :: tocopy
  END SUBROUTINE i_sub_sgrid_copy

  SUBROUTINE i_sub_sgrid_compare( sf, tocompare, is_same )
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(IN   ) :: sf
    CLASS(c_sgrid), INTENT(IN   ) :: tocompare
    LOGICAL       , INTENT(  OUT) :: is_same
  END SUBROUTINE i_sub_sgrid_compare

  FUNCTION i_fun_sgrid_find_elem( sf ,x) RESULT(iElem)
    IMPORT wp,c_sgrid
    CLASS(c_sgrid), INTENT(IN   ) :: sf
    REAL(wp)      , INTENT(IN   ) :: x
    INTEGER                       :: iElem
  END FUNCTION i_fun_sgrid_find_elem

END INTERFACE


TYPE,EXTENDS(c_sgrid) :: t_sGrid
  LOGICAL :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: nElems                   !! global number of radial elements
  INTEGER              :: nElems_str, nElems_end   !! local number of radial elements per MPI subdomain  !<<<<
  INTEGER,ALLOCATABLE  :: offset_elem(:)           !! allocated  (0:nRanks), gives range on each rank:
                                                   !!   nElems_str:nElems_end=offset_elem(rank)+1:offset_elem(myRank+1)
  INTEGER              :: grid_type                !! type of grid (stretching functions...)
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),ALLOCATABLE :: sp(:)                    !! element point positions in [0,1], size(0:nElems)
  REAL(wp),ALLOCATABLE :: ds(:)                    !! ds(i)=sp(i)-sp(i-1), size(1:nElems)

  CONTAINS
  PROCEDURE :: init          => sGrid_init
  PROCEDURE :: copy          => sGrid_copy
  PROCEDURE :: compare       => sGrid_compare
  PROCEDURE :: free          => sGrid_free
  PROCEDURE :: find_elem     => sGrid_find_elem

END TYPE t_sGrid

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize the type sgrid with number of elements
!!
!===================================================================================================================================
SUBROUTINE sGrid_init( sf, nElems_in,grid_type_in,sp_in)
! MODULES
USE MODgvec_GLobals, ONLY: PI,myRank, nRanks
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: nElems_in       !! total number of elements
  INTEGER       , INTENT(IN   ) :: grid_type_in    !! GRID_TYPE_UNIFORM, GRID_TYPE_SQRT_S, GRID_TYPE_S2, GRID_TYPE_BUMP
  REAL(wp),INTENT(IN),OPTIONAL  :: sp_in(0:nElems_in) !! inner grid point positions, first position should be 0, last should be 1.
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iElem,iRank
!===================================================================================================================================
  IF(.NOT.test_called)THEN
    SWRITE(UNIT_stdOut,'(4X,A,I6,A,I3,A)')'INIT sGrid type nElems= ',nElems_in,' grid_type= ',grid_type_in, ' ...'
  END IF

  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sGrid type!'
    CALL sf%free()
  END IF

  sf%nElems     = nElems_in
  ! MPI PSEUDO-DOMAIN DECOMPOSITION (full domain is allocated, but only part of it is actually used on the MPI task)
  ALLOCATE(sf%offset_elem(0:nRanks))
  sf%offset_elem(0)=0
  DO iRank=0,nRanks-1
    sf%offset_elem(iRank+1)=(nElems_in*(iRank+1))/nRanks
  END DO
  sf%nElems_str = sf%offset_elem(myRank  ) +1
  sf%nElems_end = sf%offset_elem(myRank+1)
  IF (sf%nElems_end-sf%nElems_str+1 .EQ. 0) THEN
    CALL abort(__STAMP__, &
         'number of MPI tasks can not be larger than nElems!')
  END IF
  sf%grid_Type  = grid_type_in
  ALLOCATE(sf%sp(0:nElems_in))
  ALLOCATE(sf%ds(1:nElems_in))

  ASSOCIATE( &
              nElems    => sf%nElems    &
            , grid_Type => sf%grid_Type )

  IF(PRESENT(sp_in))THEN
    IF(SIZE(sp_in) .NE. nElems+1) THEN
      CALL abort(__STAMP__, &
          'sGrid_init: size of sp_in does not match nElems + 1!')
    END IF
    IF(ABS(sp_in(0)).GT.EPSILON(1.0_wp) .OR. &
       ABS(sp_in(nElems)-1.0_wp).GT.EPSILON(1.0_wp)) THEN
      CALL abort(__STAMP__, &
          'sGrid_init: sp_in(0) and sp_in(nElems) should be 0 and 1!')
    END IF
    sf%sp=sp_in
  ELSE
    !uniform [0,1]
    DO iElem=0,nElems
      sf%sp(iElem)=REAL(iElem,wp)/REAL(nElems,wp)
    END DO
    SELECT CASE(grid_type)
    CASE(GRID_TYPE_UNIFORM)
      !do nothing
    CASE(GRID_TYPE_SQRT_S) !finer at the edge
      sf%sp(:)=SQRT(sf%sp(:))
    CASE(GRID_TYPE_S2)   !finer at the center
      sf%sp(:)=sf%sp(:)*sf%sp(:)
    CASE(GRID_TYPE_BUMP) !finer towards axis and edge
      sf%sp(:)=sf%sp(:)-0.05_wp*SIN(PI*2.0_wp*sf%sp(:))
    CASE(GRID_TYPE_BUMP_EDGE) ! more equidistnat at the axis and finer towards edge
      sf%sp(:)=(sf%sp(:)-0.75_wp*((sf%sp(:)-0.4_wp)**3 + 0.4_wp**3))/(1.0_wp-0.75_wp*((1.0_wp-0.4_wp)**3+0.4**3))
    CASE DEFAULT
      CALL abort(__STAMP__, &
          'given grid type does not exist')
    END SELECT
  END IF!PRESENT(sp_in)

  !compute delta s
  DO iElem=1,nElems
    sf%ds(iElem)=sf%sp(iElem)-sf%sp(iElem-1)
  END DO

  END ASSOCIATE !sf

  sf%initialized=.TRUE.

  IF(.NOT.test_called)THEN
    SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
    CALL sGrid_test(sf)
  END IF

END SUBROUTINE sGrid_init


!===================================================================================================================================
!> finalize the type sgrid
!!
!===================================================================================================================================
SUBROUTINE sGrid_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN

  sf%nElems   = -1
  sf%grid_Type= -1

  SDEALLOCATE(sf%offset_elem)
  SDEALLOCATE(sf%sp)
  SDEALLOCATE(sf%ds)

  sf%initialized=.FALSE.

END SUBROUTINE sGrid_free

!===================================================================================================================================
!> copy the type sgrid, copies sf <= tocopy ... call sf%copy(tocopy)
!!
!===================================================================================================================================
SUBROUTINE sGrid_copy( sf , tocopy)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_sgrid), INTENT(IN) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPE IS(t_sgrid)
  IF(.NOT.tocopy%initialized) THEN
    CALL abort(__STAMP__, &
        "sgrid_copy: not initialized sgrid from which to copy!")
  END IF
  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sGrid copy!'
    CALL sf%free()
  END IF
  CALL sf%init(tocopy%nElems,tocopy%grid_type,tocopy%sp)

  END SELECT !TYPE
END SUBROUTINE sGrid_copy

!===================================================================================================================================
!> compare to sf grid with input grid to see if they are the same
!!
!===================================================================================================================================
SUBROUTINE sGrid_compare( sf , tocompare,is_same)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sgrid), INTENT(IN   ) :: sf !! self
  CLASS(c_sgrid), INTENT(IN   ) :: tocompare
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  LOGICAL       , INTENT(  OUT) :: is_same   !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL :: cond(2)
!===================================================================================================================================
  SELECT TYPE(tocompare); TYPE IS(t_sgrid)
  IF(.NOT.tocompare%initialized) THEN
    CALL abort(__STAMP__, &
        "sgrid_compare: tried to compare with a not initialized sgrid!")
  END IF
  cond(1)=(sf%nElems.EQ.tocompare%nElems)
  cond(2)=(sf%grid_type.EQ.tocompare%grid_type)

  is_same=ALL(cond)

  !IF(.NOT.is_same) WRITE(*,*)'DEBUG,grid is not same... nElems ',cond(1),', grid_type', cond(2)

  END SELECT !TYPE
END SUBROUTINE sGrid_compare


!===================================================================================================================================
!> find grid cell for certain position
!!
!===================================================================================================================================
FUNCTION sGrid_find_elem( sf , x) RESULT(iElem)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sgrid), INTENT(IN   ) :: sf !! self
  REAL(wp)      , INTENT(IN   ) :: x
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER                       :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: low,high
  REAL(wp) :: xloc
!===================================================================================================================================
  iElem=1
  xloc=MIN(1.0_wp,MAX(0.0_wp,x))
  IF(xloc.LT.(sf%sp(0)+sf%ds(1))) THEN
    iElem=1
    RETURN
  END IF
  IF(xloc.GE.(sf%sp(sf%nElems)-sf%ds(sf%nElems))) THEN
    iElem=sf%nElems
    RETURN
  END IF

  SELECT CASE(sf%grid_type)
  CASE(GRID_TYPE_UNIFORM)
    iElem=CEILING(xloc*sf%nElems)
    RETURN
  CASE(GRID_TYPE_S2)   !finer at the center
    iElem=CEILING(SQRT(xloc)*sf%nElems)
    RETURN
  CASE(GRID_TYPE_SQRT_S) !finer at the edge
    iElem=CEILING((xloc**2)*sf%nElems)
    RETURN
  END SELECT

  !not efficient, bisection of sp  array is better!!
  !DO jElem=2,sf%nElems
  !  IF((xloc.GE.sf%sp(iElem-1)).AND.(xloc.LT.sf%sp(jElem)))THEN
  !    iElem=jElem
  !    EXIT
  !  END IF
  !END DO

  !bisection
  low   = 1
  high  = sf%nElems-1
  iElem = (low + high) / 2 +1
  DO WHILE ( (xloc .LT.  sf%sp(iElem-1)) .OR. (xloc .GE. sf%sp(iElem)) )
     IF (xloc .LT. sf%sp(iElem-1)) THEN
       high = iElem-1
     ELSE
       low  = iElem
     END IF
     iElem = (low + high) / 2+1
  END DO



END FUNCTION sGrid_find_elem


!===================================================================================================================================
!> test sgrid variable
!!
!===================================================================================================================================
SUBROUTINE sGrid_test( sf )
! MODULES
USE MODgvec_GLobals, ONLY: UNIT_StdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,iElem,jElem
  REAL(wp)           :: x
  CHARACTER(LEN=10)  :: fail
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  TYPE(t_sgrid)      :: testgrid
  LOGICAL            :: check
!===================================================================================================================================
  test_called=.TRUE.
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN SGRID TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT.( (ABS(sf%sp(0)).LT. realtol) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),(A,E11.3))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n =>  should be 0.0 : sp(0) = ', sf%sp(0)
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT.( (MINVAL(sf%ds).GT.realtol))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),(A,E11.3))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be >0 : MINVAL(ds)',MINVAL(sf%ds)
    END IF !TEST

    iTest=103 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT. ( (ABS(sf%sp(sf%nElems)-1.0_wp).LT.realtol) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),(A,E11.3))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n =>  should be 1.0 : sp(nElems) = ', sf%sp(sf%nElems)
    END IF !TEST

    iTest=104 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT. ( (ABS(SUM(sf%ds(:))-1.0_wp).LT.realtol) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),(A,E11.3))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n =>  should be 1.0 : SUM(ds) = ', SUM(sf%ds)
    END IF !TEST

    iTest=105 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    iElem=sf%find_elem(0.0_wp)
    IF(testdbg.OR.(.NOT.( (iElem .EQ. 1 ) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),(A,I6))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n =>   should be 1 : findelem(0.0)= ', iElem
    END IF !TEST

    iTest=106 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    iElem=sf%find_elem(1.0_wp)
    IF(testdbg.OR.(.NOT.( (iElem .EQ. sf%nElems) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),2(A,I6))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be', sf%nElems,'  :  findelem(1.0)= ', iElem
    END IF !TEST

    iTest=107 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    jElem=(sf%nElems+1)/2
    x=0.5_wp*(sf%sp(jElem-1)+sf%sp(jElem))
    iElem=sf%find_elem(x)
    IF(testdbg.OR.(.NOT.( (iElem.EQ.jElem) )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),2(A,I6))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be ',jElem,': iElem= ' , iElem
    END IF !TEST

    iTest=108 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    jElem=MIN(3,sf%nElems)
    x=sf%sp(jElem-1)+0.99_wp*sf%ds(jElem)
    iElem=sf%find_elem(x)
    IF(testdbg.OR.(.NOT.( (iElem.EQ.jElem) )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),2(A,I6))') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be ',jElem,': iElem= ' , iElem
    END IF !TEST


    !get new grid and check compare
    iTest=121 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    IF(sf%grid_type.NE.-1) THEN
      CALL testgrid%init(sf%nElems+1,sf%grid_type)
    ELSE
      CALL testgrid%init(sf%nElems+1,MERGE(1,0,(sf%grid_type.EQ.0)))
    END IF
    CALL testgrid%compare(sf,check)
    CALL testgrid%free()
    IF(check)THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),A)') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be false'
    END IF !TEST

    !get new grid and check compare
    iTest=122 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testgrid%init(sf%nElems,MERGE(1,0,(sf%grid_type.EQ.0)))
    CALL testgrid%compare(sf,check)
    CALL testgrid%free()
    IF(check)THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),A)') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be false'
    END IF !TEST

    !get new grid and check compare
    iTest=123 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    IF(sf%grid_type.NE.-1) THEN
      CALL testgrid%init(sf%nElems,sf%grid_type)
    ELSE
      CALL testgrid%init(sf%nElems,sf%grid_type,sf%sp)
    END IF
    CALL testgrid%compare(sf,check)
    CALL testgrid%free()
    IF(.NOT.check)THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SGRID TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),A)') &
      '   nElems = ', sf%nElems , ' grid_type = ', sf%grid_type , &
      '\n => should be true'
    END IF !TEST

  END IF !testlevel>=1
  test_called=.FALSE.

END SUBROUTINE sGrid_test


END MODULE MODgvec_sGrid
