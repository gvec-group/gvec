!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module ** sGrid **
!!
!! 1D grid in radial coordinate "s": Contains sgrid type definition and associated routines
!!
!===================================================================================================================================
MODULE MOD_sGrid
! MODULES
USE MOD_Globals    ,ONLY:wp
IMPLICIT NONE

PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
TYPE, ABSTRACT :: c_sgrid
  LOGICAL :: initialized
  CONTAINS
    PROCEDURE(i_sub_sgrid_init    ),DEFERRED :: init
    PROCEDURE(i_sub_sgrid_copy    ),DEFERRED :: copy
    PROCEDURE(i_sub_sgrid_free    ),DEFERRED :: free

END TYPE c_sgrid

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sgrid_init( sf , nElems_in, grid_type_in)
    IMPORT c_sgrid
    INTEGER       , INTENT(IN   ) :: nElems_in 
    INTEGER       , INTENT(IN   ) :: grid_type_in 
    CLASS(c_sgrid), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sgrid_init

  SUBROUTINE i_sub_sgrid_copy( sf, tocopy ) 
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(INOUT) :: sf
    CLASS(c_sgrid), INTENT(IN   ) :: tocopy
  END SUBROUTINE i_sub_sgrid_copy

  SUBROUTINE i_sub_sgrid_free( sf ) 
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sgrid_free

END INTERFACE
 

TYPE,EXTENDS(c_sgrid) :: t_sGrid
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: nElems                   !! total number of radial elements
  INTEGER              :: grid_type                !! type of grid (stretching functions...) 
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),ALLOCATABLE :: sp(:)                   !! element point positions in [0,1], size(0:nElems)
  REAL(wp),ALLOCATABLE :: ds(:)                    !! ds(i)=sp(i)-sp(i-1), size(1:nElems)

  CONTAINS
  PROCEDURE :: init          => sGrid_init
  PROCEDURE :: copy          => sGrid_copy
  PROCEDURE :: free          => sGrid_free

END TYPE t_sGrid


!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize the type sgrid with number of elements
!!
!===================================================================================================================================
SUBROUTINE sGrid_init( sf, nElems_in,grid_type_in)
! MODULES
USE MOD_GLobals, ONLY: PI,Unit_stdOut,abort
USE MOD_GLobals, ONLY: testlevel,nfailedMsg,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: nElems_in       !! total number of elements
  INTEGER       , INTENT(IN   ) :: grid_type_in    !! GRID_TYPE_UNIFORM, GRID_TYPE_SQRT_S, GRID_TYPE_S2, GRID_TYPE_BUMP
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iElem
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A,I6,A,I3,A)')'INIT sGrid type nElems= ',nElems_in,' grid_type= ',grid_type_in, ' ...'

  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sGrid type!'
    CALL sf%free() 
  END IF

  sf%nElems   = nElems_in
  sf%grid_Type= grid_type_in
  ALLOCATE(sf%sp(0:nElems_in))
  ALLOCATE(sf%ds(1:nElems_in))

  ASSOCIATE( &
              nElems    => sf%nElems    &
            , grid_Type => sf%grid_Type )
  
  !uniform
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
  CASE(GRID_TYPE_BUMP) !strechted towards axis and edge
    sf%sp(:)=sf%sp(:)-0.05_wp*SIN(PI*(2.0_wp*sf%sp(:)-1.0_wp))

  CASE DEFAULT
   CALL abort(__STAMP__, &
          'given grid type does not exist') 
  END SELECT
  
  !compute delta s
  DO iElem=1,nElems
    sf%ds(iElem)=sf%sp(iElem)-sf%sp(iElem-1)
  END DO
  
  END ASSOCIATE !sf
  
  sf%initialized=.TRUE.

  SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
  CALL sGrid_test(sf)

END SUBROUTINE sGrid_init

!===================================================================================================================================
!> test sgrid variable
!!
!===================================================================================================================================
SUBROUTINE sGrid_test( sf )
! MODULES
USE MOD_GLobals, ONLY: UNIT_StdOut,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iTest
!===================================================================================================================================
  IF(testlevel.LE.0) RETURN
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN TEST No.',nTestCalled,' SGRID    >>>>>>>>>'
  IF(testlevel.GT.0)THEN
    iTest=1
    IF( ABS(sf%sp(0)).GT. 1.0e-12_wp) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '!! TEST No.',nTestCalled ,': TEST ',iTest,' IN SGRID FAILED !!'
    END IF
    iTest=2
    IF( ABS(sf%sp(sf%nElems)-1.0_wp).GT.1.0e-12_wp) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '!! TEST No.',nTestCalled ,': TEST ',iTest,' IN SGRID FAILED !!'
    END IF
    iTest=3
    IF( ABS(SUM(sf%ds(:))-1.0_wp).GT.1.0e-12_wp) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '!! TEST No.',nTestCalled ,': TEST ',iTest,' IN SGRID FAILED !!'
    END IF
  END IF !testlevel>0

END SUBROUTINE sGrid_test

!===================================================================================================================================
!> copy the type sgrid, copies sf <= tocopy ... call sf%copy(tocopy)
!!
!===================================================================================================================================
SUBROUTINE sGrid_copy( sf , tocopy)
! MODULES
USE MOD_GLobals, ONLY: Unit_stdOut,abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_sgrid), INTENT(IN) :: tocopy
  CLASS(t_sgrid), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPEIS(t_sgrid)
  IF(.NOT.tocopy%initialized) THEN
    CALL abort(__STAMP__, &
        "sgrid_copy: not initialized sgrid from which to copy!")
  END IF
  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL sf%free()
  END IF
  CALL sf%init(tocopy%nElems,tocopy%grid_type) 

  END SELECT !TYPE
END SUBROUTINE sGrid_copy


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

  sf%nElems   = -1
  sf%grid_Type= -1
  
  SDEALLOCATE(sf%sp)
  SDEALLOCATE(sf%ds)
  
  sf%initialized=.FALSE.

END SUBROUTINE sGrid_free

END MODULE MOD_sGrid

