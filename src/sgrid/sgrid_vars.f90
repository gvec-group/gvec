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
!!# Module ** sGrid Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_sGrid_Vars
! MODULES
USE MOD_Globals    ,ONLY:wp
IMPLICIT NONE


PUBLIC  :: t_sgrid

PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 

TYPE, ABSTRACT :: c_sgrid
  LOGICAL :: initialized
  CONTAINS
    PROCEDURE(i_sub_sgrid_init    ),DEFERRED :: init
    PROCEDURE(i_sub_sgrid_copy    ),DEFERRED :: copy
    PROCEDURE(i_sub_sgrid_free    ),DEFERRED :: free

END TYPE c_sgrid

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sgrid_init( self , nElems, grid_type)
    IMPORT c_sgrid
    INTEGER       , INTENT(IN   ) :: nElems 
    INTEGER       , INTENT(IN   ) :: grid_type 
    CLASS(c_sgrid), INTENT(INOUT) :: self
  END SUBROUTINE i_sub_sgrid_init

  SUBROUTINE i_sub_sgrid_copy( self, tocopy ) 
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(INOUT) :: self
    CLASS(c_sgrid), INTENT(IN   ) :: tocopy
  END SUBROUTINE i_sub_sgrid_copy

  SUBROUTINE i_sub_sgrid_free( self ) 
    IMPORT c_sgrid
    CLASS(c_sgrid), INTENT(INOUT) :: self
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
SUBROUTINE sGrid_init( self, nElems,grid_type)
! MODULES
USE MOD_GLobals, ONLY: Pi,Unit_stdOut,fmt_sep,abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: nElems       !! total number of elements
  INTEGER       , INTENT(IN   ) :: grid_type    !! GRID_TYPE_UNIFORM, GRID_TYPE_SQRT_S, GRID_TYPE_S2, GRID_TYPE_BUMP
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iElem
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT sGrid type nElems= ',nElems,' grid_type= ',grid_type, ' ...'

  IF(self%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sGrid type!'
    CALL self%free() 
  END IF

  self%nElems   = nElems
  self%grid_Type= grid_type
  
  ALLOCATE(self%sp(0:self%nElems))
  ALLOCATE(self%ds( 1:self%nElems))
  ASSOCIATE(sp=>self%sp , ds=>self%ds)
  !uniform
  DO iElem=0,self%nElems
    sp(iElem)=REAL(iElem,wp)/REAL(self%nElems)
  END DO
  SELECT CASE(grid_type)
  CASE(GRID_TYPE_UNIFORM)
    !do nothing
  CASE(GRID_TYPE_SQRT_S) !finer at the edge
    sp(:)=SQRT(sp(:))
  CASE(GRID_TYPE_S2)   !finer at the center
    sp(:)=sp(:)*sp(:)
  CASE(GRID_TYPE_BUMP) !strechted towards axis and edge
    sp(:)=sp(:)-0.05_wp*SIN(Pi*(2.0_wp*sp(:)-1.0_wp))

  CASE DEFAULT
   CALL abort(__STAMP__, &
          'given grid type does not exist') 
  END SELECT
  
  !compute delta s
  DO iElem=1,self%nElems
    ds(iElem)=sp(iElem)-sp(iElem-1)
  END DO
  
  END ASSOCIATE !sp,ds
  
  self%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE sGrid_init


!===================================================================================================================================
!> copy the type sgrid, copies self <= tocopy ... call self%copy(tocopy)
!!
!===================================================================================================================================
SUBROUTINE sGrid_copy( self , tocopy)
! MODULES
USE MOD_GLobals, ONLY: Unit_stdOut,abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_sgrid), INTENT(IN) :: tocopy
  CLASS(t_sgrid), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPEIS(t_sgrid)
  IF(.NOT.tocopy%initialized) THEN
    CALL abort(__STAMP__, &
        "sgrid_copy: not initialized sgrid from which to copy!")
  END IF
  IF(self%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL self%free()
  END IF
  self%nElems    =  tocopy%nElems
  self%grid_type =  tocopy%grid_type
  ALLOCATE(self%sp(0:tocopy%nElems))
  ALLOCATE(self%ds( 1:tocopy%nElems))
  self%sp       =  tocopy%sp
  self%ds        =  tocopy%ds

  END SELECT !TYPE
END SUBROUTINE sGrid_copy


!===================================================================================================================================
!> finalize the type sgrid
!!
!===================================================================================================================================
SUBROUTINE sGrid_free( self )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sgrid), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  self%nElems   = -1
  self%grid_Type= -1
  
  SDEALLOCATE(self%sp)
  SDEALLOCATE(self%ds)
  
  self%initialized=.FALSE.

END SUBROUTINE sGrid_free

END MODULE MOD_sGrid_Vars

