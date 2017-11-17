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
!!# Module ** hmap_cylindrical **
!!
!! contains the type that points to the routines of one chosen hmap_cylindrical
!!
!===================================================================================================================================
MODULE MOD_hmap_cylindrical
! MODULES
USE MOD_Globals, ONLY:wp,Unit_stdOut,abort
USE MOD_hmap,    ONLY:t_hmap
IMPLICIT NONE

PUBLIC
 

TYPE,EXTENDS(t_hmap) :: t_hmap_cylindrical
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  PROCEDURE :: init          => hmap_cylindrical_init
  PROCEDURE :: free          => hmap_cylindrical_free
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_hmap_cylindrical

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type hmap_cylindrical with number of elements
!!
!===================================================================================================================================
SUBROUTINE hmap_cylindrical_init( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_cylindrical), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: CYLINDRICAL ...'

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_cylindrical_test(sf)

END SUBROUTINE hmap_cylindrical_init


!===================================================================================================================================
!> finalize the type hmap_cylindrical
!!
!===================================================================================================================================
SUBROUTINE hmap_cylindrical_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_hmap_cylindrical), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN

  sf%initialized=.FALSE.

END SUBROUTINE hmap_cylindrical_free

!===================================================================================================================================
!> test hmap_cylindrical 
!!
!===================================================================================================================================
SUBROUTINE hmap_cylindrical_test( sf )
USE MOD_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_cylindrical), INTENT(IN   ) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest
  REAL(wp)           :: refreal,checkreal
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_cylindrical TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.LE.1)THEN

!    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

  END IF !testlevel <1

  test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_cylindrical_test

END MODULE MOD_hmap_cylindrical

