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
!!# Module **sGrid**
!!
!! basis functions for radial coordinate: 
!! bsplines and Lagrange polynomials
!!
!===================================================================================================================================
MODULE MOD_sGrid
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitsGrid
  MODULE PROCEDURE InitsGrid
END INTERFACE

INTERFACE sGrid
  MODULE PROCEDURE sGrid
END INTERFACE

INTERFACE FinalizesGrid
  MODULE PROCEDURE FinalizesGrid
END INTERFACE

PUBLIC::InitsGrid
PUBLIC::sGrid
PUBLIC::FinalizesGrid
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitsGrid 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut
USE MOD_sGrid_Vars
USE MOD_ReadInTools,ONLY:GETLOGICAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'INIT sGrid MODULE ...'

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitsGrid


!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE sGrid()
! MODULES
USE MOD_Globals, ONLY:wp
USE MOD_sGrid_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
END SUBROUTINE sGrid 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizesGrid 
! MODULES
USE MOD_sGrid_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE FinalizesGrid

END MODULE MOD_sGrid
