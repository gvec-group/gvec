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

!===================================================================================================================================
!>
!!# Module ** Restart Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_Restart_Vars
! MODULES
USE MOD_Globals,ONLY:wp
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
LOGICAL             :: doRestart          !! true if second command line argument exists
CHARACTER(LEN=255)  :: RestartFile        !! Name of Restart files (second command line argument)

!===================================================================================================================================
END MODULE MOD_Restart_Vars

