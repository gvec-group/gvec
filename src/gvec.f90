!===================================================================================================================================
! Copyright (c) 2017-2018 Florian Hindenlang
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


!===================================================================================================================================
!> 
!!# **GVEC** Driver program 
!!
!===================================================================================================================================
PROGRAM GVEC
USE MOD_Globals    ,ONLY: fmt_sep
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
!===================================================================================================================================

WRITE(*,fmt_sep)
WRITE(*,*)
WRITE(*,'(26X,A30,26X)')'  * * * * * * * * * * * * *   '
WRITE(*,'(26X,A30,26X)')' * * *                  * * * '
WRITE(*,'(26X,A30,26X)')'* * *    G  V  E  C       * * *'
WRITE(*,'(26X,A30,26X)')' * * *                  * * * '
WRITE(*,'(26X,A30,26X)')'  * * * * * * * * * * * * *   '
WRITE(*,*)



WRITE(*,fmt_sep)
WRITE(*,'(A)') ' GVEC SUCESSFULLY FINISHED'
WRITE(*,fmt_sep)

END PROGRAM GVEC


