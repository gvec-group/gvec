!===================================================================================================================================
! Copyright (C) 2017 - 2019  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module ** Transform SFL Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_Transform_SFL_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
USE MODgvec_base,   ONLY: t_base
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER                     :: whichSFLcoord !! 
CLASS(t_base),  ALLOCATABLE :: X1sfl_base    !! container for base of variable X1 in SFL coordinates
CLASS(t_base),  ALLOCATABLE :: X2sfl_base    !! container for base of variable X2 in SFL coordinates
CLASS(t_base),  ALLOCATABLE :: GZ_base       !! container for base of variable Gzeta (transforms to BOOZER!)
CLASS(t_base),  ALLOCATABLE :: GZsfl_base    !! container for base of variable Gzeta in SFL coordinates
REAL(wp),       ALLOCATABLE :: X1sfl(:,:)    !! data (1:nBase,1:modes) of X1 in SFL coords.
REAL(wp),       ALLOCATABLE :: X2sfl(:,:)    !! data (1:nBase,1:modes) of X2 in SFL coords.
REAL(wp),       ALLOCATABLE :: GZ(:,:)       !! data (1:nBase,1:modes) of GZ in GVEC coords. (for BOOZER)
REAL(wp),       ALLOCATABLE :: GZsfl(:,:)    !! data (1:nBase,1:modes) of GZ in SFL coords.  (for BOOZER)

!===================================================================================================================================
END MODULE MODgvec_Transform_SFL_Vars

