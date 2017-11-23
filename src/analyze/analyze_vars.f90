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
!!# Module ** Analyze Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_Analyze_Vars
! MODULES
USE MOD_Globals,ONLY:wp
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER              :: which_visu             !! which data to visualize: 0: gvec, 1: vmec 
INTEGER              :: visu1D                 !! visualize 1D data (each mode). 0: off, 
                                               !! 1: R,Z,lambda  pointwise from VMEC input, (default) 
                                               !! 2: R,Z,lambda  interpolation 
                                               !! 3: R,Z  pointwise radial derivative from VMEC input, 
                                               !! 4: R,Z  radial derivative of interpolation 
                                               !! 12 : case 1 & 2, 13,23,123: with radial derivatives, combine 1,2,3,4 ascending
INTEGER              :: visu2D                 !! visualize 2D data ... 
INTEGER              :: visu3D                 !! visualize 3D data ... 
INTEGER              :: np_1d                  !! number of points for visualization in s
INTEGER              :: np_visu_bc(2)          !! number of points for visualization in theta,zeta
INTEGER              :: np_visu_planes(3)      !! number of points for visualization in s,theta,zeta
INTEGER              :: np_visu_3D(3)          !! number of points for visualization in s,theta,zeta

!===================================================================================================================================
END MODULE MOD_Analyze_Vars

