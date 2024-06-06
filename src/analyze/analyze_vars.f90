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
MODULE MODgvec_Analyze_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER              :: iAnalyze               !! counter for calls to analyze (used for output file names)
INTEGER              :: visu1D                 !! visualize 1D data (each mode). 0: off, 
                                               !! 1: R,Z,lambda  pointwise from VMEC input, (default) 
                                               !! 2: R,Z,lambda  interpolation 
                                               !! 3: R,Z  pointwise radial derivative from VMEC input, 
                                               !! 4: R,Z  radial derivative of interpolation 
                                               !! 12 : case 1 & 2, 13,23,123: with radial derivatives, combine 1,2,3,4 ascending
INTEGER              :: visu2D                 !! visualize 2D data ... 
INTEGER              :: visu3D                 !! visualize 3D data ... 
INTEGER              :: SFLout            !! input parameter: convert final state to straight-field line coordinates. =0: off, =1: PEST, =2: Boozer
INTEGER              :: SFLout_mn_max(2)   !! maximum mode number in theta and zeta. Defaults to 4*mn_max of X1_base, if set to (-1,-1)
INTEGER              :: SFLout_nrp,SFLout_mn_pts(2 )  !! number of points for SFLOut file in theta,zeta
REAL(wp),ALLOCATABLE :: SFLout_radialpos(:)     !! radial positions for output
INTEGER              :: outfileType=0          !! =1: default, vtk paraview file. =2: structured netcdf array., 12: both
INTEGER              :: np_1d                  !! number of points for visualization in s
INTEGER              :: np_visu_bc(2)          !! number of points for visualization in theta,zeta
INTEGER              :: np_visu_planes(3)      !! number of points for visualization in s,theta,zeta
INTEGER              :: np_visu_3D(3)          !! number of points for visualization in s,theta,zeta
REAL(wp)             :: visu_BC_minmax(2:3,0:1)    !! minimum and maximum in s,theta,zeta [0,1]
REAL(wp)             :: visu_planes_minmax(1:3,0:1)!! minimum and maximum in s,theta,zeta [0,1]
REAL(wp)             :: visu_3D_minmax(1:3,0:1)    !! minimum and maximum in s,theta,zeta [0,1]
LOGICAL              :: SFL_theta                  !! =T: visualize a mesh with a sstraight-field line (PEST) theta angle

!===================================================================================================================================
END MODULE MODgvec_Analyze_Vars

