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
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_MHD3D_Vars
! MODULES
USE MOD_Globals,ONLY:wp
USE MOD_sgrid_vars, ONLY: t_sgrid
USE MOD_sbase_vars, ONLY: t_sbase
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE :: t_base              
  TYPE(t_sbase)  :: s  !! container for radial basis
!  TYPE(t_fbase)  :: f  !! container for angular basis
END TYPE t_base

TYPE :: t_sol_var
  REAL(wp) ,ALLOCATABLE :: X1(:)    !! X1 variable, size (base_f%mn_mode*base_s%nBase)
  REAL(wp) ,ALLOCATABLE :: X2(:)    !! X2 variable 
  REAL(wp) ,ALLOCATABLE :: LA(:)    !! lambda variable
END TYPE t_sol_var

!-----------------------------------------------------------------------------------------------------------------------------------
! SOLUTION VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(t_base)    :: X1base  
TYPE(t_base)    :: X2base  
TYPE(t_base)    :: LAbase  

TYPE(t_sgrid)   :: grid_s  !! only one grid up to now

TYPE(t_sol_var)  :: U(-1:1)         !! solutions at levels (k-1),(k),(k+1)
TYPE(t_sol_var)  :: dUdt            !! solution update

INTEGER          :: X1_BC(2)        !! BC axis (0) and edge (1)    
INTEGER          :: X2_BC(2)        !! BC axis (0) and edge (1)    
INTEGER          :: LA_BC(2)        !! BC axis (0) and edge (1)    
INTEGER          :: nDOF_X1         !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_X2         !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_LA         !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 

!===================================================================================================================================
END MODULE MOD_MHD3D_Vars

