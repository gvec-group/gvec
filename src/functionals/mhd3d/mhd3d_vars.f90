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
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_MHD3D_Vars
! MODULES
USE MOD_Globals,ONLY: wp,Unit_stdOut,abort
USE MOD_sgrid,  ONLY: t_sgrid
USE MOD_base,   ONLY: t_base
USE MOD_Sol_Var,ONLY: t_sol_var
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------------------------------
! SOLUTION VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CLASS(t_base),ALLOCATABLE :: X1base            !! container for base of variable X1
CLASS(t_base),ALLOCATABLE :: X2base            !! container for base of variable X2
CLASS(t_base),ALLOCATABLE :: LAbase            !! container for base of variable lambda 
                             
TYPE(t_sgrid)    :: sgrid             !! only one grid up to now

TYPE(t_sol_var),ALLOCATABLE  :: U(:)  !! solutions at levels (k-1),(k),(k+1)
TYPE(t_sol_var)              :: dUdt  !! solution update

INTEGER          :: X1X2_BC(2)        !! BC axis (0) and edge (1)   for variables X1 and X2
INTEGER          :: LA_BC(2)          !! BC axis (0) and edge (1)   for variable lambda
INTEGER          :: nDOF_X1           !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_X2           !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_LA           !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 

!===================================================================================================================================

END MODULE MOD_MHD3D_Vars

