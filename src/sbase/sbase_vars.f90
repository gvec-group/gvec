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
!!# Module ** Sbase Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_Sbase_Vars
! MODULES
USE MOD_Globals    ,ONLY:wp
USE sll_m_bsplines ,ONLY: sll_c_bsplines
USE sgrid          ,ONLY: c_grid
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 

TYPE, ABSTRACT :: c_base
  LOGICAL :: initialized
  CONTAINS
    PROCEDURE(i_sub_init    ),DEFERRED :: init
    PROCEDURE(i_sub_free    ),DEFERRED :: free
    PROCEDURE(i_sub_initDOF ),DEFERRED :: initDOF
    PROCEDURE(i_fun_eval    ),DEFERRED :: eval

END TYPE c_base

ABSTRACT INTERFACE
  SUBROUTINE init( self , grid,deg,nGPelem,BC_axis,BC_edge)
    IMPORT wp,c_grid
    CLASS(c_base), INTENT(INOUT) :: self
    CLASS(c_grid), INTENT(IN   ) :: grid
    INTEGER      , INTENT(IN   ) :: deg
    INTEGER      , INTENT(IN   ) :: nGPelem 
    INTEGER      , INTENT(IN   ) :: BC_axis(BC_MODES)
    INTEGER      , INTENT(IN   ) :: BC_edge(BC_MODES)
  END SUBROUTINE init

END INTERFACE
 


TYPE,EXTENDS(c_base) :: t_sbase
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: nGP_elem                 !! number of Gauss-points per element > deg
  INTEGER              :: BC_AXIS(3)               !! Boundary condition type at axis BC for BC_MZERO,BC_MODD,BC_MEVEN
  INTEGER              :: BC_EDGE(3)               !! Boundary condition type at edge BC for BC_MZERO,BC_MODD,BC_MEVEN
                                                   !! possible values 
                                                   !! BC_OPEN      : boundary left open 
                                                   !! BC_NEUMANN   : 1st deriv fixed     
                                                   !! BC_DIRICHLET : sol fixed          (typically used at edge)
                                                   !! BC_SYMM      : all odd derivs = 0
                                                   !! BC_SYMMZERO  : all even derivs & sol = 0
                                                   !! BC_ANTISYMM  : all odd derivs & sol = 0
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: nElems                   !! total number of radial elements, taken from input grid 
  INTEGER              :: nGP                      !! total number of gausspoints
  INTEGER              :: nbase                    !! total number of degree of freedom / global basis functions
  INTEGER              :: deg                      !! input parameter: degree of Spline/polynomial 
  INTEGER ,ALLOCATABLE :: base_offset(:)           !! offset of 0:deg element local basis functions to global index of
                                                   !! degree of freedom, allocated (1:nElems). iBase = offset(iElem)+j, j=0...deg
  REAL(wp),ALLOCATABLE :: base(:,:,:)              !! basis functions, (0:deg,0:nGP,1:nElems), 
  REAL(wp),ALLOCATABLE :: base_ds(:,:,:)           !! s derivative of basis functions, (0:deg,0:nGP,1:nElems)
  REAL(wp),ALLOCATABLE :: base_dsL(:,:)            !! first derivative d/ds of basis functions at the left  boundary (1:deg+1,0:deg) 
  REAL(wp),ALLOCATABLE :: base_dsR(:,:)            !! first derivative d/ds of basis functions at the right boundary (1:deg+1,0:deg) 
  REAL(wp),ALLOCATABLE :: A_Axis(:,:,:)            !! matrix to apply boundary conditions after interpolation (direct)
  REAL(wp),ALLOCATABLE :: R_Axis(:,:,:)            !! matrix to apply boundary conditions for RHS
  REAL(wp),ALLOCATABLE :: A_Edge(:,:,:)            !! matrix to apply boundary conditions after interpolation (direct)
  REAL(wp),ALLOCATABLE :: R_Edge(:,:,:)            !! matrix to apply boundary conditions for RHS
                                                   !! size(0:deg,0:deg,3), for BC_MZERO,BC_MODD,BC_MEVEN
  CONTAINS

  PROCEDURE :: init          => init_sbase 
  PROCEDURE :: initDOF       => initDOF_sbase 
  PROCEDURE :: free          => free_sbase
  PROCEDURE :: eval          => eval_sbase 

END TYPE t_sbase
!LOGICAL              :: usethis                  !! =(whichEquilibrium>0)

!===================================================================================================================================
END MODULE MOD_Sbase_Vars

