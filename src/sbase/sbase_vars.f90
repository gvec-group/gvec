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
!!# Module ** sBase Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_sBase_Vars
! MODULES
USE MOD_Globals    ,ONLY:wp
USE sll_m_bsplines ,ONLY: sll_c_bsplines
USE MOD_sGrid_Vars ,ONLY: t_sgrid
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 

TYPE, ABSTRACT :: c_sbase
  LOGICAL :: initialized
  CONTAINS
    PROCEDURE(i_sub_sbase_init    ),DEFERRED :: init
    PROCEDURE(i_sub_sbase_copy    ),DEFERRED :: copy
    PROCEDURE(i_sub_sbase_free    ),DEFERRED :: free
!    PROCEDURE(i_sub_sbase_initDOF ),DEFERRED :: initDOF
!    PROCEDURE(i_fun_sbase_eval    ),DEFERRED :: eval

END TYPE c_sbase

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sbase_init( self ,grid,deg,continuity,nGPloc)
    IMPORT wp, c_sbase,t_sgrid
    CLASS(c_sbase), INTENT(INOUT) :: self
    CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid
    INTEGER       , INTENT(IN   )        :: deg
    INTEGER       , INTENT(IN   )        :: continuity
    INTEGER       , INTENT(IN   )        :: nGPloc 
  END SUBROUTINE i_sub_sbase_init

  SUBROUTINE i_sub_sbase_copy( self, tocopy ) 
    IMPORT c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: tocopy
    CLASS(c_sbase), INTENT(INOUT) :: self
  END SUBROUTINE i_sub_sbase_copy

  SUBROUTINE i_sub_sbase_free( self ) 
    IMPORT c_sbase
    CLASS(c_sbase), INTENT(INOUT) :: self
  END SUBROUTINE i_sub_sbase_free

END INTERFACE
 


TYPE,EXTENDS(c_sbase) :: t_sBase
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  CLASS(t_sgrid),POINTER :: grid  => NULL()        !! pointer to grid 
  INTEGER              :: deg                      !! input parameter: degree of Spline/polynomial 
  INTEGER              :: continuity               !! input parameter: full spline (=deg-1) or discontinuous (=-1)
  INTEGER              :: nGPloc                   !! number of Gauss-points per element > deg
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: nGP                      !! total number of gausspoints = nGPloc*nElems
  INTEGER              :: nbase                    !! total number of degree of freedom / global basis functions
  REAL(wp),ALLOCATABLE :: xiGPloc(:)               !! element local gauss point positions, size(0:deg)
  REAL(wp),ALLOCATABLE :: s_IP(:)                  !! position of interpolation points for initialization, size(nBase) 
  REAL(wp),ALLOCATABLE :: wGPloc(:)                !! element local gauss weights, size(0:deg)
  REAL(wp),ALLOCATABLE :: wGP(:,:)                 !! global radial integration weight size(0:deg,1:nElems)
  INTEGER ,ALLOCATABLE :: base_offset(:)           !! offset of 0:deg element local basis functions to global index of
                                                   !! degree of freedom, allocated (1:nElems). iBase = offset(iElem)+j, j=0...deg
  REAL(wp),ALLOCATABLE :: base(:,:,:)              !! basis functions, (0:deg,0:nGPloc,1:nElems), 
  REAL(wp),ALLOCATABLE :: base_ds(:,:,:)           !! s derivative of basis functions, (0:deg,0:nGPloc,1:nElems)
  REAL(wp),ALLOCATABLE :: base_dsAxis(:,:)         !! first derivative d/ds of basis functions at the left  boundary (1:deg+1,0:deg)
  REAL(wp),ALLOCATABLE :: base_dsEdge(:,:)         !! first derivative d/ds of basis functions at the right boundary (1:deg+1,0:deg)
  REAL(wp),ALLOCATABLE :: A_Axis(:,:,:)            !! matrix to apply boundary conditions after interpolation (direct)
  REAL(wp),ALLOCATABLE :: R_Axis(:,:,:)            !! matrix to apply boundary conditions for RHS
  REAL(wp),ALLOCATABLE :: A_Edge(:,:,:)            !! matrix to apply boundary conditions after interpolation (direct)
  REAL(wp),ALLOCATABLE :: R_Edge(:,:,:)            !! matrix to apply boundary conditions for RHS
                                                   !! size(0:deg,0:deg,NBC_TYPES)
                                                   !! possible Boundary conditions
                                                   !! 1=BC_TYPE_OPEN      : boundary left open 
                                                   !! 2=BC_TYPE_NEUMANN   : 1st deriv fixed     
                                                   !! 3=BC_TYPE_DIRICHLET : sol fixed          (typically used at edge)
                                                   !! 4=BC_TYPE_SYMM      : all odd derivs = 0
                                                   !! 5=BC_TYPE_SYMMZERO  : all even derivs & sol = 0
                                                   !! 6=BC_TYPE_ANTISYMM  : all odd derivs & sol = 0

  CONTAINS
  PROCEDURE :: init          => sBase_init
  PROCEDURE :: copy          => sBase_copy
  PROCEDURE :: free          => sBase_free
!  PROCEDURE :: eval          => sBase_eval
!  PROCEDURE :: initDOF       => sBase_initDOF

END TYPE t_sBase

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type sbase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE sBase_init( self, grid,deg,continuity,nGPloc)
! MODULES
USE MOD_GLobals, ONLY: Unit_stdOut,fmt_sep,abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid       !! grid information
  INTEGER       , INTENT(IN   )        :: deg        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity !! continuity: 
                                                     !! 0: disc. polynomial
                                                     !! deg-1: spline with cont. deg-1
  INTEGER       , INTENT(IN   )        :: nGPloc 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i,iElem
  REAL(wp) :: scoord,xiIP(0:deg)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT sBase type:', &
       ' degree= ',deg, &
       ' gauss points per elem = ',nGPloc, &
       ' continuity= ',continuity, ' ...'
  IF(self%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL self%free()
  END IF
  
  self%grid  => grid 
  ASSOCIATE(nElems => self%grid%nElems)
  self%deg        = deg
  self%nGPloc    = nGPloc
  self%continuity = continuity
 
  self%nGP        = self%nGPloc*nElems
  IF(self%continuity.EQ.-1)THEN !discontinuous
    self%nBase  = (self%deg+1)*self%grid%nElems
  ELSEIF(self%continuity.EQ.deg-1)THEN !bspline with full continuity 
    self%nBase  = nElems + self%deg
  ELSE
   CALL abort(__STAMP__, &
          'other spline continuities not yet implemented') 
  END IF !continuity

  CALL sbase_alloc(self)

  !TODO CALL Gauss ==> xiGPloc,wGPloc 

  DO iElem=1,nElems
    self%wGP(:,iElem)=self%wGPloc(:)*grid%ds(iElem)
  END DO 

  IF(self%continuity.EQ.-1)THEN !discontinuous
    xiIP=0.0_wp
    !TODO compute xiIP, [-1,1]
    DO iElem=1,nElems
      self%s_IP(1+(deg+1)*(iElem-1):(deg+1)*iElem)=grid%sp(iElem-1)+0.5_wp*(xiIP+1.0_wp)*grid%ds(iElem)
    END DO !iElem 
   
    !TODO eval basis basis deriv,  
    DO iElem=1,nElems
      self%base_offset(iElem)=1+(self%deg+1)*(iElem-1)
!      self%base   (:,:,iElem)=bse(:,:)
!      self%base_ds(:,:,iElem)=bsederiv(:,:)*grid%ds(iElem)
    END DO !iElem 
    !TODO eval basis deriv at boundaries  
!    self%base_dsAxis(:,:)=bsederiv(:,:)*grid%ds(1)
!    self%base_dsEdge(:,:)=bsederiv(:,:)*grid%ds(nElems)
  ELSEIF(self%continuity.EQ.deg-1)THEN !bspline with full continuity 
    DO iElem=1,nElems
      !TODO compute offset self%base_offest(iElem)
      DO i=0,self%nGPloc
        scoord=grid%sp(iElem-1)+self%xiGPloc(i)*grid%ds(iElem)
        !TODO eval basis basis deriv,  
      END DO !i=0,nGPloc
    END DO !iElem 
    !TODO eval basis deriv at boundaries  
  END IF !continuity
  END ASSOCIATE !nElems 
  self%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE sBase_init


!===================================================================================================================================
!> allocate all variables in  sbase
!!
!===================================================================================================================================
SUBROUTINE sBase_alloc( self)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  ASSOCIATE(nElems=>self%grid%nElems, nGPloc=>self%nGPloc, deg=>self%deg)
  ALLOCATE(self%xiGPloc(0:deg))
  ALLOCATE(self%wGPloc( 0:deg))
  ALLOCATE(self%wGP(    0:deg,1:nElems))
  ALLOCATE(self%A_Axis(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%A_Edge(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%R_Axis(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%R_Edge(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%base_offset(           1:nElems))
  ALLOCATE(self%base(   0:deg,0:nGPloc,1:nElems))
  ALLOCATE(self%base_ds(0:deg,0:nGPloc,1:nElems))
  ALLOCATE(self%base_dsAxis(1:deg+1,0:deg))
  ALLOCATE(self%base_dsEdge(1:deg+1,0:deg))
  ALLOCATE(self%s_IP(self%nBase))
  self%xiGPloc     =0.0_wp
  self%wGPloc      =0.0_wp            
  self%wGP         =0.0_wp            
  self%A_Axis      =0.0_wp         
  self%A_Edge      =0.0_wp         
  self%R_Axis      =0.0_wp         
  self%R_Edge      =0.0_wp         
  DO i=0,deg
    self%A_Axis(i,i,:)=1.0_wp
    self%A_Edge(i,i,:)=1.0_wp
    self%R_Axis(i,i,:)=1.0_wp
    self%R_Edge(i,i,:)=1.0_wp
  END DO
  self%base_offset =-1
  self%base        =0.0_wp           
  self%base_ds     =0.0_wp        
  self%base_dsAxis =0.0_wp    
  self%base_dsEdge =0.0_wp    
  self%s_IP=0.0_wp           
  END ASSOCIATE !nElems, nGploc, deg
END SUBROUTINE sbase_alloc


!===================================================================================================================================
!> copy the type sbase
!!
!===================================================================================================================================
SUBROUTINE sBase_copy( self , tocopy)
! MODULES
USE MOD_GLobals, ONLY: Unit_stdOut,abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_sBase), INTENT(IN   ) :: tocopy
  CLASS(t_sBase), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPEIS(t_sbase)
  IF(.NOT.tocopy%initialized) THEN
    CALL abort(__STAMP__, &
        "sBase_copy: not initialized sBase from which to copy!")
  END IF
  IF(self%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL self%free()
  END IF

  CALL self%init(tocopy%grid,tocopy%deg,tocopy%continuity,tocopy%nGPloc)

  END SELECT !TYPE
END SUBROUTINE sbase_copy


!===================================================================================================================================
!> finalize the type sBase
!!
!===================================================================================================================================
SUBROUTINE sBase_free( self )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  self%deg        =-1 
  self%continuity =-99
  self%nGPloc     =-1
  self%nGP        =-1        
  self%nbase      =-1
  NULLIFY(self%grid)
  SDEALLOCATE(self%xiGPloc)
  SDEALLOCATE(self%wGPloc)
  SDEALLOCATE(self%wGP)
  SDEALLOCATE(self%s_IP)
  SDEALLOCATE(self%base_offset)
  SDEALLOCATE(self%base)   
  SDEALLOCATE(self%base_ds)
  SDEALLOCATE(self%base_dsAxis) 
  SDEALLOCATE(self%base_dsEdge) 
  SDEALLOCATE(self%A_Axis) 
  SDEALLOCATE(self%R_Axis) 
  SDEALLOCATE(self%A_Edge) 
  SDEALLOCATE(self%R_Edge) 
  
  self%initialized=.FALSE.

END SUBROUTINE sBase_free

END MODULE MOD_sBase_Vars

