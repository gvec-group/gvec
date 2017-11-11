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
USE MOD_Globals                  ,ONLY: wp
USE sll_m_bsplines               ,ONLY: sll_c_bsplines
USE sll_m_spline_1d              ,ONLY: sll_t_spline_1d
USE sll_m_spline_interpolator_1d ,ONLY: sll_t_spline_interpolator_1d
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
  SUBROUTINE i_sub_sbase_init( self ,grid_in,deg_in,continuity_in,degGP_in)
    IMPORT wp, c_sbase,t_sgrid
    CLASS(c_sbase), INTENT(INOUT) :: self
    CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in
    INTEGER       , INTENT(IN   )        :: deg_in
    INTEGER       , INTENT(IN   )        :: continuity_in
    INTEGER       , INTENT(IN   )        :: degGP_in
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
  INTEGER              :: degGP                    !! number of Gauss-points (degGP+1) per element >= deg
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: nGP                      !! total number of gausspoints = degGP*nElems
  INTEGER              :: nbase                    !! total number of degree of freedom / global basis functions
  REAL(wp),ALLOCATABLE :: xiGPloc(:)               !! element local gauss point positions for interval [0,1], size(0:deg)
  REAL(wp),ALLOCATABLE :: wGPloc(:)                !! element local gauss weights for interval [0,1], size(0:deg)
  REAL(wp),ALLOCATABLE :: wGP(:,:)                 !! global radial integration weight size(0:deg,1:nElems)
  REAL(wp),ALLOCATABLE :: s_GP(:,:)                !! global position of gauss points  in s [0,1] , size(nGPloc,nElems) 
  REAL(wp),ALLOCATABLE :: s_IP(:)                  !! position of interpolation points for initialization, size(nBase) 
  INTEGER ,ALLOCATABLE :: base_offset(:)           !! offset of 0:deg element local basis functions to global index of
                                                   !! degree of freedom, allocated (1:nElems). iBase = offset(iElem)+j, j=0...deg
  REAL(wp),ALLOCATABLE :: baseGP(:,:,:)            !! basis functions, (0:deg,0:degGP,1:nElems), 
  REAL(wp),ALLOCATABLE :: base_dsGP(:,:,:)         !! s derivative of basis functions, (0:deg,0:degGP,1:nElems)
  REAL(wp),ALLOCATABLE :: base_dsAxis(:,:)         !! all derivatives 1..deg of all basis functions at axis (1:deg+1,0:deg)
  REAL(wp),ALLOCATABLE :: base_dsEdge(:,:)         !! all derivatives 1..deg of all basis functions at edge (1:deg+1,0:deg)
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
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl        !! contains bspline functions
  TYPE(sll_t_spline_1d)             :: spline      !! contains 1d spline functions
  TYPE(sll_t_spline_interpolator_1d):: Interpol    !! spline interpolator
  

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
SUBROUTINE sBase_init( self, grid_in,deg_in,continuity_in,degGP_in)
! MODULES
USE MOD_GLobals, ONLY: PI,Unit_stdOut,fmt_sep,abort
USE sll_m_bsplines,ONLY: sll_s_bsplines_new
USE sll_m_spline_1d                      ,ONLY: sll_t_spline_1d
USE sll_m_spline_interpolator_1d         ,ONLY: sll_t_spline_interpolator_1d
USE sll_m_boundary_condition_descriptors ,ONLY: sll_p_greville
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in       !! grid information
  INTEGER       , INTENT(IN   )        :: deg_in        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity_in !! continuity: 
                                                        !! 0: disc. polynomial
                                                        !! deg-1: spline with cont. deg-1
  INTEGER       , INTENT(IN   )        :: degGP_in 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), INTENT(INOUT) :: self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i,iGP,iElem,imin,jmin
  REAL(wp) :: xiIP(0:deg_in)
  REAL(wp) :: locbasis(0:deg_in,0:deg_in)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT sBase type:', &
       ' degree= ',deg_in, &
       ' gauss points per elem = ',degGP_in, &
       ' continuity= ',continuity_in, ' ...'
  IF(self%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL self%free()
  END IF
  
  self%grid       => grid_in
  self%deg        =  deg_in
  self%degGP      =  degGP_in
  self%continuity =  continuity_in

  ASSOCIATE(&
              nElems      =>self%grid%nElems   &
            , sp          =>self%grid%sp       &
            , ds          =>self%grid%ds       &
            , deg         =>self%deg           &
            , continuity  =>self%continuity    &
            , degGP       =>self%degGP         &
            , nGP         =>self%nGP           &
            , nBase       =>self%nBase         &
            , xiGPloc     =>self%xiGPloc       &
            , wGPloc      =>self%wGPloc        &
            , wGP         =>self%wGP           &
            , s_GP        =>self%s_GP          &
            , s_IP        =>self%s_IP          &
            , base_offset =>self%base_offset   &
            , baseGP      =>self%baseGP        &
            , base_dsGP   =>self%base_dsGP     &
            , base_dsAxis =>self%base_dsAxis   &
            , base_dsEdge =>self%base_dsEdge   &
            , bspl        =>self%bspl          &
            , spline      =>self%spline        &
            , interpol    =>self%interpol      &
            )
 
  nGP  = (degGP+1)*nElems
  IF(continuity.EQ.-1)THEN !discontinuous
    nBase  = (deg+1)*nElems
  ELSEIF(continuity.EQ.deg-1)THEN !bspline with full continuity 
    nBase  = nElems + deg
  ELSE
   CALL abort(__STAMP__, &
          'other spline continuities not yet implemented') 
  END IF !continuity

  CALL sbase_alloc(self)

             
  CALL LegendreGaussNodesAndWeights(deg,xiGPloc,wGPloc)
  ![-1,1]->[0,1]
  xiGPloc=0.5_wp*(xiGPloc+1.0_wp)
  wGPloc =0.5*wGPloc

  DO iElem=1,nElems
    wGP(:,iElem)=wGPloc(:)*ds(iElem)
  END DO 
 DO iElem=1,nElems
   s_GP(:,iElem)=sp(iElem-1)+xiGPloc(:)*ds(iElem)
 END DO !iElem 

  IF(continuity.EQ.-1)THEN !discontinuous
   
    !TODO eval basis basis deriv,  
    DO iElem=1,nElems
      base_offset(iElem)=1+(deg+1)*(iElem-1)
!      baseGP   (:,:,iElem)=bse(:,:)
!      base_dsGP(:,:,iElem)=bsederiv(:,:)*ds(iElem)
    END DO !iElem 
    !TODO eval basis deriv at boundaries  
!    base_dsAxis(:,:)=bsederiv(:,:)*ds(1)
!    base_dsEdge(:,:)=bsederiv(:,:)*ds(nElems)

    !interpolation:
    !  points are repeated at element interfaces (discontinuous)
    !  use chebychev-lobatto points for interpolation (closed form!), interval [0,1] 
    DO i=0,deg
      xiIP(i)=0.5_wp*(1.0_wp-COS(REAL(i,wp)/REAL(deg,wp)*PI))
    END DO
    ALLOCATE(self%s_IP(self%nBase)) !for spl, its allocated elsewhere...
    DO iElem=1,nElems
      s_IP(1+(deg+1)*(iElem-1):(deg+1)*iElem)=sp(iElem-1)+xiIP*ds(iElem)
    END DO !iElem 
  ELSEIF(continuity .EQ. deg-1)THEN !bspline with full continuity 
    CALL sll_s_bsplines_new(bspl ,degree=deg,periodic=.FALSE.,xmin=0.0_wp,xmax=1.0_wp,ncells=nElems,breaks=sp(:))
    !basis evaluation
    IF(bspl%nBasis.NE.nBase) STOP 'problem with bspl basis'
    DO iElem=1,nElems
      CALL bspl % eval_basis(s_GP(0,iElem),baseGP(0,0:deg,iElem),imin)
      DO iGP=1,degGP
         CALL bspl % eval_basis(s_GP(iGP,iElem),baseGP(iGP,0:deg,iElem),jmin)
         IF(jmin.NE.imin) STOP'problem, GP are not in one element!'
      END DO !iGP=0,degGP
      CALL bspl % eval_deriv(s_GP(0,iElem),base_dsGP(0,0:deg,iElem),imin)
      DO iGP=1,degGP
         CALL bspl % eval_deriv(s_GP(iGP,iElem),base_dsGP(iGP,0:deg,iElem),jmin)
         IF(jmin.NE.imin) STOP'problem, GP are not in one element!'
      END DO !iGP=0,degGP
      base_offset(iElem)=imin
    END DO !iElem=1,nElems
    !eval all basis derivatives at boundaries  
    CALL bspl % eval_basis_and_n_derivs(sp(     0),deg,locBasis,imin)
    IF(imin.NE.1) STOP'problem eval_deriv left'
    base_dsAxis=TRANSPOSE(locbasis(:,:))

    CALL bspl % eval_basis_and_n_derivs(sp(nElems),deg,locbasis,imin)
    IF(imin.NE.nBase-deg) STOP'problem eval_deriv right'
    base_dsEdge=TRANSPOSE(locbasis(:,:))

    !interpolation
    CALL Interpol%init (bspl,sll_p_greville,sll_p_greville) 
    CALL Interpol%get_interp_points ( self%s_IP ) 
    CALL spline%init( bspl ) !needed for interpolation

  END IF !continuity

  !TODO STRONG BOUNDARY CONDITIONS: A and R matrices
  
  END ASSOCIATE !self


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
  ASSOCIATE(nElems=>self%grid%nElems, degGP=>self%degGP, deg=>self%deg)
  ALLOCATE(self%xiGPloc(0:deg))
  ALLOCATE(self%wGPloc( 0:deg))
  ALLOCATE(self%wGP(    0:deg,1:nElems))
  ALLOCATE(self%base_offset(           1:nElems))
  ALLOCATE(self%baseGP(   0:deg,0:degGP,1:nElems))
  ALLOCATE(self%base_dsGP(0:deg,0:degGP,1:nElems))
  ALLOCATE(self%base_dsAxis(1:deg+1,0:deg))
  ALLOCATE(self%base_dsEdge(1:deg+1,0:deg))
  ALLOCATE(self%A_Axis(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%A_Edge(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%R_Axis(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(self%R_Edge(0:deg,0:deg,1:NBC_TYPES))
  self%xiGPloc     =0.0_wp
  self%wGPloc      =0.0_wp            
  self%wGP         =0.0_wp            
  self%s_GP        =0.0_wp            
  self%base_offset =-1
  self%baseGP      =0.0_wp           
  self%base_dsGP   =0.0_wp        
  self%base_dsAxis =0.0_wp    
  self%base_dsEdge =0.0_wp    
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
  END ASSOCIATE !nElems, degGP, deg
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

  CALL self%init(tocopy%grid,tocopy%deg,tocopy%continuity,tocopy%degGP)

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

  !pointers, classes
  NULLIFY(self%grid)
  IF(self%continuity.EQ.self%deg-1)THEN
    CALL self%interpol%free() 
    CALL self%spline%free() 
    CALL self%bspl%free() 
  END IF
  !allocatables  
  SDEALLOCATE(self%xiGPloc)
  SDEALLOCATE(self%wGPloc)
  SDEALLOCATE(self%wGP)
  SDEALLOCATE(self%s_GP)
  SDEALLOCATE(self%s_IP)
  SDEALLOCATE(self%base_offset)
  SDEALLOCATE(self%baseGP)   
  SDEALLOCATE(self%base_dsGP)
  SDEALLOCATE(self%base_dsAxis) 
  SDEALLOCATE(self%base_dsEdge) 
  SDEALLOCATE(self%A_Axis) 
  SDEALLOCATE(self%R_Axis) 
  SDEALLOCATE(self%A_Edge) 
  SDEALLOCATE(self%R_Edge) 
  
  self%continuity =-99
  self%deg        =-1 
  self%degGP      =-1
  self%nGP        =-1        
  self%nbase      =-1
  self%initialized=.FALSE.

END SUBROUTINE sBase_free

END MODULE MOD_sBase_Vars

