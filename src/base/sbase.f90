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
!!# Module ** sBase **
!!
!! 1D basis in radial coordinate "s". Contains sbase type definition and associated routines
!!
!===================================================================================================================================
MODULE MOD_sBase
! MODULES
USE MOD_Globals                  ,ONLY: wp
USE sll_m_bsplines               ,ONLY: sll_c_bsplines
USE sll_m_spline_1d              ,ONLY: sll_t_spline_1d
USE sll_m_spline_interpolator_1d ,ONLY: sll_t_spline_interpolator_1d
USE MOD_sGrid ,ONLY: c_sgrid,t_sgrid
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE, ABSTRACT :: c_sbase
  LOGICAL :: initialized
  CONTAINS
    PROCEDURE(i_sub_sbase_init    ),DEFERRED :: init
    PROCEDURE(i_sub_sbase_free    ),DEFERRED :: free
    PROCEDURE(i_sub_sbase_copy    ),DEFERRED :: copy

END TYPE c_sbase

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sbase_init( sf ,grid_in,deg_in,continuity_in,degGP_in)
    IMPORT wp, c_sbase,t_sgrid
    CLASS(c_sbase), INTENT(INOUT)        :: sf
    CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in
    INTEGER       , INTENT(IN   )        :: deg_in
    INTEGER       , INTENT(IN   )        :: continuity_in
    INTEGER       , INTENT(IN   )        :: degGP_in
  END SUBROUTINE i_sub_sbase_init

  SUBROUTINE i_sub_sbase_copy( sf, tocopy ) 
    IMPORT c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: tocopy
    CLASS(c_sbase), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sbase_copy

  SUBROUTINE i_sub_sbase_free( sf ) 
    IMPORT c_sbase
    CLASS(c_sbase), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sbase_free

END INTERFACE
 


!TYPE,EXTENDS(c_sbase) :: t_sBase
TYPE :: t_sBase
  LOGICAL :: initialized
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  CLASS(t_sgrid),POINTER :: grid  => NULL()        !! pointer to grid 
  INTEGER              :: deg                      !! input parameter: degree of Spline/polynomial 
  INTEGER              :: degGP                    !! number of Gauss-points (degGP+1) per element >= deg
  INTEGER              :: continuity               !! input parameter: full spline (=deg-1) or discontinuous (=-1)
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: nGP                      !! total number of gausspoints = degGP*nElems
  INTEGER              :: nbase                    !! total number of degree of freedom / global basis functions
  REAL(wp),ALLOCATABLE :: xiGP(:)                  !! element local gauss point positions for interval [0,1], size(0:degGP)
  REAL(wp),ALLOCATABLE :: wGPloc(:)                !! element local gauss weights for interval [0,1], size(0:degGP)
  REAL(wp),ALLOCATABLE :: wGP(:,:)                 !! global radial integration weight size(0:degGP,1:nElems)
  REAL(wp),ALLOCATABLE :: s_GP(:,:)                !! global position of gauss points  in s [0,1] , size(0:degGP,nElems) 
  REAL(wp),ALLOCATABLE :: s_IP(:)                  !! position of interpolation points for initialization, size(nBase) 
  INTEGER ,ALLOCATABLE :: base_offset(:)           !! offset of 0:deg element local basis functions to global index of
                                                   !! degree of freedom, allocated (1:nElems). iBase = offset(iElem)+j, j=0...deg
  REAL(wp),ALLOCATABLE :: baseGP(:,:,:)            !! basis functions, (0:degGP,0:deg,1:nElems), 
  REAL(wp),ALLOCATABLE :: base_dsGP(:,:,:)         !! s derivative of basis functions, (0:degGP,0:deg,1:nElems)
  REAL(wp),ALLOCATABLE :: base_dsAxis(:,:)         !! all derivatives 1..deg of all basis functions at axis (1:deg,0:deg)
  REAL(wp),ALLOCATABLE :: base_dsEdge(:,:)         !! all derivatives 1..deg of all basis functions at edge (1:deg,0:deg)
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
  PROCEDURE :: free          => sBase_free
  PROCEDURE :: copy          => sBase_copy
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
SUBROUTINE sBase_init( sf, grid_in,deg_in,continuity_in,degGP_in)
! MODULES
USE MOD_GLobals, ONLY: PI,Unit_stdOut,abort
USE MOD_Basis1D, ONLY:  LegendreGaussNodesAndWeights
USE MOD_Basis1D, ONLY:  BarycentricWeights,InitializeVandermonde,PolynomialDerivativeMatrix
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
  INTEGER       , INTENT(IN   )        :: degGP_in      !! gauss quadrature points: nGP=degGP+1 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i,iGP,iElem,imin,jmin
  REAL(wp) :: locbasis(0:deg_in,0:deg_in)
  REAL(wp) :: xiIP(0:deg_in)
  REAL(wp) :: wBaryIP(0:deg_in)
  REAL(wp) :: Vdm(   0:degGP_in,0:deg_in)
  REAL(wp) :: DmatGP(0:degGP_in,0:deg_in)
  REAL(wp) :: Dmat(  0:deg_in  ,0:deg_in)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A,3(A,I3),A)')'INIT sBase type:', &
       ' degree= ',deg_in, &
       ' gauss points per elem = ',degGP_in, &
       ' continuity= ',continuity_in, ' ...'
  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL sf%free()
  END IF
  IF(degGP_in.LT.deg_in) &
    CALL abort(__STAMP__, &
        "error in sbase: degGP must be > deg!") 
  sf%grid       => grid_in
  sf%deg        =  deg_in
  sf%degGP      =  degGP_in
  sf%continuity =  continuity_in

  ASSOCIATE(&
              nElems      => sf%grid%nElems         &
            , grid        => sf%grid                &
            , deg         => sf%deg                 &
            , degGP       => sf%degGP               &
            , continuity  => sf%continuity          &
            , nGP         => sf%nGP                 &
            , nBase       => sf%nBase               &
            )
 
  nGP  = (degGP+1)*nElems
  IF(continuity.EQ.-1)THEN !discontinuous
    nBase  = (deg+1)*nElems
  ELSEIF(continuity.EQ.deg-1)THEN !bspline with full continuity and interpolation base at boundaries 
    nBase  = nElems + deg
  ELSE
   CALL abort(__STAMP__, &
          'other spline continuities not yet implemented') 
  END IF !continuity

  CALL sbase_alloc(sf)

             
  CALL LegendreGaussNodesAndWeights(degGP,sf%xiGP,sf%wGPloc)
  ![-1,1]->[0,1]
  sf%xiGP=0.5_wp*(sf%xiGP+1.0_wp)
  sf%wGPloc =0.5*sf%wGPloc

  DO iElem=1,nElems
    sf%wGP(:,iElem)=sf%wGPloc(:)*grid%ds(iElem)
  END DO 
  DO iElem=1,nElems
    sf%s_GP(:,iElem)=grid%sp(iElem-1)+sf%xiGP(:)*grid%ds(iElem)
  END DO !iElem 

  IF(continuity.EQ.-1)THEN !discontinuous
   
    !  use chebychev-lobatto points for interpolation (closed form!), interval [-1,1] 
    DO i=0,deg
      xiIP(i)=-COS(REAL(i,wp)/REAL(deg,wp)*PI)
    END DO
    CALL BarycentricWeights(deg,xiIP,wBaryIP)
    CALL InitializeVandermonde(deg,degGP,wBaryIP,xiIP,sf%xiGP,Vdm)
    CALL PolynomialDerivativeMatrix(deg,xiIP,Dmat)
    DmatGP=MATMUL(Vdm,Dmat)
    ! eval basis and  basis derivative  
    DO iElem=1,nElems
      sf%base_offset(iElem)=1+(deg+1)*(iElem-1)
      sf%baseGP   (:,:,iElem)=Vdm(:,:)
      sf%base_dsGP(:,:,iElem)=DmatGP*(2.0_wp/grid%ds(iElem))
    END DO !iElem 
    !zero deriv: evaluation of basis functions (lagrange property!)
    sf%base_dsAxis(1      ,0)=1.0_wp
    sf%base_dsAxis(2:deg+1,0)=0.0_wp
    sf%base_dsEdge(nBase-deg:nBase-1 ,0)=0.0_wp
    sf%base_dsEdge(          nBase   ,0)=1.0_wp
    ! eval basis deriv at boundaries d/ds = d/dxi dxi/ds = 1/(0.5ds) d/dxi  
    sf%base_dsAxis(1:deg+1        ,1)=Dmat(  0,:)*(2.0_wp/grid%ds(1))
    sf%base_dsEdge(nBase-deg:nBase,1)=Dmat(deg,:)*(2.0_wp/grid%ds(nElems))
    !  higher derivatives 
    DO i=2,deg
      sf%base_dsAxis(1:deg+1        ,i)=MATMUL(TRANSPOSE(Dmat),sf%base_dsAxis(:,i-1))*(2.0_wp/grid%ds(1))
      sf%base_dsEdge(nBase-deg:nBase,i)=MATMUL(TRANSPOSE(Dmat),sf%base_dsEdge(:,i-1))*(2.0_wp/grid%ds(nElems))
    END DO
    !interpolation:
    !  points are repeated at element interfaces (discontinuous)
    ALLOCATE(sf%s_IP(sf%nBase)) !for spl, its allocated elsewhere...
    DO iElem=1,nElems
      sf%s_IP(1+(deg+1)*(iElem-1):(deg+1)*iElem)=grid%sp(iElem-1)+0.5_wp*(xiIP+1.0_wp)*grid%ds(iElem)
    END DO !iElem 
  ELSEIF(continuity .EQ. deg-1)THEN !bspline with full continuity 
    CALL sll_s_bsplines_new(sf%bspl ,degree=deg,periodic=.FALSE.,xmin=0.0_wp,xmax=1.0_wp,ncells=nElems,breaks=grid%sp(:))
    !basis evaluation
    IF(sf%bspl%nBasis.NE.nBase) STOP 'problem with bspl basis'
    DO iElem=1,nElems
      CALL sf%bspl % eval_basis(sf%s_GP(0,iElem),sf%baseGP(0,0:deg,iElem),imin)
      DO iGP=1,degGP
         CALL sf%bspl % eval_basis(sf%s_GP(iGP,iElem),sf%baseGP(iGP,0:deg,iElem),jmin)
         IF(jmin.NE.imin) STOP'problem, GP are not in one element!'
      END DO !iGP=0,degGP
      CALL sf%bspl % eval_deriv(sf%s_GP(0,iElem),sf%base_dsGP(0,0:deg,iElem),imin)
      DO iGP=1,degGP
         CALL sf%bspl % eval_deriv(sf%s_GP(iGP,iElem),sf%base_dsGP(iGP,0:deg,iElem),jmin)
         IF(jmin.NE.imin) STOP'problem, GP are not in one element!'
      END DO !iGP=0,degGP
      sf%base_offset(iElem)=imin
    END DO !iElem=1,nElems
    !eval all basis derivatives at boundaries  
    CALL sf%bspl % eval_basis_and_n_derivs(grid%sp(     0),deg,locBasis,imin) !locBasis(0:nderiv,0:deg Base)
    IF(imin.NE.1) STOP'problem eval_deriv left'
    sf%base_dsAxis(1:deg+1,0:deg) =TRANSPOSE(locbasis(:,:)) ! basis functions 1 ...deg+1

    CALL sf%bspl % eval_basis_and_n_derivs(grid%sp(nElems),deg,locbasis,imin)
    IF(imin.NE.nBase-deg) STOP'problem eval_deriv right'
    sf%base_dsEdge(nBase-deg:nBase,0:deg)=TRANSPOSE(locbasis(:,:) ) ! basis functions nBase-deg ... nbase

    !interpolation
    CALL sf%Interpol%init (sf%bspl,sll_p_greville,sll_p_greville) 
    CALL sf%Interpol%get_interp_points ( sf%s_IP ) 
    CALL sf%spline%init( sf%bspl ) !needed for interpolation

  END IF !continuity

  !TODO STRONG BOUNDARY CONDITIONS: A and R matrices
  
  END ASSOCIATE !sf


  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
  CALL sBase_test(sf)

END SUBROUTINE sBase_init


!===================================================================================================================================
!> test sbase variable
!!
!===================================================================================================================================
SUBROUTINE sBase_test( sf)
! MODULES
USE MOD_GLobals, ONLY: UNIT_stdOut,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iTest
!===================================================================================================================================
  IF(testlevel.LE.0) RETURN
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN TEST No.',nTestCalled,' SBASE    >>>>>>>>>'
  IF(testlevel.GT.0)THEN
    iTest=1
    IF( ABS(sf%s_IP(1)).GT. 1.0e-12_wp) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '!! TEST No.',nTestCalled ,': TEST ',iTest,' IN SGRID FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(2(A,I4),(A,E11.3))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity, &
        ', sIP(1)= ',sf%s_IP(1)
    END IF
    iTest=2
    IF( ABS(sf%s_IP(sf%nBase)-1.0_wp).GT. 1.0e-12_wp) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '!! TEST No.',nTestCalled ,': TEST ',iTest,' IN SGRID FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),(A,E11.3))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
        ', sIP(nBase)= ',sf%s_IP(sf%nBase)
    END IF
  ELSEIF(testlevel.GT.1)THEN
    
  END IF !testlevel>0
  
END SUBROUTINE sbase_test


!===================================================================================================================================
!> allocate all variables in  sbase
!!
!===================================================================================================================================
SUBROUTINE sBase_alloc( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  ASSOCIATE(nElems=>sf%grid%nElems, degGP=>sf%degGP, deg=>sf%deg, nBase =>sf%nBase)
  ALLOCATE(sf%xiGP(     0:degGP))
  ALLOCATE(sf%wGPloc(   0:degGP))
  ALLOCATE(sf%wGP(      0:degGP      ,1:nElems))
  ALLOCATE(sf%s_GP(     0:degGP      ,1:nElems))
  ALLOCATE(sf%baseGP(   0:degGP,0:deg,1:nElems))
  ALLOCATE(sf%base_dsGP(0:degGP,0:deg,1:nElems))
  ALLOCATE(sf%base_offset(            1:nElems))
  ALLOCATE(sf%base_dsAxis(1:deg+1        ,0:deg))
  ALLOCATE(sf%base_dsEdge(nBase-deg:nBase,0:deg))
  ALLOCATE(sf%A_Axis(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(sf%A_Edge(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(sf%R_Axis(0:deg,0:deg,1:NBC_TYPES))
  ALLOCATE(sf%R_Edge(0:deg,0:deg,1:NBC_TYPES))
  sf%xiGP        =0.0_wp
  sf%wGPloc      =0.0_wp            
  sf%wGP         =0.0_wp            
  sf%s_GP        =0.0_wp            
  sf%base_offset =-1
  sf%baseGP      =0.0_wp           
  sf%base_dsGP   =0.0_wp        
  sf%base_dsAxis =0.0_wp    
  sf%base_dsEdge =0.0_wp    
  sf%A_Axis      =0.0_wp         
  sf%A_Edge      =0.0_wp         
  sf%R_Axis      =0.0_wp         
  sf%R_Edge      =0.0_wp         
  DO i=0,deg
    sf%A_Axis(i,i,:)=1.0_wp
    sf%A_Edge(i,i,:)=1.0_wp
    sf%R_Axis(i,i,:)=1.0_wp
    sf%R_Edge(i,i,:)=1.0_wp
  END DO
  END ASSOCIATE !nElems, degGP, deg
END SUBROUTINE sbase_alloc


!===================================================================================================================================
!> finalize the type sBase
!!
!===================================================================================================================================
SUBROUTINE sBase_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  !pointers, classes
  NULLIFY(sf%grid)
  IF(sf%continuity.EQ.sf%deg-1)THEN
    CALL sf%interpol%free() 
    CALL sf%spline%free() 
    CALL sf%bspl%free() 
  END IF
  !allocatables  
  SDEALLOCATE(sf%xiGP)
  SDEALLOCATE(sf%wGPloc)
  SDEALLOCATE(sf%wGP)
  SDEALLOCATE(sf%s_GP)
  SDEALLOCATE(sf%s_IP)
  SDEALLOCATE(sf%base_offset)
  SDEALLOCATE(sf%baseGP)   
  SDEALLOCATE(sf%base_dsGP)
  SDEALLOCATE(sf%base_dsAxis) 
  SDEALLOCATE(sf%base_dsEdge) 
  SDEALLOCATE(sf%A_Axis) 
  SDEALLOCATE(sf%R_Axis) 
  SDEALLOCATE(sf%A_Edge) 
  SDEALLOCATE(sf%R_Edge) 
  
  sf%continuity =-99
  sf%deg        =-1 
  sf%degGP      =-1
  sf%nGP        =-1        
  sf%nbase      =-1
  sf%initialized=.FALSE.

END SUBROUTINE sBase_free


!===================================================================================================================================
!> copy the type sbase
!!
!===================================================================================================================================
SUBROUTINE sBase_copy( sf , tocopy)
! MODULES
USE MOD_GLobals, ONLY: Unit_stdOut,abort
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(IN   ) :: tocopy
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPEIS(t_sbase)
  IF(.NOT.tocopy%initialized) THEN
    CALL abort(__STAMP__, &
        "sBase_copy: not initialized sBase from which to copy!")
  END IF
  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase type!'
    CALL sf%free()
  END IF

  CALL sf%init(tocopy%grid,tocopy%deg,tocopy%continuity,tocopy%degGP)

  END SELECT !TYPE
END SUBROUTINE sbase_copy


END MODULE MOD_sBase

