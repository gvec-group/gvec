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
USE MOD_Globals                  ,ONLY: wp,Unit_stdOut,abort
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
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: deg                      !! input parameter: degree of Spline/polynomial 
  INTEGER              :: degGP                    !! number of Gauss-points (degGP+1) per element >= deg
  INTEGER              :: continuity               !! input parameter: full spline (=deg-1) or discontinuous (=-1)
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: nGP                      !! total number of gausspoints = degGP*nElems
  INTEGER              :: nbase                    !! total number of degree of freedom / global basis functions
  CONTAINS
    PROCEDURE(i_sub_sbase_init    ),DEFERRED :: init
    PROCEDURE(i_sub_sbase_free    ),DEFERRED :: free
    PROCEDURE(i_sub_sbase_copy    ),DEFERRED :: copy
    PROCEDURE(i_sub_sbase_eval    ),DEFERRED :: eval
    PROCEDURE(i_fun_sbase_initDOF ),DEFERRED :: initDOF

END TYPE c_sbase

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sbase_init( sf ,grid_in,degGP_in)
    IMPORT wp, c_sbase,t_sgrid
    CLASS(c_sbase), INTENT(INOUT)        :: sf
    CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in
    INTEGER       , INTENT(IN   )        :: degGP_in
  END SUBROUTINE i_sub_sbase_init

  SUBROUTINE i_sub_sbase_free( sf ) 
    IMPORT c_sbase
    CLASS(c_sbase), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sbase_free

  SUBROUTINE i_sub_sbase_copy( sf, tocopy ) 
    IMPORT c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: tocopy
    CLASS(c_sbase), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sbase_copy

  SUBROUTINE i_sub_sBase_eval( sf , x, deriv,iElem,base_x)
    IMPORT wp,c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: sf !! self
    REAL(wp)      , INTENT(IN   ) :: x
    INTEGER       , INTENT(IN   ) :: deriv 
    INTEGER       , INTENT(  OUT) :: iElem
    REAL(wp)      , INTENT(  OUT) :: base_x(:)
  END SUBROUTINE i_sub_sbase_eval

  FUNCTION i_fun_sbase_initDOF( sf, g_IP ) RESULT(DOFs) 
    IMPORT wp,c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: sf
    REAL(wp)      , INTENT(IN   ) :: g_IP(1:sf%nBase)
    REAL(wp)                      :: DOFs(1:sf%nBase)
  END FUNCTION i_fun_sbase_initDOF
END INTERFACE
 


TYPE,EXTENDS(c_sbase) :: t_sBase
  LOGICAL :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  CLASS(t_sgrid),POINTER :: grid  => NULL()        !! pointer to grid 
  !---------------------------------------------------------------------------------------------------------------------------------
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
  INTEGER ,ALLOCATABLE :: nDOF_BC(:)               !! number of boudnary dofs involved in bc of BC_TYPE, size(NBC_TYPES)
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
  PROCEDURE :: free          => sBase_free
  PROCEDURE :: copy          => sBase_copy
  PROCEDURE :: eval          => sBase_eval
  PROCEDURE :: initDOF       => sBase_initDOF

END TYPE t_sBase

TYPE,EXTENDS(t_sbase) :: t_sBase_disc
  REAL(wp),ALLOCATABLE :: xiIP(:)                  !! element local interpolation points (continuity =-1)
  REAL(wp),ALLOCATABLE :: wbaryIP(:)               !! element local interpolation points (continuity =-1)
  REAL(wp),ALLOCATABLE :: DmatIP(:,:)              !! element local interpolation points (continuity =-1)
END TYPE t_sBase_disc

TYPE,EXTENDS(t_sbase) :: t_sBase_spl
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl        !! contains bspline functions
  TYPE(sll_t_spline_1d)             :: spline      !! contains 1d spline functions
  TYPE(sll_t_spline_interpolator_1d):: Interpol    !! spline interpolator
END TYPE t_sBase_spl

LOGICAL  :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type sbase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE sBase_new( sf,deg_in,continuity_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   )        :: deg_in        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity_in !! continuity: 
                                                        !! 0: disc. polynomial
                                                        !! deg-1: spline with cont. deg-1
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(continuity_in.EQ.-1)THEN
    ALLOCATE(t_sbase_disc :: sf)
  ELSEIF(continuity_in.EQ.deg_in-1)THEN
    ALLOCATE(t_sbase_spl :: sf)
  ELSE
    CALL abort(__STAMP__,&
        " error in sbase new: continuity only full (deg-1) or discontinuous (-1) !") 
  END IF
  sf%deg        =deg_in
  sf%continuity =continuity_in
END SUBROUTINE sbase_new

!===================================================================================================================================
!> initialize the type sbase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE sBase_init( sf, grid_in,degGP_in)
! MODULES
USE MOD_GLobals, ONLY: PI
USE MOD_LinAlg , ONLY: INV
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
  INTEGER       , INTENT(IN   )        :: degGP_in      !! gauss quadrature points: nGP=degGP+1 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i,iGP,iElem,imin,jmin
  REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: locbasis,Vdm,DmatGP
  INTEGER  :: iBC,j,diri,even 
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A,3(A,I3),A)')'INIT sBase type:', &
       ' degree= ',sf%deg, &
       ' gauss points per elem = ',degGP_in, &
       ' continuity= ',sf%continuity, ' ...'
  IF(sf%initialized) THEN
    CALL abort(__STAMP__, &
        "Trying to reinit sbase!") 
  END IF
  IF(degGP_in.LT.sf%deg) &
    CALL abort(__STAMP__, &
        "error in sbase: degGP must be > deg!") 
  SELECT TYPE(sf)
  TYPEIS(t_sbase_disc)
    IF(sf%continuity.NE.-1) &
      CALL abort(__STAMP__, &
          "error in sbase init: type is disc but continuity is not -1, mabye sbase_new was not called before!") 
  TYPEIS(t_sbase_spl)
    IF(sf%continuity.NE.sf%deg-1) &
      CALL abort(__STAMP__, &
          "error in sbase init: type is spl but continuity is not deg-1, mabye sbase_new was not called before!") 
  END SELECT !Type
  sf%grid       => grid_in
  sf%degGP      =  degGP_in

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

  SELECT TYPE(sf)
  TYPEIS(t_sbase_disc)   
    ALLOCATE(Vdm(   0:degGP,0:deg))
    ALLOCATE(DmatGP(0:degGP,0:deg))
    !  use chebychev-lobatto points for interpolation (closed form!), interval [-1,1] 
    DO i=0,deg
      sf%xiIP(i)=-COS(REAL(i,wp)/REAL(deg,wp)*PI)
    END DO
    CALL BarycentricWeights(deg,sf%xiIP,sf%wBaryIP)
    CALL InitializeVandermonde(deg,degGP,sf%wBaryIP,sf%xiIP,sf%xiGP,Vdm)
    CALL PolynomialDerivativeMatrix(deg,sf%xiIP,sf%DmatIP)
    DmatGP=MATMUL(Vdm,sf%DmatIP)
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
    sf%base_dsAxis(1:deg+1        ,1)=sf%DmatIP(  0,:)*(2.0_wp/grid%ds(1))
    sf%base_dsEdge(nBase-deg:nBase,1)=sf%DmatIP(deg,:)*(2.0_wp/grid%ds(nElems))
    !  higher derivatives 
    DO i=2,deg
      sf%base_dsAxis(1:deg+1        ,i)=MATMUL(TRANSPOSE(sf%DmatIP),sf%base_dsAxis(:,i-1))*(2.0_wp/grid%ds(1))
      sf%base_dsEdge(nBase-deg:nBase,i)=MATMUL(TRANSPOSE(sf%DmatIP),sf%base_dsEdge(:,i-1))*(2.0_wp/grid%ds(nElems))
    END DO
    !interpolation:
    !  points are repeated at element interfaces (discontinuous)
    ALLOCATE(sf%s_IP(sf%nBase)) !for spl, its allocated elsewhere...
    DO iElem=1,nElems
      sf%s_IP(1+(deg+1)*(iElem-1):(deg+1)*iElem)=grid%sp(iElem-1)+0.5_wp*(sf%xiIP+1.0_wp)*grid%ds(iElem)
    END DO !iElem 
    DEALLOCATE(Vdm,DmatGP)
  TYPEIS(t_sbase_spl)   
    ALLOCATE(locbasis(0:deg,0:deg))
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

    DEALLOCATE(locbasis)
  END SELECT !TYPE is t_sbase_disc/spl

  !TODO STRONG BOUNDARY CONDITIONS: A and R matrices

   
  DO iBC=1,NBC_TYPES
    ASSOCIATE(nD=>sf%nDOF_BC(iBC))
    SELECT CASE(iBC)
    CASE(BC_TYPE_OPEN)     ; nD =0
    CASE(BC_TYPE_NEUMANN)  ; nD =1          ;   diri=0 ;  even=0
    CASE(BC_TYPE_DIRICHLET); nD =1          ;   diri=1 ;  even=0 !even not used
    CASE(BC_TYPE_SYMM)     ; nD =(deg+1)/2  ;   diri=0 ;  even=0
    CASE(BC_TYPE_SYMMZERO) ; nD =(deg+1)/2+1;   diri=1 ;  even=0
    CASE(BC_TYPE_ANTISYMM) ; nD =deg/2+1    ;   diri=1 ;  even=1
    END SELECT !iBC 
    !A and R are already initialized unit matrices!!
    IF(nD.GT.0)THEN
      DO i=diri+1,nD
        j=2*(i-diri)-(1-even) !even=0 odd derivs, even=1 even derivatives
        sf%A_Axis(        i,:,iBC)=sf%base_dsAxis(:,j)/sf%base_dsAxis(i,j)
        sf%A_Edge(nBase+1-i,:,iBC)=sf%base_dsEdge(:,j)/sf%base_dsEdge(nBase+1-i,j)
      END DO
      !invert BC part
      sf%R_Axis(1:nD,1:nD,iBC)=INV(sf%A_Axis(1:nD,1:nD,iBC))
      sf%R_Axis(:,:,iBC)=TRANSPOSE(MATMUL(sf%R_Axis(:,:,iBC),sf%A_Axis(:,:,iBC)))
      DO i=1,deg+1; DO j=1,i-1
        sf%R_Axis(i,j,iBC)=-sf%R_Axis(i,j,iBC)
      END DO; END DO
      sf%R_Edge(nBase-nD+1:nBase,nBase-nD+1:nBase,iBC)=INV(sf%A_Edge(nBase-nD+1:nBase,nBase-nD+1:nBase,iBC))
      sf%R_Edge(:,:,iBC)=TRANSPOSE(MATMUL(sf%R_Edge(:,:,iBC),sf%A_Edge(:,:,iBC)))
      DO i=nBase-deg,nBase; DO j=i+1,nBase
        sf%R_Edge(i,j,iBC)=-sf%R_Edge(i,j,iBC)
      END DO; END DO
    END IF
    END ASSOCIATE !nD=>nDOF_BC(iBC)
  END DO!iBC=1,NBC_TYPES
  END ASSOCIATE !sf


  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
  IF(.NOT.test_called) CALL sBase_test(sf)

END SUBROUTINE sBase_init


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
  ALLOCATE(sf%nDOF_BC(            1:NBC_TYPES))
  ALLOCATE(sf%A_Axis(1:deg+1,1:deg+1,1:NBC_TYPES))
  ALLOCATE(sf%R_Axis(1:deg+1,1:deg+1,1:NBC_TYPES))
  ALLOCATE(sf%A_Edge(nBase-deg:nBase,nBase-deg:nBase,1:NBC_TYPES))
  ALLOCATE(sf%R_Edge(nBase-deg:nBase,nBase-deg:nBase,1:NBC_TYPES))
  sf%xiGP        =0.0_wp
  sf%wGPloc      =0.0_wp            
  sf%wGP         =0.0_wp            
  sf%s_GP        =0.0_wp            
  sf%base_offset =-1
  sf%baseGP      =0.0_wp           
  sf%base_dsGP   =0.0_wp        
  sf%base_dsAxis =0.0_wp    
  sf%base_dsEdge =0.0_wp    
  sf%nDOF_BC     =0
  sf%A_Axis      =0.0_wp         
  sf%A_Edge      =0.0_wp         
  sf%R_Axis      =0.0_wp         
  sf%R_Edge      =0.0_wp         
  DO i=0,deg
    sf%A_Axis(1+i,1+i,:)=1.0_wp
    sf%R_Axis(1+i,1+i,:)=1.0_wp
    sf%A_Edge(nBase-deg+i,nBase-deg+i,:)=1.0_wp
    sf%R_Edge(nBase-deg+i,nBase-deg+i,:)=1.0_wp
  END DO
  SELECT TYPE(sf)
  TYPEIS(t_sbase_disc)
    ALLOCATE(sf%xiIP(   0:deg))
    ALLOCATE(sf%wbaryIP(0:deg))
    ALLOCATE(sf%DmatIP( 0:deg,0:deg))
    sf%xiIP   =0.0_wp
    sf%wbaryIP=0.0_wp
    sf%DmatIP =0.0_wp
  END SELECT !TYPE is t_sbase_disc
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
  IF(.NOT.sf%initialized) RETURN
  !pointers, classes
  NULLIFY(sf%grid)
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
  SDEALLOCATE(sf%nDOF_BC)
  SDEALLOCATE(sf%A_Axis) 
  SDEALLOCATE(sf%R_Axis) 
  SDEALLOCATE(sf%A_Edge) 
  SDEALLOCATE(sf%R_Edge) 
  SELECT TYPE (sf) 
  TYPEIS(t_sbase_spl)
    CALL sf%interpol%free() 
    CALL sf%spline%free() 
    CALL sf%bspl%free() 
  TYPEIS(t_sbase_disc)
    SDEALLOCATE(sf%xiIP)
    SDEALLOCATE(sf%wbaryIP)
    SDEALLOCATE(sf%DmatIP)
  END SELECT !TYPE
  
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
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_sBase), INTENT(IN   ) :: tocopy
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
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of sBase in copy!'
    CALL sf%free()
  END IF
  sf%deg=tocopy%deg
  sf%continuity=tocopy%continuity
  CALL sf%init(tocopy%grid,tocopy%degGP)

  END SELECT !TYPE
END SUBROUTINE sbase_copy


!===================================================================================================================================
!> evaluate sbase at position x [0,1], NOT EFFICIENT!!
!!
!===================================================================================================================================
SUBROUTINE sBase_eval( sf , x,deriv,iElem,base_x)
! MODULES
USE MOD_Basis1D, ONLY:LagrangeInterpolationPolys
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  REAL(wp)      , INTENT(IN   ) :: x     !! position [0,1] where to evaluate
  INTEGER       , INTENT(IN   ) :: deriv !! 0: evaluation,1: 1st derivative
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER       , INTENT(  OUT) :: iElem
  REAL(wp)      , INTENT(  OUT) :: base_x(:)  !! all basis functions (0:deg) evaluated
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: ideriv
  REAL(wp):: xiloc,baseloc(0:sf%deg,0:sf%deg)
!===================================================================================================================================

SELECT TYPE(sf)
TYPEIS(t_sbase_disc)
  iElem=sf%grid%find_elem(x)
  xiloc  =(x-sf%grid%sp(iElem-1))*2.0_wp/sf%grid%ds(iElem)-1.0_wp !in [-1,1]

  IF(deriv.EQ.0)THEN
    CALL LagrangeInterpolationPolys(xiloc,sf%deg,sf%xiIP,sf%wBaryIP,base_x(:))
  ELSE
    CALL LagrangeInterpolationPolys(xiloc,sf%deg,sf%xiIP,sf%wBaryIP,baseloc(:,0))
    DO ideriv=1,deriv
      baseloc(:,ideriv)=MATMUL(TRANSPOSE(sf%DmatIP),baseloc(:,ideriv-1))*(2.0_wp/sf%grid%ds(iElem))
    END DO
    base_x=baseloc(:,deriv)
  END IF!deriv
TYPEIS(t_sbase_spl)
  IF(deriv.EQ.0)THEN
    CALL sf%bspl%eval_basis(x,base_x(:),iElem)
  ELSE
    CALL sf%bspl%eval_basis_and_n_derivs(x,deriv,baseloc,iElem)
    base_x=baseloc(:,deriv)
  END IF

CLASS DEFAULT
  CALL abort(__STAMP__, &
    "this type of continuity not implemented!")
END SELECT !TYPE

END SUBROUTINE sbase_eval


!===================================================================================================================================
!>  take values interpolated at sf%s_IP positions and give back the degrees of freedom 
!!
!===================================================================================================================================
FUNCTION sBase_initDOF( sf , g_IP) RESULT(DOFs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  REAL(wp)      , INTENT(IN   ) :: g_IP(sf%nBase)  !!  interpolation values at s_IP positions [0,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: DOFs(sf%nBase)  !! result of interpolation 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT TYPE(sf)
TYPEIS(t_sbase_disc)
  DOFs(:)=g_IP
TYPEIS(t_sbase_spl)
  CALL sf%interpol%compute_interpolant( sf%spline, g_IP )
  DOFs(:)=sf%spline%bcoef(:) !somewhat not perfect, since interpolant saves result to bcoef of spline 
CLASS DEFAULT
  CALL abort(__STAMP__, &
    "this type of continuity not implemented!")
END SELECT !TYPE

END FUNCTION sbase_initDOF

!===================================================================================================================================
!> apply strong boundary conditions at axis and edge 
!!
!===================================================================================================================================
SUBROUTINE sBase_applyBCtoDOF(sf ,DOFs,BC_Type,BC_val)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  INTEGER       , INTENT(IN   ) :: BC_Type(2)           !! bc type on axis (1) and edge (2)
  REAL(wp)      , INTENT(IN   ) :: BC_Val(2)           !! for dirichlet BC : value
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)      , INTENT(INOUT) :: DOFs(sf%nBase)  !! DOFs with boundary conditions applied 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp):: rhs(0:sf%deg)
!===================================================================================================================================
  SELECT CASE(BC_Type(1)) !AXIS
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    DOFs(1)=BC_Val(1)
  CASE(BC_TYPE_NEUMANN,BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM)
    !DOFs(1:sf%deg+1)=
    STOP' only open and dirichlet BC are yet implemented'
  END SELECT !BCtype(1) !AXIS

  SELECT CASE(BC_Type(2)) !Edge
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    DOFs(sf%nBase)=BC_Val(2)
  CASE(BC_TYPE_NEUMANN,BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM)
    !DOFs(sf%nBase-deg:sf%nBase)=
    STOP' only open and dirichlet BC are yet implemented'
  END SELECT !BCtype(2) !Edge

END SUBROUTINE sbase_applyBCtoDOF

!===================================================================================================================================
!> test sbase variable
!!
!===================================================================================================================================
SUBROUTINE sBase_test( sf)
! MODULES
USE MOD_GLobals, ONLY: testdbg,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iTest,iElem,jElem
  REAL(wp):: x,base_x(0:sf%deg) 
!===================================================================================================================================
  test_called=.TRUE.
  IF(testlevel.LE.0) RETURN
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN SBASE TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.LE.1)THEN
    iTest=101
    IF(testdbg.OR.(.NOT.( ABS(sf%s_IP(1)).LT. 1.0e-12_wp))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(2(A,I4),(A,E11.3))') &
      '   degree = ',sf%deg, '   continuity = ',sf%continuity, &
      '\n =>  should be 0.0 : sIP(1)= ',sf%s_IP(1)
    END IF
    iTest=102
    IF(testdbg.OR.(.NOT.( ABS(sf%s_IP(sf%nBase)-1.0_wp).LT. 1.0e-12_wp))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),(A,E11.3))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
      '\n =>  should be 1.0 : sIP(nBase)= ',sf%s_IP(sf%nBase)
    END IF
    iTest=103
    IF(testdbg.OR.(.NOT.( ABS(SUM(sf%wGP)-1.0_wp).LT. 1.0e-12_wp))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),(A,E11.3))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
      '\n => should be 1.0 : SUM(wGP)= ',SUM(sf%wGP)
    END IF
    iTest=104 
    CALL sf%eval(0.0_wp,0,iElem,base_x) 
    IF(testdbg.OR.(.NOT.((ABS(base_x(0)-1.).LT.1.0e-12).AND. &
                         ( SUM(ABS(base_x(1:sf%deg))).LT. 1.0e-12_wp).AND. &
                         (iElem.EQ.1) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),A,I6,*(E11.3)))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
      '\n => should be 1 and (1,0,0...) : base_x(x=0)= ',iElem,base_x
    END IF
    iTest=105
    CALL sf%eval(1.0_wp,0,iElem,base_x) 
    IF(testdbg.OR.(.NOT.((ABS(base_x(sf%deg)-1.).LT.1.0e-12).AND. &
                         ( SUM(ABS(base_x(0:sf%deg-1))).LT. 1.0e-12_wp).AND. &
                         (iElem.EQ.sf%grid%nElems)))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),A,*(E11.3)))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
      '\n => should be ...,0,0,1 : base_x(x=0)= ',base_x
    END IF
    iTest=106
    jElem=sf%grid%nElems/2
    x=0.5_wp*(sf%grid%sp(jElem-1)+sf%grid%sp(jElem))
    CALL sf%eval(x,0,iElem,base_x) 
    IF(testdbg.OR.(.NOT.(iElem.EQ.jElem)))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),2(A,I6))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
      '\n => should be ',jElem,': iElem= ' , iElem
    END IF
    iTest=107
    jElem=MIN(2,sf%grid%nElems)
    x=sf%grid%sp(jElem-1)+0.01_wp*sf%grid%ds(jElem)
    CALL sf%eval(x ,0,iElem,base_x) 
    IF(testdbg.OR.(.NOT.(iElem.EQ.jElem)))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,' FAILED !!'
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(3(A,I4),2(A,I6))') &
      '   degree = ',sf%deg, &
      '   continuity = ',sf%continuity,  ', nBase= ',sf%nBase, &
      '\n => should be ',jElem,': iElem= ' , iElem
    END IF
  END IF !testlevel<1
  !IF(testlevel.LE.2)THEN
  !  
  !END IF !testlevel<2
  test_called=.FALSE.
   
END SUBROUTINE sbase_test


END MODULE MOD_sBase

