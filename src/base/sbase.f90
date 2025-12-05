!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** sBase **
!!
!! 1D basis in radial coordinate "s". Contains sbase type definition and associated routines
!!
!===================================================================================================================================
MODULE MODgvec_sBase
! MODULES
USE MODgvec_Globals                  ,ONLY: wp,Unit_stdOut,abort,MPIRoot
USE sll_m_bsplines               ,ONLY: sll_c_bsplines
USE sll_m_spline_interpolator_1d ,ONLY: sll_t_spline_interpolator_1d
USE sll_m_spline_matrix          ,ONLY: sll_c_spline_matrix
USE MODgvec_sGrid ,ONLY: c_sgrid,t_sgrid
IMPLICIT NONE

PRIVATE
PUBLIC t_sbase,sbase_new

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
  INTEGER              :: nGP                      !! global number of gausspoints = (degGP+1)*nElems
  INTEGER              :: nGP_str, nGP_end         !! local number of gausspoints = (degGP+1)*nElems per MPI subdomain
  INTEGER              :: nbase                    !! total number of degree of freedom / global basis functions
  CLASS(sll_c_spline_matrix),ALLOCATABLE :: mass
  CONTAINS
    PROCEDURE(i_sub_sbase_init          ),DEFERRED :: init
    PROCEDURE(i_sub_sbase_free          ),DEFERRED :: free
    PROCEDURE(i_sub_sbase_copy          ),DEFERRED :: copy
    PROCEDURE(i_sub_sbase_compare       ),DEFERRED :: compare
    PROCEDURE(i_sub_sbase_change_base   ),DEFERRED :: change_base
    PROCEDURE(i_sub_sbase_eval          ),DEFERRED :: eval
    PROCEDURE(i_fun_sbase_evalDOF_s     ),DEFERRED :: evalDOF_s
    PROCEDURE(i_fun_sbase_evalDOF2D_s   ),DEFERRED :: evalDOF2D_s
    PROCEDURE(i_fun_sbase_evalDOF_base  ),DEFERRED :: evalDOF_base
    PROCEDURE(i_fun_sbase_evalDOF_GP    ),DEFERRED :: evalDOF_GP
    PROCEDURE(i_fun_sbase_initDOF       ),DEFERRED :: initDOF
    PROCEDURE(i_sub_sbase_applyBCtoDOF  ),DEFERRED :: applyBCtoDOF
    PROCEDURE(i_sub_sbase_applyBCtoRHS  ),DEFERRED :: applyBCtoRHS

END TYPE c_sbase

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sbase_init( sf ,deg_in,continuity_in,grid_in,degGP_in)
    IMPORT wp, c_sbase,t_sgrid
    CLASS(c_sbase), INTENT(INOUT)        :: sf
    INTEGER       , INTENT(IN   )        :: deg_in
    INTEGER       , INTENT(IN   )        :: continuity_in
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

  SUBROUTINE i_sub_sbase_compare( sf, tocompare, is_same, cond_out )
    IMPORT c_sbase
    CLASS(c_sbase)  , INTENT(IN   ) :: sf
    CLASS(c_sbase)  , INTENT(IN   ) :: tocompare
    LOGICAL,OPTIONAL, INTENT(  OUT) :: is_same
    LOGICAL,OPTIONAL, INTENT(  OUT) :: cond_out(:)
  END SUBROUTINE i_sub_sbase_compare

  SUBROUTINE i_sub_sbase_change_base( sf, old_sBase, iterDim, old_data, sf_data )
    IMPORT c_sbase,wp
    CLASS(c_sbase) , INTENT(IN   ) :: sf
    CLASS(c_sbase) , INTENT(IN   ) :: old_sBase
    INTEGER        , INTENT(IN   ) :: iterDim
    REAL(wp)       , INTENT(IN   ) :: old_data(:,:)
    REAL(wp)       , INTENT(  OUT) :: sf_data(:,:)
  END SUBROUTINE i_sub_sbase_change_base

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
    REAL(wp)      , INTENT(IN   ) :: g_IP(:)
    REAL(wp)                      :: DOFs(1:sf%nBase)
  END FUNCTION i_fun_sbase_initDOF

  FUNCTION i_fun_sbase_evalDOF_s( sf, x,deriv,DOFs ) RESULT(y)
    IMPORT wp,c_sbase
  CLASS(c_sbase), INTENT(IN   ) :: sf
  REAL(wp)      , INTENT(IN   ) :: x
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)
  REAL(wp)                      :: y
  END FUNCTION i_fun_sbase_evalDOF_s

  FUNCTION i_fun_sbase_evalDOF2D_s( sf, x,nd,deriv,DOFs ) RESULT(y)
    IMPORT wp,c_sbase
  CLASS(c_sbase), INTENT(IN   ) :: sf
  REAL(wp)      , INTENT(IN   ) :: x
  INTEGER       , INTENT(IN   ) :: nd
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%nBase,1:nd)
  REAL(wp)                      :: y(1:nd)
  END FUNCTION i_fun_sbase_evalDOF2D_s

  FUNCTION i_fun_sbase_evalDOF_base( sf, iElem,base_x,DOFs ) RESULT(y)
    IMPORT wp,c_sbase
  CLASS(c_sbase), INTENT(IN   ) :: sf
  INTEGER       , INTENT(IN   ) :: iElem
  REAL(wp)      , INTENT(IN   ) :: base_x(0:sf%deg)
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)
  REAL(wp)                      :: y  !out
  END FUNCTION i_fun_sbase_evalDOF_base

  FUNCTION i_fun_sbase_evalDOF_GP( sf,deriv,DOFs ) RESULT(y_GP)
    IMPORT wp,c_sbase
  CLASS(c_sbase), INTENT(IN   ) :: sf
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)
  REAL(wp)                      :: y_GP(sf%nGP)
  END FUNCTION i_fun_sbase_evalDOF_GP

  SUBROUTINE i_sub_sBase_applyBCtoDOF( sf ,DOFs,BC_Type,BC_val)
    IMPORT wp,c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: sf
    INTEGER       , INTENT(IN   ) :: BC_Type(2)
    REAL(wp)      , INTENT(IN   ) :: BC_Val(2)
    REAL(wp)      , INTENT(INOUT) :: DOFs(1:sf%nBase)
  END SUBROUTINE i_sub_sBase_applyBCtoDOF

  SUBROUTINE i_sub_sBase_applyBCtoRHS( sf ,RHS,BC_Type)
    IMPORT wp,c_sbase
    CLASS(c_sbase), INTENT(IN   ) :: sf
    INTEGER       , INTENT(IN   ) :: BC_Type(2)
    REAL(wp)      , INTENT(INOUT) :: RHS(1:sf%nBase)
  END SUBROUTINE i_sub_sBase_applyBCtoRHS

END INTERFACE



TYPE,EXTENDS(c_sbase) :: t_sBase
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  CLASS(t_sgrid),POINTER :: grid                   !! pointer to grid
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),ALLOCATABLE :: xi_GP(:)                  !! element local gauss point positions for interval [-1,1], size(0:degGP)
  REAL(wp),ALLOCATABLE :: w_GPloc(:)                !! element local gauss weights for interval [-1,1], size(0:degGP)
  REAL(wp),ALLOCATABLE :: w_GP(:)                   !! global radial integration weight size((degGP+1)*nElems)
  REAL(wp),ALLOCATABLE :: s_GP(:)                  !! global position of gauss points  in s [0,1] , size((degGP+1)*nElems)
  REAL(wp),ALLOCATABLE :: s_IP(:)                  !! position of interpolation points for initialization, size(nBase)
  INTEGER ,ALLOCATABLE :: base_offset(:)           !! offset of 0:deg element local basis functions to global index of
                                                   !! degree of freedom, allocated (1:nElems). iBase = offset(iElem)+j, j=0...deg
  REAL(wp),ALLOCATABLE :: base_GP(:,:,:)            !! basis functions, (0:degGP,0:deg,1:nElems),
  REAL(wp),ALLOCATABLE :: base_ds_GP(:,:,:)         !! s derivative of basis functions, (0:degGP,0:deg,1:nElems)
  REAL(wp),ALLOCATABLE :: base_dsAxis(:,:)         !! all derivatives 1..deg of all basis functions at axis size(1:deg+1,0:deg)
  REAL(wp),ALLOCATABLE :: base_dsEdge(:,:)         !! all derivatives 1..deg of all basis functions at edge size(nBase-deg:nBase,0:deg)
  INTEGER ,ALLOCATABLE :: nDOF_BC(:)               !! number of boundary dofs involved in bc of BC_TYPE, size(-(deg+1):NBC_TYPES)
  REAL(wp),ALLOCATABLE :: A_Axis(:,:,:)            !! matrix to apply boundary conditions after interpolation (direct)
  REAL(wp),ALLOCATABLE :: invA_Axis(:,:,:)         !! inverse of A_Axis
  REAL(wp),ALLOCATABLE :: R_Axis(:,:,:)            !! matrix to apply boundary conditions for RHS (testfunction)
                                                   !! size(1:deg+1,1:deg+1,-(deg+1):NBC_TYPES)
  REAL(wp),ALLOCATABLE :: AR_Axis(:,:,:)           !! matrix to apply boundary conditions via LGM with identity
                                                   !! size(1:deg+1,1:deg+1,-(deg+1):NBC_TYPES)
  REAL(wp),ALLOCATABLE :: A_Edge(:,:,:)            !! matrix to apply boundary conditions after interpolation (direct)
  REAL(wp),ALLOCATABLE :: invA_Edge(:,:,:)         !! inverse of A_Edge
  REAL(wp),ALLOCATABLE :: R_Edge(:,:,:)            !! matrix to apply boundary conditions for RHS
                                                   !! size(nBase-deg:nBase,nBase-deg:nBase,NBC_TYPES)
  REAL(wp),ALLOCATABLE :: AR_Edge(:,:,:)           !! matrix to apply boundary conditions via LGM with identity
                                                   !! size(nBase-deg:nBase,nBase-deg:nBase,NBC_TYPES)
                                                   !! possible Boundary conditions
                                                   !! 1=BC_TYPE_OPEN      : boundary left open
                                                   !! 2=BC_TYPE_NEUMANN   : 1st deriv fixed
                                                   !! 3=BC_TYPE_DIRICHLET : sol fixed          (typically used at edge)
                                                   !! 4=BC_TYPE_SYMM      : all odd derivs = 0
                                                   !! 5=BC_TYPE_SYMMZERO  : all even derivs & sol = 0
                                                   !! 6=BC_TYPE_ANTISYMM  : all odd derivs & sol = 0
                                                   !! <=0: m-dependent BC:
                                                   !!  0, m=0:  BC_TYPE_SYMM
                                                   !! -1, m=1:  BC_TYPE_ANTISYMM
                                                   !! -2, m=2:  BC_TYPE_SYMMZERO
                                                   !! ... odd m<=deg: ANTISYMM + all derivatives 1,...,m-1 =0
                                                   !! ...even m<=deg: SYMMZERO + all derivatives 1,...,m-1 =0

  CONTAINS

  PROCEDURE :: init          => sBase_init
  PROCEDURE :: free          => sBase_free
  PROCEDURE :: copy          => sBase_copy
  PROCEDURE :: compare       => sBase_compare
  PROCEDURE :: change_base   => sBase_change_base
  PROCEDURE :: eval          => sBase_eval
  PROCEDURE :: evalDOF_s     => sBase_evalDOF_s
  PROCEDURE :: evalDOF2D_s   => sBase_evalDOF2D_s
  PROCEDURE :: evalDOF_base  => sBase_evalDOF_base
  PROCEDURE :: evalDOF_GP    => sBase_evalDOF_GP
  PROCEDURE :: initDOF       => sBase_initDOF
  PROCEDURE :: applyBCtoDOF  => sBase_applyBCtoDOF_LGM
  PROCEDURE :: applyBCtoRHS  => sBase_applyBCtoRHS

END TYPE t_sBase

TYPE,EXTENDS(t_sbase) :: t_sBase_disc
  REAL(wp),ALLOCATABLE :: xiIP(:)                  !! element local interpolation points in [-1,1] (continuity =-1) size(0:deg)
  REAL(wp),ALLOCATABLE :: wbaryIP(:)               !! barycentric weights for xiIP size(0:deg)
  REAL(wp),ALLOCATABLE :: DmatIP(:,:)              !! Nodal derivative matrix size(0:deg,0:deg)
END TYPE t_sBase_disc

TYPE,EXTENDS(t_sbase) :: t_sBase_spl
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl        !! contains bspline functions
  TYPE(sll_t_spline_interpolator_1d):: Interpol    !! spline interpolator
END TYPE t_sBase_spl

LOGICAL, PRIVATE  :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type sbase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE sBase_new(sbase_in,deg_in,continuity_in,grid_in,degGP_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   )        :: deg_in        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity_in !! continuity:
                                                        !! 0: disc. polynomial
                                                        !! deg-1: spline with cont. deg-1
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in       !! grid information
  INTEGER       , INTENT(IN   )        :: degGP_in      !! gauss quadrature points: nGP=degGP+1
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), ALLOCATABLE,INTENT(INOUT)        :: sbase_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(continuity_in.EQ.-1)THEN
    ALLOCATE(t_sbase_disc :: sbase_in)
  ELSEIF(continuity_in.EQ.deg_in-1)THEN
    ALLOCATE(t_sbase_spl :: sbase_in)
  ELSE
    CALL abort(__STAMP__,&
        " error in sbase new: continuity only full (deg-1) or discontinuous (-1) !",&
        TypeInfo="InvalidParameterError")
  END IF

  CALL sbase_in%init(deg_in,continuity_in,grid_in,degGP_in)

END SUBROUTINE sbase_new

!===================================================================================================================================
!> initialize the type sbase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE sBase_init( sf,deg_in,continuity_in,grid_in,degGP_in)
! MODULES
USE MODgvec_GLobals,   ONLY: PI
USE MODgvec_LinAlg ,   ONLY: INV,SOLVEMAT,getLU
USE MODgvec_Basis1D,   ONLY:  LegendreGaussNodesAndWeights
USE MODgvec_Basis1D,   ONLY:  BarycentricWeights,InitializeVandermonde,MthPolynomialDerivativeMatrix
USE sll_m_bsplines,ONLY: sll_s_bsplines_new
USE sll_m_spline_interpolator_1d         ,ONLY: sll_t_spline_interpolator_1d
USE sll_m_boundary_condition_descriptors ,ONLY: sll_p_greville
USE sll_m_spline_matrix                  ,ONLY: sll_s_spline_matrix_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   )        :: deg_in        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity_in !! continuity:
                                                        !! 0: disc. polynomial
                                                        !! deg-1: spline with cont. deg-1
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in       !! grid information
  INTEGER       , INTENT(IN   )        :: degGP_in      !! gauss quadrature points: nGP=degGP+1 per elements
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sbase), INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i,m,iGP,iElem,imin,jmin,nElems,deg,degGP,nBase
  INTEGER  :: iBC,j,diri,odd_even,nD,BC_derivs(deg_in+1),drv
  REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: locbasis,VdmGP
  REAL(wp),DIMENSION(1:deg_in+1,1:deg_in+1) :: PLmat,Umat
!===================================================================================================================================
  IF(.NOT.test_called) THEN
    SWRITE(UNIT_stdOut,'(4X,A,3(A,I3),A)')'INIT sBase type:', &
         ' degree= ',deg_in, &
         ' gauss points per elem = ',degGP_in, &
         ' continuity= ',continuity_in, ' ...'
  END IF
  IF(sf%initialized) THEN
    CALL abort(__STAMP__, &
        "Trying to reinit sbase!")
  END IF
  IF(degGP_in.LT.deg_in) &
    CALL abort(__STAMP__, &
        "error in sbase: degGP must be > deg!",&
        TypeInfo="InvalidParameterError")
  SELECT TYPE(sf)
  TYPE IS(t_sbase_disc)
    IF(continuity_in.NE.-1) &
      CALL abort(__STAMP__, &
          "error in sbase init: type is disc but continuity is not -1, maybe sbase_new was not called before!")
  TYPE IS(t_sbase_spl)
    IF(continuity_in.NE.deg_in-1) &
      CALL abort(__STAMP__, &
          "error in sbase init: type is spl but continuity is not deg-1, mabye sbase_new was not called before!")
  CLASS DEFAULT
      CALL abort(__STAMP__, &
          "error in sbase init: type is neither disc or spl!")
  END SELECT !Type
  sf%deg        =deg_in
  sf%continuity =continuity_in
  sf%grid       => grid_in
  sf%degGP      =  degGP_in

  nElems      = sf%grid%nElems
  deg         = sf%deg
  degGP       = sf%degGP

  sf%nGP  = (degGP+1)*nElems
  sf%nGP_str = (degGP+1)*(sf%grid%nElems_str-1)+1 !< for MPI
  sf%nGP_end = (degGP+1)*sf%grid%nElems_end

  IF(sf%continuity.EQ.-1)THEN !discontinuous
    sf%nBase  = (deg+1)*nElems
  ELSEIF(sf%continuity.EQ.deg-1)THEN !bspline with full continuity and interpolation base at boundaries
    sf%nBase  = nElems + deg
  ELSE
   CALL abort(__STAMP__, &
          'other spline continuities not yet implemented')
  END IF !continuity

  nBase = sf%nBase
  CALL sbase_alloc(sf)


  CALL LegendreGaussNodesAndWeights(degGP,sf%xi_GP,sf%w_GPloc) ![-1,1] !!!

  DO iElem=1,nElems
    sf%w_GP(1+(degGP+1)*(iElem-1):(degGP+1)*iElem)=                    0.5_wp *sf%w_GPloc(:)      *sf%grid%ds(iElem)
    sf%s_GP(1+(degGP+1)*(iElem-1):(degGP+1)*iElem)=sf%grid%sp(iElem-1)+0.5_wp*(sf%xi_GP(:)+1.0_wp)*sf%grid%ds(iElem)
  END DO !iElem

  SELECT TYPE(sf)
  TYPE IS(t_sbase_disc)
    IF(deg.EQ.0)THEN
      sf%xiIP(0)=0.0_wp
      CALL BarycentricWeights(deg,sf%xiIP,sf%wBaryIP)
      sf%DmatIP=0.0_wp
      DO iElem=1,nElems
        sf%base_offset(iElem)=1+(deg+1)*(iElem-1)
        sf%base_GP   (0:degGP,0,iElem)=1.0_wp
        sf%base_ds_GP(0:degGP,0,iElem)=0.0_wp
      END DO !iElem
      !zero deriv: evaluation of basis functions
      sf%base_dsAxis(0,1      )=1.0_wp
      sf%base_dsEdge(0,nBase  )=1.0_wp
    ELSE
      ALLOCATE(VdmGP( 0:degGP,0:deg))
      !  use chebychev-lobatto points for interpolation (closed form!), interval [-1,1]
      DO i=0,deg
        sf%xiIP(i)=-COS(REAL(i,wp)/REAL(deg,wp)*PI)
      END DO
      CALL BarycentricWeights(deg,sf%xiIP,sf%wBaryIP)
      CALL InitializeVandermonde(deg,degGP,sf%wBaryIP,sf%xiIP,sf%xi_GP,VdmGP)
      CALL MthPolynomialDerivativeMatrix(deg,sf%xiIP,1,sf%DmatIP)
      ! eval basis and  basis derivative
      DO iElem=1,nElems
        sf%base_offset(iElem)=1+(deg+1)*(iElem-1)
        sf%base_GP   (0:degGP,0:deg,iElem)=VdmGP(:,:)
        sf%base_ds_GP(0:degGP,0:deg,iElem)=MATMUL(VdmGP,sf%DmatIP)*(2.0_wp/sf%grid%ds(iElem))
      END DO !iElem
      !zero deriv: evaluation of basis functions (lagrange property!)
      sf%base_dsAxis(0,1      )=1.0_wp
      sf%base_dsAxis(0,2:deg+1)=0.0_wp
      sf%base_dsEdge(0,nBase-deg:nBase-1)=0.0_wp
      sf%base_dsEdge(0,          nBase  )=1.0_wp
      ! eval basis deriv at boundaries d/ds = d/dxi dxi/ds = 1/(0.5ds) d/dxi
      sf%base_dsAxis(1,1:deg+1        )=sf%DmatIP(  0,:)*(2.0_wp/sf%grid%ds(1))
      sf%base_dsEdge(1,nBase-deg:nBase)=sf%DmatIP(deg,:)*(2.0_wp/sf%grid%ds(nElems))
      !  higher derivatives
      DO i=2,deg
        sf%base_dsAxis(i,1:deg+1        )=MATMUL(TRANSPOSE(sf%DmatIP),sf%base_dsAxis(i-1,:))*(2.0_wp/sf%grid%ds(1))
        sf%base_dsEdge(i,nBase-deg:nBase)=MATMUL(TRANSPOSE(sf%DmatIP),sf%base_dsEdge(i-1,:))*(2.0_wp/sf%grid%ds(nElems))
      END DO
      DEALLOCATE(VdmGP)
    END IF !deg=0
    !interpolation:
    !  points are repeated at element interfaces (discontinuous)
    ALLOCATE(sf%s_IP(nBase)) !for spl, its allocated elsewhere...
    DO iElem=1,nElems
      sf%s_IP(1+(deg+1)*(iElem-1):(deg+1)*iElem)=sf%grid%sp(iElem-1)+0.5_wp*(sf%xiIP+1.0_wp)*sf%grid%ds(iElem)
    END DO !iElem
    sf%s_IP(1)=0.0_wp
    sf%s_IP(nBase)=1.0_wp
  TYPE IS(t_sbase_spl)
    ALLOCATE(locbasis(0:deg,0:deg))
    CALL sll_s_bsplines_new(sf%bspl ,degree=deg,periodic=.FALSE., &
                            xmin=sf%grid%sp(0),xmax=sf%grid%sp(nElems),&
                            ncells=nElems,breaks=sf%grid%sp)
    !basis evaluation

    IF(sf%bspl%nBasis.NE.nBase) CALL abort(__STAMP__, &
                                           'problem with bspl basis')
    DO iElem=1,nElems
      j=1+(degGP+1)*(iElem-1)
      CALL sf%bspl % eval_basis(sf%s_GP(j),sf%base_GP(0,0:deg,iElem),imin)
      IF(imin.EQ.-1) CALL abort(__STAMP__,&
                                'problem, element in eval_basis not found!')
      DO iGP=1,degGP
        CALL sf%bspl % eval_basis(sf%s_GP(j+iGP),sf%base_GP(iGP,0:deg,iElem),jmin)
        IF(jmin.EQ.-1) CALL abort(__STAMP__,&
                                 'problem, element in eval_basis not found!')
        IF(jmin.NE.imin) CALL abort(__STAMP__,&
                                    'problem, GP are not in one element!')
      END DO !iGP=0,degGP
      CALL sf%bspl % eval_deriv(sf%s_GP(j),sf%base_ds_GP(0,0:deg,iElem),imin)
      IF(imin.EQ.-1) CALL abort(__STAMP__,&
                                'problem, element in eval_basis not found!')
      DO iGP=1,degGP
        CALL sf%bspl % eval_deriv(sf%s_GP(j+iGP),sf%base_ds_GP(iGP,0:deg,iElem),jmin)
        IF(jmin.EQ.-1) CALL abort(__STAMP__,&
                                  'problem, element in eval_basis not found!')
        IF(jmin.NE.imin) CALL abort(__STAMP__,&
                                    'problem, GP are not in one element!')
      END DO !iGP=0,degGP
      sf%base_offset(iElem)=imin
    END DO !iElem=1,nElems
    !eval all basis derivatives at boundaries
    CALL sf%bspl % eval_basis_and_n_derivs(sf%grid%sp(     0),deg,locBasis,imin) !locBasis(0:nderiv,0:deg Base)
    IF(imin.NE.1) CALL abort(__STAMP__,&
                             'problem eval_deriv left')
    sf%base_dsAxis(0:deg,1:deg+1) =locbasis(:,:) ! basis functions 1 ...deg+1

    CALL sf%bspl % eval_basis_and_n_derivs(sf%grid%sp(nElems),deg,locbasis,imin)
    IF(imin.NE.nBase-deg) CALL abort(__STAMP__,&
                                     'problem eval_deriv right')
    sf%base_dsEdge(0:deg,nBase-deg:nBase)=locbasis(:,:) ! basis functions nBase-deg ... nbase

    !interpolation
    CALL sf%Interpol%init (sf%bspl,sll_p_greville,sll_p_greville)
    CALL sf%Interpol%get_interp_points ( sf%s_IP )

    DEALLOCATE(locbasis)
  END SELECT !TYPE is t_sbase_disc/spl


  !mass matrix
  CALL sll_s_spline_matrix_new(sf%mass , "banded",nBase,sf%deg,sf%deg)
  DO iElem=1,nElems
    jmin=sf%base_offset(iElem)
    DO i=0,deg
      DO j=0,deg
        CALL sf%mass%add_element(jmin+i,jmin+j, &
             (SUM(sf%w_GP(1+(degGP+1)*(iElem-1):(degGP+1)*iElem)*sf%base_GP(:,i,iElem)*sf%base_GP(:,j,iElem))))
      END DO !j=0,deg
    END DO !i=0,deg
  END DO !iElem=1,nElems
  CALL sf%mass % factorize()


  !STRONG BOUNDARY CONDITIONS: A and R matrices

  DO iBC=1,NBC_TYPES
    SELECT CASE(iBC)      !nDOF involved:    Dirichlet?  odd(0),even(1)
    CASE(BC_TYPE_OPEN)     ; nD =0                                   ! do nothing
    CASE(BC_TYPE_NEUMANN)  ; nD =1          ;   diri=0 ;  odd_even=0 ! first derivative=0
    CASE(BC_TYPE_DIRICHLET); nD =1          ;   diri=1 ;  odd_even=0 ! odd_even not used, first DOF=BC_val
    CASE(BC_TYPE_SYMM)     ; nD =(deg+1)/2  ;   diri=0 ;  odd_even=0 ! open,    all derivatives (2*k-1)=0, k=1,...(deg+1)/2
    CASE(BC_TYPE_SYMMZERO) ; nD =1+(deg+1)/2;   diri=1 ;  odd_even=0 ! dirichlet=0+ derivatives (2*k-1)=0, k=1,...(deg+1)/2
    CASE(BC_TYPE_ANTISYMM) ; nD =1+deg/2    ;   diri=1 ;  odd_even=1 ! dirichlet=0+ derivatives (2*k  )=0  k=1,... deg/2
    END SELECT !iBC
    sf%nDOF_BC(iBC)=nD
    !A and R are already initialized as unit matrices!!
    IF((nD.GT.0).AND.(deg.GT.0))THEN
      IF(diri.EQ.1) BC_derivs(1)=0  !A(1,:) is already unit matrix
      DO i=diri+1,nD
         BC_derivs(i)=2*(i-diri)-(1-odd_even)
      END DO
      DO i=diri+1,nD
        drv=BC_derivs(i)
        sf%A_Axis(        i,:,iBC)=sf%base_dsAxis(drv,:)/sf%base_dsAxis(drv,i) !normalized with diagonal entry
        sf%A_Edge(nBase+1-i,:,iBC)=sf%base_dsEdge(drv,:)/sf%base_dsEdge(drv,nBase+1-i)
      END DO

      !PRECONDITIONING OF AMAT, FORWARD SWEEP:

      !lower left triangular part to zero, and normalize diagonal
      !! DO j=1,nD
      !!   DO i=j+1,nD
      !!     sf%A_axis(i,:,iBC)=sf%A_axis(i,:,iBC)-sf%A_axis(j,:,iBC)*(sf%A_axis(i,j,iBC)/sf%A_axis(j,j,iBC))
      !!   END DO
      !! END DO
      !! DO i=diri+1,nD
      !!   sf%A_axis(i,:,iBC)=sf%A_axis(i,:,iBC)/sf%A_axis(i,i,iBC)
      !! END DO

      ! A*x=0 , [A_s A_r]*[x1 x2] = 0  A_s*x1=-A_r*x2, U^T (PL)^T x1 = -A_r*x2 -> NEW A_s= (PL)^T, A_r = U^{-T} A_r
      CALL getLU(nD,TRANSPOSE(sf%A_axis(1:nD,1:nD,iBC)),PLmat(1:nD,1:nD),Umat(1:nD,1:nD))
      sf%A_axis(1:nD,1:nD      ,iBC)=TRANSPOSE(PLmat(1:nD,1:nD))
      sf%A_axis(1:nD,nD+1:deg+1,iBC)=SOLVEMAT(TRANSPOSE(Umat(1:nD,1:nD)),sf%A_axis(1:nD,nD+1:deg+1,iBC))

      !upper right triangular part to zero, and normalize diagonal
      DO j=1,nD
        DO i=j+1,nD
          sf%A_Edge(nBase+1-i,:,iBC)= sf%A_Edge(nBase+1-i,:,iBC)  &
                                     -sf%A_Edge(nBase+1-j,:,iBC)*(sf%A_Edge(nBase+1-i,nBase+1-j,iBC)/sf%A_Edge(nBase+1-j,nBase+1-j,iBC))
        END DO
      END DO
      DO i=diri+1,nD
        sf%A_Edge(nBase+1-i,:,iBC)=sf%A_Edge(nBase+1-i,:,iBC)/sf%A_Edge(nBase+1-i,nBase+1-i,iBC)
      END DO

      !! NOT CLEAR HOW TO PRECONDITION A_EDGE WITH LU...



      !invert BC part (A_s is ( nD x nD ) and invertible):
      !     |  A_s   A_r |        |  0                    0 |
      !  A= |            |  , R=  |                         |
      !     |   0     I  |        | -(A_s^{-1}A_r)^T     I  |
      !
      sf%R_axis(   1:nD   ,:   ,iBC)=0.0_wp
      sf%R_axis(nD+1:deg+1,1:nD,iBC)=-TRANSPOSE(SOLVEMAT(sf%A_axis(1:nD,1:nD,iBC),sf%A_axis(1:nD,nD+1:deg+1,iBC)))

      sf%AR_axis(      1:nD,:,iBC)=sf%A_axis(      1:nD,:,iBC)
      sf%AR_axis(nD+1:deg+1,:,iBC)=sf%R_axis(nD+1:deg+1,:,iBC)
      !     |   I     0  |        |  I    -(B_s^{-1}B_r)^T  |
      !  B= |            |  , R=  |                         |
      !     |  B_r   B_s |        |  0                    0 |
      !
      sf%R_edge(nBase-deg:nBase-nD,nBase+1-nD:nBase,iBC)=-TRANSPOSE(SOLVEMAT(sf%A_edge(nBase+1-nD:nBase,nBase+1-nD:nBase  ,iBC), &
                                                                             sf%A_edge(nBase+1-nD:nBase,nBase-deg:nBase-nD,iBC)))
      sf%R_edge(nBase+1-nD:nBase,:,iBC)=0.0_wp


      sf%AR_edge(nBase-deg:nBase-nD,:,iBC)= sf%R_edge(nBase-deg:nBase-nD,:,iBC)
      sf%AR_edge(nBase+1-nD:nBase  ,:,iBC)= sf%A_edge(nBase+1-nD:nBase  ,:,iBC)


      !! sf%R_Axis(1:nD,1:nD,iBC)=INV(sf%A_Axis(1:nD,1:nD,iBC))
      !! sf%R_Axis(:,:,iBC)=TRANSPOSE(MATMUL(sf%R_Axis(:,:,iBC),sf%A_Axis(:,:,iBC)))
      !! DO i=1,deg+1; DO j=1,i-1
      !!   sf%R_Axis(i,j,iBC)=-sf%R_Axis(i,j,iBC)
      !! END DO; END DO
      !! sf%R_Edge(nBase-nD+1:nBase,nBase-nD+1:nBase,iBC)=INV(sf%A_Edge(nBase-nD+1:nBase,nBase-nD+1:nBase,iBC))
      !! sf%R_Edge(:,:,iBC)=TRANSPOSE(MATMUL(sf%R_Edge(:,:,iBC),sf%A_Edge(:,:,iBC)))
      !! DO i=nBase-deg,nBase; DO j=i+1,nBase
      !!   sf%R_Edge(i,j,iBC)=-sf%R_Edge(i,j,iBC)
      !! END DO; END DO
      !prepare for applyBC
      sf%invA_axis(:,:,iBC)=INV(sf%A_axis(:,:,iBC))
      sf%invA_edge(:,:,iBC)=INV(sf%A_edge(:,:,iBC))
      !automatically set rows 1:nD to zero for R matrices (no contribution from these DOF)
      !! sf%R_axis(         1:nD   ,:,iBC)=0.0_wp
      !! sf%R_edge(nBase-nD+1:nBase,:,iBC)=0.0_wp
    END IF
  END DO!iBC=1,NBC_TYPES

  ! mode number dependent boundary condition at axis, m  <=> -iBC
  ! iBC= 0, m=0:  BC_TYPE_SYMM
  ! iBC=-1, m=1:  BC_TYPE_ANTISYMM
  ! iBC=-2, m=2:  BC_TYPE_SYMMZERO
  ! ... odd m<=deg: ANTISYMM + all derivatives 1,...,m-1 =0
  ! ...even m<=deg: SYMMZERO + all derivatives 1,...,m-1 =0
  !m=0
  sf%nDOF_BC(   0) = sf%nDOF_BC(   BC_TYPE_SYMM)
  sf%A_Axis(:,:,0) = sf%A_Axis(:,:,BC_TYPE_SYMM)

  !  odd m <=deg
  DO m=1,deg+1,2
    iBC=-m
    !BC_TYPE_ANTISYMM: nD =1+deg/2    ;   diri=1 ;  odd_even=1 ! dirichlet=0+ derivatives (2*k  )=0  k=1,... deg/2
    !and also derivatives (2*k-1)=0, k=1,...,(m-1)/2  (2*k-1 < m-1)
    nD=1+deg/2 +(m-1)/2
    sf%nDOF_BC(iBC)=nD
    IF(nD.EQ.deg+1) CYCLE ! All deg+1 DOFs at the boundary just set to zero (A_axis is already identity matrix), if deg even and m=deg+1
    i=1
    BC_derivs(1)=0 !diri
    DO drv=1,(m-1)
      i=i+1
      BC_derivs(i)=drv
    END DO
    DO drv=m+1,deg,2
      i=i+1
      BC_derivs(i)=drv
    END DO
    IF(i.NE.nD) CALL abort(__STAMP__,&
                           "DEBUG odd m, i NE nD")
    DO i=2,nD
      drv=BC_derivs(i)
      sf%A_Axis(i,:,iBC)=sf%base_dsAxis(drv,:)/sf%base_dsAxis(drv,i) !normalized with diagonal entry
    END DO
  END DO
  !  even m <=deg
  DO m=2,deg+1,2
    iBC=-m
    !BC_TYPE_SYMMZERO: nD =1+(deg+1)/2;   diri=1 ;  odd_even=0 ! dirichlet=0+ derivatives (2*k-1)=0, k=1,...(deg+1)/2
    !and also derivatives (2*k)=0, k=1,...,(m-2)/2  (2*k < m-1)
    nD=1+(deg+1)/2 +(m-2)/2
    sf%nDOF_BC(iBC)=nD
    IF(nD.EQ.deg+1) CYCLE ! All deg+1 DOFs at the boundary just set to zero (A_axis is already identity matrix), if deg odd and m=deg+1
    i=1
    BC_derivs(1)=0 !diri
    DO drv=1,(m-1)
      i=i+1
      BC_derivs(i)=drv
    END DO
    DO drv=m+1,deg,2
      i=i+1
      BC_derivs(i)=drv
    END DO
    IF(i.NE.nD) CALL abort(__STAMP__,&
                           "DEBUG even m,i NE nD")
    DO i=2,nD
      drv=BC_derivs(i)
      sf%A_Axis(i,:,iBC)=sf%base_dsAxis(drv,:)/sf%base_dsAxis(drv,i) !normalized with diagonal entry
    END DO
  END DO

  DO m=0,deg+1
    iBC=-m
    nD=sf%nDOF_BC(iBC)
    IF(nD.LT.deg+1)THEN
      !precondition A via LU decomp. of A_s^T
      ! A*x=0 , [A_s A_r]*[x1 x2] = 0  A_s*x1=-A_r*x2, U^T (PL)^T x1 = -A_r*x2 -> NEW A_s= (PL)^T, A_r = U^{-T} A_r
      CALL getLU(nD,TRANSPOSE(sf%A_axis(1:nD,1:nD,iBC)),PLmat(1:nD,1:nD),Umat(1:nD,1:nD))
      sf%A_axis(1:nD,1:nD      ,iBC)=TRANSPOSE(PLmat(1:nD,1:nD))
      sf%A_axis(1:nD,nD+1:deg+1,iBC)=SOLVEMAT(TRANSPOSE(Umat(1:nD,1:nD)),sf%A_axis(1:nD,nD+1:deg+1,iBC))

      !invert BC part (A_s is ( nD x nD ) and invertible):
      !     |  A_s   A_r |        |  0                    0 |
      !  A= |            |  , R=  |                         |
      !     |   0     I  |        | -(A_s^{-1}A_r)^T     I  |
      !
      sf%R_axis(   1:nD   ,:   ,iBC)=0.0_wp
      sf%R_axis(nD+1:deg+1,1:nD,iBC)=-TRANSPOSE(SOLVEMAT(sf%A_axis(1:nD,1:nD,iBC),sf%A_axis(1:nD,nD+1:deg+1,iBC)))

      sf%AR_axis(      1:nD,:,iBC)=sf%A_axis(      1:nD,:,iBC)
      sf%AR_axis(nD+1:deg+1,:,iBC)=sf%R_axis(nD+1:deg+1,:,iBC)
      !invert BC part
      !!sf%R_Axis(1:nD,1:nD,iBC)=INV(sf%A_Axis(1:nD,1:nD,iBC))
      !!sf%R_Axis( :  , :  ,iBC)=TRANSPOSE(MATMUL(sf%R_Axis(:,:,iBC),sf%A_Axis(:,:,iBC)))
      !!DO i=1,deg+1; DO j=1,i-1
      !!  sf%R_Axis(i,j,iBC)=-sf%R_Axis(i,j,iBC)
      !!END DO; END DO
      !prepare for applyBC
      sf%invA_axis(:,:,iBC)=INV(sf%A_axis(:,:,iBC))
      !automatically set rows 1:nD to zero for R matrices (no contribution from these DOF)
      !!sf%R_axis(1:nD,:,iBC)=0.0_wp
    ELSE
      !nD=deg+1, A_axis is simply the identity!
      sf%R_axis(   :,:,iBC)=0.0_wp
      sf%AR_axis(  :,:,iBC)=sf%A_axis(:,:,iBC)
      sf%invA_axis(:,:,iBC)=sf%A_axis(:,:,iBC)
    END IF
  END DO!m=0,deg+1


  sf%initialized=.TRUE.
  IF(.NOT.test_called) THEN
    SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
    CALL sBase_test(sf)
  END IF

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
  ALLOCATE(sf%xi_GP(     0:degGP))
  ALLOCATE(sf%w_GPloc(   0:degGP))
  ALLOCATE(sf%w_GP((degGP+1)*nElems))
  ALLOCATE(sf%s_GP((degGP+1)*nElems))
  ALLOCATE(sf%base_GP(   0:degGP,0:deg,1:nElems))
  ALLOCATE(sf%base_ds_GP(0:degGP,0:deg,1:nElems))
  ALLOCATE(sf%base_offset(            1:nElems))
  ALLOCATE(sf%base_dsAxis(0:deg,1:deg+1        ))
  ALLOCATE(sf%base_dsEdge(0:deg,nBase-deg:nBase))
  ALLOCATE(sf%nDOF_BC(                  -(deg+1):NBC_TYPES))
  ALLOCATE(sf%A_Axis(   1:deg+1,1:deg+1,-(deg+1):NBC_TYPES))
  ALLOCATE(sf%invA_Axis(1:deg+1,1:deg+1,-(deg+1):NBC_TYPES))
  ALLOCATE(sf%R_Axis(   1:deg+1,1:deg+1,-(deg+1):NBC_TYPES))
  ALLOCATE(sf%AR_Axis(  1:deg+1,1:deg+1,-(deg+1):NBC_TYPES))
  ALLOCATE(sf%A_Edge(   nBase-deg:nBase,nBase-deg:nBase,1:NBC_TYPES))
  ALLOCATE(sf%invA_Edge(nBase-deg:nBase,nBase-deg:nBase,1:NBC_TYPES))
  ALLOCATE(sf%R_Edge(   nBase-deg:nBase,nBase-deg:nBase,1:NBC_TYPES))
  ALLOCATE(sf%AR_Edge(  nBase-deg:nBase,nBase-deg:nBase,1:NBC_TYPES))
  sf%xi_GP        =0.0_wp
  sf%w_GPloc      =0.0_wp
  sf%w_GP         =0.0_wp
  sf%s_GP        =0.0_wp
  sf%base_offset =-1
  sf%base_GP      =0.0_wp
  sf%base_ds_GP   =0.0_wp
  sf%base_dsAxis =0.0_wp
  sf%base_dsEdge =0.0_wp
  sf%nDOF_BC     =0
  sf%A_Axis      =0.0_wp
  sf%invA_Axis   =0.0_wp
  sf%A_Edge      =0.0_wp
  sf%invA_Edge   =0.0_wp
  sf%R_Axis      =0.0_wp
  sf%R_Edge      =0.0_wp
  DO i=0,deg
    sf%A_Axis(   1+i,1+i,:)=1.0_wp
    sf%invA_Axis(1+i,1+i,:)=1.0_wp
    sf%R_Axis(   1+i,1+i,:)=1.0_wp
    sf%A_Edge(   nBase-deg+i,nBase-deg+i,:)=1.0_wp
    sf%invA_Edge(nBase-deg+i,nBase-deg+i,:)=1.0_wp
    sf%R_Edge(   nBase-deg+i,nBase-deg+i,:)=1.0_wp
  END DO
  SELECT TYPE(sf)
  TYPE IS(t_sbase_disc)
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
  SDEALLOCATE(sf%xi_GP)
  SDEALLOCATE(sf%w_GPloc)
  SDEALLOCATE(sf%w_GP)
  SDEALLOCATE(sf%s_GP)
  SDEALLOCATE(sf%s_IP)
  SDEALLOCATE(sf%base_offset)
  SDEALLOCATE(sf%base_GP)
  SDEALLOCATE(sf%base_ds_GP)
  SDEALLOCATE(sf%base_dsAxis)
  SDEALLOCATE(sf%base_dsEdge)
  SDEALLOCATE(sf%nDOF_BC)
  SDEALLOCATE(sf%A_Axis)
  SDEALLOCATE(sf%invA_Axis)
  SDEALLOCATE(sf%R_Axis)
  SDEALLOCATE(sf%AR_Axis)
  SDEALLOCATE(sf%A_Edge)
  SDEALLOCATE(sf%invA_Edge)
  SDEALLOCATE(sf%R_Edge)
  SDEALLOCATE(sf%AR_Edge)
  SELECT TYPE (sf)
  TYPE IS(t_sbase_spl)
    CALL sf%interpol%free()
    CALL sf%bspl%free()
  TYPE IS(t_sbase_disc)
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
!> copy  onto sf <-- tocopy
!!
!===================================================================================================================================
SUBROUTINE sBase_copy( sf , tocopy)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_sBase), INTENT(IN   ) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); CLASS IS(t_sbase)
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
  CALL sf%init(tocopy%deg,tocopy%continuity,tocopy%grid,tocopy%degGP)

  END SELECT !TYPE
END SUBROUTINE sbase_copy


!===================================================================================================================================
!> compare sf with input sbase
!!
!===================================================================================================================================
SUBROUTINE sBase_compare( sf , tocompare,is_same,cond_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sBase),  INTENT(IN   ) :: sf !! self
  CLASS(c_sBase),  INTENT(IN   ) :: tocompare
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  LOGICAL,OPTIONAL,INTENT(  OUT) :: is_same
  LOGICAL,OPTIONAL,INTENT(  OUT) :: cond_out(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL  :: cond(5)
!===================================================================================================================================
  SELECT TYPE(tocompare); CLASS IS(t_sbase)
  IF(.NOT.tocompare%initialized) THEN
    CALL abort(__STAMP__, &
        "sBase_compare: tried to compare to non-initialized sBase!")
  END IF

  CALL sf%grid%compare(tocompare%grid,cond(1))
  IF(cond(1)) THEN
    !same grid, compare basis
    cond(2)= (sf%nbase      .EQ. tocompare%nbase     )
    cond(3)= (sf%deg        .EQ. tocompare%deg       )
    cond(4)= (sf%continuity .EQ. tocompare%continuity)
  ELSE
    cond(2:4)=.FALSE.
  END IF

  cond(5)=.FALSE.
  SELECT TYPE(tocompare)
  TYPE IS(t_sbase_disc)
    SELECT TYPE(sf)
    TYPE IS(t_sbase_disc)
      cond(5)=.TRUE.
    END SELECT
  TYPE IS(t_sbase_spl)
    SELECT TYPE(sf)
    TYPE IS(t_sbase_spl)
     cond(5)=.TRUE.
    END SELECT
  END SELECT !TYPE(tocompare)

  IF(PRESENT(is_same)) is_same=ALL(cond)
  !IF(.NOT.ALL(cond)) WRITE(*,*)'DEBUG, not all cond in sbase for compare', cond

  IF(PRESENT(cond_out)) cond_out(1:5)=cond

  END SELECT !TYPE(tocompare)

END SUBROUTINE sbase_compare


!===================================================================================================================================
!> change data from old_sBase to self.
!! using interpolations of the old data at the new interpolation points
!!
!===================================================================================================================================
SUBROUTINE sBase_change_base( sf,old_sBase,iterDim,old_data,sf_data)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sBase),  INTENT(IN   ) :: sf !! self
  CLASS(c_sBase),  INTENT(IN   ) :: old_sBase       !! base of old_data
  INTEGER         ,INTENT(IN   ) :: iterDim        !! iterate on first or second dimension or old_data/sf_data
  REAL(wp)        ,INTENT(IN   ) :: old_data(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)        ,INTENT(  OUT) :: sf_data(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL              :: sameBase
  INTEGER              :: iIP,iter,iterSize
  REAL(wp),ALLOCATABLE :: eval_old_sbase(:,:),gIP(:,:)
  INTEGER ,ALLOCATABLE :: Elem_old_sbase(:)
!===================================================================================================================================
  SELECT TYPE(old_sBase); CLASS IS(t_sBase)
  IF(.NOT.old_sBase%initialized) THEN
    CALL abort(__STAMP__, &
        "sBase_chamge_base: tried to change base with non-initialized sBase!")
  END IF
  IF((iterDim.LT.1).OR.(iterDim.GT.2))THEN
    CALL abort(__STAMP__, &
        "sBase_chamge_base: iterDim can only be 1 or 2!")
  END IF
  iterSize=SIZE(old_data,iterDim)
  IF(iterSize.NE.SIZE(sf_data,iterDim)) THEN
    CALL abort(__STAMP__, &
        "sBase_chamge_base: iteration dimenion of old_data and sf_data have to be the same!")
  END IF
  IF(SIZE(old_data,3-iterDim).NE.old_sBase%nbase) THEN
    CALL abort(__STAMP__, &
        "sBase_chamge_base: old_data size does not match old_sBase!")
  END IF
  IF(SIZE( sf_data,3-iterDim).NE.      sf%nbase) THEN
    CALL abort(__STAMP__, &
        "sBase_chamge_base: sf_data size does not match sf sBase!")
  END IF

  CALL sf%compare(old_sBase,is_same=sameBase)

  IF(sameBase)THEN
   !same base
   sf_data=old_data
  ELSE
    !actually change base
    ALLOCATE(Elem_old_sbase(1:sf%nBase),eval_old_sbase(0:old_sbase%deg,1:sf%nBase))
    ALLOCATE(gIP(1:sf%nBase,iterSize))
    DO iIP=1,sf%nBase
      CALL old_sbase%eval(sf%S_IP(iIP),0,Elem_old_sbase(iIP),eval_old_sbase(:,iIP))
    END DO
    SELECT CASE(iterDim)
    CASE(1)
      DO iIP=1,sf%nBase
        DO iter=1,iterSize
          gIP(iIP,iter)=old_sbase%evalDOF_base(Elem_old_sbase(iIP),eval_old_sbase(:,iIP),old_data(iter,:))
        END DO
      END DO
      DO iter=1,iterSize
        sf_data(iter,:)=sf%initDOF(gIP(:,iter))
      END DO
    CASE(2)
      DO iIP=1,sf%nBase
        DO iter=1,iterSize
          gIP(iIP,iter)=old_sbase%evalDOF_base(Elem_old_sbase(iIP),eval_old_sbase(:,iIP),old_data(:,iter))
        END DO
      END DO
      DO iter=1,iterSize
        sf_data(:,iter)=sf%initDOF(gIP(:,iter))
      END DO
    END SELECT !iterDim

    DEALLOCATE(Elem_old_sBase,eval_old_sbase,gIP)
  END IF !same base

  END SELECT !TYPE
END SUBROUTINE sBase_change_base


!===================================================================================================================================
!> evaluate sbase at position x [0,1], NOT EFFICIENT!!
!!
!===================================================================================================================================
SUBROUTINE sBase_eval( sf , x,deriv,iElem,base_x)
! MODULES
USE MODgvec_Basis1D, ONLY:LagrangeInterpolationPolys
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
TYPE IS(t_sbase_disc)
  iElem=sf%grid%find_elem(x)
  xiloc  =(x-sf%grid%sp(iElem-1))*2.0_wp/sf%grid%ds(iElem)-1.0_wp !in [-1,1]

  IF(deriv.EQ.0)THEN
    CALL LagrangeInterpolationPolys(xiloc,sf%deg,sf%xiIP,sf%wBaryIP,base_x(:))
  ELSEIF(deriv.GT.sf%deg) THEN
    base_x=0.0_wp
  ELSE
    CALL LagrangeInterpolationPolys(xiloc,sf%deg,sf%xiIP,sf%wBaryIP,baseloc(:,0))
    DO ideriv=1,deriv
      baseloc(:,ideriv)=MATMUL(TRANSPOSE(sf%DmatIP),baseloc(:,ideriv-1))*(2.0_wp/sf%grid%ds(iElem))
    END DO
    base_x=baseloc(:,deriv)
  END IF!deriv
TYPE IS(t_sbase_spl)
  IF(deriv.EQ.0)THEN
    CALL sf%bspl%eval_basis(x,base_x(:),iElem)
  ELSEIF(deriv.GT.sf%deg) THEN
    iElem=sf%grid%find_elem(x)
    base_x=0.0_wp
  ELSE
    CALL sf%bspl%eval_basis_and_n_derivs(x,deriv,baseloc(0:deriv,:),iElem)
    base_x=baseloc(deriv,:)
  END IF

  IF(iElem.EQ.-1)CALL abort(__STAMP__,&
                            'PROBLEM, iElem not found in spline eval (sbase_eval)...')
CLASS DEFAULT
  CALL abort(__STAMP__, &
    "this type of continuity not implemented!")
END SELECT !TYPE

END SUBROUTINE sbase_eval

!===================================================================================================================================
!> simply evaluate function or derivative at point x
!!
!===================================================================================================================================
FUNCTION sBase_evalDOF_s(sf,x,deriv,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf     !! self
  REAL(wp)      , INTENT(IN   ) :: x      !! point positions in [0,1]
  INTEGER       , INTENT(IN   ) :: deriv  !! derivative (=0: solution)
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all degrees of freedom
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iElem
  REAL(wp)                      :: base_x(0:sf%deg)
!===================================================================================================================================
IF(SIZE(DOFs,1).NE.sf%nBase) CALL abort(__STAMP__, &
             'nDOF not correct when calling sBase_evalDOF_s')
  CALL sf%eval(x,deriv,iElem,base_x)
  y=sf%evalDOF_base(iElem,base_x,DOFs)

END FUNCTION sbase_evalDOF_s


!===================================================================================================================================
!> simply evaluate function with a base or base derivative evaluated at a point and its corresponding iElem
!! use together with sBase_eval(x) => iElem,base
!!
!===================================================================================================================================
FUNCTION sBase_evalDOF_base(sf ,iElem,base_x,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  INTEGER       , INTENT(IN   ) :: iElem  !! element where evaluation point
  REAL(wp)      , INTENT(IN   ) :: base_x(0:sf%deg)  !! evaluation of base or its derivative in element  iElem at a point position
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! degrees of freedom
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(SIZE(DOFs,1).NE.sf%nBase) CALL abort(__STAMP__, &
             'nDOF not correct when calling sBase_evalDOF_base')
  ASSOCIATE(j=>sf%base_offset(iElem))
  y = SUM(DOFs(j:j+sf%deg)*base_x)
  END ASSOCIATE

END FUNCTION sbase_evalDOF_base

!===================================================================================================================================
!> simply evaluate function or derivative at point x, for multiple DOF vectors
!!
!===================================================================================================================================
FUNCTION sBase_evalDOF2D_s(sf,x,nd,deriv,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf     !! self
  REAL(wp)      , INTENT(IN   ) :: x      !! point positions in [0,1]
  INTEGER       , INTENT(IN   ) :: nd     !! number of DOF vectors
  INTEGER       , INTENT(IN   ) :: deriv  !! derivative (=0: solution)
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%nBase,1:nd)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y(1:nd)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iElem,j
  REAL(wp)                      :: base_x(0:sf%deg)
!===================================================================================================================================
  __PERFON('eval_dof2d_s')
  CALL sf%eval(x,deriv,iElem,base_x)
  j=sf%base_offset(iElem)
  !y(1:nd) =MATMUL(base_x(0:sf%deg),DOFs(j:j+sf%deg,1:nd))
  __MATVEC_T(y,DOFs(j:j+sf%deg,:),base_x(0:sf%deg))

  __PERFOFF('eval_dof2d_s')
END FUNCTION sbase_evalDOF2D_s

!===================================================================================================================================
!> evaluate all degrees of freedom at all Gauss Points (deriv=0 solution, deriv=1 first derivative d/ds)
!!
!===================================================================================================================================
FUNCTION sBase_evalDOF_GP(sf,deriv,DOFs) RESULT(y_GP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: deriv  !! only 0 or 1
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all degrees of freedom
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y_GP(1:sf%nGP) ! will be be 1D array on input/output
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iElem,j,k,deg,degGP,nelems
!===================================================================================================================================
  IF(SIZE(DOFs,1).NE.sf%nBase) CALL abort(__STAMP__, &
               'nDOF not correct when calling sBase_evalDOF_GP')
  deg=sf%deg; degGP=sf%degGP; nElems=sf%grid%nElems
  SELECT CASE(deriv)
  CASE(0)
!    k=1
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) &
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iElem,j,k)  &
!$OMP   SHARED(sf,y_GP,DOFs,deg,degGP,nElems)
     DO iElem=1,nElems
      j=sf%base_offset(iElem)
      k=(iElem-1)*(degGP+1)+1
      y_GP(k:k+degGP)=MATMUL(sf%base_GP(   :,:,iElem),DOFs(j:j+deg))
      !!__MATVEC_N(y_GP(k:k+degGP),sf%base_GP(:,:,iElem),DOFs(j:j+deg))
    END DO
!$OMP END PARALLEL DO
  CASE(DERIV_S)
!    k=1
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) &
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iElem,j,k)  &
!$OMP   SHARED(sf,y_GP,DOFs,deg,degGP,nElems)
    DO iElem=1,nElems
      j=sf%base_offset(iElem)
      k=(iElem-1)*(degGP+1)+1
      y_GP(k:k+degGP)=MATMUL(sf%base_ds_GP(:,:,iElem),DOFs(j:j+deg))
      !!__MATVEC_N(y_GP(k:k+degGP),sf%base_ds_GP(:,:,iElem),DOFs(j:j+deg))
    END DO
!$OMP END PARALLEL DO
  CASE DEFAULT
    CALL abort(__STAMP__, &
       'called evalDOF_GP: deriv must be 0 or DERIV_S!' )
  END SELECT !deriv
END FUNCTION sbase_evalDOF_GP

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
  REAL(wp)      , INTENT(IN   ) :: g_IP(:)  !!  interpolation values at s_IP positions [0,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: DOFs(1:sf%nBase)  !! result of interpolation
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(SIZE(g_IP,1).NE.sf%nBase) CALL abort(__STAMP__, &
               'nDOF not correct when calling sBase_initDOF')
  SELECT TYPE(sf)
  TYPE IS(t_sbase_disc)
    DOFs(:)=g_IP
  TYPE IS(t_sbase_spl)
    CALL sf%interpol%compute_interpolant( DOFs(:), g_IP )
  CLASS DEFAULT
    CALL abort(__STAMP__, &
      "this type of continuity not implemented!")
  END SELECT !TYPE

END FUNCTION sbase_initDOF

!===================================================================================================================================
!> apply strong boundary conditions at axis and edge
!! Not used anymore, WAS FOUND TO BE NOT STABLE (OSCILLATORY) IN GENERAL, especially for BCs with only derivatives prescribed...
!! new LGM method below is used now!
!!
!===================================================================================================================================
SUBROUTINE sBase_applyBCtoDOF_STRONG(sf ,DOFs,BC_Type,BC_val)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  INTEGER       , INTENT(IN   ) :: BC_Type(2)           !! bc type on axis (1) and edge (2)
  REAL(wp)      , INTENT(IN   ) :: BC_Val(2)           !! for dirichlet BC : value
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)      , INTENT(INOUT) :: DOFs(1:sf%nBase)  !! DOFs with boundary conditions applied
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp):: raxis(1:sf%deg+1),redge(sf%nBase-sf%deg:sf%nbase)
!===================================================================================================================================
  ASSOCIATE( tBCaxis => BC_TYPE(BC_AXIS)             &
            ,tBCedge => BC_TYPE(BC_EDGE)             &
            ,nDaxis  => sf%nDOF_BC(BC_Type(BC_AXIS)) & !number of dofs involved in BC axis
            ,nDedge  => sf%nDOF_BC(BC_Type(BC_EDGE)) & !number of dofs involved in BC edge
            ,nB      => sf%nBase                     &
            ,deg     => sf%deg                       )
  SELECT CASE(tBCaxis)
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    DOFs(1)=BC_Val(BC_AXIS)
  CASE DEFAULT !BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM
    raxis(1:nDaxis)      =0.0_wp
    raxis(nDaxis+1:deg+1)= DOFs(nDaxis+1:deg+1)
    DOFs(1:deg+1)= MATMUL(sf%invA_axis(:,:,tBCaxis),raxis(:))
  END SELECT !tBCaxis

  SELECT CASE(tBCedge)
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    DOFs(nB)=BC_Val(BC_EDGE)
  CASE DEFAULT !BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM
    redge(nB-deg:nB-nDedge) = DOFs(nB-deg:nB-nDedge)
    redge(nB-nDedge+1:nB)   = 0.0_wp
    DOFs(nB-deg:nB)=MATMUL(sf%invA_edge(:,:,tBCedge),redge(:))
  END SELECT !tBCedge
  END ASSOCIATE

END SUBROUTINE sbase_applyBCtoDOF_STRONG

!===================================================================================================================================
!> apply boundary conditions at axis and edge, via solving the Lagrange multiplier problem:
!! x_new=x_old & A*x_new = d
!!
!!  | 0    A  | |lambda|     |  d   |
!!  |         | |      |  =  |      |
!!  | A^T  I  | | xnew |     | xold |
!!
!!
!===================================================================================================================================
SUBROUTINE sBase_applyBCtoDOF_LGM(sf ,DOFs,BC_Type,BC_val)
! MODULES
USE MODgvec_linalg, ONLY: SOLVE
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  INTEGER       , INTENT(IN   ) :: BC_Type(2)           !! bc type on axis (1) and edge (2)
  REAL(wp)      , INTENT(IN   ) :: BC_Val(2)           !! for dirichlet BC : value
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)      , INTENT(INOUT) :: DOFs(1:sf%nBase)  !! DOFs with boundary conditions applied
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp):: r(0:sf%deg)
!===================================================================================================================================
  ASSOCIATE( tBCaxis => BC_TYPE(BC_AXIS)             &
            ,tBCedge => BC_TYPE(BC_EDGE)             &
            ,nDaxis  => sf%nDOF_BC(BC_Type(BC_AXIS)) & !number of dofs involved in BC axis
            ,nDedge  => sf%nDOF_BC(BC_Type(BC_EDGE)) & !number of dofs involved in BC edge
            ,nB      => sf%nBase                     &
            ,deg     => sf%deg                       )
  SELECT CASE(tBCaxis)
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    DOFs(1)=BC_Val(BC_AXIS)
  CASE DEFAULT !BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM,m-dependent
    r(:) = DOFs(1:deg+1)
    r(:) = MATMUL(sf%R_axis(:,:,tBCaxis),r(:))
    DOFs(1:deg+1)= SOLVE(sf%AR_axis(:,:,tBCaxis),r(:))
  END SELECT !tBCaxis

  SELECT CASE(tBCedge)
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    DOFs(nB)=BC_Val(BC_EDGE)
  CASE DEFAULT !BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM
    r(:) = DOFs(nB-deg:nB)
    r(:) = MATMUL(sf%R_edge(:,:,tBCedge),r(:))
    DOFs(nB-deg:nB)= SOLVE(sf%AR_edge(:,:,tBCedge),r(:))
  END SELECT !tBCedge
  END ASSOCIATE

END SUBROUTINE sbase_applyBCtoDOF_LGM


!===================================================================================================================================
!> apply strong boundary conditions at axis and edge for solution update
!!
!===================================================================================================================================
SUBROUTINE sBase_applyBCtoRHS(sf ,RHS,BC_Type)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sbase), INTENT(IN   ) :: sf    !! self
  INTEGER       , INTENT(IN   ) :: BC_Type(2)           !! bc type on axis (1) and edge (2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)      , INTENT(INOUT) :: RHS(sf%nBase)  !! DOFs with boundary conditions applied
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp):: r(0:sf%deg)
!===================================================================================================================================
  ASSOCIATE( tBCaxis => BC_Type(BC_AXIS)             &
            ,tBCedge => BC_Type(BC_EDGE)             &
!            ,nDaxis  => sf%nDOF_BC(BC_Type(BC_AXIS)) & !number of dofs involved in BC axis
!            ,nDedge  => sf%nDOF_BC(BC_Type(BC_EDGE)) & !number of dofs involved in BC edge
            ,nB      => sf%nBase                     &
            ,deg     => sf%deg                       )
  SELECT CASE(tBCaxis)
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    RHS(1)= 0.0_wp
  CASE DEFAULT !BC_TYPE_NEUMANN,BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM
    r(:)         = RHS(1:deg+1)
    RHS(1:deg+1)= MATMUL(sf%R_axis(1:deg+1,1:deg+1,tBCaxis),r(:))
    !RHS(1:nDaxis)= 0.0_wp
  END SELECT !tBCaxis

  SELECT CASE(tBCedge)
  CASE(BC_TYPE_OPEN)
    !do noting
  CASE(BC_TYPE_DIRICHLET)
    RHS(nB)= 0.0_wp
  CASE DEFAULT !BC_TYPE_NEUMANN,BC_TYPE_SYMM,BC_TYPE_SYMMZERO,BC_TYPE_ANTISYMM
    r(:)                   = RHS(nB-deg:nB)
    RHS(nB-deg:nB)= MATMUL(sf%R_edge(nB-deg:nB,nB-deg:nB,tBCedge),r(:))
    !RHS(nB-nDedge:nB)= 0.0_wp
  END SELECT !tBCedge
  END ASSOCIATE

END SUBROUTINE sbase_applyBCtoRHS


!===================================================================================================================================
!> test sbase variable
!!
!===================================================================================================================================
SUBROUTINE sBase_test( sf)
! MODULES
USE MODgvec_GLobals, ONLY: testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: i,iTest,iElem,jElem,BC_Type(2),deg,degGP,nElems,cont,nBase
  REAL(wp)           :: x,y,y2,dy,dy2
  REAL(wp)           :: y_BC(0:sf%deg),y2_BC(0:sf%deg),base_x(0:sf%deg)
  REAL(wp)           :: g_IP(1:sf%nBase),dofs(1:sf%nBase)
  REAL(wp)           :: y_GP(1:sf%nGP)
  REAL(wp)           :: dof_GP(1:sf%nGP),BC_Val(2)
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
  CLASS(t_sbase),ALLOCATABLE :: testsbase
  LOGICAL            :: check(5)
  REAL(wp),ALLOCATABLE :: oldDOF(:,:),newDOF(:,:)
!===================================================================================================================================
  test_called=.TRUE.
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  deg=sf%deg
  degGP=sf%degGP
  cont=sf%continuity
  nBase=sf%nbase
  nElems=sf%grid%nElems
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN SBASE TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT.( ABS(sf%s_IP(1)).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,I4),(A,E11.3))') &
      '   degree = ',deg, '   continuity = ',cont, &
      '\n =>  should be 0.0 : sIP(1)= ',sf%s_IP(1)
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT.( ABS(sf%s_IP(nBase)-1.0_wp).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),(A,E11.3))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n =>  should be 1.0 : sIP(nBase)= ',sf%s_IP(nBase)
    END IF !TEST

    iTest=103 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    IF(testdbg.OR.(.NOT.( ABS(SUM(sf%w_GP)-1.0_wp).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),(A,E11.3))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be 1.0 : SUM(w_GP)= ',SUM(sf%w_GP)
    END IF !TEST

    iTest=104 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    CALL sf%eval(0.0_wp,0,iElem,base_x)
    IF(testdbg.OR.(.NOT.((ABS(base_x(0)-1.).LT.realtol).AND. &
                         ( SUM(ABS(base_x(1:deg))).LT. realtol).AND. &
                         (iElem.EQ.1) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,I6,*(E11.3))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be 1 and (1,0,0...) : base_x(x=0)= ',iElem,base_x
    END IF !TEST

    iTest=105 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    CALL sf%eval(1.0_wp,0,iElem,base_x)
    IF(testdbg.OR.(.NOT.((ABS(base_x(deg)-1.)         .LT. realtol ).AND. &
                         (MAXVAL(ABS(base_x(0:deg-1))).LT. realtol ).AND. &
                         (iElem                       .EQ. nElems     )      ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E11.3))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be ...,0,0,1 : base_x(x=1)= ',base_x
    END IF !TEST

    iTest=106 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    jElem=(nElems+1)/2
    x=0.5_wp*(sf%grid%sp(jElem-1)+sf%grid%sp(jElem))
    CALL sf%eval(x,0,iElem,base_x)
    IF(testdbg.OR.(.NOT.(iElem.EQ.jElem)))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,I6))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be ',jElem,': iElem= ' , iElem
    END IF !TEST

    iTest=107 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    jElem=MIN(2,nElems)
    x=sf%grid%sp(jElem-1)+0.01_wp*sf%grid%ds(jElem)
    CALL sf%eval(x ,0,iElem,base_x)
    IF(testdbg.OR.(.NOT.(iElem.EQ.jElem)))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,I6))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be ',jElem,': iElem= ' , iElem
    END IF !TEST

    !check compare
    iTest=111 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL sbase_new(testsbase,sf%deg,sf%continuity,sf%grid, sf%degGP)
    CALL testsbase%compare(sf,is_same=check(1))

    IF(.NOT.check(1))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),(A))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be true'
    END IF !TEST

    !check compare
    iTest=112 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testsbase%compare(sf, cond_out=check(1:5))

    IF(.NOT.ALL(check(:)))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),(A,5L))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be all true', check(:)
    END IF !TEST
    CALL testsbase%free()
    DEALLOCATE(testsbase)

    !check compare
    iTest=113 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL sbase_new(testsbase,sf%deg+1,MERGE(-1,sf%deg,(sf%continuity.EQ.-1)),sf%grid, sf%degGP)
    CALL testsbase%compare(sf,is_same=check(1))

    IF(check(1))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),(A))') &
      '   degree = ',deg, &
      '   continuity = ',cont,  ', nBase= ',nBase, &
      '\n => should be false'
    END IF !TEST
    CALL testsbase%free()
    DEALLOCATE(testsbase)

    !check change_base, execution (aborts if not correct)
    iTest=121 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL sbase_new(testsbase,sf%deg,sf%continuity,sf%grid, sf%degGP)
    ALLOCATE(oldDOF(1:sf%nBase,2),newDOF(1:testsBase%nBase,2))
    oldDOF(:,1)=1.0_wp
    oldDOF(:,2)=2.0_wp
    CALL testsbase%change_base(sf, 2, oldDOF,newDOF)
    y  =  MAXVAL(ABS(newDOF(:,1)-oldDOF(:,1)))
    y2 =  MINVAL(ABS(newDOF(:,2)-oldDOF(:,2)))

    IF(testdbg.OR.(.NOT.((y  .LT.realtol ).AND. &
                         (y2 .LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be zero ...',y  &
      ,'\n => should be zero ...',y2
    END IF !TEST
    DEALLOCATE(oldDOF,newDOF)
    CALL testsbase%free()
    DEALLOCATE(testsbase)

    !check change_base,only execution (aborts if not correct)
    iTest=122 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL sbase_new(testsbase,sf%deg+1,MERGE(-1,sf%deg,(sf%continuity.EQ.-1)),sf%grid, sf%degGP)
    ALLOCATE(oldDOF(2,1:sf%nBase),newDOF(2,1:testsBase%nBase))
    oldDOF(1,:)=-1.0_wp
    oldDOF(2,:)=-2.0_wp
    CALL testsbase%change_base(sf, 1, oldDOF,newDOF)
    y =  MINVAL(ABS(newDOF(1,:)))-MINVAL(ABS(oldDOF(1,:)))
    y2=  MINVAL(ABS(newDOF(2,:)))-MINVAL(ABS(oldDOF(2,:)))
    DEALLOCATE(oldDOF,newDOF)
    CALL testsbase%free()
    DEALLOCATE(testsbase)
    IF(testdbg.OR.(.NOT.((y  .LT.realtol ).AND. &
                         (y2 .LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be zero ...',y &
      ,'\n => should be zero ...',y2
    END IF !TEST

  END IF !testlevel>=1
  IF(testlevel.GE.2)THEN
    jElem=MAX(1,nElems/2)
    x=sf%grid%sp(jElem-1)
    IF((cont.EQ.-1).AND.(nElems.GT.1))THEN
      g_IP = testf(sf%s_IP(:))
      g_IP(1:(deg+1)*jElem) =g_IP(1:(deg+1)*jElem)-0.113_wp !discontinuous between jElem and jElem+1
    ELSE
      g_IP = testf(sf%s_IP)
    END IF !TEST

    dofs(1:sf%nBase)=sf%initDOF(g_IP)

    iTest=201 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !CHECK boundary values
    y=testf(0.0_wp)
    IF((cont.EQ.-1).AND.(nElems.GT.1))THEN
      y=y-0.113_wp
    END IF
    y2=testf(1.0_wp)
    IF(testdbg.OR.(.NOT.((ABS( dofs(1)-y      ).LT.realtol ).AND.&
                         (ABS( dofs(nBase)-y2 ).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y  ,' : dofs(1)     = ' , dofs(1)     &
      ,'\n => should be ',y2 ,' : dofs(nBase) = ' , dofs(nBase)
    END IF !TEST

    !check dsAxis and dsEdge, basis evaluation
    iTest=202 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    y=SUM(sf%base_dsAxis(0,:)*dofs(1:deg+1))
    IF((cont.EQ.-1).AND.(nElems.GT.1))THEN
      y=y+0.113_wp
    END IF
    y2=SUM(sf%base_dsEdge(0,:)*dofs(nbase-deg:nBase))

    IF(testdbg.OR.(.NOT.((ABS( y  - testf(0.0_wp) ).LT.realtol ).AND.&
                         (ABS( y2 - testf(1.0_wp) ).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',testf(0.0_wp) ,' : dofs(y=0)     = ' , y  &
      ,'\n => should be ',testf(1.0_wp) ,' : dofs(y=1)     = ' , y2
    END IF !TEST

    !check dsAxis and dsEdge, first derivative
    iTest=203 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    y =sf%evalDOF_s(0.0_wp,1,dofs)
    y2=sf%evalDOF_s(1.0_wp,1,dofs)
    IF(deg.GE.1)THEN
      dy =SUM(sf%base_dsAxis(1,:)*dofs(        1:deg+1))
      dy2=SUM(sf%base_dsEdge(1,:)*dofs(nbase-deg:nBase))

      IF(testdbg.OR.(.NOT.((ABS( dy  - testf_dx(0.0_wp) ).LT.realtol/sf%grid%ds(     1) ).AND.&
                           (ABS( dy2 - testf_dx(1.0_wp) ).LT.realtol/sf%grid%ds(nElems) ).AND.&
                           (ABS( y   - dy               ).LT.realtol/sf%grid%ds(     1) ).AND.&
                           (ABS( y2  - dy2              ).LT.realtol/sf%grid%ds(nElems) ))))THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
        '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),8(A,E11.3))') &
         '   degree = ',deg &
        ,'   continuity = ',cont,  ', nBase= ',nBase  &
        ,'\n => should be ',testf_dx(0.0_wp) ,' : dofs_ds(y=0)     = ' , dy  &
        ,'\n => should be ',testf_dx(1.0_wp) ,' : dofs_ds(y=1)     = ' , dy2 &
        ,'\n => should be ',dy                ,' : dofs_ds(y=0)     = ' , y  &
        ,'\n => should be ',dy2               ,' : dofs_ds(y=1)     = ' , y2
      END IF !TEST
    ELSE
      IF(testdbg.OR.(.NOT.((ABS( y   ).LT.realtol ).AND.&
                           (ABS( y2  ).LT.realtol ))))THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
        '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
         '   degree = ',deg &
        ,'   continuity = ',cont,  ', nBase= ',nBase  &
        ,'\n => should be ',0.0_wp ,' : dofs_ds^2(y=0)     = ' , dy  &
        ,'\n => should be ',0.0_wp ,' : dofs_ds^2(y=1)     = ' , dy2
      END IF !TEST
    END IF

    !check dsAxis and dsEdge, second derivative
    iTest=204 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    dy =sf%evalDOF_s(0.0_wp,2,dofs)
    dy2=sf%evalDOF_s(1.0_wp,2,dofs)
    IF(deg.GE.2)THEN
      y =SUM(sf%base_dsAxis(2,:)*dofs(1:deg+1))
      y2=SUM(sf%base_dsEdge(2,:)*dofs(nbase-deg:nBase))
      IF(testdbg.OR.(.NOT.((ABS( dy  - testf_dxdx(0.0_wp) ).LT.realtol/(sf%grid%ds(     1)**2) ).AND.&
                           (ABS( dy2 - testf_dxdx(1.0_wp) ).LT.realtol/(sf%grid%ds(nElems)**2) ).AND.&
                           (ABS( y   - dy                 ).LT.realtol/(sf%grid%ds(     1)**2) ).AND.&
                           (ABS( y2  - dy2                ).LT.realtol/(sf%grid%ds(nElems)**2) ))))THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
        '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),8(A,E11.3))') &
         '   degree = ',deg &
        ,'   continuity = ',cont,  ', nBase= ',nBase  &
        ,'\n => should be ',testf_dxdx(0.0_wp) ,' : dofs_ds^2(y=0)     = ' , dy  &
        ,'\n => should be ',testf_dxdx(1.0_wp) ,' : dofs_ds^2(y=1)     = ' , dy2  &
        ,'\n => should be ',dy                ,' : dofs_ds(y=0)     = ' , y  &
        ,'\n => should be ',dy2               ,' : dofs_ds(y=1)     = ' , y2
      END IF !TEST
    ELSE
      IF(testdbg.OR.(.NOT.((ABS( dy   ).LT.realtol ).AND.&
                           (ABS( dy2  ).LT.realtol ))))THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
        '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
         '   degree = ',deg &
        ,'   continuity = ',cont,  ', nBase= ',nBase  &
        ,'\n => should be ',0.0_wp ,' : dofs_ds^2(y=0)     = ' , dy  &
        ,'\n => should be ',0.0_wp ,' : dofs_ds^2(y=1)     = ' , dy2
      END IF !TEST
    END IF

    IF(nElems.GE.2)THEN

    !check eval, evalDOF_base and evalDOF_s, function evaluation
    iTest=211 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    x=sf%grid%sp(jElem)+0.9503_wp*sf%grid%ds(jElem+1) !element jElem+1
    CALL sf%eval(x ,0,iElem,base_x)
    y=sf%evalDOF_base(iElem,base_x,dofs(:))
    y2=sf%evalDOF_s(x,0,dofs(:))
    IF(testdbg.OR.(.NOT.((ABS( y-testf(x)   ).LT.realtol ).AND.&
                         (ABS( y-y2         ).LT.realtol ).AND.&
                         (iElem              .EQ.jElem+1 )     )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,I6),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',jElem,': iElem= ' , iElem &
      ,'\n => should be ',testf(x),' : y = ' , y     &
      ,'\n => should be ',y       ,' : y2= ' , y2
    END IF !TEST
    END IF ! nElems.GE.2


      !check eval, evalDOF_base and evalDOF_s, function evaluation, in jElem
      iTest=212 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    x=sf%grid%sp(jElem-1)+ 0.7353_wp*sf%grid%ds(jElem) !in element jElem
      CALL sf%eval(x ,0,iElem,base_x)
      y=sf%evalDOF_base(iElem,base_x,dofs(:))
      y2=sf%evalDOF_s(x,0,dofs(:))
    IF((cont.EQ.-1).AND.(nElems.GT.1))THEN
        y=y+0.113_wp
        y2=y2+0.113_wp
      END IF
      IF(testdbg.OR.(.NOT.((ABS( y-testf(x)   ).LT.realtol ).AND.&
                           (ABS( y-y2         ).LT.realtol ).AND.&
                           (iElem              .EQ.jElem   )     )))THEN
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
        '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
        nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,I6),4(A,E11.3))') &
         '   degree = ',deg &
        ,'   continuity = ',cont,  ', nBase= ',nBase  &
        ,'\n => should be ',jElem,': iElem= ' , iElem &
        ,'\n => should be ',testf(x),' : y = ' , y     &
        ,'\n => should be ',y       ,' : y2= ' , y2
      END IF !test

    IF(nElems.GE.2)THEN
    !test first derivative
    iTest=213 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    x=sf%grid%sp(jElem)+0.64303_wp*sf%grid%ds(jElem+1) !in elem jElem+1
    CALL sf%eval(x ,1,iElem,base_x)
    dy=sf%evalDOF_base(iElem,base_x,dofs(:))
    dy2=sf%evalDOF_s(x,1,dofs(:))
    IF(testdbg.OR.(.NOT.((ABS(dy-testf_dx(x)).LT.realtol/sf%grid%ds(jElem+1) ).AND.&
                         (ABS(dy-dy2        ).LT.realtol/sf%grid%ds(jElem+1) ).AND.&
                         (iElem              .EQ.jElem+1 ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,I6),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',jElem,': iElem= ' , iElem &
      ,'\n => should be ',testf_dx(x),': dy = ' , dy  &
      ,'\n => should be ',y          ,': dy2= ' , dy
    END IF !test

    !test second derivative
    iTest=214 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    x=sf%grid%sp(jElem)+0.17313_wp*sf%grid%ds(jElem+1) !in elem jElem+1
    CALL sf%eval(x ,2,iElem,base_x)
    dy=sf%evalDOF_base(iElem,base_x,dofs(:)) !second derivative
    dy2=sf%evalDOF_s(x,2,dofs(:))
    IF(testdbg.OR.(.NOT.((ABS(dy-testf_dxdx(x)).LT.realtol/(sf%grid%ds(jElem+1)**2) ).AND.&
                         (ABS(dy-dy2          ).LT.realtol/(sf%grid%ds(jElem+1)**2) ).AND.&
                         (iElem                .EQ.jElem+1 ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,I6),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',jElem,': iElem= ' , iElem &
      ,'\n => should be ',testf_dxdx(x),': dy = ' , dy  &
      ,'\n => should be ',y          ,': dy2= ' , dy2
    END IF !test
    END IF ! nElems.GE.2

    !check evalDOF_GP
    iTest=221 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    y_GP = testf(sf%s_GP)
    IF((cont.EQ.-1).AND.(nElems.GT.1)) y_GP(1:(degGP+1)*jElem)=y_GP(1:(degGP+1)*jElem)-0.113_wp
    dof_GP = sf%evalDOF_GP(0,dofs)

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_GP-dof_GP)).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y_GP(1)        ,':     dof_GP(1) = ' , dof_GP(1)  &
      ,'\n => should be ',MAXVAL(ABS(y_GP)) ,': MAX(|dof_GP|) = ' , MAXVAL(ABS(dof_GP))
    END IF !test

    !check integration, full domain
    iTest=222 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    y  = testf_int(0.0_wp,1.0_wp)
    IF((cont.EQ.-1).AND.(nElems.GT.1)) y=y-0.113_wp*(sf%grid%sp(jElem)-sf%grid%sp(0))
    y2 = SUM(dof_GP*sf%w_GP)

    IF(testdbg.OR.(.NOT.( (ABS(y-y2).LT.realtol ) )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y   ,': SUM(dof_GP*wGP) = ' , y2
    END IF !test

    !check integration, last element
    iTest=223 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    y  = testf_int(sf%grid%sp(nElems-1),1.0_wp)
    y2 = SUM(  dof_GP(1+(degGP+1)*(nElems-1):(degGP+1)*nElems) &
             *sf%w_GP(1+(degGP+1)*(nElems-1):(degGP+1)*nElems))

    IF(testdbg.OR.(.NOT.( (ABS(y-y2).LT.realtol ) )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),2(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y   ,': SUM(dof_GP*wGP)|_lastElem = ' , y2
    END IF !test

    !check evalDOF_GP, first derivative
    iTest=224 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    y_GP = testf_dx(sf%s_GP)
    dof_GP = sf%evalDOF_GP(DERIV_S,dofs)

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_GP-dof_GP)).LT.realtol/MINVAL(sf%grid%ds) ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y_GP(1)        ,':     dof_ds_GP(1) = ' , dof_GP(1)  &
      ,'\n => should be ',MAXVAL(ABS(y_GP)) ,': MAX(|dof_ds_GP|) = ' , MAXVAL(ABS(dof_GP))
    END IF !test

    !apply BC
    iTest=231 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_DIRICHLET
    BC_Type(2)=BC_TYPE_OPEN
    BC_Val(1:2)=(/0.96_wp,0.63_wp/)
    CALL sf%applyBCtoDOF(g_IP,BC_Type,BC_val)

    y_BC=g_IP(1:deg+1)

    y2_BC(    0)=0.96_wp
    y2_BC(1:deg)=dofs(2:deg+1)

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC-y2_BC)).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y_BC(0)        ,':  y(BCaxis)        = ' , y2_BC(0)  &
      ,'\n => should be ',MAXVAL(ABS(y_BC)) ,': MAX(|y(BC_axis)|) = ' , MAXVAL(ABS(y2_BC))
    END IF !test

    iTest=232 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_OPEN
    BC_Type(2)=BC_TYPE_DIRICHLET
    CALL sf%applyBCtoDOF(g_IP,BC_type,BC_Val)

    y_BC=g_IP(nBase-deg:nBase)

    y2_BC(0:deg-1)=dofs(nBase-deg:nBase-1)
    y2_BC(    deg)=0.63_wp

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC-y2_BC)).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),4(A,E11.3))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,'\n => should be ',y_BC(deg)      ,':  y(BCedge)       = ' , y2_BC(deg)  &
      ,'\n => should be ',MAXVAL(ABS(y_BC)) ,': MAX(|y(BCedge)|) = ' , MAXVAL(ABS(y2_BC))
    END IF !test

    iTest=233 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_NEUMANN
    BC_Type(2)=BC_TYPE_OPEN
    BC_Val(1:2)=(/0.96_wp,0.63_wp/)
    CALL sf%applyBCtoDOF(g_IP,BC_Type,BC_val)

    y_BC=MATMUL(sf%base_dsAxis(:,:),g_IP(1:deg+1))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(1)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((ABS(y_BC(1)).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y''... : dy/ds BC_NEUMANN axis = " , y_BC(:)
    END IF !test

    iTest=234 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_OPEN
    BC_Type(2)=BC_TYPE_NEUMANN
    CALL sf%applyBCtoDOF(g_IP,BC_type,BC_Val)

    y_BC=MATMUL(sf%base_dsEdge(:,:),g_IP(nBase-deg:nBase))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(nElems)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((ABS(y_BC(1)).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y'',... : dy/ds BC_NEUMANN edge = " ,  y_BC(:)
    END IF !test

    iTest=235 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_SYMM
    BC_Type(2)=BC_TYPE_OPEN
    BC_Val(1:2)=(/0.96_wp,0.63_wp/)
    CALL sf%applyBCtoDOF(g_IP,BC_Type,BC_val)

    y_BC=MATMUL(sf%base_dsAxis(:,:),g_IP(1:deg+1))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(1)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(1:deg:2))).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y''... : dy/ds BC_SYMM axis = " , y_BC(:)
    END IF !test

    iTest=236 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_OPEN
    BC_Type(2)=BC_TYPE_SYMM
    CALL sf%applyBCtoDOF(g_IP,BC_type,BC_Val)

    y_BC=MATMUL(sf%base_dsEdge(:,:),g_IP(nBase-deg:nBase))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(nElems)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(1:deg:2))).LT.realtol ))))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y'',... : dy/ds BC_SYMM edge = " ,  y_BC(:)
    END IF !test

    iTest=237 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_SYMMZERO
    BC_Type(2)=BC_TYPE_OPEN
    BC_Val(1:2)=(/0.96_wp,0.63_wp/)
    CALL sf%applyBCtoDOF(g_IP,BC_Type,BC_val)

    y_BC=MATMUL(sf%base_dsAxis(:,:),g_IP(1:deg+1))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(1)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(1:deg:2))).LT.realtol ).AND. &
                         (ABS(y_BC(0)              ).LT.realtol )     )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y''... : dy/ds BC_SYMMZERO axis = " , y_BC(:)
    END IF !test

    iTest=238 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_OPEN
    BC_Type(2)=BC_TYPE_SYMMZERO
    CALL sf%applyBCtoDOF(g_IP,BC_type,BC_Val)

    y_BC=MATMUL(sf%base_dsEdge(:,:),g_IP(nBase-deg:nBase))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(nElems)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(1:deg:2))).LT.realtol ).AND. &
                         (ABS(y_BC(0)              ).LT.realtol )     )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y'',... : dy/ds BC_SYMMZERO edge = " ,  y_BC(:)
    END IF !test

    iTest=239 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_ANTISYMM
    BC_Type(2)=BC_TYPE_OPEN
    BC_Val(1:2)=(/0.96_wp,0.63_wp/)
    CALL sf%applyBCtoDOF(g_IP,BC_Type,BC_val)

    y_BC=MATMUL(sf%base_dsAxis(:,:),g_IP(1:deg+1))
    DO i=1,deg
      y_BC(i)=y_BC(i)*(sf%grid%ds(1)**i)/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(2:deg:2))).LT.realtol ).AND. &
                         (ABS(y_BC(0)              ).LT.realtol )     )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y''... : dy/ds BC_ANTISYMM axis = " , y_BC(:)
    END IF !test

    iTest=240 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=BC_TYPE_OPEN
    BC_Type(2)=BC_TYPE_ANTISYMM
    CALL sf%applyBCtoDOF(g_IP,BC_type,BC_Val)

    y_BC=MATMUL(sf%base_dsEdge(:,:),g_IP(nBase-deg:nBase))
    DO i=1,deg
      y_BC(i)=y_BC(i)*sf%grid%ds(nElems)**i/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(2:deg:2))).LT.realtol ).AND. &
                         (ABS(y_BC(0)              ).LT.realtol )     )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y,y',y'',... : dy/ds BC_ANTISYMM edge = " ,  y_BC(:)
    END IF !test

    iTest=241 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP = dofs !
    BC_Type(1)=-deg  ! TEST  m-dependent BC_AXIS...iBC=-m,  m=deg <=> function and all derivatives up to deg-1 are zero
    BC_Type(2)=BC_TYPE_OPEN
    CALL sf%applyBCtoDOF(g_IP,BC_type,BC_Val)

    y_BC=MATMUL(sf%base_dsAxis(:,:),g_IP(1:deg+1))
    DO i=1,deg
      y_BC(i)=y_BC(i)*(sf%grid%ds(1)**i)/REAL(i*i,wp)
    END DO

    IF(testdbg.OR.(.NOT.((MAXVAL(ABS(y_BC(1:deg-1))).LT.realtol ).AND. &
                         (ABS(y_BC(0)              ).LT.realtol )     )))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(3(A,I4),A,*(E13.5))') &
       '   degree = ',deg &
      ,'   continuity = ',cont,  ', nBase= ',nBase  &
      ,"\n => should be y, dy^k/ds, k=1,deg-1 at axis = " , y_BC(:)
    END IF !test

  END IF !testlevel>=2

  test_called=.FALSE.

  CONTAINS

    !FUNCTION WITH CONTINUITY deg-1 at sp(jElem) (discontinuity must be treated on an array level, see above)
    ELEMENTAL FUNCTION testf(x)
      REAL(wp),INTENT(IN) :: x
      REAL(wp)          :: testf
      IF(sf%deg.EQ.0) THEN
        IF(x.GT.sf%grid%sp(jElem))THEN
          testf=0.21_wp
        ELSE
          testf=0.05_wp
        END IF
      ELSE
        testf=(1.13_wp-0.3_wp*x)**sf%deg &
               + MAX(0.0_wp,0.21_wp*(x-sf%grid%sp(jElem)))**sf%deg
      END IF
    END FUNCTION testf

    ELEMENTAL FUNCTION testf_dx(x)
      REAL(wp),INTENT(IN) :: x
      REAL(wp)          :: testf_dx
      IF(sf%deg.EQ.0)THEN
        testf_dx=0.0_wp
      ELSEIF(sf%deg.EQ.1)THEN
        IF(x.GT.sf%grid%sp(jElem))THEN
          testf_dx=-0.3_wp*REAL(sf%deg,wp)+ 0.21_wp*REAL(sf%deg,wp)
        ELSE
          testf_dx=-0.3_wp*REAL(sf%deg,wp)
        END IF
      ELSE
        testf_dx=-0.3_wp*REAL(sf%deg,wp)*(1.13_wp-0.3_wp*x)**(sf%deg-1)  &
                 + 0.21_wp*REAL(sf%deg,wp)*MAX(0.0_wp,0.21_wp*(x-sf%grid%sp(jElem)))**(sf%deg-1)
      END IF
    END FUNCTION testf_dx

    ELEMENTAL FUNCTION testf_dxdx(x)
      REAL(wp),INTENT(IN) :: x
      REAL(wp)          :: testf_dxdx
      IF(sf%deg.LE.1)THEN
        testf_dxdx=0.0_wp
      ELSEIF(sf%deg.EQ.2)THEN
        IF(x.GT.sf%grid%sp(jElem))THEN
          testf_dxdx=  (0.3_wp)**2*REAL(sf%deg*(sf%deg-1),wp) &
                     + (0.21_wp)**2*REAL(sf%deg*(sf%deg-1),wp)
        ELSE
          testf_dxdx=(0.3_wp)**2*REAL(sf%deg*(sf%deg-1),wp)
        END IF
      ELSE
        testf_dxdx=(0.3_wp)**2*REAL(sf%deg*(sf%deg-1),wp)*(1.13_wp-0.3_wp*x)**(sf%deg-2) &
                   + (0.21_wp)**2*REAL(sf%deg*(sf%deg-1),wp)*MAX(0.0_wp,0.21_wp*(x-sf%grid%sp(jElem)))**(sf%deg-2)
      END IF
    END FUNCTION testf_dxdx

    FUNCTION testf_int(a,b)
      REAL(wp),INTENT(IN) :: a,b !!a<=b
      REAL(wp)          :: testf_int
      IF(a.GT.sf%grid%sp(jElem))THEN
        IF(sf%deg.EQ.0)THEN
          testf_int=0.21_wp*(b-a)
        ELSE
          testf_int=   1.0_wp/(-0.3_wp*REAL(sf%deg+1,wp))*( (1.13_wp-0.3_wp*b)**(sf%deg+1) &
                                                           -(1.13_wp-0.3_wp*a)**(sf%deg+1)) &
                     + 1.0_wp/(0.21_wp*REAL(sf%deg+1,wp))*( (0.21_wp*(b-sf%grid%sp(jElem)))**(sf%deg+1) &
                                                           -(0.21_wp*(a-sf%grid%sp(jElem)))**(sf%deg+1) )
        END IF
      ELSEIF(b.GE.sf%grid%sp(jElem))THEN
        IF(sf%deg.EQ.0)THEN
          testf_int=0.21_wp*(b-sf%grid%sp(jElem))+0.05_wp*(sf%grid%sp(jElem)-a)
        ELSE
          testf_int=   1.0_wp/(-0.3_wp*REAL(sf%deg+1,wp))*( (1.13_wp-0.3_wp*b)**(sf%deg+1) &
                                                           -(1.13_wp-0.3_wp*a)**(sf%deg+1)) &
                     + 1.0_wp/(0.21_wp*REAL(sf%deg+1,wp))*( (0.21_wp*(b-sf%grid%sp(jElem)))**(sf%deg+1) )
        END IF
      ELSE
        IF(sf%deg.EQ.0)THEN
          testf_int=0.05_wp*(b-a)
        ELSE
          testf_int=   1.0_wp/(-0.3_wp*REAL(sf%deg+1,wp))*( (1.13_wp-0.3_wp*b)**(sf%deg+1) &
                                                           -(1.13_wp-0.3_wp*a)**(sf%deg+1))
        END IF
      END IF
    END FUNCTION testf_int

END SUBROUTINE sbase_test


END MODULE MODgvec_sBase
