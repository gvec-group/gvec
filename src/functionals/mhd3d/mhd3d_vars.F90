!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_MHD3D_Vars
! MODULES
USE MODgvec_Globals,ONLY: PI,wp,Unit_stdOut,abort
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base
USE MODgvec_Sol_Var_MHD3D,ONLY: t_sol_var_MHD3D
USE MODgvec_hmap
USE MODgvec_rProfile_base, ONLY: c_rProfile
USE MODgvec_boundaryFromFile , ONLY: t_boundaryFromFile

IMPLICIT NONE
PUBLIC


!-----------------------------------------------------------------------------------------------------------------------------------
! derived type variables
!-----------------------------------------------------------------------------------------------------------------------------------
CLASS(t_base),  ALLOCATABLE,TARGET :: X1_base   !! container for base of variable X1
CLASS(t_base),  ALLOCATABLE,TARGET :: X2_base   !! container for base of variable X2
CLASS(t_base),  ALLOCATABLE,TARGET :: LA_base   !! container for base of variable lambda

TYPE(t_sgrid)               :: sgrid     !! only one grid up to now

TYPE(t_sol_var_MHD3D),ALLOCATABLE :: U(:)      !! solutions at levels (k-1),(k),(k+1)
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: V(:)      !! 'velocity' in minimizer
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: F(:)      !! force
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: P(:)      !! temporary for update
INTEGER                     :: nDOF_X1   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes
INTEGER                     :: nDOF_X2   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes
INTEGER                     :: nDOF_LA   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes
INTEGER,ALLOCATABLE         :: X1_BC_type(:,:) !! X1 var: BC type for axis and edge for each mode (1:2,1:modes) (1=axis,2=edge)
INTEGER,ALLOCATABLE         :: X2_BC_type(:,:) !! X2 var: BC type for axis and edge for each mode (1:2,1:modes) (1=axis,2=edge)
INTEGER,ALLOCATABLE         :: LA_BC_type(:,:) !! LA var: BC type for axis and edge for each mode (1:2,1:modes) (1=axis,2=edge)
#ifdef PP_WHICH_HMAP
TYPE(PP_T_HMAP),  ALLOCATABLE :: hmap      !! type containing subroutines for evaluating the map h (Omega_p x S^1) --> Omega
TYPE(PP_T_HMAP_AUXVAR),ALLOCATABLE :: hmap_auxvar(:) !! auxiliary variables for hmap
#else
CLASS(PP_T_HMAP),  ALLOCATABLE :: hmap      !! type containing subroutines for evaluating the map h (Omega_p x S^1) --> Omega
CLASS(PP_T_HMAP_AUXVAR),ALLOCATABLE :: hmap_auxvar(:) !! auxiliary variables for hmap
#endif
!===================================================================================================================================
INTEGER              :: which_init      !! select initialization. 0: only using input parameter, 1: using a VMEC equilibrium
INTEGER              :: which_hmap
LOGICAL              :: init_fromBCOnly !! default=TRUE, for VMEC only, if set false: initial mapping is interpolated for s=0..1
LOGICAL              :: init_with_profile_pressure !! default=FALSE, if True, overwrite profile from VMEC ini using  profile from parameterfile
LOGICAL              :: init_with_profile_iota     !! default=FALSE, if True, overwrite profile from VMEC ini using  profile from parameterfile
LOGICAL              :: init_average_axis !! default=FALSE, if true, use outer boundary to estimate axis position (center of closed line)
LOGICAL              :: boundary_perturb !! default=FALSE, if true, mapping is perturbed with a given modal perturbation of the boundary (X1pert_b,X2pert_b)
INTEGER              :: boundary_perturb_type
INTEGER, PARAMETER   :: BLEND_LEGACY=0, BLEND_COSM=1 !! types of blending functions for boundary_perturb_type
REAL(wp)             :: boundary_perturb_depth !! depth of boundary perturbation
REAL(wp)             :: average_axis_move(2) !! used if init_average_axis=True to additionally move axis in X1,X2
INTEGER              :: init_BC         !! active if init_fromBC_only=T: -1: keep vmec axis and boundary (default), 0: keep vmec boundary, overwrite axis, 1: keep vmec axis, overwrite boundary, 2: overwrite axis and boundary
INTEGER              :: getBoundaryFromFile !! -1: off, 1: read from specific netcdf file
LOGICAL              :: init_LA         !! false: lambda=0 at initialization, true: lambda is computed from initial mapping
INTEGER              :: PrecondType     !! -1: off: 1: ..
! input parameters for minimization
INTEGER              :: MinimizerType   !! which mimimizer to use: 0: gradient descent (default) , 10: accelerated gradient descent
INTEGER              :: maxIter         !! maximum iteration count for minimization
INTEGER              :: outputIter      !! number of iterations after which output files are written
INTEGER              :: logIter         !! number of iterations after which a screen log is written
INTEGER              :: nlogScreen      !! number of log outputs after a screen output is written
REAL(wp)             :: minimize_tol    !! absolute tolerance for minimization of functional
REAL(wp)             :: start_dt        !! starting time step, is adapted during iteration
REAL(wp)             :: dW_allowed      !! for minimizer, accept step if dW<dW_allowed*W_MHD(iter=0) default +10e-10
LOGICAL              :: DoCheckDistance !! TRUE: check distance between solutions of two log output states (default: false)
LOGICAL              :: DoCheckAxis     !! TRUE: check axis position (default: true)
! input parameters for functional
REAL(wp)             :: Phi_edge        !! toroidal flux at the last flux surface of the domain

!constants
REAL(wp)             :: mu_0            !! permeability
REAL(wp)             :: gamm            !! isentropic exponent, if gamma /= 0 pres ~ mass profile
REAL(wp)             :: sgammM1         !! =1/(gamm-1)


REAL(wp),ALLOCATABLE :: X1_b(:)         !! fourier modes of the edge boundary for X1
REAL(wp),ALLOCATABLE :: X2_b(:)         !! fourier modes of the edge boundary for X2
REAL(wp),ALLOCATABLE :: LA_b(:)         !! fourier modes of the edge boundary for LA
REAL(wp),ALLOCATABLE :: X1_a(:)         !! fourier modes of the axis boundary for X1
REAL(wp),ALLOCATABLE :: X2_a(:)         !! fourier modes of the axis boundary for X2
REAL(wp),ALLOCATABLE :: X1pert_b(:)     !! fourier modes of the boundary perturbation for X1 (if boundary_perturb=T)
REAL(wp),ALLOCATABLE :: X2pert_b(:)     !! fourier modes of the boundary perturbation for X2 (if boundary_perturb=T)
CLASS(t_boundaryFromFile),ALLOCATABLE:: BFF  !! class for reading a boundary from file

CLASS(c_rProfile), ALLOCATABLE   :: iota_profile
CLASS(c_rProfile), ALLOCATABLE   :: pres_profile
CLASS(c_rProfile), ALLOCATABLE   :: Phi_profile
CLASS(c_rProfile), ALLOCATABLE   :: chi_profile
!===================================================================================================================================

END MODULE MODgvec_MHD3D_Vars
