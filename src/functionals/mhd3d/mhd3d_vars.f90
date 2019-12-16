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
MODULE MODgvec_MHD3D_Vars
! MODULES
USE MODgvec_Globals,ONLY: PI,wp,Unit_stdOut,abort
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base
USE MODgvec_fbase,  ONLY: t_fbase
USE MODgvec_Sol_Var_MHD3D,ONLY: t_sol_var_MHD3D
USE MODgvec_c_hmap, ONLY: c_hmap
IMPLICIT NONE
PUBLIC


!-----------------------------------------------------------------------------------------------------------------------------------
! derived type variables
!-----------------------------------------------------------------------------------------------------------------------------------
CLASS(t_base),  ALLOCATABLE :: X1_base   !! container for base of variable X1
CLASS(t_base),  ALLOCATABLE :: X2_base   !! container for base of variable X2
CLASS(t_base),  ALLOCATABLE :: LA_base   !! container for base of variable lambda 
                             
TYPE(t_sgrid)               :: sgrid     !! only one grid up to now
                                         
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: U(:)      !! solutions at levels (k-1),(k),(k+1)
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: F(:)      !! force
TYPE(t_sol_var_MHD3D),ALLOCATABLE :: P(:)      !! temporary for update
INTEGER                     :: nDOF_X1   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER                     :: nDOF_X2   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER                     :: nDOF_LA   !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER,ALLOCATABLE         :: X1_BC_type(:,:) !! X1 var: BC type for axis and edge for each mode (1:2,1:modes) (1=axis,2=edge) 
INTEGER,ALLOCATABLE         :: X2_BC_type(:,:) !! X2 var: BC type for axis and edge for each mode (1:2,1:modes) (1=axis,2=edge) 
INTEGER,ALLOCATABLE         :: LA_BC_type(:,:) !! LA var: BC type for axis and edge for each mode (1:2,1:modes) (1=axis,2=edge) 
                                         
                                         
CLASS(c_hmap),  ALLOCATABLE :: hmap      !! type containing subroutines for evaluating the map h (Omega_p x S^1) --> Omega

!===================================================================================================================================
INTEGER              :: which_init      !! select initialization. 0: only using input parameter, 1: using a VMEC equilibrium
INTEGER              :: which_hmap
LOGICAL              :: init_fromBCOnly !! default=TRUE, for VMEC only, if set false: initial mapping is interpolated for s=0..1
LOGICAL              :: init_average_axis !! default=FALSE, if true, use outer boundary to estimate axis position (center of closed line)
REAL(wp)             :: average_axis_move(2) !! used if init_average_axis=True to additionally move axis in X1,X2   
INTEGER              :: init_BC         !! active if init_fromBC_only=T: =0: keep vmec axis and boundary, =1: overwrite axis, =2: overwrite boundary, =3: overwrite  axis and boundary 
LOGICAL              :: init_LA         !! false: lambda=0 at initialization, true: lambda is computed from initial mapping
INTEGER              :: PrecondType     !! -1: off: 1: .. 
! input parameters for minimization
INTEGER              :: MinimizerType   !! which mimimizer to use: 0: steepest descent (default) , 1: LBFGS
INTEGER              :: maxIter         !! maximum iteration count for minimization 
INTEGER              :: outputIter      !! number of iterations after which output files are written
INTEGER              :: logIter         !! number of iterations after which a screen log is written
REAL(wp)             :: minimize_tol    !! absolute tolerance for minimization of functional
REAL(wp)             :: start_dt        !! starting time step, is adapted during iteration
! input parameters for functional
REAL(wp)             :: Phi_edge        !! toroidal flux at the last flux surface of the domain
INTEGER              :: n_pres_coefs    !! number of polynomial coeffients for mass profile
INTEGER              :: n_iota_coefs    !! number of polynomial coeffients for iota profile
REAL(wp),ALLOCATABLE :: pres_coefs(:)   !! polynomial coefficients of the mass profile
REAL(wp),ALLOCATABLE :: iota_coefs(:)   !! polynomial coefficients of the iota profile
!constants
REAL(wp)             :: mu_0            !! permeability
REAL(wp)             :: gamm            !! isentropic exponent, if gamma /= 0 pres ~ mass profile
REAL(wp)             :: sgammM1         !! =1/(gamm-1)

                                        
REAL(wp),ALLOCATABLE :: X1_b(:)         !! fourier modes of the edge boundary for X1
REAL(wp),ALLOCATABLE :: X2_b(:)         !! fourier modes of the edge boundary for X2
REAL(wp),ALLOCATABLE :: LA_b(:)         !! fourier modes of the edge boundary for LA
REAL(wp),ALLOCATABLE :: X1_a(:)         !! fourier modes of the axis boundary for X1
REAL(wp),ALLOCATABLE :: X2_a(:)         !! fourier modes of the axis boundary for X2


!===================================================================================================================================

END MODULE MODgvec_MHD3D_Vars

