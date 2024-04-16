!================================================================================================================================!
! Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!================================================================================================================================!
! PyGVEC postprocessing
!
! This module contains the Python interface for postprocessing the GVEC output, starting with a parameter file and a statefile.
!================================================================================================================================!

#include "defines.h"

MODULE MODpygvec_post

USE MODgvec_c_functional, ONLY: t_functional

IMPLICIT NONE

PUBLIC

CLASS(t_functional), ALLOCATABLE :: functional

! TYPE :: t_evaluations
!   REAL, ALLOCATABLE, DIMENSION(  :,:,:) :: X1, X2
!   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: position
! END TYPE t_evaluations

CONTAINS

!================================================================================================================================!
! a simple scalar function to test the wrapper
SUBROUTINE Init(parameterfile)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Analyze,        ONLY: InitAnalyze
  USE MODgvec_Output,         ONLY: InitOutput
  USE MODgvec_Restart,        ONLY: InitRestart
  USE MODgvec_ReadInTools,    ONLY: FillStrings,GETLOGICAL,GETINT,IgnoredStrings 
  USE MODgvec_Functional,     ONLY: InitFunctional
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=*) :: parameterfile
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: which_functional
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A)') "GVEC POST ! GVEC POST ! GVEC POST ! GVEC POST"
  WRITE(Unit_stdOut,'(132("="))')

  ! read parameter file
  CALL FillStrings(parameterfile)
  
  ! initialization phase
  CALL InitRestart()
  CALL InitOutput()
  CALL InitAnalyze()
  
  ! initialize the functional
  which_functional = GETINT('which_functional',Proposal=1)
  CALL InitFunctional(functional,which_functional)
  
  ! print the ignored parameters
  CALL IgnoredStrings()
END SUBROUTINE init

!================================================================================================================================!
SUBROUTINE ReadState(statefile)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Output_Vars,    ONLY: outputLevel,ProjectName
  USE MODgvec_ReadState_Vars, ONLY: fileID_r,outputLevel_r
  USE MODgvec_Restart,        ONLY: RestartFromState
  USE MODgvec_MHD3D_Vars,     ONLY: U
  USE MODgvec_MHD3D_visu,     ONLY: Get_SFL_theta
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=255) :: statefile
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  ProjectName='POST_'//TRIM(StateFile(1:INDEX(StateFile,'_State_')-1))
  CALL RestartFromState(StateFile,U(0))
  ! CALL EvalForce(U(0),.TRUE.,2, F(0))
  ! CALL Analyze(FileID_r)
  outputLevel=outputLevel_r
END SUBROUTINE ReadState

!================================================================================================================================!
! Handle the selection of the functional and derivatives, based on the selection strings
SUBROUTINE evaluate_base_select(var, sel_deriv_s, sel_deriv_f, base, solution_dofs, seli_deriv_s, seli_deriv_f)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,U
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var                 ! selection string: which variable to evaluate
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_s         ! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_f         ! selection string: which derivative to evaluate for the fourier series
  CLASS(t_base), POINTER, INTENT(OUT) :: base         ! pointer to the base object (X1, X2, LA)
  REAL, POINTER, INTENT(OUT) :: solution_dofs(:,:)    ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  INTEGER, INTENT(OUT) :: seli_deriv_s, seli_deriv_f  ! integer values for the derivative selection
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(var)
    CASE('X1')
      base => X1_base
      solution_dofs => U(0)%X1
    CASE('X2')
      base => X2_base
      solution_dofs => U(0)%X2
    CASE('LA')
      base => LA_base
      solution_dofs => U(0)%LA
    CASE DEFAULT
      WRITE(*,*) 'ERROR: variable', var, 'not recognized'
      STOP
  END SELECT
  SELECT CASE(TRIM(sel_deriv_s))
    CASE('s')
      seli_deriv_s = DERIV_S
    CASE('ss')
      seli_deriv_s = DERIV_S_S
    CASE DEFAULT
      seli_deriv_s = 0
  END SELECT
  SELECT CASE(TRIM(sel_deriv_f))
    CASE('t')
      seli_deriv_f = DERIV_THET
    CASE('tt')
      seli_deriv_f = DERIV_THET_THET
    CASE('z')
      seli_deriv_f = DERIV_ZETA
    CASE('zz')
      seli_deriv_f = DERIV_ZETA_ZETA
    CASE("tz")
      seli_deriv_f = DERIV_THET_ZETA
    CASE DEFAULT
      seli_deriv_f = 0
  END SELECT
END SUBROUTINE evaluate_base_select

!================================================================================================================================!
! Evaluate the basis for a list of (theta, zeta) positions on all flux surfaces given by s
SUBROUTINE evaluate_base_list_tz(n_s, n_tz, s, thetazeta, var, sel_deriv_s, sel_deriv_f, result)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s, n_tz              ! number of evaluation points
  REAL, INTENT(IN) :: s(n_s), thetazeta(2,n_tz) ! evaluation points
  CHARACTER(LEN=2) :: var                       ! selection string: which variable to evaluate
  CHARACTER(LEN=2) :: sel_deriv_s               ! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2) :: sel_deriv_f               ! selection string: which derivative to evaluate for the fourier series
  REAL, INTENT(OUT) :: result(n_s,n_tz)         ! output array
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i_s, i_t, i_z                ! loop variables
  INTEGER :: seli_deriv_s, seli_deriv_f   ! integer values for the derivative selection
  CLASS(t_base), POINTER :: base          ! pointer to the base object (X1, X2, LA)
  REAL, POINTER :: solution_dofs(:,:)     ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  REAL, ALLOCATABLE :: fourier_dofs(:)    ! DOFs for the fourier series, calculated from the spline
  REAL, ALLOCATABLE :: intermediate(:)    ! intermediate result array before reshaping
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL evaluate_base_select(var, sel_deriv_s, sel_deriv_f, base, solution_dofs, seli_deriv_s, seli_deriv_f)

  DO i_s=1,n_s
    ! evaluate spline to get the fourier dofs
    fourier_dofs = base%s%evalDOF2D_s(s(i_s),base%f%modes,seli_deriv_s,solution_dofs)
    result(i_s,:) = base%f%evalDOF_xn(n_tz, thetazeta, seli_deriv_f, fourier_dofs)
  END DO
  DEALLOCATE(fourier_dofs)
END SUBROUTINE evaluate_base_list_tz

!================================================================================================================================!
! Evaluate the basis with a tensorproduct for the given 1D (s, theta, zeta) values
SUBROUTINE evaluate_base_tens(s, theta, zeta, var, sel_deriv_s, sel_deriv_f, result)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  REAL, INTENT(IN) :: s(:), theta(:), zeta(:)   ! evaluation points to construct a mesh
  CHARACTER(LEN=2) :: var                       ! selection string: which variable to evaluate
  CHARACTER(LEN=2) :: sel_deriv_s               ! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2) :: sel_deriv_f               ! selection string: which derivative to evaluate for the fourier series
  REAL, INTENT(OUT) :: result(:,:,:)            ! output array
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i_s, i_t, i_z                ! loop variables
  INTEGER :: seli_deriv_s, seli_deriv_f   ! integer values for the derivative selection
  CLASS(t_base), POINTER :: base          ! pointer to the base object (X1, X2, LA)
  REAL, POINTER :: solution_dofs(:,:)     ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  REAL, ALLOCATABLE :: fourier_dofs(:)    ! DOFs for the fourier series, calculated from the spline
  REAL, ALLOCATABLE :: intermediate(:)    ! intermediate result array before reshaping
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL evaluate_base_select(var, sel_deriv_s, sel_deriv_f, base, solution_dofs, seli_deriv_s, seli_deriv_f)

  DO i_s=1,SIZE(s)
    ! evaluate spline to get the fourier dofs
    fourier_dofs = base%s%evalDOF2D_s(s(i_s), base%f%modes, seli_deriv_s, solution_dofs(:,:))
    ! use the tensorproduct for theta and zeta
    intermediate = base%f%evalDOF_xn_tens(SIZE(theta), SIZE(zeta), theta, zeta, seli_deriv_f, fourier_dofs)
    result(i_s,:,:) = RESHAPE(intermediate, (/SIZE(theta), SIZE(zeta)/))
  END DO
  DEALLOCATE(intermediate)
  DEALLOCATE(fourier_dofs)
END SUBROUTINE evaluate_base_tens

!================================================================================================================================!
SUBROUTINE Finalize()
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Analyze,        ONLY: FinalizeAnalyze
  USE MODgvec_Output,         ONLY: FinalizeOutput
  USE MODgvec_Restart,        ONLY: FinalizeRestart
  USE MODgvec_Functional,     ONLY: FinalizeFunctional
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL FinalizeFunctional(functional)
  DEALLOCATE(functional)
  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  CALL FinalizeRestart()
  
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A)') "GVEC POST FINISHED !"
  WRITE(Unit_stdOut,'(132("="))')
END SUBROUTINE Finalize

END MODULE MODpygvec_post