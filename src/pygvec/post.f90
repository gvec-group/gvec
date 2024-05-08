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
LOGICAL :: initialized = .FALSE.

CONTAINS

!================================================================================================================================!
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

  initialized = .TRUE.
END SUBROUTINE Init

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
  CHARACTER(LEN=2), INTENT(IN) :: var           ! selection string: which variable to evaluate
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_s   ! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_f   ! selection string: which derivative to evaluate for the fourier series
  REAL, INTENT(OUT) :: result(n_s,n_tz)         ! output array
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i_s                          ! loop variables
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
  INTEGER :: i_s                          ! loop variables
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
! Evaluate the basis with a tensorproduct for the given 1D (s, theta, zeta) values
SUBROUTINE evaluate_base_tens_all(n_s, n_t, n_z, s, theta, zeta, Qsel, Q, dQ_ds, dQ_dthet, dQ_dzeta, &
                                  dQ_dss, dQ_dst, dQ_dsz, dQ_dtt, dQ_dtz, dQ_dzz)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,U
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s, n_t, n_z                                            ! number of evaluation points
  REAL, INTENT(IN) :: s(n_s), theta(n_t), zeta(n_z)                               ! evaluation points to construct a mesh
  CHARACTER(LEN=2) :: Qsel                                                        ! selection string: which variable to evaluate
  REAL, INTENT(OUT), DIMENSION(n_s,n_t,n_z) :: Q, dQ_ds, dQ_dthet, dQ_dzeta       ! reference space position & derivatives
  REAL, INTENT(OUT), DIMENSION(n_s,n_t,n_z) :: dQ_dss, dQ_dst, dQ_dsz, dQ_dtt, dQ_dtz, dQ_dzz
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i                                                         ! loop variables
  CLASS(t_base), POINTER :: base                                       ! pointer to the base object (X1, X2, LA)
  REAL, POINTER :: solution_dofs(:,:)                                  ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  REAL, ALLOCATABLE, DIMENSION(:) :: Q_dofs, dQ_ds_dofs, dQ_dss_dofs   ! DOFs for the fourier series
  REAL, ALLOCATABLE :: intermediate(:)                                 ! intermediate result array before reshaping
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(Qsel)
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
      WRITE(*,*) 'ERROR: variable', Qsel, 'not recognized'
      STOP
  END SELECT
  DO i=1,n_s
    ! evaluate spline to get the fourier dofs
    Q_dofs = base%s%evalDOF2D_s(s(i), base%f%modes, 0, solution_dofs(:,:))
    dQ_ds_dofs = base%s%evalDOF2D_s(s(i), base%f%modes, DERIV_S, solution_dofs(:,:))
    dQ_dss_dofs = base%s%evalDOF2D_s(s(i), base%f%modes, DERIV_S_S, solution_dofs(:,:))
    ! use the tensorproduct for theta and zeta
    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, 0, Q_dofs)
    Q(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, 0, dQ_ds_dofs)
    dQ_ds(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_THET, Q_dofs)
    dQ_dthet(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_ZETA, Q_dofs)
    dQ_dzeta(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, 0, dQ_dss_dofs)
    dQ_dss(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_THET, dQ_ds_dofs)
    dQ_dst(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_ZETA, dQ_ds_dofs)
    dQ_dsz(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_THET_THET, Q_dofs)
    dQ_dtt(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_THET_ZETA, Q_dofs)
    dQ_dtz(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))

    intermediate = base%f%evalDOF_xn_tens(n_t, n_z, theta, zeta, DERIV_ZETA_ZETA, Q_dofs)
    dQ_dzz(i,:,:) = RESHAPE(intermediate, (/n_t, n_z/))
  END DO
  DEALLOCATE(intermediate)
  DEALLOCATE(Q_dofs, dQ_ds_dofs)
END SUBROUTINE evaluate_base_tens_all

!================================================================================================================================!
! Evaluate the mapping from reference to physical space (hmap)
SUBROUTINE evaluate_hmap(n, X1, X2, zeta, dX1_ds, dX2_ds, dX1_dthet, dX2_dthet, dX1_dzeta, dX2_dzeta, coord, e_s, e_thet, e_zeta)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                                      ! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta, dX1_ds, dX2_ds                ! reference space position & derivatives
  REAL, INTENT(IN), DIMENSION(n) :: dX1_dthet, dX2_dthet, dX1_dzeta, dX2_dzeta  ! reference space derivatives
  REAL, INTENT(OUT), DIMENSION(3,n) :: coord, e_s, e_thet, e_zeta               ! real space position and basis vectors
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i              ! loop variable
  REAL, DIMENSION(3) :: q   ! position in reference space
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  DO i=1,n
    q = (/X1(i),X2(i),zeta(i)/)
    coord(:,i)  = hmap%eval(q)
    e_s(:,i)    = hmap%eval_dxdq(q,(/dX1_ds(i),dX2_ds(i),0.0/))
    e_thet(:,i) = hmap%eval_dxdq(q,(/dX1_dthet(i),dX2_dthet(i),0.0/))
    e_zeta(:,i) = hmap%eval_dxdq(q,(/dX1_dzeta(i),dX2_dzeta(i),1.0/))
  END DO
END SUBROUTINE

!================================================================================================================================!
! Evaluate the mapping from reference to physical space (hmap) without logical coordinates
SUBROUTINE evaluate_hmap_only(n, X1, X2, zeta, pos, e_X1, e_X2, e_zeta3)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                      ! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta                ! reference space position
  REAL, INTENT(OUT), DIMENSION(3,n) :: pos, e_X1, e_X2, e_zeta3 ! real space position and reference tangent basis vectors
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i              ! loop variable
  REAL, DIMENSION(3) :: q   ! position in reference space
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  DO i=1,n
    q = (/X1(i),X2(i),zeta(i)/)
    pos(:,i)     = hmap%eval(q)
    e_X1(:,i)    = hmap%eval_dxdq(q,(/1.0, 0.0, 0.0/))
    e_X2(:,i)    = hmap%eval_dxdq(q,(/0.0, 1.0, 0.0/))
    e_zeta3(:,i) = hmap%eval_dxdq(q,(/0.0, 0.0, 1.0/))
  END DO
END SUBROUTINE

!================================================================================================================================!
! evaluate components of the metric tensor and their derivatives
SUBROUTINE evaluate_metric(n, X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz, &
                           dX1_dss, dX2_dss, dX1_dst, dX2_dst, dX1_dsz, dX2_dsz, &
                           dX1_dtt, dX2_dtt, dX1_dtz, dX2_dtz, dX1_dzz, dX2_dzz, &
                           g_ss, g_st, g_sz, g_tt, g_tz, g_zz, &
                           dg_ss_ds, dg_st_ds, dg_sz_ds, dg_tt_ds, dg_tz_ds, dg_zz_ds, &
                           dg_ss_dt, dg_st_dt, dg_sz_dt, dg_tt_dt, dg_tz_dt, dg_zz_dt, &
                           dg_ss_dz, dg_st_dz, dg_sz_dz, dg_tt_dz, dg_tz_dz, dg_zz_dz)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                                                        ! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz  ! reference coordinates
  REAL, INTENT(IN), DIMENSION(n) :: dX1_dss, dX2_dss, dX1_dst, dX2_dst, dX1_dsz, dX2_dsz          ! and their derivatives
  REAL, INTENT(IN), DIMENSION(n) :: dX1_dtt, dX2_dtt, dX1_dtz, dX2_dtz, dX1_dzz, dX2_dzz
  REAL, INTENT(OUT), DIMENSION(n) :: g_ss, g_st, g_sz, g_tt, g_tz, g_zz                           ! metric coefficients
  REAL, INTENT(OUT), DIMENSION(n) :: dg_ss_ds, dg_st_ds, dg_sz_ds, dg_tt_ds, dg_tz_ds, dg_zz_ds   ! derivatives of the m. coef.
  REAL, INTENT(OUT), DIMENSION(n) :: dg_ss_dt, dg_st_dt, dg_sz_dt, dg_tt_dt, dg_tz_dt, dg_zz_dt   ! derivatives of the m. coef.
  REAL, INTENT(OUT), DIMENSION(n) :: dg_ss_dz, dg_st_dz, dg_sz_dz, dg_tt_dz, dg_tz_dz, dg_zz_dz   ! derivatives of the m. coef.
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i                                                                ! loop variable
  REAL, DIMENSION(3) :: q, q_s, q_t, q_z, q_ss, q_st, q_sz, q_tt, q_tz, q_zz  ! position in reference space
  REAL :: g_ss_dq1, g_ss_dq2, g_st_dq1, g_st_dq2, g_sz_dq1, g_sz_dq2
  REAL :: g_tt_dq1, g_tt_dq2, g_tz_dq1, g_tz_dq2, g_zz_dq1, g_zz_dq2
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  DO i=1,n
    q   = (/X1(i), X2(i), zeta(i)/)
    q_s = (/dX1_ds(i), dX2_ds(i), 0.0/)
    q_t = (/dX1_dt(i), dX2_dt(i), 0.0/)
    q_z = (/dX1_dz(i), dX2_dz(i), 1.0/)
    q_ss = (/dX1_dss(i), dX2_dss(i), 0.0/)
    q_st = (/dX1_dst(i), dX2_dst(i), 0.0/)
    q_sz = (/dX1_dsz(i), dX2_dsz(i), 0.0/)
    q_tt = (/dX1_dtt(i), dX2_dtt(i), 0.0/)
    q_tz = (/dX1_dtz(i), dX2_dtz(i), 0.0/)
    q_zz = (/dX1_dzz(i), dX2_dzz(i), 0.0/)

    g_ss(i) = hmap%eval_gij(q_s, q, q_s)
    g_ss_dq1 = hmap%eval_gij_dq1(q_s, q, q_s)
    g_ss_dq2 = hmap%eval_gij_dq2(q_s, q, q_s)
    dg_ss_ds(i) = 2 * hmap%eval_gij(q_ss, q, q_s) + q_s(1) * g_st_dq1 + q_s(2) * g_st_dq2
    dg_ss_dt(i) = 2 * hmap%eval_gij(q_st, q, q_s) + q_t(1) * g_st_dq1 + q_t(2) * g_st_dq2
    dg_ss_dz(i) = 2 * hmap%eval_gij(q_sz, q, q_s) + q_z(1) * g_st_dq1 + q_z(2) * g_st_dq2
    g_st(i) = hmap%eval_gij(q_s, q, q_t)
    g_st_dq1 = hmap%eval_gij_dq1(q_s, q, q_t)
    g_st_dq2 = hmap%eval_gij_dq2(q_s, q, q_t)
    dg_st_ds(i) = hmap%eval_gij(q_ss, q, q_t) + hmap%eval_gij(q_s, q, q_st) + q_s(1) * g_st_dq1 + q_s(2) * g_st_dq2
    dg_st_dt(i) = hmap%eval_gij(q_st, q, q_t) + hmap%eval_gij(q_s, q, q_tt) + q_t(1) * g_st_dq1 + q_t(2) * g_st_dq2
    dg_st_dz(i) = hmap%eval_gij(q_sz, q, q_t) + hmap%eval_gij(q_s, q, q_tz) + q_z(1) * g_st_dq1 + q_z(2) * g_st_dq2
    g_sz(i) = hmap%eval_gij(q_s, q, q_z)
    g_sz_dq1 = hmap%eval_gij_dq1(q_s, q, q_z)
    g_sz_dq2 = hmap%eval_gij_dq2(q_s, q, q_z)
    dg_sz_ds(i) = hmap%eval_gij(q_ss, q, q_z) + hmap%eval_gij(q_s, q, q_sz) + q_s(1) * g_sz_dq1 + q_s(2) * g_sz_dq2
    dg_sz_dt(i) = hmap%eval_gij(q_st, q, q_z) + hmap%eval_gij(q_s, q, q_tz) + q_t(1) * g_sz_dq1 + q_t(2) * g_sz_dq2
    dg_sz_dz(i) = hmap%eval_gij(q_sz, q, q_z) + hmap%eval_gij(q_s, q, q_zz) + q_z(1) * g_sz_dq1 + q_z(2) * g_sz_dq2
    g_tt(i) = hmap%eval_gij(q_t, q, q_t)
    g_tt_dq1 = hmap%eval_gij_dq1(q_t, q, q_t)
    g_tt_dq2 = hmap%eval_gij_dq2(q_t, q, q_t)
    dg_tt_ds(i) = 2 * hmap%eval_gij(q_st, q, q_t) + q_s(1) * g_tt_dq1 + q_s(2) * g_tt_dq2
    dg_tt_dt(i) = 2 * hmap%eval_gij(q_tt, q, q_t) + q_t(1) * g_tt_dq1 + q_t(2) * g_tt_dq2
    dg_tt_dz(i) = 2 * hmap%eval_gij(q_tz, q, q_t) + q_z(1) * g_tt_dq1 + q_z(2) * g_tt_dq2
    g_tz(i) = hmap%eval_gij(q_t, q, q_z)
    g_tz_dq1 = hmap%eval_gij_dq1(q_t, q, q_z)
    g_tz_dq2 = hmap%eval_gij_dq2(q_t, q, q_z)
    dg_tz_ds(i) = hmap%eval_gij(q_st, q, q_z) + hmap%eval_gij(q_t, q, q_sz) + q_s(1) * g_tz_dq1 + q_s(2) * g_tz_dq2
    dg_tz_dt(i) = hmap%eval_gij(q_tt, q, q_z) + hmap%eval_gij(q_t, q, q_tz) + q_t(1) * g_tz_dq1 + q_t(2) * g_tz_dq2
    dg_tz_dz(i) = hmap%eval_gij(q_tz, q, q_z) + hmap%eval_gij(q_t, q, q_zz) + q_z(1) * g_tz_dq1 + q_z(2) * g_tz_dq2
    g_zz(i) = hmap%eval_gij(q_z, q, q_z)
    g_zz_dq1 = hmap%eval_gij_dq1(q_z, q, q_z)
    g_zz_dq2 = hmap%eval_gij_dq2(q_z, q, q_z)
    dg_zz_ds(i) = 2 * hmap%eval_gij(q_sz, q, q_z) + q_s(1) * g_zz_dq1 + q_s(2) * g_zz_dq2
    dg_zz_dt(i) = 2 * hmap%eval_gij(q_tz, q, q_z) + q_t(1) * g_zz_dq1 + q_t(2) * g_zz_dq2
    dg_zz_dz(i) = 2 * hmap%eval_gij(q_zz, q, q_z) + q_z(1) * g_zz_dq1 + q_z(2) * g_zz_dq2
  END DO
END SUBROUTINE

!================================================================================================================================!
! evaluate the jacobian determinant and its derivatives
SUBROUTINE evaluate_jacobian(n, X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz, Jh, dJh_ds, dJh_dt, dJh_dz)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                                                        ! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz  ! reference coordinates
  REAL, INTENT(OUT), DIMENSION(n) :: Jh, dJh_ds, dJh_dt, dJh_dz                                   ! jacobian det. and derivatives
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i              ! loop variable
  REAL, DIMENSION(3) :: q   ! position in reference space
  REAL :: Jh_dq1, Jh_dq2
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  DO i=1,n
    q   = (/X1(i), X2(i), zeta(i)/)
    Jh(i)  = hmap%eval_Jh(q)
    Jh_dq1 = hmap%eval_Jh_dq1(q)
    Jh_dq2 = hmap%eval_Jh_dq2(q)
    dJh_ds(i) = dX1_ds(i) * Jh_dq1 + dX2_ds(i) * Jh_dq2
    dJh_dt(i) = dX1_dt(i) * Jh_dq1 + dX2_dt(i) * Jh_dq2
    dJh_dz(i) = dX1_dz(i) * Jh_dq1 + dX2_dz(i) * Jh_dq2
  END DO
END SUBROUTINE

!================================================================================================================================!
SUBROUTINE evaluate_profile(n_s, s, var, result)
  ! MODULES
  USE MODgvec_MHD3D_profiles, ONLY: Eval_iota, Eval_iota_Prime, Eval_pres, Eval_p_prime, Eval_mass, Eval_chi, Eval_chiPrime, &
                                    Eval_Phi, Eval_PhiPrime, Eval_Phi_TwoPrime, Eval_PhiNorm, Eval_PhiNormPrime
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s                  ! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n_s) :: s       ! radial evaluation points
  CHARACTER(LEN=*), INTENT(IN) :: var         ! selection string: which profile to evaluate
  REAL, INTENT(OUT), DIMENSION(n_s) :: result ! values of the profile
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i  ! loop variable
  PROCEDURE(Eval_iota), POINTER :: eval_profile
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(TRIM(var))
    CASE("iota")
      eval_profile => Eval_iota
    CASE("iota_prime")
      eval_profile => Eval_iota_Prime
    CASE("p")
      eval_profile => Eval_pres
    CASE("p_prime")
      eval_profile => Eval_p_prime
    CASE("mass")
      eval_profile => Eval_mass
    CASE("chi")
      eval_profile => Eval_chi
    CASE("chi_prime")
      eval_profile => Eval_chiPrime
    CASE("Phi")
      eval_profile => Eval_Phi
    CASE("Phi_prime")
      eval_profile => Eval_PhiPrime
    CASE("Phi_2prime")
      eval_profile => Eval_Phi_TwoPrime
    CASE("PhiNorm")
      eval_profile => Eval_PhiNorm
    CASE("PhiNorm_prime")
      eval_profile => Eval_PhiNormPrime
    CASE DEFAULT
      WRITE(*,*) 'ERROR: variable', var, 'not recognized'
      STOP
  END SELECT

  DO i = 1,n_s
    result(i) = eval_profile(s(i))
  END DO
END SUBROUTINE evaluate_profile

!================================================================================================================================!
SUBROUTINE Finalize()
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_Analyze,        ONLY: FinalizeAnalyze
  USE MODgvec_Output,         ONLY: FinalizeOutput
  USE MODgvec_Restart,        ONLY: FinalizeRestart
  USE MODgvec_Functional,     ONLY: FinalizeFunctional
  USE MODgvec_ReadInTools,    ONLY: FinalizeReadIn
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL FinalizeFunctional(functional)
  DEALLOCATE(functional)
  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  CALL FinalizeRestart()
  CALL FinalizeReadIn()
  initialized = .FALSE.

  WRITE(Unit_stdOut,'(132("="))')
  WRITE(UNIT_stdOut,'(A)') "GVEC POST FINISHED !"
  WRITE(Unit_stdOut,'(132("="))')
END SUBROUTINE Finalize

END MODULE MODpygvec_post
