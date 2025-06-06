!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

MODULE MODgvec_py_post

USE MODgvec_c_functional, ONLY: t_functional
USE MODgvec_base,         ONLY: t_base

IMPLICIT NONE
PUBLIC

CLASS(t_functional), ALLOCATABLE :: functional
LOGICAL :: initialized = .FALSE.
INTEGER :: nfp = 0

CONTAINS

!================================================================================================================================!
SUBROUTINE Init(parameterfile)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut
  USE MODgvec_MHD3D_vars,     ONLY: X1_base
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
  CALL InitOutput()
  CALL InitAnalyze()

  ! initialize the functional
  which_functional = GETINT('which_functional',Proposal=1)
  CALL InitFunctional(functional,which_functional)

  ! print the ignored parameters
  CALL IgnoredStrings()

  ! additional global variables
  nfp = X1_base%f%nfp
  initialized = .TRUE.
END SUBROUTINE Init

!================================================================================================================================!
SUBROUTINE InitSolution()
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL functional%InitSolution()
END SUBROUTINE InitSolution

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
  outputLevel=outputLevel_r
END SUBROUTINE ReadState

!================================================================================================================================!
!> Handle the selection of the base, based on the selection string
!================================================================================================================================!
SUBROUTINE select_base_dofs(var, base, dofs)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut,abort
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,U
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var          !! selection string: which variable to evaluate
  CLASS(t_base), POINTER, INTENT(OUT) :: base  !! pointer to the base object (X1, X2, LA)
  REAL, POINTER, INTENT(OUT) :: dofs(:,:)      !! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(var)
    CASE('X1')
      base => X1_base
      dofs => U(0)%X1
    CASE('X2')
      base => X2_base
      dofs => U(0)%X2
    CASE('LA')
      base => LA_base
      dofs => U(0)%LA
    CASE DEFAULT
      CALL abort(__STAMP__, &
      'ERROR: variable "'//TRIM(var)//'" not recognized')
  END SELECT
END SUBROUTINE select_base_dofs

!================================================================================================================================!
!> Handle the selection of the base, based on the selection string
!================================================================================================================================!
SUBROUTINE select_base(var, base)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut,abort
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var                 !! selection string: which variable to evaluate
  CLASS(t_base), POINTER, INTENT(OUT) :: base         !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(var)
    CASE('X1')
      base => X1_base
    CASE('X2')
      base => X2_base
    CASE('LA')
      base => LA_base
    CASE DEFAULT
      CALL abort(__STAMP__, &
      'ERROR: variable "'//TRIM(var)//'" not recognized')
  END SELECT
END SUBROUTINE select_base

!================================================================================================================================!
!> Handle the selection of the functional and derivatives, based on the selection strings
!================================================================================================================================!
SUBROUTINE evaluate_base_select(var, sel_deriv_s, sel_deriv_f, base, solution_dofs, seli_deriv_s, seli_deriv_f)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var                 !! selection string: which variable to evaluate
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_s         !! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_f         !! selection string: which derivative to evaluate for the fourier series
  CLASS(t_base), POINTER, INTENT(OUT) :: base         !! pointer to the base object (X1, X2, LA)
  REAL, POINTER, INTENT(OUT) :: solution_dofs(:,:)    !! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  INTEGER, INTENT(OUT) :: seli_deriv_s, seli_deriv_f  !! integer values for the derivative selection
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base_dofs(var, base, solution_dofs)
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
SUBROUTINE get_integration_points_num(var, n_s, n_t, n_z)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut,abort
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var           !! selection string: which variable to evaluate
  INTEGER, INTENT(OUT) :: n_s, n_t, n_z         !! number of integration points
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base                !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  n_s = base%s%nGP
  IF (base%s%nGP /= SIZE(base%s%s_GP)) THEN
    CALL abort(__STAMP__, &
    'ERROR: number of integration points does not match the size of the array')
  END IF
  n_t = base%f%mn_nyq(1)
  n_z = base%f%mn_nyq(2)
END SUBROUTINE get_integration_points_num

!================================================================================================================================!
!> Retrieve the integration points and weights (gauss points for radial integration)
!================================================================================================================================!
SUBROUTINE get_integration_points(var, s_GP, s_w, t_w, z_w)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var           !! selection string: which variable to evaluate
  REAL, DIMENSION(:), INTENT(OUT) :: s_GP, s_w  !! output arrays for the gauss points and weights
  REAL, INTENT(OUT) :: t_w, z_w                 !! output array for the fourier interpolation weights (equidistant points)
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base                !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  s_GP = base%s%s_GP
  s_w = base%s%w_GP
  t_w = base%f%d_thet
  z_w = base%f%d_zeta
END SUBROUTINE get_integration_points

!================================================================================================================================!
SUBROUTINE get_modes(var, modes)
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var   !! selection string: which variable to evaluate
  INTEGER, INTENT(OUT) :: modes         !! total number of modes in basis (depends if only sin/cos or sin & cos are used)
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base        !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  modes = base%f%modes
END SUBROUTINE get_modes

!================================================================================================================================!
SUBROUTINE get_mn_max(var, m_max, n_max)
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var   !! selection string: which variable to evaluate
  INTEGER, INTENT(OUT) :: m_max, n_max  !! maximum number of poloidal, toroidal modes
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base        !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  m_max = base%f%mn_max(1)
  n_max = base%f%mn_max(2)
END SUBROUTINE get_mn_max

!================================================================================================================================!
SUBROUTINE get_mn_IP(var, mn_IP)
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var   !! selection string: which variable to evaluate
  INTEGER, INTENT(OUT) :: mn_IP         !! =mn_nyq(1)*mn_nyq(2)
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base        !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  mn_IP = base%f%mn_IP
END SUBROUTINE get_mn_IP

!================================================================================================================================!
SUBROUTINE get_s_nBase(var, s_nbase)
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var   !! selection string: which variable to evaluate
  INTEGER, INTENT(OUT) :: s_nbase       !! total number of degree of freedom / global basis functions
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base        !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  s_nbase = base%s%nbase
END SUBROUTINE get_s_nBase

!================================================================================================================================!
SUBROUTINE get_s_IP(var, s_IP)
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=2), INTENT(IN) :: var       !! selection string: which variable to evaluate
  REAL, DIMENSION(:), INTENT(OUT) :: s_IP   !! position of interpolation points for initialization, size(nBase)
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  CLASS(t_base), POINTER :: base            !! pointer to the base object (X1, X2, LA)
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base(var, base)
  s_IP = base%s%s_IP
END SUBROUTINE get_s_IP

!================================================================================================================================!
!> Evaluate the basis for a list of (theta, zeta) positions on all flux surfaces given by s
!================================================================================================================================!`
SUBROUTINE evaluate_base_list_tz(n_s, n_tz, s, thetazeta, var, sel_deriv_s, sel_deriv_f, result)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s, n_tz              !! number of evaluation points
  REAL, INTENT(IN) :: s(n_s), thetazeta(2,n_tz) !! evaluation points
  CHARACTER(LEN=2), INTENT(IN) :: var           !! selection string: which variable to evaluate
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_s   !! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2), INTENT(IN) :: sel_deriv_f   !! selection string: which derivative to evaluate for the fourier series
  REAL, INTENT(OUT) :: result(n_s,n_tz)         !! output array
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i_s                          ! loop variables
  INTEGER :: seli_deriv_s, seli_deriv_f   ! integer values for the derivative selection
  CLASS(t_base), POINTER :: base          ! pointer to the base object (X1, X2, LA)
  REAL, POINTER :: solution_dofs(:,:)     ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  REAL, ALLOCATABLE :: fourier_dofs(:)    ! DOFs for the fourier series, calculated from the spline
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
!> Evaluate the basis and all derivatives for a list of (theta, zeta) positions on all flux surfaces given by s
!================================================================================================================================!
SUBROUTINE evaluate_base_list_tz_all(n_s, n_tz, s, thetazeta, Qsel, Q, dQ_ds, dQ_dthet, dQ_dzeta, &
                                     dQ_dss, dQ_dst, dQ_dsz, dQ_dtt, dQ_dtz, dQ_dzz)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s, n_tz                !! number of evaluation points
  REAL, INTENT(IN) :: s(n_s), thetazeta(2,n_tz)   !! evaluation points
  CHARACTER(LEN=2), INTENT(IN) :: Qsel            !! selection string: which variable to evaluate
  REAL, INTENT(OUT), DIMENSION(n_s,n_tz) :: Q, &  !! reference space position and derivatives
    dQ_ds, dQ_dthet, dQ_dzeta, dQ_dss, dQ_dst, dQ_dsz, dQ_dtt, dQ_dtz, dQ_dzz
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i                                                        ! loop variables
  INTEGER :: seli_deriv_s, seli_deriv_f                               ! integer values for the derivative selection
  CLASS(t_base), POINTER :: base                                      ! pointer to the base object (X1, X2, LA)
  REAL, POINTER :: solution_dofs(:,:)                                 ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  REAL, ALLOCATABLE, DIMENSION(:) :: Q_dofs, dQ_ds_dofs, dQ_dss_dofs  ! DOFs for the fourier series
  REAL, ALLOCATABLE :: intermediate(:)                                ! intermediate result array before reshaping
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base_dofs(Qsel, base, solution_dofs)
  DO i=1,n_s
    ! evaluate spline to get the fourier dofs
    Q_dofs = base%s%evalDOF2D_s(s(i), base%f%modes, 0, solution_dofs(:,:))
    dQ_ds_dofs = base%s%evalDOF2D_s(s(i), base%f%modes, DERIV_S, solution_dofs(:,:))
    dQ_dss_dofs = base%s%evalDOF2D_s(s(i), base%f%modes, DERIV_S_S, solution_dofs(:,:))
    ! use the tensorproduct for theta and zeta
    Q(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, 0, Q_dofs)
    dQ_ds(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, 0, dQ_ds_dofs)
    dQ_dthet(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_THET, Q_dofs)
    dQ_dzeta(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_ZETA, Q_dofs)
    dQ_dss(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, 0, dQ_dss_dofs)
    dQ_dst(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_THET, dQ_ds_dofs)
    dQ_dsz(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_ZETA, dQ_ds_dofs)
    dQ_dtt(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_THET_THET, Q_dofs)
    dQ_dtz(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_THET_ZETA, Q_dofs)
    dQ_dzz(i,:) = base%f%evalDOF_xn(n_tz, thetazeta, DERIV_ZETA_ZETA, Q_dofs)
  END DO
  DEALLOCATE(Q_dofs, dQ_ds_dofs, dQ_dss_dofs)
END SUBROUTINE evaluate_base_list_tz_all

!================================================================================================================================!
!> Evaluate the basis with a tensorproduct for the given 1D (s, theta, zeta) values
!================================================================================================================================!
SUBROUTINE evaluate_base_tens(s, theta, zeta, var, sel_deriv_s, sel_deriv_f, result)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  REAL, INTENT(IN) :: s(:), theta(:), zeta(:)   !! evaluation points to construct a mesh
  CHARACTER(LEN=2) :: var                       !! selection string: which variable to evaluate
  CHARACTER(LEN=2) :: sel_deriv_s               !! selection string: which derivative to evaluate for the spline
  CHARACTER(LEN=2) :: sel_deriv_f               !! selection string: which derivative to evaluate for the fourier series
  REAL, INTENT(OUT) :: result(:,:,:)            !! output array
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
!> Evaluate the basis with a tensorproduct for the given 1D (s, theta, zeta) values
!================================================================================================================================!
SUBROUTINE evaluate_base_tens_all(n_s, n_t, n_z, s, theta, zeta, Qsel, Q, dQ_ds, dQ_dthet, dQ_dzeta, &
                                  dQ_dss, dQ_dst, dQ_dsz, dQ_dtt, dQ_dtz, dQ_dzz)
  ! MODULES
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s, n_t, n_z                                            !! number of evaluation points
  REAL, INTENT(IN) :: s(n_s), theta(n_t), zeta(n_z)                               !! evaluation points to construct a mesh
  CHARACTER(LEN=2) :: Qsel                                                        !! selection string: which variable to evaluate
  REAL, INTENT(OUT), DIMENSION(n_s,n_t,n_z) :: Q, dQ_ds, dQ_dthet, dQ_dzeta       !! reference space position & derivatives
  REAL, INTENT(OUT), DIMENSION(n_s,n_t,n_z) :: dQ_dss, dQ_dst, dQ_dsz, dQ_dtt, dQ_dtz, dQ_dzz
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i                                                         ! loop variables
  CLASS(t_base), POINTER :: base                                       ! pointer to the base object (X1, X2, LA)
  REAL, POINTER :: solution_dofs(:,:)                                  ! pointer to the solution dofs (U(0)%X1, U(0)%X2, U(0)%LA)
  REAL, ALLOCATABLE, DIMENSION(:) :: Q_dofs, dQ_ds_dofs, dQ_dss_dofs   ! DOFs for the fourier series
  REAL, ALLOCATABLE :: intermediate(:)                                 ! intermediate result array before reshaping
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL select_base_dofs(Qsel, base, solution_dofs)
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
  DEALLOCATE(Q_dofs, dQ_ds_dofs, dQ_dss_dofs)
END SUBROUTINE evaluate_base_tens_all

!================================================================================================================================!
!> Evaluate the mapping from reference to physical space (hmap)
!================================================================================================================================!
SUBROUTINE evaluate_hmap(n, X1, X2, zeta, dX1_ds, dX2_ds, dX1_dthet, dX2_dthet, dX1_dzeta, dX2_dzeta, coord, e_s, e_thet, e_zeta)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                                      !! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta, dX1_ds, dX2_ds                !! reference space position & derivatives
  REAL, INTENT(IN), DIMENSION(n) :: dX1_dthet, dX2_dthet, dX1_dzeta, dX2_dzeta  !! reference space derivatives
  REAL, INTENT(OUT), DIMENSION(3,n) :: coord, e_s, e_thet, e_zeta               !! real space position and basis vectors
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
!> Evaluate the mapping from reference to physical space (hmap) without logical coordinates
!================================================================================================================================!
SUBROUTINE evaluate_hmap_only(n, X1, X2, zeta, pos, e_X1, e_X2, e_zeta3)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                      !! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta                !! reference space position
  REAL, INTENT(OUT), DIMENSION(3,n) :: pos, e_X1, e_X2, e_zeta3 !! real space position and reference tangent basis vectors
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
!> evaluate components of the metric tensor and their derivatives
!================================================================================================================================!
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
  INTEGER, INTENT(IN) :: n                                                                        !! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz  !! reference coordinates
  REAL, INTENT(IN), DIMENSION(n) :: dX1_dss, dX2_dss, dX1_dst, dX2_dst, dX1_dsz, dX2_dsz          !! and their derivatives
  REAL, INTENT(IN), DIMENSION(n) :: dX1_dtt, dX2_dtt, dX1_dtz, dX2_dtz, dX1_dzz, dX2_dzz
  REAL, INTENT(OUT), DIMENSION(n) :: g_ss, g_st, g_sz, g_tt, g_tz, g_zz                           !! metric coefficients
  REAL, INTENT(OUT), DIMENSION(n) :: dg_ss_ds, dg_st_ds, dg_sz_ds, dg_tt_ds, dg_tz_ds, dg_zz_ds   !! derivatives of the m. coef.
  REAL, INTENT(OUT), DIMENSION(n) :: dg_ss_dt, dg_st_dt, dg_sz_dt, dg_tt_dt, dg_tz_dt, dg_zz_dt   !! derivatives of the m. coef.
  REAL, INTENT(OUT), DIMENSION(n) :: dg_ss_dz, dg_st_dz, dg_sz_dz, dg_tt_dz, dg_tz_dz, dg_zz_dz   !! derivatives of the m. coef.
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
!> evaluate the jacobian determinant and its derivatives
!================================================================================================================================!
SUBROUTINE evaluate_jacobian(n, X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz, Jh, dJh_ds, dJh_dt, dJh_dz)
  ! MODULES
  USE MODgvec_MHD3D_vars,     ONLY: hmap
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n                                                                        !! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n) :: X1, X2, zeta, dX1_ds, dX2_ds, dX1_dt, dX2_dt, dX1_dz, dX2_dz  !! reference coordinates
  REAL, INTENT(OUT), DIMENSION(n) :: Jh, dJh_ds, dJh_dt, dJh_dz                                   !! jacobian det. and derivatives
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
!> evaluate iota/pressure profile and its derivatives with respect to rho2=rho^2
!================================================================================================================================!
SUBROUTINE evaluate_rho2_profile(n_s, rho2, deriv, var, result)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut,abort
  USE MODgvec_rProfile_base,  ONLY: c_rProfile
  USE MODgvec_MHD3D_Vars,     ONLY: iota_profile, pres_profile, chi_profile, Phi_profile
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s                  !! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n_s) :: rho2    !! radial evaluation points (in rho^2)
  INTEGER, INTENT(IN) :: deriv                !! order of the derivative in rho^2
  CHARACTER(LEN=*), INTENT(IN) :: var         !! selection string: which profile to evaluate
  REAL, INTENT(OUT), DIMENSION(n_s) :: result !! values of the profile
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i  ! loop variable
  CLASS(c_rProfile), ALLOCATABLE :: input_profile
  SELECT CASE(TRIM(var))
  CASE("iota")
    input_profile = iota_profile
  CASE("p")
    input_profile = pres_profile
  CASE("chi")
    input_profile = chi_profile
  CASE("Phi")
    input_profile = Phi_profile
  CASE DEFAULT
    CALL abort(__STAMP__, &
    'ERROR: variable "'//TRIM(var)//'" not recognized')
  END SELECT
  DO i = 1,n_s
    result(i) = input_profile%eval_at_rho2(rho2(i), deriv)
  END DO
  SDEALLOCATE(input_profile)
END SUBROUTINE evaluate_rho2_profile


!================================================================================================================================!
!> evaluate 1D-radial profiles and their derivatives with respect to rho
!================================================================================================================================!
SUBROUTINE evaluate_profile(n_s, s, deriv, var, result)
  ! MODULES
  USE MODgvec_Globals,        ONLY: Unit_stdOut,abort
  USE MODgvec_rProfile_base,  ONLY: c_rProfile
  USE MODgvec_MHD3D_Vars,     ONLY: which_init, init_with_profile_iota, init_with_profile_pressure, iota_profile, pres_profile, &
                                    Phi_profile, chi_profile
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: n_s                  !! number of evaluation points
  REAL, INTENT(IN), DIMENSION(n_s) :: s       !! radial evaluation points
  INTEGER, INTENT(IN) :: deriv                !! order of the derivative in rho
  CHARACTER(LEN=*), INTENT(IN) :: var         !! selection string: which profile to evaluate
  REAL, INTENT(OUT), DIMENSION(n_s) :: result !! values of the profile
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i  ! loop variable
  CLASS(c_rProfile), ALLOCATABLE :: input_profile
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(TRIM(var))
    CASE("iota")
      input_profile = iota_profile
    CASE("p")
      input_profile = pres_profile
    CASE("chi")
      input_profile = chi_profile
    CASE("Phi")
      input_profile = Phi_profile
    CASE DEFAULT
      CALL abort(__STAMP__, &
      'ERROR: variable "'//TRIM(var)//'" not recognized')
  END SELECT

  DO i = 1,n_s
    result(i) = input_profile%eval_at_rho(s(i),deriv)
  END DO
  SDEALLOCATE(input_profile)
END SUBROUTINE evaluate_profile

!================================================================================================================================!
!> initialize a SFL-Boozer object, with some parameters taken from the state (globals)
!> Note: as of v0.2.16 (Thanks Christopher Albert) f90wrap supports ALLOCATABLEs in the return value
!================================================================================================================================!
FUNCTION init_boozer(mn_max, mn_nyq, sin_cos, nrho, rho_pos, relambda) RESULT(sfl_boozer)
  ! MODULES
  USE MODgvec_SFL_Boozer,   ONLY: t_sfl_boozer, sfl_boozer_new
  USE MODgvec_MHD3D_vars,     ONLY: hmap, iota_profile, Phi_profile
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  TYPE(t_sfl_boozer), ALLOCATABLE :: sfl_boozer                 !! SFL-Boozer object
  INTEGER, INTENT(IN) :: mn_max(2), mn_nyq(2), nrho             !! parameters for the Boozer object
  CHARACTER(LEN=8), INTENT(IN)   :: sin_cos      !! can be either only sine: " _sin_" only cosine: " _cos_" or full: "_sin_cos_"
  REAL, INTENT(IN), DIMENSION(nrho) :: rho_pos                  !! radial positions
  LOGICAL, INTENT(IN) :: relambda                               !! recompute lambda flag
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i
  REAL, DIMENSION(nrho) :: iota, phiPrime                       ! iota and phiPrime
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  DO i = 1,nrho
    iota(i) = iota_profile%eval_at_rho(rho_pos(i))
    phiPrime(i) = Phi_profile%eval_at_rho(rho_pos(i),deriv=1)
  END DO
  ! ALLOCATE is called within sfl_boozer_new
  CALL sfl_boozer_new(sfl_boozer, mn_max, mn_nyq, nfp, sin_cos, hmap, nrho, rho_pos, iota, phiPrime, relambda)
END FUNCTION init_boozer

SUBROUTINE get_boozer(sfl_boozer)
  ! MODULES
  USE MODgvec_SFL_Boozer,   ONLY: t_sfl_boozer
  USE MODgvec_MHD3D_vars,     ONLY: X1_base, X2_base, LA_base, U
  USE MODgvec_base,           ONLY: t_base
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  TYPE(t_sfl_boozer), INTENT(INOUT) :: sfl_boozer  ! SFL-Boozer object
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  CALL sfl_boozer%get_boozer(X1_base, X2_base, LA_base, U(0)%X1, U(0)%X2, U(0)%LA)
END SUBROUTINE get_boozer

!================================================================================================================================!
!> Evaluate LA or NU and all derivatives for a list of (theta, zeta) positions on all flux surfaces given by s
!================================================================================================================================!
SUBROUTINE evaluate_boozer_list_tz_all(sfl_boozer, n_s, n_tz, irho, thetazeta, Qsel, Q, dQ_dthet, dQ_dzeta, &
                                       dQ_dtt, dQ_dtz, dQ_dzz)
  ! MODULES
  USE MODgvec_Globals,      ONLY: Unit_stdOut,abort
  USE MODgvec_fbase,        ONLY: t_fbase
  USE MODgvec_SFL_Boozer,   ONLY: t_sfl_boozer
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  TYPE(t_sfl_boozer), INTENT(IN), TARGET :: sfl_boozer  ! SFL-Boozer object
  INTEGER, INTENT(IN) :: n_s, n_tz                      ! number of evaluation points
  INTEGER, INTENT(IN) :: irho(1:n_s)                      ! indices of the flux surfaces
  REAL, INTENT(IN) :: thetazeta(2,n_tz)                 ! evaluation points
  CHARACTER(LEN=2), INTENT(IN) :: Qsel                  ! selection string: which variable to evaluate
  REAL, INTENT(OUT), DIMENSION(1:n_s,n_tz) :: Q, &        ! reference space position and derivatives
    dQ_dthet, dQ_dzeta, dQ_dtt, dQ_dtz, dQ_dzz
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: i                                                        ! loop variables
  INTEGER :: seli_deriv_s, seli_deriv_f                               ! integer values for the derivative selection
  CLASS(t_fbase), POINTER :: base                                     ! pointer to the base object (LA, NU)
  REAL, POINTER :: dofs(:,:)                                          ! pointer to the solution dofs
  REAL, ALLOCATABLE, DIMENSION(:) :: Q_dofs                           ! DOFs for the fourier series
  REAL, ALLOCATABLE :: intermediate(:)                                ! intermediate result array before reshaping
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  SELECT CASE(Qsel)
    CASE('LA')
      base => sfl_boozer%nu_fbase
      dofs => sfl_boozer%lambda
    CASE('NU')
      base => sfl_boozer%nu_fbase
      dofs => sfl_boozer%nu
    CASE DEFAULT
      WRITE(UNIT_stdout,*) 'ERROR: variable', Qsel, 'not recognized'
      CALL abort(__STAMP__,"")
  END SELECT
  DO i=1,n_s
    ! there is no radial basis in the Boozer object, so we only select the correct radial position
    Q_dofs = dofs(:, irho(i)+1)
    ! use the tensorproduct for theta and zeta
    Q(i,:) = base%evalDOF_xn(n_tz, thetazeta, 0, Q_dofs)
    dQ_dthet(i,:) = base%evalDOF_xn(n_tz, thetazeta, DERIV_THET, Q_dofs)
    dQ_dzeta(i,:) = base%evalDOF_xn(n_tz, thetazeta, DERIV_ZETA, Q_dofs)
    dQ_dtt(i,:) = base%evalDOF_xn(n_tz, thetazeta, DERIV_THET_THET, Q_dofs)
    dQ_dtz(i,:) = base%evalDOF_xn(n_tz, thetazeta, DERIV_THET_ZETA, Q_dofs)
    dQ_dzz(i,:) = base%evalDOF_xn(n_tz, thetazeta, DERIV_ZETA_ZETA, Q_dofs)
  END DO
  DEALLOCATE(Q_dofs)
END SUBROUTINE evaluate_boozer_list_tz_all

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
  FLUSH(Unit_stdOut)
END SUBROUTINE Finalize

END MODULE MODgvec_py_post
