!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_LBFGS
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_LBFGS
    ! USE WITH CARE HARDLY TESTED!
    ! Limited memory BFGS (L-BFGS) does not store a full Hessian but rather just some defined number of updates

    INTEGER :: m_updates ! number of saved updates (should be lower than n_dofs)
    REAL(wp), ALLOCATABLE :: y_n(:,:) ! gradient steps => g_{n+1} - g_{n}; the first entry in the second dim is the last update step
    REAL(wp), ALLOCATABLE :: s_n(:,:) ! dof steps => dofs_{n+1} - dofs_{n}; the first entry in the second dim is the last update step
    REAL(wp), ALLOCATABLE :: rho_n(:) ! y_n^T s_n
    REAL(wp), ALLOCATABLE :: alpha(:) ! last m stepsizes
    REAL(wp), ALLOCATABLE :: gradf_old(:) ! gradient of step k

    CONTAINS
    PROCEDURE :: get_update => LBFGS_get_update
    PROCEDURE :: reset      => LBFGS_reset
    FINAL     :: LBFGS_free

  END TYPE t_free_b_minimizer_LBFGS

  INTERFACE t_free_b_minimizer_LBFGS
      MODULE PROCEDURE LBFGS_new
  END INTERFACE t_free_b_minimizer_LBFGS

  CONTAINS

  FUNCTION LBFGS_new(n_dofs, dt, max_dt, m_updates) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)  :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN) :: dt ! initial / default time step
    REAL(wp), INTENT(IN) :: max_dt
    INTEGER, INTENT(IN)  :: m_updates  ! number of remembered tau_hat

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_LBFGS) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    sf%dt_initial = dt

    sf%m_updates = m_updates
    sf%dt = sf%dt_initial
    sf%max_dt = max_dt

    ALLOCATE(sf%y_n(sf%n_dofs,sf%m_updates))
    ALLOCATE(sf%s_n(sf%n_dofs,sf%m_updates))
    ALLOCATE(sf%rho_n(sf%m_updates))
    ALLOCATE(sf%alpha(sf%m_updates))
    ALLOCATE(sf%gradf_old(1:sf%n_dofs))

    sf%gradf_old = 0.0_wp
    sf%y_n = 0.0_wp
    sf%s_n = 0.0_wp
    sf%rho_n = 0.0_wp
    sf%n_updates = 0
  END FUNCTION LBFGS_new

  SUBROUTINE LBFGS_reset(sf)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_LBFGS), INTENT(INOUT) :: sf
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%y_n = 0.0_wp
    sf%s_n = 0.0_wp
    sf%rho_n = 0.0_wp
    sf%n_updates = 0
  END SUBROUTINE LBFGS_reset

  FUNCTION LBFGS_get_update(sf, force_in) RESULT(push_out)
    USE MODgvec_MHD3D_Vars, ONLY: Precond_b_X1, Precond_b_X2, X1_base, X2_base
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_LBFGS), INTENT(INOUT) :: sf
    REAL(wp)                       , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES
    INTEGER :: i, m, n_saved_updates
    REAL(wp):: q(sf%n_dofs), H0(sf%n_dofs,sf%n_dofs), alpha_i(sf%m_updates), beta, r(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    q = force_in
    alpha_i = 0.0_wp
    ! identity matrix
    H0 = 0.0_wp
    r = 0.0_wp
    n_saved_updates = MIN(sf%m_updates,sf%n_updates)

    IF (sf%n_updates == 0) THEN
      ! use diag. preconditioner as initial guess for the hessian
      DO i=1,X1_base%f%modes
        H0(i,i) = 1 / Precond_b_X1(i)
      END DO

      DO i=1,X2_base%f%modes
        H0(X1_base%f%modes+i,X1_base%f%modes+i) = 1 / Precond_b_X2(i)
      END DO

    ELSE
      DO i=1,sf%n_dofs
        H0(i,i) = 1.0_wp
      END DO
      H0 = DOT_PRODUCT(sf%s_n(:,1),sf%y_n(:,1)) / DOT_PRODUCT(sf%y_n(:,1),sf%y_n(:,1))*H0
    END IF

    DO i=1,n_saved_updates
      alpha_i(i) = sf%rho_n(i)*DOT_PRODUCT(sf%s_n(:,i),q)
      q = q - alpha_i(i)*sf%y_n(:,i)
    END DO

    r = MATMUL(H0,q)

    DO m=1,n_saved_updates
      i = n_saved_updates + 1 - m
      beta = sf%rho_n(i)*DOT_PRODUCT(sf%y_n(:,i),r)
      r = r + sf%s_n(:,i)*(alpha_i(i) - beta)
    END DO

    ! -H_k (nabla*f_k)
    push_out = -sf%dt*r
    sf%n_updates = sf%n_updates + 1

    ! update history
    sf%s_n(:,2:) = sf%s_n(:,1:sf%m_updates-1)
    sf%y_n(:,2:) = sf%y_n(:,1:sf%m_updates-1)
    sf%rho_n(2:) = sf%rho_n(1:sf%m_updates-1)

    sf%s_n(:,1) = r
    sf%y_n(:,1) = force_in - sf%gradf_old
    sf%rho_n(1) = 1 / DOT_PRODUCT(sf%y_n(:,1),sf%s_n(:,1))
  END FUNCTION LBFGS_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE LBFGS_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_LBFGS), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%y_n)
    SDEALLOCATE(sf%s_n)
    SDEALLOCATE(sf%rho_n)
    SDEALLOCATE(sf%alpha)
    SDEALLOCATE(sf%gradf_old)
  END SUBROUTINE LBFGS_free

END MODULE MODgvec_free_b_minimizer_LBFGS
