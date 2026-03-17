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
    INTEGER :: m_updates ! number of saved updates (should be lower than n_dofs)
    REAL(wp), ALLOCATABLE :: y_n(:,:) ! gradient steps => g_{n+1} - g_{n}; the first entry in the second dim is the last update step
    REAL(wp), ALLOCATABLE :: s_n(:,:) ! dof steps => dofs_{n+1} - dofs_{n}; the first entry in the second dim is the last update step
    REAL(wp), ALLOCATABLE :: rho_n(:) ! y_n^T s_n
    REAL(wp), ALLOCATABLE :: alpha(:) ! last m stepsizes

    CONTAINS
    PROCEDURE :: get_update => LBFGS_get_update
    PROCEDURE :: reset      => LBFGS_reset
    FINAL     :: LBFGS_free

  END TYPE t_free_b_minimizer_LBFGS

  INTERFACE t_free_b_minimizer_LBFGS
      MODULE PROCEDURE LBFGS_new
  END INTERFACE t_free_b_minimizer_LBFGS

  CONTAINS

  FUNCTION LBFGS_new(n_dofs, dt, m_updates) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)            :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN), OPTIONAL :: dt ! initial / default time step
    INTEGER, INTENT(IN) , OPTIONAL :: m_updates  ! number of remembered tau_hat

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_LBFGS) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    IF (PRESENT(dt)) THEN
        sf%dt_initial = dt
    ELSE
        sf%dt_initial = 1e-3_wp
    END IF

    IF (PRESENT(m_updates)) THEN
        sf%m_updates = m_updates
    ELSE
        sf%m_updates = 10
    END IF

    sf%dt =  sf%dt_initial

    ALLOCATE(sf%y_n(sf%n_dofs,sf%m_updates))
    ALLOCATE(sf%s_n(sf%n_dofs,sf%m_updates))
    ALLOCATE(sf%rho_n(sf%m_updates))
    ALLOCATE(sf%alpha(sf%m_updates))

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

  SUBROUTINE set_sk_yk(sf, s_k, y_k)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_LBFGS), INTENT(INOUT) :: sf
    REAL(wp), INTENT(IN) :: s_k(sf%n_dofs), y_k(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%s_n(:,2:) = sf%s_n(:,1:sf%m_updates-1)
    sf%y_n(:,2:) = sf%y_n(:,1:sf%m_updates-1)

    sf%s_n(:,1) = s_k
    sf%y_n(:,1) = y_k
  END SUBROUTINE set_sk_yk

  FUNCTION LBFGS_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_LBFGS), INTENT(INOUT) :: sf
    REAL(wp)                       , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES
    INTEGER :: i, m
    REAL(wp):: q(sf%n_dofs), H0(sf%n_dofs,sf%n_dofs), alpha_i(sf%m_updates), beta
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    q = force_in
    alpha_i = 0.0_wp
    ! identity matrix
    H0 = 0.0_wp
    DO i=1,sf%n_dofs
        H0(i,i) = 1.0_wp
    END DO

    IF(sf%n_updates>1) THEN
        H0 = DOT_PRODUCT(sf%s_n(:,1),sf%y_n(:,1)) / DOT_PRODUCT(sf%y_n(:,1),sf%y_n(:,1))
        ! first loop
        DO i=1,MIN(sf%m_updates,sf%n_updates)
            alpha_i(i) = sf%rho_n(i)*DOT_PRODUCT(sf%s_n(:,i),q)
            q = q - alpha_i(i)*sf%y_n(:,i)
        END DO
        push_out = MATMUL(H0,q)
        DO i=1,MIN(sf%m_updates,sf%n_updates)
            beta = sf%rho_n(i)*DOT_PRODUCT(sf%y_n(:,i),push_out)
            push_out = push_out + sf%s_n(:,i)*(alpha_i(i) - beta)
        END DO
        ! -H_k (nabla*f_k)
        push_out = -sf%dt*push_out
        sf%n_updates = sf%n_updates + 1
    ELSE
        ! -(nabla*f_k)
        push_out = -sf%dt*q
        sf%n_updates = sf%n_updates + 1
    END IF
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
    SDEALLOCATE(sf%y_n)
    SDEALLOCATE(sf%s_n)
    SDEALLOCATE(sf%rho_n)
    SDEALLOCATE(sf%alpha)
  END SUBROUTINE LBFGS_free

END MODULE MODgvec_free_b_minimizer_LBFGS
