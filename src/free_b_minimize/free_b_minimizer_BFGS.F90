!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_BFGS
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_BFGS

    REAL(wp), ALLOCATABLE :: s_k(:) ! last update
    REAL(wp), ALLOCATABLE :: y_k(:) ! difference of gradients
    REAL(wp), ALLOCATABLE :: p_k (:) ! search direction => -H_k grad f(x_k)
    REAL(wp), ALLOCATABLE :: gradf_old(:) ! gradient of step k
    REAL(wp), ALLOCATABLE :: H_k(:,:) ! approximate inverse Hessian
    REAL(wp), ALLOCATABLE :: I(:,:)
    INTEGER  :: n_grad_steps

    CONTAINS
    PROCEDURE :: get_update => BFGS_get_update
    PROCEDURE :: reset      => BFGS_reset
    FINAL     :: BFGS_free

  END TYPE t_free_b_minimizer_BFGS

  INTERFACE t_free_b_minimizer_BFGS
      MODULE PROCEDURE BFGS_new
  END INTERFACE t_free_b_minimizer_BFGS

  CONTAINS

  FUNCTION BFGS_new(n_dofs, dt, n_grad_steps, max_dt, mn_X1X2) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)            :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN), OPTIONAL :: dt ! initial / default time step
    REAL(wp), INTENT(IN), OPTIONAL :: max_dt
    INTEGER, INTENT(IN), OPTIONAL  :: n_grad_steps
    REAL(wp), INTENT(IN) :: mn_X1X2(2,n_dofs)

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_BFGS) :: sf
    ! LOCAL VARIABLES -------------------------!
    INTEGER :: i
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    IF (PRESENT(n_grad_steps)) THEN
        sf%n_grad_steps = n_grad_steps
    ELSE
        sf%n_grad_steps = 0
    END IF

    IF (PRESENT(dt)) THEN
        sf%dt_initial = dt
    ELSE
        sf%dt_initial = 1e-3_wp
    END IF

    IF (PRESENT(max_dt)) THEN
        sf%max_dt = max_dt
    ELSE
        sf%max_dt = 0.1_wp
    END IF

    ALLOCATE(sf%s_k(1:sf%n_dofs))
    ALLOCATE(sf%y_k(1:sf%n_dofs))
    ALLOCATE(sf%p_k(1:sf%n_dofs))
    ALLOCATE(sf%gradf_old(1:sf%n_dofs))
    ALLOCATE(sf%H_k(1:sf%n_dofs, 1:sf%n_dofs))
    ALLOCATE(sf%I(1:sf%n_dofs, 1:sf%n_dofs))

    sf%s_k = 0.0_wp
    sf%y_k = 0.0_wp
    sf%p_k = 0.0_wp
    sf%gradf_old = 0.0_wp
    sf%H_k = 0.0_wp
    sf%I = 0.0_wp
    DO i=1,sf%n_dofs
        sf%H_k(i,i) = 1.0_wp
        sf%I(i,i) = 1.0_wp
    END DO

    sf%n_updates = 0
  END FUNCTION BFGS_new

  SUBROUTINE BFGS_reset(sf)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_BFGS), INTENT(INOUT) :: sf
    ! LOCAL VARIABLES
    INTEGER :: i
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%s_k = 0.0_wp
    sf%y_k = 0.0_wp
    sf%p_k = 0.0_wp
    sf%gradf_old = 0.0_wp
    sf%H_k = 0.0_wp
    DO i=1,sf%n_dofs
        sf%H_k(i,i) = 1.0_wp
    END DO
    sf%n_updates = 0
  END SUBROUTINE BFGS_reset

  FUNCTION OUTER_PRODUCT(A,B) result(AB)
    REAL(wp), INTENT(IN) :: A(:),B(:)
    REAL(wp), ALLOCATABLE :: AB(:,:)
    INTEGER :: nA,nB
    nA=SIZE(A)
    nB=SIZE(B)
    ALLOCATE(AB(nA,nB))
    AB = SPREAD(source = A, dim = 2, ncopies = nB) * &
        SPREAD(source = B, dim = 1, ncopies = nA)
  END FUNCTION OUTER_PRODUCT

  FUNCTION BFGS_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_BFGS), INTENT(INOUT) :: sf
    REAL(wp)                      , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES
    REAL(wp) :: dot_sy, norm_p_k, p_k(sf%n_dofs), g_k(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    g_k = force_in

    IF (sf%n_updates.LE.sf%n_grad_steps) THEN
        p_k = -g_k
    ELSE
        sf%y_k = g_k - sf%gradf_old
        dot_sy = DOT_PRODUCT(sf%s_k,sf%y_k)
        sf%H_k = OUTER_PRODUCT(sf%s_k,sf%s_k) / dot_sy + MATMUL(sf%I-OUTER_PRODUCT(sf%s_k,sf%y_k) / dot_sy,MATMUL(sf%H_k, sf%I-OUTER_PRODUCT(sf%y_k,sf%s_k) / dot_sy))
        p_k = -MATMUL(sf%H_k,g_k)
    END IF
    sf%gradf_old = g_k

    push_out = sf%dt*p_k
    sf%s_k = push_out
    sf%n_updates = sf%n_updates+1
    push_out =  push_out
  END FUNCTION BFGS_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE BFGS_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_BFGS), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%s_k)
    SDEALLOCATE(sf%y_k)
    SDEALLOCATE(sf%p_k)
    SDEALLOCATE(sf%gradf_old)
    SDEALLOCATE(sf%H_k)
  END SUBROUTINE BFGS_free

END MODULE MODgvec_free_b_minimizer_BFGS
