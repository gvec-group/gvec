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
    ! Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm is a quasi Newton method that iteratively builds an approximation of the inverse Hessian H_k.
    ! Typically one would perform line-search to find a valid step-size. We do not do this since we cannot evaluate the action.
    ! Thus, we opt to use the adaptive BFGS from Yu et al., which estimates the step size from the Hessian B_k.
    ! Has 3 hyper parameters:
    ! the initial step-size, the maximum allowed step size, the number of pure gradient descents steps done before BFGS kicks in

    REAL(wp), ALLOCATABLE :: s_k(:) ! last update
    REAL(wp), ALLOCATABLE :: y_k(:) ! difference of gradients
    REAL(wp), ALLOCATABLE :: p_k (:) ! search direction => -H_k grad f(x_k)
    REAL(wp), ALLOCATABLE :: gradf_old(:) ! gradient of step k
    REAL(wp), ALLOCATABLE :: H_k(:,:) ! approximate inverse Hessian
    REAL(wp), ALLOCATABLE :: B_k(:,:) ! approximate Hessian
    REAL(wp), ALLOCATABLE :: I(:,:)
    INTEGER  :: n_grad_steps
    REAL(wp) :: M ! stepsize controll parameter / self-concordance

    CONTAINS
    PROCEDURE :: get_update => BFGS_get_update
    PROCEDURE :: reset      => BFGS_reset
    FINAL     :: BFGS_free

  END TYPE t_free_b_minimizer_BFGS

  INTERFACE t_free_b_minimizer_BFGS
      MODULE PROCEDURE BFGS_new
  END INTERFACE t_free_b_minimizer_BFGS

  CONTAINS

  FUNCTION BFGS_new(n_dofs, dt, n_grad_steps, max_dt, M) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)             :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN)            :: dt ! initial / default time step
    REAL(wp), INTENT(IN)            :: max_dt
    INTEGER, INTENT(IN), OPTIONAL   :: n_grad_steps
    REAL(wp), INTENT(IN), OPTIONAL  :: M

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

    sf%dt_initial = dt
    sf%dt = dt
    sf%max_dt = max_dt

    IF (PRESENT(M)) THEN
        sf%M = M
    ELSE
        sf%M = 1e-3_wp
    END IF


    ALLOCATE(sf%s_k(1:sf%n_dofs))
    ALLOCATE(sf%y_k(1:sf%n_dofs))
    ALLOCATE(sf%p_k(1:sf%n_dofs))
    ALLOCATE(sf%gradf_old(1:sf%n_dofs))
    ALLOCATE(sf%H_k(1:sf%n_dofs, 1:sf%n_dofs))
    ALLOCATE(sf%B_k(1:sf%n_dofs, 1:sf%n_dofs))
    ALLOCATE(sf%I(1:sf%n_dofs, 1:sf%n_dofs))

    sf%s_k = 0.0_wp
    sf%y_k = 0.0_wp
    sf%p_k = 0.0_wp
    sf%gradf_old = 0.0_wp
    sf%H_k = 0.0_wp
    sf%B_k = 0.0_wp
    sf%I = 0.0_wp
    DO i=1,sf%n_dofs
        sf%H_k(i,i) = 1.0_wp
        sf%B_k(i,i) = 1.0_wp
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
    USE MODgvec_MHD3D_Vars, ONLY: Precond_b_X1, Precond_b_X2, X1_base, X2_base
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_BFGS), INTENT(INOUT) :: sf
    REAL(wp)                      , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES
    REAL(wp) :: dot_sy, norm_p_k, p_k(sf%n_dofs), g_k(sf%n_dofs), dot_gp
    INTEGER :: i
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    g_k = force_in

    IF (sf%n_updates.LE.sf%n_grad_steps) THEN
        p_k = -g_k
        DO i=1,X1_base%f%modes
            sf%B_k(i,i) = Precond_b_X1(i)
            sf%H_k(i,i) = 1 / Precond_b_X1(i)
        END DO

        DO i=1,X2_base%f%modes
            sf%B_k(X1_base%f%modes+i,X1_base%f%modes+i) = Precond_b_X2(i)
            sf%H_k(X1_base%f%modes+i,X1_base%f%modes+i) = 1 / Precond_b_X2(i)
        END DO
    ELSE
        ! update s_k with the actual used dt
        sf%s_k = sf%dt*sf%s_k
        sf%y_k = g_k - sf%gradf_old
        dot_sy = DOT_PRODUCT(sf%s_k,sf%y_k)

        ! only update valid stepsizes
        IF (dot_sy>0.0_wp) THEN
          sf%B_k = sf%B_k - MATMUL(sf%B_k,OUTER_PRODUCT(sf%s_k,MATMUL(sf%s_k,sf%B_k))) / (DOT_PRODUCT(sf%s_k,MATMUL(sf%B_k,sf%s_k))) &
            + OUTER_PRODUCT(sf%y_k,sf%y_k) / dot_sy
          sf%H_k = OUTER_PRODUCT(sf%s_k,sf%s_k) / dot_sy + MATMUL(sf%I-OUTER_PRODUCT(sf%s_k,sf%y_k) / dot_sy,MATMUL(sf%H_k, sf%I-OUTER_PRODUCT(sf%y_k,sf%s_k) / dot_sy))
        END IF

        p_k = -MATMUL(sf%H_k,g_k)
        dot_gp = DOT_PRODUCT(g_k,p_k)
        norm_p_k = SQRT(DOT_PRODUCT(p_k,MATMUL(sf%B_k,p_k)))
        sf%dt = MIN(-dot_gp / (norm_p_k*(norm_p_k-dot_gp)),sf%max_dt)
    END IF
    sf%gradf_old = g_k

    push_out = sf%dt*p_k
    ! update after adapting dt outside
    sf%s_k = p_k
    sf%n_updates = sf%n_updates+1
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
    SDEALLOCATE(sf%B_k)
  END SUBROUTINE BFGS_free

END MODULE MODgvec_free_b_minimizer_BFGS
