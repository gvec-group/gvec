!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_SR1
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_SR1
    ! Symmetric-rank-1 method is a quasi newton method similar to BFGS but does not guarantee a positive definite matrix
    ! It might skip updates of the Hessian to retain stability
    ! Has 3 hyper parameters:
    ! the initial step-size, the maximum allowed step size, the number of pure gradient descents steps done before SR1 kicks in

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
    PROCEDURE :: get_update => SR1_get_update
    PROCEDURE :: reset      => SR1_reset
    FINAL     :: SR1_free

  END TYPE t_free_b_minimizer_SR1

  INTERFACE t_free_b_minimizer_SR1
      MODULE PROCEDURE SR1_new
  END INTERFACE t_free_b_minimizer_SR1

  CONTAINS

  FUNCTION SR1_new(n_dofs, dt, n_grad_steps, max_dt) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)            :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN)           :: dt ! initial / default time step
    REAL(wp), INTENT(IN)           :: max_dt
    INTEGER, INTENT(IN), OPTIONAL  :: n_grad_steps

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_SR1) :: sf
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
    sf%I = 0.0_wp
    DO i=1,sf%n_dofs
        sf%H_k(i,i) = 1.0_wp
        sf%I(i,i) = 1.0_wp
        sf%B_k(i,i) = 1.0_wp
    END DO

    sf%n_updates = 0
  END FUNCTION SR1_new

  SUBROUTINE SR1_reset(sf)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_SR1), INTENT(INOUT) :: sf
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
        sf%B_k(i,i) = 1.0_wp
    END DO
    sf%n_updates = 0
  END SUBROUTINE SR1_reset

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

  FUNCTION SR1_get_update(sf, force_in) RESULT(push_out)
    USE MODgvec_MHD3D_Vars, ONLY: Precond_b_X1, Precond_b_X2, X1_base, X2_base
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_SR1), INTENT(INOUT) :: sf
    REAL(wp)                     , INTENT(IN)    :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES
    REAL(wp) :: p_k(sf%n_dofs), g_k(sf%n_dofs), sHy(sf%n_dofs), sHyy, yBs(sf%n_dofs), yBss
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
        sHy = sf%s_k-MATMUL(sf%H_k,sf%y_k)
        sHyy = DOT_PRODUCT(sHy,sf%y_k)
        yBs = sf%y_k-MATMUL(sf%B_k,sf%s_k)
        yBss = DOT_PRODUCT(sf%s_k, yBs)
        IF ((ABS(yBss).GE.1e-8_wp*NORM2(sf%s_k)*NORM2(yBs)).AND.(ABS(yBss)>1e-12_wp)) THEN
          sf%H_k = sf%H_k + OUTER_PRODUCT(sHy,sHy) / sHyy
          sf%B_k = sf%B_k + OUTER_PRODUCT(yBs,yBs) / yBss
        END IF
        p_k = -MATMUL(sf%H_k,g_k)
    END IF

    sf%gradf_old = g_k

    push_out = sf%dt*p_k
    ! update after adapting dt outside
    sf%s_k = p_k
    sf%n_updates = sf%n_updates+1
  END FUNCTION SR1_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE SR1_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_SR1), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%s_k)
    SDEALLOCATE(sf%y_k)
    SDEALLOCATE(sf%p_k)
    SDEALLOCATE(sf%gradf_old)
    SDEALLOCATE(sf%H_k)
    SDEALLOCATE(sf%B_k)
  END SUBROUTINE SR1_free

END MODULE MODgvec_free_b_minimizer_SR1
