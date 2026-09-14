!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_aGD
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_aGD
    ! Accelerated gradient descent. Same as the fixed boundary minimizer, using a momentum and damping term.
    ! One controllable hyper parameter, the step size dt.

    REAL(wp), ALLOCATABLE :: velocity(:) ! velocity -> gradient average
    REAL(wp)              :: force_hist ! force from the previous update |F_{n-1}|^2
    ! damping parameters
    INTEGER               :: n_tau_hist
    REAL(wp)              :: tau_bar ! tau_bar = 0.5*<tau_hist>
    REAL(wp), ALLOCATABLE :: tau_hist(:) ! history of last 10 tau_hat=min(|ln(|F|_n^2)| / force_hist, 0.15), initialized to 0.15

    CONTAINS
    PROCEDURE :: get_update => aGD_get_update
    PROCEDURE :: reset      => aGD_reset
    FINAL     :: aGD_free

  END TYPE t_free_b_minimizer_aGD

  INTERFACE t_free_b_minimizer_aGD
      MODULE PROCEDURE aGD_new
  END INTERFACE t_free_b_minimizer_aGD

  CONTAINS

  FUNCTION aGD_new(n_dofs, dt, n_tau_hist, max_dt) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)            :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN)           :: dt ! initial / default time step
    REAL(wp), INTENT(IN)           :: max_dt
    INTEGER, INTENT(IN) , OPTIONAL :: n_tau_hist  ! number of remembered tau_hat

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_aGD) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    sf%dt_initial = dt
    sf%dt =  sf%dt_initial
    sf%max_dt = max_dt


    IF (PRESENT(n_tau_hist)) THEN
        sf%n_tau_hist = n_tau_hist
    ELSE
        sf%n_tau_hist = 10
    END IF


    ALLOCATE(sf%velocity(sf%n_dofs))
    ALLOCATE(sf%tau_hist(sf%n_tau_hist))

    sf%force_hist = 1e-8_wp
    sf%velocity = 0.0_wp
    sf%tau_bar = 0.15_wp
    sf%tau_hist = 0.15_wp
    sf%n_updates = 0
  END FUNCTION aGD_new

  SUBROUTINE aGD_reset(sf)
    ! INPUT VARIABLES -------------------------!
    ! OUTPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_aGD), INTENT(INOUT) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%force_hist = 1e-8_wp
    sf%velocity = 0.0_wp
    sf%tau_bar = 0.15_wp
    sf%tau_hist = 0.15_wp
    sf%n_updates = 0
  END SUBROUTINE aGD_reset

  FUNCTION aGD_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_aGD), INTENT(INOUT) :: sf
    REAL(wp)                           , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES
    REAL(wp) :: force_norm
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    force_norm = SUM(force_in**2)

    sf%tau_hist(1:sf%n_tau_hist-1) = sf%tau_hist(2:sf%n_tau_hist)
    sf%tau_hist(sf%n_tau_hist) = MIN(ABS(LOG(force_norm / sf%force_hist)), 0.15_wp)
    sf%tau_bar = 0.5*SUM(sf%tau_hist) / sf%n_tau_hist

    sf%velocity = 1 / (1 + sf%tau_bar) * ((1-sf%tau_bar)*sf%velocity + sf%dt*force_in)
    push_out = -sf%dt*sf%velocity

    sf%force_hist = force_norm
    sf%n_updates =  sf%n_updates + 1
  END FUNCTION aGD_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE aGD_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_aGD), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%velocity)
    SDEALLOCATE(sf%tau_hist)
  END SUBROUTINE aGD_free

END MODULE MODgvec_free_b_minimizer_aGD
