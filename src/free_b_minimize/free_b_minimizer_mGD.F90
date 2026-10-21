!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_mGD
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_mGD
    ! momentum gradient descent also includes some history of the gradients via a momentum term.
    ! Has 2 hyper parameters:
    ! step size dt, and momentum weighting factor beta.
    ! lower betas discard older gradients faster.

    REAL(wp), ALLOCATABLE :: velocity(:) ! velocity -> gradient average
    REAL(wp)              :: beta ! decay rate
    CONTAINS
    PROCEDURE :: get_update => mGD_get_update
    PROCEDURE :: reset      => mGD_reset
    FINAL     :: mGD_free

  END TYPE t_free_b_minimizer_mGD

  INTERFACE t_free_b_minimizer_mGD
      MODULE PROCEDURE mGD_new
  END INTERFACE t_free_b_minimizer_mGD

  CONTAINS

  FUNCTION mGD_new(n_dofs, beta, dt, max_dt) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)  :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN) :: beta
    REAL(wp), INTENT(IN) :: dt ! initial / default time step
    REAL(wp), INTENT(IN) :: max_dt

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_mGD) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    sf%beta  = beta

    sf%dt_initial = dt
    sf%dt = sf%dt_initial
    sf%max_dt = max_dt

    ALLOCATE(sf%velocity(sf%n_dofs))
    sf%velocity = 0.0_wp
    sf%n_updates = 0
  END FUNCTION mGD_new

  SUBROUTINE mGD_reset(sf)
    ! INPUT VARIABLES -------------------------!
    ! OUTPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_mGD), INTENT(INOUT) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%velocity = 0.0_wp
    sf%n_updates = 0
  END SUBROUTINE mGD_reset

  FUNCTION mGD_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_mGD), INTENT(INOUT) :: sf
    REAL(wp)                           , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%velocity = sf%beta*sf%velocity + (1-sf%beta)*force_in
    push_out = -sf%dt*sf%velocity
    sf%n_updates =  sf%n_updates + 1
  END FUNCTION mGD_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE mGD_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_mGD), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%velocity)
  END SUBROUTINE mGD_free

END MODULE MODgvec_free_b_minimizer_mGD
