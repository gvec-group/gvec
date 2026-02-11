!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

!===================================================================================================================================
!>
!!# Module ** polyProfile **
!!
!! Defines a 1-D profile in rho^2 via a power polynomial.
!===================================================================================================================================
MODULE MODgvec_free_b_minimizer_moment_GD
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_moment_GD
    REAL(wp), ALLOCATABLE :: velocity(:) ! velocity -> gradient average
    REAL(wp)              :: beta ! decay rate
    CONTAINS
    PROCEDURE :: get_update => moment_GD_get_update
    PROCEDURE :: reset      => moment_GD_reset
    FINAL     :: moment_GD_free

  END TYPE t_free_b_minimizer_moment_GD

  INTERFACE t_free_b_minimizer_moment_GD
      MODULE PROCEDURE moment_GD_new
  END INTERFACE t_free_b_minimizer_moment_GD

  CONTAINS

  FUNCTION moment_GD_new(n_dofs, beta, dt) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)            :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN)           :: beta ! initial / default time step
    REAL(wp), INTENT(IN), OPTIONAL :: dt ! initial / default time step

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_moment_GD) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    sf%beta  = beta
    IF (PRESENT(dt)) THEN
        sf%dt_initial = dt
    ELSE
        sf%dt_initial = 1e-3_wp
    END IF
    sf%dt =  sf%dt_initial

    ALLOCATE(sf%velocity(sf%n_dofs))
    sf%velocity = 0.0_wp
    sf%n_updates = 0
  END FUNCTION moment_GD_new

  SUBROUTINE moment_GD_reset(sf)
    ! INPUT VARIABLES -------------------------!
    ! OUTPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_moment_GD), INTENT(INOUT) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%velocity = 0.0_wp
    sf%n_updates = 0
  END SUBROUTINE moment_GD_reset

  FUNCTION moment_GD_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_moment_GD), INTENT(INOUT) :: sf
    REAL(wp)                           , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%velocity = sf%beta*sf%velocity + (1-sf%beta)*force_in
    push_out = -sf%dt*sf%velocity
    sf%n_updates =  sf%n_updates + 1
  END FUNCTION moment_GD_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE moment_GD_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_moment_GD), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%velocity)
  END SUBROUTINE moment_GD_free

END MODULE MODgvec_free_b_minimizer_moment_GD
