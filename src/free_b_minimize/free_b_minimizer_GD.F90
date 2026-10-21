!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_GD
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_GD
    ! Basic gradient descent for a fixed given step size
    ! One hyper-parameter, the step size

    CONTAINS
    PROCEDURE :: get_update => GD_get_update
    PROCEDURE :: reset      => GD_reset
    FINAL     :: GD_free

  END TYPE t_free_b_minimizer_GD

  INTERFACE t_free_b_minimizer_GD
      MODULE PROCEDURE GD_new
  END INTERFACE t_free_b_minimizer_GD

  CONTAINS

  FUNCTION GD_new(n_dofs, dt, max_dt) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER,  INTENT(IN)  :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN)  :: dt ! initial / default time step
    REAL(wp), INTENT(IN)  :: max_dt ! initial / default time step
    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_GD) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs

    sf%dt_initial = dt
    sf%max_dt = max_dt
    sf%dt = sf%dt_initial

    sf%n_updates = 0
  END FUNCTION GD_new

  SUBROUTINE GD_reset(sf)
    ! INPUT VARIABLES -------------------------!
    ! OUTPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_GD), INTENT(INOUT) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%n_updates = 0
  END SUBROUTINE GD_reset

  FUNCTION GD_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_GD), INTENT(INOUT) :: sf
    REAL(wp)                    , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    push_out = -sf%dt*force_in
    sf%n_updates = sf%n_updates+1
  END FUNCTION GD_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE GD_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_GD), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
  END SUBROUTINE GD_free

END MODULE MODgvec_free_b_minimizer_GD
