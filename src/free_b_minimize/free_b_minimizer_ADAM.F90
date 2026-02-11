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
MODULE MODgvec_free_b_minimizer_ADAM
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp
  USE MODgvec_free_b_minimizer_base, ONLY: c_free_b_minimizer
  IMPLICIT NONE

  PUBLIC

  TYPE, EXTENDS(c_free_b_minimizer) :: t_free_b_minimizer_ADAM
    REAL(wp), ALLOCATABLE :: velocity(:)
    REAL(wp), ALLOCATABLE :: momentum(:)

    REAL(wp)              :: beta1 ! decay rate momentum
    REAL(wp)              :: beta2 ! decay rate velocity
    CONTAINS
    PROCEDURE :: get_update => ADAM_get_update
    PROCEDURE :: reset      => ADAM_reset
    FINAL     :: ADAM_free

  END TYPE t_free_b_minimizer_ADAM

  INTERFACE t_free_b_minimizer_ADAM
      MODULE PROCEDURE ADAM_new
  END INTERFACE t_free_b_minimizer_ADAM

  CONTAINS

  FUNCTION ADAM_new(n_dofs, beta1, beta2, dt) RESULT(sf)
    ! INPUT VARIABLES -------------------------!
    INTEGER, INTENT(IN)            :: n_dofs  ! number of degrees of freedom
    REAL(wp), INTENT(IN)           :: beta1 ! decay rate momentum
    REAL(wp), INTENT(IN)           :: beta2 ! decay rate velocity
    REAL(wp), INTENT(IN), OPTIONAL :: dt ! initial / default time step

    ! OUTPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_ADAM) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%n_dofs = n_dofs
    sf%beta1  = beta1
    sf%beta2 = beta2
    IF (PRESENT(dt)) THEN
        sf%dt_initial = dt
    ELSE
        sf%dt_initial = 1e-3_wp
    END IF
    sf%dt =  sf%dt_initial

    ALLOCATE(sf%momentum(sf%n_dofs))
    ALLOCATE(sf%velocity(sf%n_dofs))

    sf%momentum = 0.0_wp
    sf%velocity = 0.0_wp
    sf%n_updates = 0
  END FUNCTION ADAM_new

  SUBROUTINE ADAM_reset(sf)
    ! INPUT VARIABLES -------------------------!
    ! OUTPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_ADAM), INTENT(INOUT) :: sf
    ! LOCAL VARIABLES -------------------------!
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    sf%dt = sf%dt_initial
    sf%momentum = 0.0_wp
    sf%velocity = 0.0_wp
    sf%n_updates = 0
  END SUBROUTINE ADAM_reset

  FUNCTION ADAM_get_update(sf, force_in) RESULT(push_out)
    ! INPUT VARIABLES -------------------------!
    CLASS(t_free_b_minimizer_ADAM), INTENT(INOUT) :: sf
    REAL(wp)                      , INTENT(IN) :: force_in(sf%n_dofs)
    ! OUTPUT VARIABLES -------------------------!
    REAL(wp) :: push_out(sf%n_dofs)
    ! LOCAL VARIABLES -------------------------!
    REAL(wp) :: m_hat(sf%n_dofs)
    REAL(wp) :: v_hat(sf%n_dofs)
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    push_out = 0.0_wp
    sf%momentum = sf%beta1*sf%momentum + (1-sf%beta1)*force_in
    sf%velocity = sf%beta2*sf%velocity + (1-sf%beta2)*force_in**2

    ! NOTE: free_boundary_updates starts at 0, hence the +1
    m_hat = sf%momentum / (1-sf%beta1**(sf%n_updates+1))
    v_hat = sf%velocity / (1-sf%beta2**(sf%n_updates+1))

    push_out = -sf%dt*m_hat/(SQRT(v_hat)+1e-12_wp)

    sf%n_updates = sf%n_updates + 1
  END FUNCTION ADAM_get_update

  !===================================================================================================================================
  !> finalize the type gradient descent
  !!
  !===================================================================================================================================
  SUBROUTINE ADAM_free(sf)
    ! INPUT VARIABLES -------------------------!
    TYPE(t_free_b_minimizer_ADAM), INTENT(INOUT)  :: sf !! self
    ! CODE --------------------------------------------------------------------------------------------------------------------------!
    SDEALLOCATE(sf%velocity)
    SDEALLOCATE(sf%momentum)
  END SUBROUTINE ADAM_free

END MODULE MODgvec_free_b_minimizer_ADAM
