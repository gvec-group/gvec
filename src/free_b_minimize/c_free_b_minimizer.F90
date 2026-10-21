!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

MODULE MODgvec_free_b_minimizer_base
  ! MODULES
  USE MODgvec_Globals ,ONLY: wp, abort
  IMPLICIT NONE

  PUBLIC

  TYPE, ABSTRACT :: c_free_b_minimizer
      REAL(wp) :: dt_initial ! initial time-step
      REAL(wp) :: dt ! current time step
      REAL(wp) :: max_dt ! limit the maximum time step, especially relevant for quasi-newton methods
      INTEGER  :: n_dofs ! number of degrees of freedom
      INTEGER  :: n_updates ! number of calls to get_update
      contains

      PROCEDURE(i_fun_get_update), DEFERRED :: get_update
      PROCEDURE(i_fun_reset),      DEFERRED :: reset
  end type c_free_b_minimizer

  ABSTRACT INTERFACE

      FUNCTION i_fun_get_update( sf, force_in) RESULT(push_out)
          IMPORT c_free_b_minimizer
          IMPORT wp
          CLASS(c_free_b_minimizer), INTENT(INOUT)  :: sf
          REAL(wp)                 , INTENT(IN)  :: force_in(sf%n_dofs)
          REAL(wp)                               :: push_out(sf%n_dofs)
      END FUNCTION i_fun_get_update

      SUBROUTINE i_fun_reset(sf)
          IMPORT c_free_b_minimizer
          CLASS(c_free_b_minimizer), INTENT(INOUT)   :: sf
      END SUBROUTINE i_fun_reset
  END INTERFACE

END MODULE MODgvec_free_b_minimizer_base
