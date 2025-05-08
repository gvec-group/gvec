!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** HMAP new **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_hmap
! MODULES
USE MODgvec_c_hmap    , ONLY: c_hmap,c_hmap_auxvar
IMPLICIT NONE
PUBLIC


CONTAINS

!===================================================================================================================================
!> initialize the type hmap, also readin parameters here if necessary
!!
!===================================================================================================================================
SUBROUTINE hmap_new( sf, which_hmap,hmap_in)
! MODULES
USE MODgvec_Globals   , ONLY: abort
USE MODgvec_hmap_RZ   , ONLY: t_hmap_RZ
USE MODgvec_hmap_cyl  , ONLY: t_hmap_cyl
USE MODgvec_hmap_knot , ONLY: t_hmap_knot
USE MODgvec_hmap_frenet,ONLY: t_hmap_frenet
USE MODgvec_hmap_axisNB,ONLY: t_hmap_axisNB
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: which_hmap         !! input number of field periods
  CLASS(c_hmap), INTENT(IN),OPTIONAL :: hmap_in       !! if present, copy this hmap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_hmap),ALLOCATABLE,INTENT(INOUT) :: sf !! self
!===================================================================================================================================
  IF(.NOT. PRESENT(hmap_in))THEN
    SELECT CASE(which_hmap)
    CASE(1)
      sf=t_hmap_RZ()
    !CASE(2)
    !  ALLOCATE(t_hmap_RphiZ :: sf)
    CASE(3)
      sf=t_hmap_cyl()
    CASE(10)
      sf=t_hmap_knot()
    CASE(20)
      sf=t_hmap_frenet()
    CASE(21)
      sf=t_hmap_axisNB()
    CASE DEFAULT
      CALL abort(__STAMP__, &
           "this hmap choice does not exist  !")
    END SELECT
    sf%which_hmap=which_hmap
  ELSE
    IF(which_hmap.NE.hmap_in%which_hmap) CALL abort(__STAMP__, &
       "hmap_in does not coincide with requested hmap in hmap_new")
    ALLOCATE(sf,source=hmap_in)
  END IF

END SUBROUTINE hmap_new


END MODULE MODgvec_hmap
