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
USE MODgvec_c_hmap   , ONLY: c_hmap,c_hmap_auxvar
USE MODgvec_hmap_RZ   , ONLY: t_hmap_RZ    ,t_hmap_RZ_auxvar
USE MODgvec_hmap_cyl  , ONLY: t_hmap_cyl   ,t_hmap_cyl_auxvar
USE MODgvec_hmap_knot , ONLY: t_hmap_knot  ,t_hmap_knot_auxvar
USE MODgvec_hmap_frenet,ONLY: t_hmap_frenet,t_hmap_frenet_auxvar
USE MODgvec_hmap_axisNB,ONLY: t_hmap_axisNB,t_hmap_axisNB_auxvar

IMPLICIT NONE
PUBLIC


CONTAINS

#ifdef PP_WHICH_HMAP
!===================================================================================================================================
!> initialize the type hmap, also readin parameters here if necessary
!!
!===================================================================================================================================
SUBROUTINE hmap_new( sf, which_hmap,hmap_in)
! MODULES
USE MODgvec_Globals   , ONLY: abort,Unit_stdOut
USE PP_MOD_HMAP   , ONLY: PP_T_HMAP,PP_T_HMAP_AUXVAR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: which_hmap         !! input number of field periods
  TYPE(PP_T_HMAP), INTENT(IN),OPTIONAL :: hmap_in       !! if present, copy this hmap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(PP_T_HMAP),ALLOCATABLE,INTENT(INOUT) :: sf !! self
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A,I4,A,I4)')'INIT PRECOMPILED HMAP, with PP_WHICH_HMAP = ',PP_WHICH_HMAP,' ! input which_hmap=',which_hmap
  IF(.NOT. PRESENT(hmap_in))THEN
    SELECT CASE(which_hmap)
    CASE(PP_WHICH_HMAP)
      sf=PP_T_HMAP()
    CASE DEFAULT
      CALL abort(__STAMP__, &
           "FIXED HMAP TO PP_WHICH_HMAP AT COMPILE TIME,  hmap choice is therefore not compatible  !")
    END SELECT
    sf%which_hmap=which_hmap
  ELSE
    ALLOCATE(sf,source=hmap_in)
  END IF

END SUBROUTINE hmap_new

!===================================================================================================================================
!> initialize the  hmap auxiliary variables, depends on hmap type
!!
!===================================================================================================================================
SUBROUTINE hmap_new_auxvar(hmap,zeta,xv)
! MODULES
USE MODgvec_Globals   , ONLY: abort,wp
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  TYPE(PP_T_HMAP), INTENT(IN) :: hmap
  REAL(wp)     , INTENT(IN) :: zeta(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(PP_T_HMAP_AUXVAR),ALLOCATABLE,INTENT(INOUT) :: xv(:) !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,nzeta
!===================================================================================================================================
  nzeta=SIZE(zeta)
  SELECT CASE(hmap%which_hmap)
  CASE(PP_WHICH_HMAP)
    ALLOCATE(PP_T_HMAP_AUXVAR :: xv(nzeta))
    !$OMP PARALLEL DO &
    !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
    DO i=1,nzeta
      xv(i)= PP_T_HMAP_AUXVAR(hmap,zeta(i))
    END DO !i
    !$OMP END PARALLEL DO
  CASE DEFAULT
    CALL abort(__STAMP__, &
           "FIXED HMAP TO PP_WHICH_HMAP AT COMPILE TIME,  which_hmap choice is therefore not compatible  !")
  END SELECT


END SUBROUTINE hmap_new_auxvar


#else /*PP_WHICH_HMAP not defined*/
!===================================================================================================================================
!> initialize the type hmap, also readin parameters here if necessary
!!
!===================================================================================================================================
SUBROUTINE hmap_new( sf, which_hmap,hmap_in)
! MODULES
USE MODgvec_Globals   , ONLY: abort
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

!===================================================================================================================================
!> initialize the  hmap auxiliary variables, depends on hmap type
!!
!===================================================================================================================================
SUBROUTINE hmap_new_auxvar(hmap,zeta,xv)
! MODULES
USE MODgvec_Globals   , ONLY: abort,wp
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_hmap), INTENT(IN) :: hmap
  REAL(wp)     , INTENT(IN) :: zeta(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_hmap_auxvar),ALLOCATABLE,INTENT(INOUT) :: xv(:) !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,nzeta
!===================================================================================================================================
  nzeta=SIZE(zeta)
  SELECT TYPE(hmap)
  CLASS IS(t_hmap_RZ)
    ALLOCATE(t_hmap_RZ_auxvar :: xv(nzeta))
    SELECT TYPE(xv)
    TYPE IS(t_hmap_RZ_auxvar)
      !$OMP PARALLEL DO &
      !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
      DO i=1,nzeta
        xv(i)= t_hmap_RZ_auxvar(hmap,zeta(i))
      END DO !i
      !$OMP END PARALLEL DO
    END SELECT !TYPE(xv) 
  CLASS IS(t_hmap_cyl)
    ALLOCATE(t_hmap_cyl_auxvar :: xv(nzeta))
    SELECT TYPE(xv)
    TYPE IS(t_hmap_cyl_auxvar)
      !$OMP PARALLEL DO &
      !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
      DO i=1,nzeta
        xv(i)= t_hmap_cyl_auxvar(hmap,zeta(i))
      END DO !i
      !$OMP END PARALLEL DO
    END SELECT !TYPE(xv) 
  CLASS IS(t_hmap_knot)
    ALLOCATE(t_hmap_knot_auxvar :: xv(nzeta))
    SELECT TYPE(xv)
    TYPE IS(t_hmap_knot_auxvar)
      !$OMP PARALLEL DO &
      !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
      DO i=1,nzeta
        xv(i)= t_hmap_knot_auxvar(hmap,zeta(i))
      END DO !i
      !$OMP END PARALLEL DO
    END SELECT !TYPE(xv) 
  CLASS IS(t_hmap_frenet)
    ALLOCATE(t_hmap_frenet_auxvar :: xv(nzeta))
    SELECT TYPE(xv)
    TYPE IS(t_hmap_frenet_auxvar)
      !$OMP PARALLEL DO &
      !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
      DO i=1,nzeta
        xv(i)= t_hmap_frenet_auxvar(hmap,zeta(i))
      END DO !i
      !$OMP END PARALLEL DO
    END SELECT !TYPE(xv) 
  CLASS IS(t_hmap_axisNB)
    ALLOCATE(t_hmap_axisNB_auxvar :: xv(nzeta))
    SELECT TYPE(xv)
    TYPE IS(t_hmap_axisNB_auxvar)
      !$OMP PARALLEL DO &
      !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
      DO i=1,nzeta
        xv(i)= t_hmap_axisNB_auxvar(hmap,zeta(i))
      END DO !i
      !$OMP END PARALLEL DO
    END SELECT !TYPE(xv) 
  CLASS DEFAULT
    CALL abort(__STAMP__, &
          "hmap_new_auxvar: this hmap class is not implemented  !")
  END SELECT

END SUBROUTINE hmap_new_auxvar
#endif /*PP_WHICH_HMAP defined*/

END MODULE MODgvec_hmap
