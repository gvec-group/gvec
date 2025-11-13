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
USE MODgvec_c_hmap,     ONLY: c_hmap       ,c_hmap_auxvar
#ifdef PP_WHICH_HMAP
#  if PP_WHICH_HMAP == 1
USE MODgvec_hmap_RZ,    ONLY: t_hmap_RZ    ,t_hmap_RZ_auxvar
#  elif PP_WHICH_HMAP == 3
USE MODgvec_hmap_cyl,   ONLY: t_hmap_cyl   ,t_hmap_cyl_auxvar
#  elif PP_WHICH_HMAP == 10
USE MODgvec_hmap_knot,  ONLY: t_hmap_knot  ,t_hmap_knot_auxvar
#  elif PP_WHICH_HMAP == 20
USE MODgvec_hmap_frenet,ONLY: t_hmap_frenet,t_hmap_frenet_auxvar
#  elif PP_WHICH_HMAP == 21
USE MODgvec_hmap_axisNB,ONLY: t_hmap_axisNB,t_hmap_axisNB_auxvar
#  endif
#else
USE MODgvec_hmap_RZ,    ONLY: t_hmap_RZ    ,t_hmap_RZ_auxvar
USE MODgvec_hmap_cyl,   ONLY: t_hmap_cyl   ,t_hmap_cyl_auxvar
USE MODgvec_hmap_knot,  ONLY: t_hmap_knot  ,t_hmap_knot_auxvar
USE MODgvec_hmap_frenet,ONLY: t_hmap_frenet,t_hmap_frenet_auxvar
USE MODgvec_hmap_axisNB,ONLY: t_hmap_axisNB,t_hmap_axisNB_auxvar
#endif /*defined(PP_WHICH_HMAP)*/

IMPLICIT NONE
PUBLIC


CONTAINS


!===================================================================================================================================
!> initialize the type hmap, also readin parameters here if necessary
!!
!===================================================================================================================================
SUBROUTINE hmap_new( sf, which_hmap,hmap_in)
! MODULES
USE MODgvec_Globals   , ONLY: abort,enter_subregion,exit_subregion
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: which_hmap         !! input number of field periods
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP), INTENT(IN),OPTIONAL :: hmap_in       !! if present, copy this hmap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(PP_T_HMAP),ALLOCATABLE,INTENT(INOUT) :: sf !! self
#else
  CLASS(c_hmap), INTENT(IN),OPTIONAL :: hmap_in       !! if present, copy this hmap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_hmap),ALLOCATABLE,INTENT(INOUT) :: sf !! self
#endif /*defined(PP_WHICH_HMAP)*/
!===================================================================================================================================
  CALL enter_subregion("hmap")
  IF(.NOT. PRESENT(hmap_in))THEN
    SELECT CASE(which_hmap)
#ifdef PP_WHICH_HMAP
    CASE(PP_WHICH_HMAP)
      sf=PP_T_HMAP()
    CASE DEFAULT
      CALL abort(__STAMP__, &
           "FIXED HMAP TO PP_WHICH_HMAP AT COMPILE TIME,  hmap choice is therefore not compatible  !")
#else
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
#endif /*defined(PP_WHICH_HAP)*/
    END SELECT
    sf%which_hmap=which_hmap
  ELSE
    IF(which_hmap.NE.hmap_in%which_hmap) CALL abort(__STAMP__, &
       "hmap_in does not coincide with requested hmap in hmap_new")
    ALLOCATE(sf,source=hmap_in)
  END IF
  CALL exit_subregion("hmap")
END SUBROUTINE hmap_new

!===================================================================================================================================
!> initialize the  hmap auxiliary variables, depends on hmap type
!!
!===================================================================================================================================
SUBROUTINE hmap_new_auxvar(hmap,zeta,xv,do_2nd_der)
! MODULES
USE MODgvec_Globals   , ONLY: abort,wp
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP),INTENT(IN) :: hmap
#else
  CLASS(c_hmap),  INTENT(IN) :: hmap
#endif
  REAL(wp)     ,  INTENT(IN) :: zeta(:)
  LOGICAL      ,  INTENT(IN) :: do_2nd_der
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP_AUXVAR),ALLOCATABLE,INTENT(INOUT) :: xv(:) !! self
#else
  CLASS(c_hmap_auxvar),ALLOCATABLE,INTENT(INOUT) :: xv(:) !! self
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,nzeta
!===================================================================================================================================
  nzeta=SIZE(zeta)
#ifdef PP_WHICH_HMAP
  IF(hmap%which_hmap .NE. PP_WHICH_HMAP) CALL abort(__STAMP__, &
           "FIXED HMAP TO PP_WHICH_HMAP AT COMPILE TIME,  which_hmap choice is therefore not compatible  !")
  ALLOCATE(PP_T_HMAP_AUXVAR :: xv(nzeta))
  !$OMP PARALLEL DO &
  !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,nzeta
    xv(i)= PP_T_HMAP_AUXVAR(hmap,zeta(i),do_2nd_der)
  END DO !i
  !$OMP END PARALLEL DO

#else

  SELECT TYPE(hmap)
  CLASS IS(t_hmap_RZ)
    ALLOCATE(t_hmap_RZ_auxvar :: xv(nzeta))
    SELECT TYPE(xv)
    TYPE IS(t_hmap_RZ_auxvar)
      !$OMP PARALLEL DO &
      !$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
      DO i=1,nzeta
        xv(i)= t_hmap_RZ_auxvar(hmap,zeta(i),do_2nd_der)
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
        xv(i)= t_hmap_cyl_auxvar(hmap,zeta(i),do_2nd_der)
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
        xv(i)= t_hmap_knot_auxvar(hmap,zeta(i),do_2nd_der)
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
        xv(i)= t_hmap_frenet_auxvar(hmap,zeta(i),do_2nd_der)
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
        xv(i)= t_hmap_axisNB_auxvar(hmap,zeta(i),do_2nd_der)
      END DO !i
      !$OMP END PARALLEL DO
    END SELECT !TYPE(xv)
  CLASS DEFAULT
    CALL abort(__STAMP__, &
          "hmap_new_auxvar: this hmap class is not implemented  !")
  END SELECT
#endif /*defined(PP_WHICH_HAP)*/

END SUBROUTINE hmap_new_auxvar

END MODULE MODgvec_hmap
