!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **write_modes**
!!
!! Analyze and output equilibrium data
!!
!===================================================================================================================================
MODULE MODgvec_write_modes
! MODULES
USE MODgvec_Globals, ONLY:wp,abort
IMPLICIT NONE
PRIVATE

INTERFACE write_modes
  MODULE PROCEDURE write_modes
END INTERFACE

PUBLIC::write_modes
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> write modes prepared above
!!
!===================================================================================================================================
SUBROUTINE write_modes(fname,vname,nval,modes,xm,xn,coord,rho_first,values_in,VarNames_in)
! MODULES
USE MODgvec_Output_CSV, ONLY:WriteDataToCSV
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=*),INTENT(IN)   :: fname
  CHARACTER(LEN=*),INTENT(IN)   :: vname
  INTEGER         ,INTENT(IN)   :: modes
  INTEGER         ,INTENT(IN)   :: xm(1:modes)
  INTEGER         ,INTENT(IN)   :: xn(1:modes)
  REAL(wp),INTENT(IN)           :: coord(:)
  REAL(wp),INTENT(IN)           :: rho_first
  REAL(wp),INTENT(INOUT)        :: values_in(:,:)
  CHARACTER(LEN=*),INTENT(INOUT):: VarNames_in(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER         ,INTENT(INOUT):: nval
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                    :: i,iMode
  REAL(wp)                   :: minmaxval(2)
  REAL(wp) ,ALLOCATABLE      :: max_loc_val(:)
  CHARACTER(LEN=100),ALLOCATABLE :: varnames_max(:)
!===================================================================================================================================

  minmaxval(1)=MINVAL(values_in(nVal-modes:nVal,:))
  minmaxval(2)=MAXVAL(values_in(nVal-modes:nVal,:))

!  DO iMode=1,modes
!    nVal=nVal+1
!    WRITE(VarNames_in(nVal),'(A)')TRIM(VarNames_in(nVal-modes))//'_norm'
!    values_in(nVal,:)=values_in(nVal-modes,:)/(MAXVAL(ABS(values_in(nVal-modes,:)))+1.0E-12)
!  END DO

!  nVal=nVal+2
!  Varnames_in(nVal-1)=TRIM(vname)//', m= odd, n= 000'
!  Varnames_in(nVal)=  TRIM(vname)//', m=even, n= 000'
!  values_in(nVal-1:nVal,:)=0.
!  DO iMode=1,modes
!    IF((xn(iMode)).EQ.0)THEN
!      IF(MOD((xm(iMode)),2).NE.0)THEN
!        values_in(nVal-1,:)= values_in(nVal-1,:)+values_in(nVal-2-2*modes+iMode,:)
!      ELSE
!        values_in(nVal,:)= values_in(nVal,:)+values_in(nVal-2-2*modes+iMode,:)
!      END IF
!    END IF !n=0
!  END DO

  CALL WriteDataToCSV(VarNames_in(1:nVal),Values_in(1:nVal,:), TRIM(fname))

!  ALLOCATE(max_loc_val(nVal),Varnames_max(nVal))
!  DO i=1,nVal
!    max_loc_val(i)=coord(MAXLOC(ABS(values_in(i,:)),1))
!    Varnames_max(i)=TRIM(VarNames_in(i))//'_maxloc'
!  END DO
!  CALL WriteDataToCSV(VarNames_max(:) ,RESHAPE(max_loc_val(:),(/nval,1/)) &
!                             ,TRIM(fname) &
!                             ,append_in=.TRUE.,vfmt_in='E15.5')
!  DO i=1,nVal
!    max_loc_val(i)=      MAXVAL(ABS(values_in(i,:)))+1.0E-12
!    Varnames_max(i)=TRIM(VarNames_in(i))//'_maxval'
!  END DO
!  CALL WriteDataToCSV(VarNames_max(:) ,RESHAPE(max_loc_val(:),(/nval,1/)) &
!                             ,TRIM(fname) &
!                             ,append_in=.TRUE.,vfmt_in='E15.5')
!  DEALLOCATE(max_loc_val,Varnames_max)
!  !write position of first flux surface
!  CALL WriteDataToCSV((/'rhoFirst'/) ,RESHAPE((/rho_First/),(/1,1/)) &
!                             ,TRIM(fname) &
!                             ,append_in=.TRUE.)
!  !write position of first flux surface
!  CALL WriteDataToCSV((/'minval_total','maxval_total'/) ,RESHAPE(minmaxval,(/2,1/)) &
!                             ,TRIM(fname) &
!                             ,append_in=.TRUE.)

END SUBROUTINE write_modes


END MODULE MODgvec_write_modes
