!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **MHD3D**
!!
!! CONTAINS INITIALIZATION OF MHD 3D Energy functional that will be minimized
!!
!===================================================================================================================================
MODULE MOD_Functional
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitFunctional
  MODULE PROCEDURE InitMHD3D
END INTERFACE

INTERFACE FinalizeFunctional
  MODULE PROCEDURE FinalizeMHD3D
END INTERFACE

PUBLIC::InitFunctional
PUBLIC::FinalizeFunctional
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitMHD3D 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MOD_MHD3D_Vars
USE MOD_sgrid, ONLY: t_sgrid
USE MOD_sbase, ONLY: t_sbase
USE MOD_ReadInTools,ONLY:GETSTR,GETINT,GETINTARRAY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i,nElems
INTEGER          :: grid_type
INTEGER          :: X1_deg,X1_cont,X1_mn_max(2)
CHARACTER(LEN=8) :: X1_sincos
INTEGER          :: X2_deg,X2_cont,X2_mn_max(2)
CHARACTER(LEN=8) :: X2_sincos
INTEGER          :: LA_deg,LA_cont,LA_mn_max(2)
CHARACTER(LEN=8) :: LA_sincos
INTEGER          :: degGP,mn_nyq(2),fac_nyq,nfp
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'

nElems   =GETINT("sgrid_nElems","10")
grid_type=GETINT("sgrid_grid_type","0")

!mandatory global input parameters
degGP   = GETINT("degGP","4")
fac_nyq  = GETINT("fac_nyq","4")
nfp     = GETINT("nfp","1")

X1_BC   = GETINTARRAY(   "X1_BC",2,"0 1")
X1_deg     = GETINT(     "X1_deg","3")
X1_cont    = GETINT(     "X1_continuity","2")
X1_mn_max  = GETINTARRAY("X1_mn_max",2,"2 0")
X1_sincos  = GETSTR(     "X1_sincos","_COS_")  !_SIN_,_COS_,_SINCOS_


X2_BC   = GETINTARRAY(   "X2_BC",2,"0 1")
X2_deg     = GETINT(     "X2_deg","3")
X2_cont    = GETINT(     "X2_continuity","2")
X2_mn_max  = GETINTARRAY("X2_mn_max",2,"2 0")
X2_sincos  = GETSTR(     "X2_sincos","_SIN_")


LA_BC   = GETINTARRAY(   "LA_BC",2,"0 0")
LA_deg     = GETINT(     "LA_deg","3")
LA_cont    = GETINT(     "LA_continuity","-1")
LA_mn_max  = GETINTARRAY("LA_mn_max",2,"2 0")
LA_sincos  = GETSTR(     "LA_sincos","_SIN_")

mn_nyq(1)=MAX(1,MAXVAL(fac_nyq*(/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/)))
mn_nyq(2)=MAX(1,MAXVAL(fac_nyq*(/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/)))
SWRITE(*,'(A,I4,A,I6," , ",I6)')'fac_nyq = ', fac_nyq,' ==> interpolation points mIP,nIP',mn_nyq(:)

ALLOCATE(t_sgrid :: sgrid)
CALL sgrid%init(nElems,grid_type)

ALLOCATE(t_base :: X1base)
ALLOCATE(t_base :: X2base)
ALLOCATE(t_base :: LAbase)
ALLOCATE(t_sbase :: X1base%s)
ALLOCATE(t_sbase :: X2base%s)
ALLOCATE(t_sbase :: LAbase%s)
CALL X1base%s%init(sgrid,X1_deg,X1_cont,degGP)
CALL X2base%s%init(sgrid,X2_deg,X2_cont,degGP)
CALL LAbase%s%init(sgrid,LA_deg,LA_cont,degGP)

!ALLOCATE(t_fbase :: X1base%f)
!ALLOCATE(t_fbase :: X2base%f)
!ALLOCATE(t_fbase :: LAbase%f)
!TODO CALL X1base%f%init(nfp,X1_mn_max,mn_nyq,X1_sincos)
!TODO CALL X1base%f%init(nfp,X1_mn_max,mn_nyq,X1_sincos)
!TODO CALL X1base%f%init(nfp,X1_mn_max,mn_nyq,X1_sincos)

nDOF_X1 = X1base%s%nBase* 1 ! TODO X1base%f%mn_modes
nDOF_X2 = X2base%s%nBase* 1 ! TODO X2base%f%mn_modes
nDOF_LA = LAbase%s%nBase* 1 ! TODO LAbase%f%mn_modes

ALLOCATE(t_sol_var :: U(-1:1))
DO i=-1,1
  CALL U(i)%init((/nDOF_X1,nDOF_X2,nDOF_LA/))
END DO
ALLOCATE(t_sol_var :: dUdt)
CALL dudt%init((/nDOF_X1,nDOF_X2,nDOF_LA/))

SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitMHD3D


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeMHD3D 
! MODULES
USE MOD_MHD3D_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE FinalizeMHD3D

END MODULE MOD_Functional
