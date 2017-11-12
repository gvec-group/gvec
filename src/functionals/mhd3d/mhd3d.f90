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
!!
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
INTEGER          :: deg,cont,degGP
INTEGER          :: mn_max(2),mn_nyq(2),nfp
CHARACTER(LEN=8) :: sincos
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'

nElems   =GETINT("sgrid_nElems","10")
grid_type=GETINT("sgrid_grid_type","0")

CALL grid_s%init(nElems,grid_type)

!mandatory global input parameters
degGP   = GETINT("degGP","4")
mn_nyq  = GETINTARRAY("mn_nyq",2,"4 1")
nfp     = GETINT("nfp","1")

X1_BC   = GETINTARRAY("X1_BC",2,"0 1")
deg     = GETINT(     "X1_sbase_deg","3")
cont    = GETINT(     "X1_sbase_continuity","2")
mn_max  = GETINTARRAY("X1_fbase_mn_max",2,"2 0")
sincos  = GETSTR(     "X1_fbase_sincos","_COS_")  !_SIN_,_COS_,_SINCOS_

CALL X1base%s%init(grid_s,deg,cont,degGP)
!TODO CALL X1base%f%init(nfp,mn_max,nyq_fac,sincos)

X2_BC   = GETINTARRAY("X2_BC",2,"0 1")
deg     = GETINT(     "X2_sbase_deg","3")
cont    = GETINT(     "X2_sbase_continuity","2")
mn_max  = GETINTARRAY("X2_fbase_mn_max",2,"2 0")
sincos  = GETSTR(     "X2_fbase_sincos","_SIN_")

CALL X2base%s%init(grid_s,deg,cont,degGP)
!TODO CALL X2base%f%init(nfp,mn_max,nyq_fac,sincos)

X2_BC   = GETINTARRAY("LA_BC",2,"0 0")
deg     = GETINT(     "LA_sbase_deg","3")
cont    = GETINT(     "LA_sbase_continuity","-1")
mn_max  = GETINTARRAY("LA_fbase_mn_max",2,"2 0")
sincos  = GETSTR(     "LA_fbase_sincos","_SIN_")

CALL LAbase%s%init(grid_s,deg,cont,degGP)
!TODO CALL LAbase%f%init(nfp,mn_max,nyq_fac,sincos)

nDOF_X1 = X1base%s%nBase* 1 ! TODO X1base%f%mn_modes
nDOF_X2 = X2base%s%nBase* 1 ! TODO X2base%f%mn_modes
nDOF_LA = LAbase%s%nBase* 1 ! TODO LAbase%f%mn_modes

DO i=-1,1
  ALLOCATE(U(i)%X1(nDOF_X1))
  ALLOCATE(U(i)%X2(nDOF_X2))
  ALLOCATE(U(i)%LA(nDOF_LA))
END DO

ALLOCATE(dUdt%X1(nDOF_X1))
ALLOCATE(dUdt%X2(nDOF_X2))
ALLOCATE(dUdt%LA(nDOF_LA))

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
