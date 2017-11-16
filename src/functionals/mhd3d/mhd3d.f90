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
USE MOD_sbase, ONLY: t_sbase,sbase_new
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
INTEGER          :: X1X2_deg,X1X2_cont,X1X2_mn_max(2)
INTEGER          :: LA_deg,LA_cont,LA_mn_max(2)
CHARACTER(LEN=8) :: X1_sincos
CHARACTER(LEN=8) :: X2_sincos
CHARACTER(LEN=8) :: LA_sincos
CHARACTER(LEN=8) :: defstr
INTEGER          :: degGP,mn_nyq(2),fac_nyq,nfp
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'

nElems   =GETINT("sgrid_nElems","10")
grid_type=GETINT("sgrid_grid_type","0")

!mandatory global input parameters
degGP   = GETINT("degGP","4")
fac_nyq  = GETINT("fac_nyq","4")
nfp     = GETINT("nfp","1")

X1X2_BC   = GETINTARRAY(   "X1X2_BC",2,"0 1")

X1X2_deg     = GETINT(     "X1X2_deg","3")
WRITE(defStr,'(I4)') X1X2_deg-1
X1X2_cont    = GETINT(     "X1X2_continuity",defStr)
X1X2_mn_max  = GETINTARRAY("X1X2_mn_max",2,"2 0")
X1_sincos    = GETSTR(     "X1_sincos","_COS_")  !_SIN_,_COS_,_SINCOS_
X2_sincos    = GETSTR(     "X2_sincos","_SIN_")


LA_BC   = GETINTARRAY(   "LA_BC",2,"0 0")
LA_deg     = GETINT(     "LA_deg","3")
LA_cont    = GETINT(     "LA_continuity","-1")
LA_mn_max  = GETINTARRAY("LA_mn_max",2,"2 0")
LA_sincos  = GETSTR(     "LA_sincos","_SIN_")

mn_nyq(1)=MAX(1,fac_nyq*MAX(X1X2_mn_max(1),LA_mn_max(1)))
mn_nyq(2)=MAX(1,fac_nyq*MAX(X1X2_mn_max(2),LA_mn_max(2)))

SWRITE(UNIT_stdOut,*)
SWRITE(UNIT_stdOut,'(A,I4,A,I6," , ",I6,A)')'    fac_nyq = ', fac_nyq,'  ==> interpolation points mn_nyq=( ',mn_nyq(:),' )'
SWRITE(UNIT_stdOut,*)

CALL sgrid%init(nElems,grid_type)

CALL sbase_new(X1base%s,X1X2_deg,X1X2_cont)
CALL sbase_new(X2base%s,X1X2_deg,X1X2_cont)
CALL sbase_new(LAbase%s,LA_deg,LA_cont)


CALL X1base%s%init(sgrid,degGP)
CALL X2base%s%copy(X1base%s)
CALL LAbase%s%init(sgrid,degGP)

!TODO CALL X1base%f%init(nfp,X1_mn_max,mn_nyq,X1_sincos)
!TODO CALL X1base%f%init(nfp,X1_mn_max,mn_nyq,X1_sincos)
!TODO CALL X1base%f%init(nfp,X1_mn_max,mn_nyq,X1_sincos)

nDOF_X1 = X1base%s%nBase* 1 ! TODO X1base%f%mn_modes
nDOF_X2 = X2base%s%nBase* 1 ! TODO X2base%f%mn_modes
nDOF_LA = LAbase%s%nBase* 1 ! TODO LAbase%f%mn_modes

ALLOCATE(U(-1:1))
CALL U(1)%init((/nDOF_X1,nDOF_X2,nDOF_LA/))
DO i=-1,0
  CALL U(i)%copy(U(1))
END DO
CALL U(1)%AXBY(0.4_wp ,U(-1),-0.25_wp,U(0))
CALL dUdt%copy(U(1))

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
