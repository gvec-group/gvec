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
MODULE MOD_MHD3D
! MODULES
USE MOD_Globals, ONLY:wp,UNIT_stdOut,fmt_sep
USE MOD_c_functional,   ONLY: t_functional
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE,EXTENDS(t_functional) :: t_functional_mhd3d
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL :: initialized
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS
    PROCEDURE :: init => InitMHD3D
    PROCEDURE :: free => FinalizeMHD3D
END TYPE t_functional_mhd3d

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitMHD3D(sf) 
! MODULES
USE MOD_Globals,ONLY:PI
USE MOD_MHD3D_Vars
USE MOD_sgrid, ONLY: t_sgrid
USE MOD_base,  ONLY: t_base,base_new
USE MOD_ReadInTools,ONLY:GETSTR,GETINT,GETINTARRAY,GETREAL,GETREALALLOCARRAY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i,nElems
INTEGER          :: grid_type
INTEGER          :: X1X2_deg,X1X2_cont,X1X2_mn_max(2)
INTEGER          :: LA_deg,LA_cont,LA_mn_max(2)
CHARACTER(LEN=8) :: X1_sin_cos
CHARACTER(LEN=8) :: X2_sin_cos
CHARACTER(LEN=8) :: LA_sin_cos
CHARACTER(LEN=8) :: defstr
INTEGER          :: degGP,mn_nyq(2),fac_nyq
INTEGER          :: which_hmap 
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'
  
  nElems   =GETINT("sgrid_nElems","10")
  grid_type=GETINT("sgrid_grid_type","0")
  
  !mandatory global input parameters
  degGP   = GETINT( "degGP","4")
  fac_nyq = GETINT( "fac_nyq","4")
  nfp     = GETINT( "nfp","1")
  
  !constants
  gamm    = GETREAL("GAMMA","0.")
  
  mu_0    = 4.0e-07_wp*PI
  
  which_init = GETINT("which_init","0")
  
  !hmap
  SELECT CASE(which_init)
  CASE(0)
    which_hmap=GETINT("which_hmap","1")
    CALL GETREALALLOCARRAY("iota_coefs",iota_coefs,n_iota_coefs,"1.0 -0.5") !a+b*s+c*s^2...
    CALL GETREALALLOCARRAY("mass_coefs",mass_coefs,n_mass_coefs,"1.1 0.2 0.1") !a+b*s+c*s^2...
   
!  CASE(1) !VMEC init
!    which_hmap=1 !hmap_RZ
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "the choice for which_init input parameter does not exist!")
  END SELECT
  CALL hmap_new(hmap,which_hmap)
  
  X1X2_BC   = GETINTARRAY(   "X1X2_BC",2,"0 1")
  
  X1X2_deg     = GETINT(     "X1X2_deg","3")
  WRITE(defStr,'(I4)') X1X2_deg-1
  X1X2_cont    = GETINT(     "X1X2_continuity",defStr)
  X1X2_mn_max  = GETINTARRAY("X1X2_mn_max",2,"2 0")

  SELECT CASE(which_init)
  CASE(0)
    X1_sin_cos    = GETSTR(     "X1_sin_cos","_cos_")  !_SIN_,_COS_,_sin_cos_
    X2_sin_cos    = GETSTR(     "X2_sin_cos","_sin_")
!  CASE(1) !VMEC
  END SELECT !which_init
  
  
  LA_BC   = GETINTARRAY(   "LA_BC",2,"0 0")
  LA_deg     = GETINT(     "LA_deg","3")
  LA_cont    = GETINT(     "LA_continuity","-1")
  LA_mn_max  = GETINTARRAY("LA_mn_max",2,"2 0")
  SELECT CASE(which_init)
  CASE(0)
    LA_sin_cos  = GETSTR(     "LA_sin_cos","_sin_")
!  CASE(1) !VMEC
  END SELECT !which_init
  
  mn_nyq(1)=MAX(1,fac_nyq*MAX(X1X2_mn_max(1),LA_mn_max(1)))
  mn_nyq(2)=MAX(1,fac_nyq*MAX(X1X2_mn_max(2),LA_mn_max(2)))
  
  SWRITE(UNIT_stdOut,*)
  SWRITE(UNIT_stdOut,'(A,I4,A,I6," , ",I6,A)')'    fac_nyq = ', fac_nyq,'  ==> interpolation points mn_nyq=( ',mn_nyq(:),' )'
  SWRITE(UNIT_stdOut,*)

  !INITIALIZE GRID  
  CALL sgrid%init(nElems,grid_type)

  !INITIALIZE BASE        !sbase parameter                 !fbase parameter               ...exclude_mn_zero
  CALL base_new(X1base  , X1X2_deg,X1X2_cont,sgrid,degGP , X1X2_mn_max,mn_nyq,nfp,X1_sin_cos,.FALSE.)
  CALL base_new(X2base  , X1X2_deg,X1X2_cont,sgrid,degGP , X1X2_mn_max,mn_nyq,nfp,X2_sin_cos,.FALSE.)
  CALL base_new(LAbase  ,   LA_deg,  LA_cont,sgrid,degGP ,   LA_mn_max,mn_nyq,nfp,LA_sin_cos,.TRUE. )
  
  nDOF_X1 = X1base%s%nBase* X1base%f%modes
  nDOF_X2 = X2base%s%nBase* X2base%f%modes
  nDOF_LA = LAbase%s%nBase* LAbase%f%modes
  
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
SUBROUTINE FinalizeMHD3D(sf) 
! MODULES
USE MOD_MHD3D_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!===================================================================================================================================
  CALL X1base%free()
  CALL X2base%free()
  CALL LAbase%free()
  
  DO i=-1,1
    CALL U(i)%free()
  END DO
  CALL dUdt%free()
  CALL sgrid%free()

END SUBROUTINE FinalizeMHD3D

END MODULE MOD_MHD3D
