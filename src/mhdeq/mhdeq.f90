!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2015  Prof. Claus-Dieter Munz (github.com/flexi-framework/hopr)
!
! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
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
!!# Module **MHDEQ**
!!
!! For initialization of an MHD equilibrium. 
!!
!===================================================================================================================================
MODULE MOD_MHDEQ
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitMHDEQ 
  MODULE PROCEDURE InitMHDEQ
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToMHDEQ
!  MODULE PROCEDURE MapToMHDEQ
!END INTERFACE

PUBLIC::InitMHDEQ
PUBLIC::MapToMHDEQ
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE InitMHDEQ 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,abort
USE MOD_ReadInTools,ONLY:GETINT,GETREALARRAY
USE MOD_MHDEQ_Vars
USE MOD_VMEC, ONLY:InitVMEC
USE MOD_Solov, ONLY:InitSolov
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'INIT MHD EQUILIBRIUM INPUT ...'
whichEquilibrium    = GETINT('whichEquilibrium','0')   
IF(WhichEquilibrium.EQ.0) THEN 
  WRITE(UNIT_stdOut,'(A)')'... NOTHING TO BE DONE'
  RETURN
END IF
SELECT CASE(whichEquilibrium)
CASE(1)
  useMHDEQ=.TRUE.
  WRITE(*,*)'Using VMEC as equilibrium solution...'
  CALL InitVMEC()
CASE(2)
  useMHDEQ=.TRUE.
  WRITE(*,*)'Using Soloviev as equilibrium solution...'
  CALL InitSolov()
CASE DEFAULT
  WRITE(*,*)'WARNING: No Equilibrium solution for which Equilibrium= ', whichEquilibrium
  STOP
END SELECT
  !density coefficients of the polynomial coefficients: rho_1+rho_2*x + rho_3*x^2 ...
  nRhoCoefs=GETINT("nRhoCoefs","0")
  IF(nRhoCoefs.GT.0)THEN
    RhoFluxVar=GETINT("RhoFluxVar") ! dependant variable: =0: psinorm (tor. flux), =1:chinorm (pol. flux)
    ALLOCATE(RhoCoefs(nRhoCoefs))
    RhoCoefs=GETREALARRAY("RhoCoefs",nRhoCoefs)
  END IF
  InputCoordSys=GETINT("MHDEQ_inputCoordSys","0")
  ! =0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
  ! =1: x_in(1:3) are (r,zeta,theta) coordinates r= [0;1], zeta= [0;1], theta=[0;1]

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitMHDEQ


!===================================================================================================================================
!> Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from VMEC data. 
!! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!!
!===================================================================================================================================
SUBROUTINE MapToMHDEQ(nTotal,x_in,x_out,MHDEQdata)
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars, ONLY: nVarMHDEQ,whichEquilibrium,InputCoordSys
USE MOD_VMEC, ONLY:MapToVMEC
USE MOD_Solov, ONLY:MapToSolov
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,  INTENT(IN) :: nTotal         !! total number of points
REAL(wp), INTENT(IN) :: x_in(3,nTotal) !! input coordinates represent a cylinder: 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT) :: x_out(3,nTotal) !! mapped x,y,z coordinates with vmec data
REAL(wp),INTENT(OUT) :: MHDEQdata(nVarMHDEQ,nTotal) !! output data  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(whichEquilibrium)
CASE(1)
  CALL MapToVMEC(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
CASE(2)
  CALL MapToSolov(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
END SELECT
END SUBROUTINE MapToMHDEQ 

END MODULE MOD_MHDEQ
