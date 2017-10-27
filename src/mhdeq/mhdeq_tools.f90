!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **MHD tools**
!!
!! tools for the MHDEQ module 
!!
!===================================================================================================================================
MODULE MOD_MHDEQ_Tools
! MODULES
USE MOD_Globals,ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE Eval1DPoly 
  MODULE PROCEDURE Eval1DPoly
END INTERFACE

PUBLIC::Eval1DPoly
!===================================================================================================================================

CONTAINS
 

!===================================================================================================================================
!> evalute monomial polynomial c_1+c_2*x+c_3*x^2 ...
!!
!===================================================================================================================================
FUNCTION Eval1DPoly(nCoefs,Coefs,x)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,  INTENT(IN)  :: nCoefs                   !! number of coefficients 
REAL(wp), INTENT(IN)  :: Coefs(nCoefs)            !! coefficients
REAL(wp), INTENT(IN)  :: x                        !! evaluation position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)              :: Eval1DPoly
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!===================================================================================================================================
Eval1DPoly=0.
DO i=nCoefs,1,-1
  Eval1DPoly=Eval1DPoly*x+Coefs(i)
END DO

END FUNCTION Eval1DPoly

END MODULE MOD_MHDEQ_Tools
