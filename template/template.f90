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
!!# Module **TEMPLATE**
!!
!!
!!
!===================================================================================================================================
MODULE MOD_TEMPLATE
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PRIVATE

INTERFACE InitTEMPLATE
  MODULE PROCEDURE InitTEMPLATE
END INTERFACE

INTERFACE TEMPLATE
  MODULE PROCEDURE TEMPLATE
END INTERFACE

INTERFACE FinalizeTEMPLATE
  MODULE PROCEDURE FinalizeTEMPLATE
END INTERFACE

PUBLIC::InitTEMPLATE
PUBLIC::TEMPLATE
PUBLIC::FinalizeTEMPLATE
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitTEMPLATE 
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MOD_TEMPLATE_Vars
USE MOD_ReadInTools,ONLY:GETLOGICAL
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
SWRITE(UNIT_stdOut,'(A)')'INIT TEMPLATE ...'
useThis    = GETLOGICAL('useThis','F')   

SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitTEMPLATE


!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE TEMPLATE()
! MODULES
USE MOD_Globals, ONLY:wp
USE MOD_TEMPLATE_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
END SUBROUTINE TEMPLATE 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeTEMPLATE 
! MODULES
USE MOD_TEMPLATE_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE FinalizeTEMPLATE

END MODULE MOD_TEMPLATE
