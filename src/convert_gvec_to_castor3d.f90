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
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!> 
!!# **GVEC_TO_CASTOR3D**  converter program 
!!
!===================================================================================================================================
PROGRAM CONVERT_GVEC_TO_CASTOR3D
USE MODgvec_Globals
USE MODgvec_gvec_to_castor3d
!$ USE omp_lib
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
REAL(wp)                :: StartTime,EndTime
!===================================================================================================================================
  call perfinit()
  call perfon ('main')

  CALL CPU_TIME(StartTime)
!$ StartTime=OMP_GET_WTIME()
    
  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(5(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - CONVERT GVEC => CASTOR3D  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '&
,'  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
  WRITE(Unit_stdOut,'(132("="))')
  CALL GET_CLA_gvec_to_castor3d()
  
  !initialization phase
  CALL Init_gvec_to_castor3d()
 
  CALL gvec_to_castor3d_writeToFile()
  CALL Finalize_gvec_to_castor3d()

  CALL CPU_TIME(EndTime)
!$ EndTime=OMP_GET_WTIME()
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' CONVERT GVEC TO CASTOR3D FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)

  call perfoff('main')
  call perfout('main')
END PROGRAM CONVERT_GVEC_TO_CASTOR3D


