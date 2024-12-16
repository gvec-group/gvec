!================================================================================================================================!
! Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!================================================================================================================================!
! PyGVEC postprocessing
!
! This module contains the Python interface for postprocessing the GVEC output, starting with a parameter file and a statefile.
!================================================================================================================================!

#include "defines.h"

MODULE MODgvec_py_run

IMPLICIT NONE
PUBLIC

CONTAINS

!================================================================================================================================!
SUBROUTINE start_rungvec(parameterfile,restartfile_in,comm_in)
  ! MODULES
  USE MODgvec_Globals, ONLY: Unit_stdOut
  USE MODgvec_MPI    ,ONLY  : par_Init,par_finalize
  USE MODgvec_rungvec, ONLY: rungvec
  ! INPUT/OUTPUT VARIABLES ------------------------------------------------------------------------------------------------------!
  CHARACTER(LEN=*),INTENT(IN) :: parameterfile
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: restartfile_in
  INTEGER,INTENT(IN),OPTIONAL :: comm_in 
  ! LOCAL VARIABLES -------------------------------------------------------------------------------------------------------------!
  INTEGER :: comm
  ! CODE ------------------------------------------------------------------------------------------------------------------------!
  IF(.NOT.PRESENT(comm_in))THEN 
    CALL par_init() !USE MPI_COMM_WORLD
    IF(PRESENT(restartfile_in))THEN
      CALL rungvec(parameterfile,restartfile_in=restartfile_in)
    ELSE                                        
      CALL rungvec(parameterfile)
    END IF
    CALL par_finalize()
  ELSE
    IF(PRESENT(restartfile_in))THEN
      CALL rungvec(parameterfile,restartfile_in=restartfile_in,comm_in=comm_in)
    ELSE                                        
      CALL rungvec(parameterfile,comm_in=comm_in)
    END IF
  END IF
END SUBROUTINE start_rungvec

END MODULE MODgvec_py_run
