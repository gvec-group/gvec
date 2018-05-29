!===================================================================================================================================
! Copyright (C) 2018  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2018  Maurice Maurer <maurice_maurer@gmx.de>
! Copyright (C) 2018  Alejandro Banon Navarro <abanonna@ipp.mpg.de>
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
!!# **GVEC TO GENE** Driver program 
!!
!===================================================================================================================================
PROGRAM TEST_GVEC_TO_GENE
USE MODgvec_Globals
USE MODgvec_gvec_to_gene
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: nArgs
CHARACTER(LEN=255)      :: filename 
REAL(wp)                :: StartTime,EndTime
REAL(wp)                :: Fa,minor_r 
INTEGER                 :: n0_global
!===================================================================================================================================
  CALL CPU_TIME(StartTime)
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.GE.1)THEN
    CALL GET_COMMAND_ARGUMENT(1,filename)
  ELSE
    STOP 'GVEC_TO_GENE: gvec filename not given, usage: "./executable gvec_file.dat"'
  END IF
    
  
  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(5(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - -  GVEC ==> GENE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - -  GVEC ==> GENE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - -  GVEC ==> GENE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '
  WRITE(Unit_stdOut,'(132("="))')
  
  !initialization phase
  CALL Init_gvec_to_gene(filename)
 
  CALL gvec_to_gene_scalars(Fa,minor_r,n0_global)
  WRITE(*,*)'Fa',Fa
  WRITE(*,*)'minor_r',minor_r
  WRITE(*,*)'n0_global',n0_global

  CALL Finalize_gvec_to_gene()

  CALL CPU_TIME(EndTime)
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' GVEC_TO_GENE FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)

END PROGRAM TEST_GVEC_TO_GENE


