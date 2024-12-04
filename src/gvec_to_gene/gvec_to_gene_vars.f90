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
!=================================================================================================================================

!===================================================================================================================================
!>
!!# Module ** gvec_to_gene Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_gene_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
USE MODgvec_transform_sfl     ,ONLY: t_transform_sfl
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER       :: SFLcoord            !! =0: 'old way' of PEST with newton iteration, =1: PEST, =2: Boozer
INTEGER       :: factorSFL           !! factor of the SFL coordinate mode numbers over the number of GVEC modes in X1/X2/LA
TYPE(t_transform_sfl),ALLOCATABLE :: trafoSFL
                                         
!===================================================================================================================================
END MODULE MODgvec_gvec_to_gene_Vars

