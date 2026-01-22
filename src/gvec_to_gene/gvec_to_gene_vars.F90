!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================

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
