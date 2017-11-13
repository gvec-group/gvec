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
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MOD_MHD3D_Vars
! MODULES
USE MOD_Globals,ONLY:wp
USE MOD_sgrid, ONLY: c_sgrid,t_sgrid
USE MOD_sbase, ONLY: c_sbase,t_sbase
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE :: t_base              
  CLASS(t_sbase),ALLOCATABLE  :: s  !! container for radial basis
!  CLASS(t_fbase),ALLOCATABLE  :: f  !! container for angular basis
END TYPE t_base


TYPE,ABSTRACT :: c_sol_var
  LOGICAL               :: initialized
  CONTAINS
  PROCEDURE(i_sub_sol_var_init  ),DEFERRED :: init
  PROCEDURE(i_sub_sol_var_free  ),DEFERRED :: free
  PROCEDURE(i_sub_sol_var_copy  ),DEFERRED :: copy
END TYPE c_sol_var

ABSTRACT INTERFACE
  SUBROUTINE i_sub_sol_var_init( sf ,size_sf)
    IMPORT c_sol_var
    CLASS(c_sol_var), INTENT(INOUT) :: sf
    INTEGER         , INTENT(IN   ) :: size_sf(3)
  END SUBROUTINE i_sub_sol_var_init

  SUBROUTINE i_sub_sol_var_free( sf ) 
    IMPORT c_sol_var
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var_free

  SUBROUTINE i_sub_sol_var_copy( sf, tocopy ) 
    IMPORT c_sol_var
    CLASS(c_sol_var), INTENT(IN   ) :: tocopy
    CLASS(c_sol_var), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_sol_var_copy

END INTERFACE


TYPE,EXTENDS(c_sol_var) :: t_sol_var
  REAL(wp) ,ALLOCATABLE :: X1(:)    !! X1 variable, size (base_f%mn_mode*base_s%nBase)
  REAL(wp) ,ALLOCATABLE :: X2(:)    !! X2 variable 
  REAL(wp) ,ALLOCATABLE :: LA(:)    !! lambda variable
  CONTAINS

  PROCEDURE        :: init   => sol_var_init
  PROCEDURE        :: free   => sol_var_free
  PROCEDURE        :: copy   => sol_var_copy
!  PROCEDURE        :: norm_2 => sol_var_norm_2  
!  PROCEDURE,NOPASS :: AXBY   => sol_var_AXBY
END TYPE t_sol_var

!-----------------------------------------------------------------------------------------------------------------------------------
! SOLUTION VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CLASS(t_base),ALLOCATABLE     :: X1base  
CLASS(t_base),ALLOCATABLE     :: X2base  
CLASS(t_base),ALLOCATABLE     :: LAbase  
                             
CLASS(t_sgrid),ALLOCATABLE    :: sgrid  !! only one grid up to now

CLASS(t_sol_var),ALLOCATABLE  :: U(:)         !! solutions at levels (k-1),(k),(k+1)
CLASS(t_sol_var),ALLOCATABLE  :: dUdt            !! solution update

INTEGER          :: X1_BC(2)        !! BC axis (0) and edge (1)    
INTEGER          :: X2_BC(2)        !! BC axis (0) and edge (1)    
INTEGER          :: LA_BC(2)        !! BC axis (0) and edge (1)    
INTEGER          :: nDOF_X1         !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_X2         !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 
INTEGER          :: nDOF_LA         !! total number of degrees of freedom, sBase%nBase * fbase%mn_modes 

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> allocate t_base type 
!!
!===================================================================================================================================
SUBROUTINE base_new( base)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base),ALLOCATABLE,INTENT(INOUT) :: base
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(t_base :: base)
  ALLOCATE(t_sbase :: base%s)
!  ALLOCATE(t_fbase :: base%f)

END SUBROUTINE base_new

!===================================================================================================================================
!> initialize (=allocate) sf of type t_sol_var 
!!
!===================================================================================================================================
SUBROUTINE sol_var_init( sf,size_sf)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var), INTENT(INOUT) :: sf  !!sf
  INTEGER         , INTENT(IN   ) :: size_sf(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(sf%initialized)THEN
    WRITE(*,*)'WARNING: REALLOCATION OF SOL_VAR!!'
    CALL sf%free()
  END IF
  ALLOCATE(sf%X1(size_sf(1)))
  ALLOCATE(sf%X2(size_sf(2)))
  ALLOCATE(sf%LA(size_sf(3)))
  sf%initialized=.TRUE.

END SUBROUTINE sol_var_init

!===================================================================================================================================
!> free (=deallocate) sf of type t_sol_var 
!!
!===================================================================================================================================
SUBROUTINE sol_var_free( sf )
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var), INTENT(INOUT) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized)RETURN
  DEALLOCATE(sf%X1)
  DEALLOCATE(sf%X2)
  DEALLOCATE(sf%LA)
  sf%initialized=.FALSE.
END SUBROUTINE sol_var_free

!===================================================================================================================================
!> copy tocopy  => sf
!!
!===================================================================================================================================
SUBROUTINE sol_var_copy( sf ,tocopy) 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_sol_var), INTENT(IN   ) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var), INTENT(INOUT) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPEIS(t_sol_var)
  IF(.NOT.sf%initialized)STOP 'copying to sol_var that is not initialized'
  IF(.NOT.tocopy%initialized)STOP 'copying from sol_var that is not initialized'
  sf%X1=tocopy%X1
  sf%X2=tocopy%X2
  sf%LA=tocopy%LA

  END SELECT !TYPE
END SUBROUTINE sol_var_copy

!===================================================================================================================================
!> |X|^2, where X is of type t_var_sol, so three values are returned: |X1|^2,|X2|^2,|LA|^2 
!!
!===================================================================================================================================
FUNCTION sol_var_norm_2( sf ) RESULT(norm_2)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var), INTENT(IN   ) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: norm_2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized)STOP 'taking the norm of a sol_var that is not initialized'
  norm_2(1)=SUM(sf%X1*sf%X1)
  norm_2(2)=SUM(sf%X2*sf%X2)
  norm_2(3)=SUM(sf%LA*sf%LA)
END FUNCTION sol_var_norm_2

!===================================================================================================================================
!> res=a*X+b*Y , where X,Y,res are of type t_var_sol
!!
!===================================================================================================================================
SUBROUTINE sol_var_AXBY(res,a,X,b,Y)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  TYPE(t_sol_var), INTENT(IN   ) :: X
  TYPE(t_sol_var), INTENT(IN   ) :: Y
  REAL(wp)       , INTENT(IN   )  :: a,b 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_sol_var), INTENT(  OUT) :: res
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.X%initialized) STOP 'AXBY: X not initialized'
  IF(.NOT.Y%initialized) STOP 'AXBY: Y not initialized'
  res%X1 = a*X%X1 + b*Y%X1
  res%X2 = a*X%X2 + b*Y%X2
  res%LA = a*X%LA + b*Y%LA
END SUBROUTINE sol_var_AXBY

END MODULE MOD_MHD3D_Vars

