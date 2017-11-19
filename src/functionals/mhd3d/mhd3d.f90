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
USE MOD_mhdeq_Vars, ONLY: whichInitEquilibrium
USE MOD_sgrid,      ONLY: t_sgrid
USE MOD_fbase,      ONLY: t_fbase,fbase_new
USE MOD_base,       ONLY: t_base,base_new
USE MOD_ReadInTools,ONLY: GETSTR,GETINT,GETINTARRAY,GETREAL,GETREALALLOCARRAY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i,iMode,nElems
INTEGER          :: grid_type
INTEGER          :: X1X2_deg,X1X2_cont
INTEGER          :: X1_mn_max(2),X2_mn_max(2)
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
  
  which_init = whichInitEquilibrium ! GETINT("which_init","0")
  
  !hmap
  which_hmap=GETINT("which_hmap","1")
  CALL GETREALALLOCARRAY("iota_coefs",iota_coefs,n_iota_coefs,"1.0 -0.5") !a+b*s+c*s^2...
  CALL GETREALALLOCARRAY("mass_coefs",mass_coefs,n_mass_coefs,"1.1 0.2 0.1") !a+b*s+c*s^2...
  SELECT CASE(which_init)
  CASE(0)
    
  CASE(1) !VMEC init
    which_hmap=1 !hmap_RZ
  END SELECT !which_init

  CALL hmap_new(hmap,which_hmap)
  
  
  X1X2_deg     = GETINT(     "X1X2_deg")
  WRITE(defStr,'(I4)') X1X2_deg-1
  X1X2_cont    = GETINT(     "X1X2_continuity",defStr)
  X1X2_BC      = GETINTARRAY("X1X2_BC",2,"0 1")
  X1_mn_max    = GETINTARRAY("X1_mn_max",2,"2 0")
  X2_mn_max    = GETINTARRAY("X2_mn_max",2,"2 0")
  X1_sin_cos   = GETSTR(     "X1_sin_cos","_cos_")  !_sin_,_cos_,_sin_cos_
  X2_sin_cos   = GETSTR(     "X2_sin_cos","_sin_")
  
  
  LA_deg     = GETINT(     "LA_deg")
  LA_cont    = GETINT(     "LA_continuity","-1")
  LA_BC      = GETINTARRAY("LA_BC",2,"0 0")
  LA_mn_max  = GETINTARRAY("LA_mn_max",2,"2 0")
  LA_sin_cos = GETSTR(     "LA_sin_cos","_sin_")
  
  mn_nyq(1)=MAX(1,fac_nyq*MAXVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/)))
  mn_nyq(2)=MAX(1,fac_nyq*MAXVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/)))
  
  SWRITE(UNIT_stdOut,*)
  SWRITE(UNIT_stdOut,'(A,I4,A,I6," , ",I6,A)')'    fac_nyq = ', fac_nyq,'  ==> interpolation points mn_nyq=( ',mn_nyq(:),' )'
  SWRITE(UNIT_stdOut,*)

  !INITIALIZE GRID  
  CALL sgrid%init(nElems,grid_type)

  !INITIALIZE BASE        !sbase parameter                 !fbase parameter               ...exclude_mn_zero
  CALL base_new(X1_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X1_mn_max,mn_nyq,nfp,X1_sin_cos,.FALSE.)
  CALL base_new(X2_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X2_mn_max,mn_nyq,nfp,X2_sin_cos,.FALSE.)
  CALL base_new(LA_base  ,   LA_deg,  LA_cont,sgrid,degGP , LA_mn_max,mn_nyq,nfp,LA_sin_cos,.TRUE. )

  CALL fbase_new(X1_b_base , X1_mn_max,mn_nyq,nfp,X1_sin_cos,.FALSE.)
  CALL fbase_new(X2_b_base , X2_mn_max,mn_nyq,nfp,X2_sin_cos,.FALSE.)
  CALL fbase_new(LA_b_base , LA_mn_max,mn_nyq,nfp,LA_sin_cos,.TRUE. )

  ALLOCATE(X1_b(1:X1_b_base%modes) )
  ALLOCATE(X2_b(1:X2_b_base%modes) )
  ALLOCATE(LA_b(1:LA_b_base%modes) )
  X1_b=0.0_wp
  X2_b=0.0_wp
  LA_b=0.0_wp
  
  nDOF_X1 = X1_base%s%nBase* X1_base%f%modes
  nDOF_X2 = X2_base%s%nBase* X2_base%f%modes
  nDOF_LA = LA_base%s%nBase* LA_base%f%modes
  
  ALLOCATE(U(-1:1))
  CALL U(1)%init((/nDOF_X1,nDOF_X2,nDOF_LA/))
  DO i=-1,0
    CALL U(i)%copy(U(1))
  END DO
  CALL U(1)%AXBY(0.4_wp ,U(-1),-0.25_wp,U(0))
  CALL dUdt%copy(U(1))


  !auxiliary variables
  ASSOCIATE( nGP   => (degGP+1)*sgrid%nElems,&
             mn_IP => mn_nyq(1)*mn_nyq(2)    ) !same for all variables

  ALLOCATE(mass_GP(         nGP) )
  ALLOCATE(iota_GP(         nGP) )
  ALLOCATE(PhiPrime_GP(     nGP) )
  ALLOCATE(Vprime_GP(       nGP) )
  ALLOCATE(J_h(       mn_IP,nGP) )
  ALLOCATE(J_p(       mn_IP,nGP) )
  ALLOCATE(dX1_ds(    mn_IP,nGP) )
  ALLOCATE(dX2_ds(    mn_IP,nGP) )
  ALLOCATE(dX1_dthet( mn_IP,nGP) )
  ALLOCATE(dX2_dthet( mn_IP,nGP) )
  ALLOCATE(dLA_dthet( mn_IP,nGP) )
  ALLOCATE(dX1_dzeta( mn_IP,nGP) )
  ALLOCATE(dX2_dzeta( mn_IP,nGP) )
  ALLOCATE(dLA_dzeta( mn_IP,nGP) )
  ALLOCATE(b_a (    2,mn_IP,nGP) )
  ALLOCATE(g_ab(  2,2,mn_IP,nGP) )
 
  END ASSOCIATE !mn_IP,nGP

  SELECT CASE(which_init)
  CASE(0)
    WRITE(UNIT_stdOut,'(4X,A)')'... read boundary data for X1:'
    !READ boudnary values from input file
    ASSOCIATE(modes=>X1_b_base%modes,sin_range=>X1_b_base%sin_range,cos_range=>X1_b_base%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X1_b(iMode)=get_iMode('X1_b_sin',X1_b_base%Xmn(:,iMode))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1_b(iMode)=get_iMode('X1_b_cos',X1_b_base%Xmn(:,iMode))
    END DO !iMode
    END ASSOCIATE
    WRITE(UNIT_stdOut,'(4X,A)')'... read boundary data for X2:'
    ASSOCIATE(modes=>X2_b_base%modes,sin_range=>X2_b_base%sin_range,cos_range=>X2_b_base%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X2_b(iMode)=get_iMode('X2_b_sin',X2_b_base%Xmn(:,iMode))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2_b(iMode)=get_iMode('X2_b_cos',X2_b_base%Xmn(:,iMode))
    END DO !iMode
    END ASSOCIATE
  CASE(1) !VMEC
!    X1_mn_max  = MAX( 
!    X2_mn_max  = MAX( 
!    X1_sin_cos   = ? 
!    X2_sin_cos   = ?
!  LA_mn_max  = 
!  LA_sin_cos = ? 
  END SELECT !which_init
  CALL InitializeSolution(U(-1))
  
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)
  
  CONTAINS 
  FUNCTION get_iMode(varname_in,mn)
    IMPLICIT NONE
    !-------------------------------------------
    !input/output
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
    INTEGER         ,INTENT(IN) :: mn(2)
    REAL(wp)                    :: get_iMode
    !local variables
    CHARACTER(LEN=100) :: varstr
    !-------------------------------------------
    WRITE(varstr,'(A,"("I4,";",I4,")")')TRIM(varname_in),mn(:)
    varstr=delete_spaces(varstr)
    get_iMode=GETREAL(TRIM(varstr),"0.0")
   
  END FUNCTION get_iMode

  FUNCTION delete_spaces(str_in)
    USE ISO_VARYING_STRING,ONLY:VARYING_STRING,VAR_STR,CHAR,replace
    IMPLICIT NONE
    !-------------------------------------------
    !input/output
    CHARACTER(LEN=*),INTENT(IN) :: str_in
    CHARACTER(LEN=LEN(str_in))  :: delete_spaces
    TYPE(VARYING_STRING) :: tmpstr
    !-------------------------------------------
    tmpstr=VAR_STR(str_in);delete_spaces=CHAR(Replace(tmpstr," ","",Every=.true.))
   
  END FUNCTION delete_spaces

END SUBROUTINE InitMHD3D

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE InitializeSolution(U0) 
! MODULES
USE MOD_MHD3D_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_sol_var), INTENT(INOUT) :: U0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(which_init)
CASE(0)
CASE(1)
END SELECT 

END SUBROUTINE InitializeSolution

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
  CALL X1_base%free()
  CALL X2_base%free()
  CALL LA_base%free()
  CALL X1_b_base%free()
  CALL X2_b_base%free()
  CALL LA_b_base%free()
  
  DO i=-1,1
    CALL U(i)%free()
  END DO
  CALL dUdt%free()
  CALL sgrid%free()

  SDEALLOCATE(mass_GP     )
  SDEALLOCATE(iota_GP     )
  SDEALLOCATE(PhiPrime_GP )
  SDEALLOCATE(Vprime_GP   )
  SDEALLOCATE(J_h         )
  SDEALLOCATE(J_p         )
  SDEALLOCATE(dX1_ds      )
  SDEALLOCATE(dX2_ds      )
  SDEALLOCATE(dX1_dthet   )
  SDEALLOCATE(dX2_dthet   )
  SDEALLOCATE(dLA_dthet   )
  SDEALLOCATE(dX1_dzeta   )
  SDEALLOCATE(dX2_dzeta   )
  SDEALLOCATE(dLA_dzeta   )
  SDEALLOCATE(b_a         )
  SDEALLOCATE(g_ab        )

END SUBROUTINE FinalizeMHD3D

END MODULE MOD_MHD3D
