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
USE MOD_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep
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
USE MOD_VMEC_Readin,ONLY: nfp
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
INTEGER          :: nfp_loc,which_hmap 
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'
  
  nElems   =GETINT("sgrid_nElems","10")
  grid_type=GETINT("sgrid_grid_type","0")
  
  !mandatory global input parameters
  degGP   = GETINT( "degGP","4")
  fac_nyq = GETINT( "fac_nyq","4")
  
  !constants
  gamm    = GETREAL("GAMMA","0.")
  
  mu_0    = 4.0e-07_wp*PI
  
  which_init = whichInitEquilibrium ! GETINT("which_init","0")
  
  SELECT CASE(which_init)
  CASE(0)
    nfp_loc  = GETINT( "nfp","1")
    !hmap
    which_hmap=GETINT("which_hmap","1")
    CALL GETREALALLOCARRAY("iota_coefs",iota_coefs,n_iota_coefs,"1.1 0.1") !a+b*s+c*s^2...
    CALL GETREALALLOCARRAY("mass_coefs",mass_coefs,n_mass_coefs,"1.0 -0.9") !a+b*s+c*s^2...
  CASE(1) !VMEC init
    nfp_loc = nfp
    !hmap
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
  CALL base_new(X1_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X1_mn_max,mn_nyq,nfp_loc,X1_sin_cos,.FALSE.)
  CALL base_new(X2_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X2_mn_max,mn_nyq,nfp_loc,X2_sin_cos,.FALSE.)
  CALL base_new(LA_base  ,   LA_deg,  LA_cont,sgrid,degGP , LA_mn_max,mn_nyq,nfp_loc,LA_sin_cos,.TRUE. )

  nDOF_X1 = X1_base%s%nBase* X1_base%f%modes
  nDOF_X2 = X2_base%s%nBase* X2_base%f%modes
  nDOF_LA = LA_base%s%nBase* LA_base%f%modes

  ALLOCATE(X1_b(1:X1_base%f%modes) )
  ALLOCATE(X2_b(1:X2_base%f%modes) )
  ALLOCATE(LA_b(1:LA_base%f%modes) )
  ALLOCATE(X1_a(1:X1_base%f%modes) )
  ALLOCATE(X2_a(1:X2_base%f%modes) )
  X1_b=0.0_wp
  X2_b=0.0_wp
  LA_b=0.0_wp
  X1_a=0.0_wp
  X2_a=0.0_wp

  SELECT CASE(which_init)
  CASE(0)
    !READ boudnary values from input file
    WRITE(UNIT_stdOut,'(4X,A)')'... read axis boundary data for X1:'
    ASSOCIATE(modes=>X1_base%f%modes,sin_range=>X1_base%f%sin_range,cos_range=>X1_base%f%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X1_a(iMode)=get_iMode('X1_a_sin',X1_base%f%Xmn(:,iMode))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1_a(iMode)=get_iMode('X1_a_cos',X1_base%f%Xmn(:,iMode))
    END DO !iMode
    WRITE(UNIT_stdOut,'(4X,A)')'... read edge boundary data for X1:'
    DO iMode=sin_range(1)+1,sin_range(2)
      X1_b(iMode)=get_iMode('X1_b_sin',X1_base%f%Xmn(:,iMode))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1_b(iMode)=get_iMode('X1_b_cos',X1_base%f%Xmn(:,iMode))
    END DO !iMode
    END ASSOCIATE
    WRITE(UNIT_stdOut,'(4X,A)')'... read axis boundary data for X2:'
    ASSOCIATE(modes=>X2_base%f%modes,sin_range=>X2_base%f%sin_range,cos_range=>X2_base%f%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X2_a(iMode)=get_iMode('X2_a_sin',X2_base%f%Xmn(:,iMode))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2_a(iMode)=get_iMode('X2_a_cos',X2_base%f%Xmn(:,iMode))
    END DO !iMode
    WRITE(UNIT_stdOut,'(4X,A)')'... read edge boundary data for X2:'
    DO iMode=sin_range(1)+1,sin_range(2)
      X2_b(iMode)=get_iMode('X2_b_sin',X2_base%f%Xmn(:,iMode))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2_b(iMode)=get_iMode('X2_b_cos',X2_base%f%Xmn(:,iMode))
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

  ALLOCATE(U(-1:1))
  CALL U(1)%init((/X1_base%s%nbase,X2_base%s%nbase,LA_base%s%nBase,   &
                   X1_base%f%modes,X2_base%f%modes,LA_base%f%modes/)  )
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

  CALL InitializeSolution(U(-1) )
  CALL U(0)%set_to(U(-1))
  
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
    WRITE(varstr,'(A,"("I4,";",I4,")")')TRIM(varname_in),mn(1),mn(2)/nfp_loc
    varstr=delete_spaces(varstr)         !quiet on default=0.0
    get_iMode=GETREAL(TRIM(varstr),"0.0",.TRUE.)
   
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
!> evaluate iota at position s
!!
!===================================================================================================================================
FUNCTION Eval_iota(spos)
! MODULES
USE MOD_Globals, ONLY:EVAL1DPOLY
USE MOD_MHD3D_Vars,ONLY: which_init,n_iota_coefs,iota_coefs
USE MOD_VMEC , ONLY:VMEC_EvalSpl
USE MOD_VMEC_vars , ONLY:iota_spl
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos(:) !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_iota(size(spos,1))
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
SELECT CASE(which_init)
CASE(0)
  DO i=1,size(spos,1)
    eval_iota(i)=Eval1DPoly(n_iota_coefs,iota_coefs,spos(i))
  END DO
CASE(1)
!  eval_iota=VMEC_EvalSpl(1,spos,chi_Spl)/VMEC_EvalSpl(1,spos,Phi_spl) 
  eval_iota=VMEC_EvalSpl(1,spos,iota_Spl)
END SELECT
END FUNCTION Eval_iota

!===================================================================================================================================
!> evaluate pressure profile at position s
!!
!===================================================================================================================================
FUNCTION Eval_pres(spos)
! MODULES
USE MOD_Globals, ONLY:EVAL1DPOLY
USE MOD_MHD3D_Vars,ONLY: which_init,gamm,n_mass_coefs,mass_coefs
USE MOD_VMEC , ONLY:VMEC_EvalSpl
USE MOD_VMEC_vars , ONLY:pres_spl
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos(:) !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_pres(size(spos,1))
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
SELECT CASE(which_init)
CASE(0)
  IF(ABS(gamm).LT.1.0E-12)THEN
    DO i=1,size(spos,1)
      eval_pres(i)=Eval1DPoly(n_mass_coefs,mass_coefs,spos(i)) 
    END DO
  ELSE
    !TODO still need to divide by V'^gamma
    CALL abort(__STAMP__, &
        ' gamma>0: pressure profile not yet implemented!')
  END IF
CASE(1)
  eval_pres=VMEC_EvalSpl(0,spos,pres_Spl)
END SELECT
END FUNCTION Eval_pres

!===================================================================================================================================
!> Initialize the solution with the given boundary condition 
!!
!===================================================================================================================================
SUBROUTINE InitializeSolution(U0)! X1_in,X2_in,LA_in) 
! MODULES
USE MOD_Globals, ONLY:EVAL1DPOLY
USE MOD_MHD3D_Vars
USE MOD_sol_var_MHD3D, ONLY:t_sol_var_MHD3D
USE MOD_lambda_solve,  ONLY:lambda_solve
USE MOD_VMEC_Vars,ONLY:Rmnc_spl,Rmns_spl,Zmnc_spl,Zmns_spl
USE MOD_VMEC_Readin,ONLY:lasym
USE MOD_VMEC     ,ONLY:VMEC_EvalSplMode
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CLASS(t_sol_var_MHD3D),INTENT(INOUT) :: U0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: is,iMode
INTEGER  :: BC_type(2)
REAL(wp) :: BC_val(2)
REAL(wp) :: spos,iota_s 
REAL(wp) :: X1_gIP(1:X1_base%s%nBase)
REAL(wp) :: X2_gIP(1:X2_base%s%nBase)
REAL(wp) :: LA_gIP(1:LA_base%s%nBase,1:LA_base%f%modes)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)') "INTIALIZE SOLUTION..."
  SELECT CASE(which_init)
  CASE(0)
    ASSOCIATE(s_IP        => X1_base%s%s_IP, &
              modes        =>X1_base%f%modes, &
              zero_odd_even=>X1_base%f%zero_odd_even)
    DO imode=1,modes
      SELECT CASE(zero_odd_even(iMode))
      CASE(MN_ZERO,M_ZERO)
        X1_gIP(:)=(1.0_wp-(s_IP(:)**2))*X1_a(iMode)+(s_IP(:)**2)*X1_b(iMode)  ! meet edge and axis, ~(1-s^2)
      CASE(M_ODD)
        X1_gIP(:)=s_IP(:)*X1_b(iMode)      ! first odd mode ~s
      CASE(M_EVEN)
        X1_gIP(:)=(1.0_wp-(s_IP(:)**2))*X1_a(iMode)+(s_IP(:)**2)*X1_b(iMode)   !even mode ~s^2
      END SELECT !X1(:,iMode) zero odd even
      U0%X1(:,iMode)=X1_base%s%initDOF( X1_gIP(:) )
    END DO 
    END ASSOCIATE

    ASSOCIATE(s_IP        => X2_base%s%s_IP, &
              modes        =>X2_base%f%modes, &
              zero_odd_even=>X2_base%f%zero_odd_even)
    DO imode=1,modes
      SELECT CASE(zero_odd_even(iMode))
      CASE(MN_ZERO,M_ZERO)
        X2_gIP(:)=(1.0_wp-(s_IP(:)**2))*X2_a(iMode)+(s_IP(:)**2)*X2_b(iMode) ! meet edge and axis, ~(1-s^2)
      CASE(M_ODD)
        X2_gIP(:)=s_IP(:)*X2_b(iMode)      ! first odd mode ~s
      CASE(M_EVEN)
        X2_gIP(:)=(s_IP(:)**2)*X2_b(iMode) !even mode ~s^2
      END SELECT !X1(:,iMode) zero odd even
      U0%X2(:,iMode)=X2_base%s%initDOF( X2_gIP(:))
    END DO 
    END ASSOCIATE

  CASE(1) !VMEC
    ASSOCIATE(s_IP         => X1_base%s%s_IP, &
              sin_range    => X1_base%f%sin_range,&
              cos_range    => X1_base%f%cos_range )
    DO imode=cos_range(1)+1,cos_range(2)
      X1_gIP(:)       =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,s_IP,Rmnc_Spl)
      U0%X1(:,iMode)  =X1_base%s%initDOF( X1_gIP(:) )
    END DO 
    IF(lasym)THEN
      DO imode=sin_range(1)+1,sin_range(2)
        X1_gIP(:)     =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,s_IP,Rmns_Spl)
        U0%X1(:,iMode)=X1_base%s%initDOF( X1_gIP(:) )
      END DO 
    END IF !lasym
    END ASSOCIATE !X1

    ASSOCIATE(s_IP        => X2_base%s%s_IP, &
              sin_range    =>X2_base%f%sin_range,&
              cos_range    =>X2_base%f%cos_range )
    DO imode=sin_range(1)+1,sin_range(2)
      X2_gIP(:)       =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,s_IP,Zmns_Spl)
      U0%X2(:,iMode)  =X2_base%s%initDOF( X2_gIP(:) )
    END DO 
    IF(lasym)THEN
      DO imode=cos_range(1)+1,cos_range(2)
        X2_gIP(:)     =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,s_IP,Zmnc_Spl)
        U0%X2(:,iMode)=X2_base%s%initDOF( X2_gIP(:) )
      END DO 
    END IF !lasym
    END ASSOCIATE !X2
    !get boundary conditions
    DO imode=1,X1_base%f%modes
      X1_a(iMode)=U0%X1(              1,iMode)
      X1_b(iMode)=U0%X1(X1_base%s%nBase,iMode)
    END DO !iMode
    DO imode=1,X2_base%f%modes
      X2_a(iMode)=U0%X2(              1,iMode)
      X2_b(iMode)=U0%X2(X2_base%s%nBase,iMode)
    END DO !iMode
  END SELECT !which_init
 
  !apply strong boundary conditions
  ASSOCIATE(modes        =>X1_base%f%modes, &
            zero_odd_even=>X1_base%f%zero_odd_even)
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      BC_type=(/BC_TYPE_SYMM    ,BC_TYPE_DIRICHLET/)
      BC_val =(/ X1_a(iMode)    ,      X1_b(iMode)/)
    CASE(M_ODD)
      BC_type=(/BC_TYPE_ANTISYMM,BC_TYPE_DIRICHLET/)
      BC_val =(/          0.0_wp,      X1_b(iMode)/)
    CASE(M_EVEN)
      BC_type=(/BC_TYPE_SYMMZERO,BC_TYPE_DIRICHLET/)
      BC_val =(/          0.0_wp,      X1_b(iMode)/)
    END SELECT !X1(:,iMode) zero odd even
    CALL X1_base%s%applyBCtoDOF(U0%X1(:,iMode),BC_type,BC_val)
  END DO 
  END ASSOCIATE !X1

  ASSOCIATE(modes        =>X2_base%f%modes, &
            zero_odd_even=>X2_base%f%zero_odd_even)
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      BC_type=(/BC_TYPE_SYMM    ,BC_TYPE_DIRICHLET/)
      BC_val =(/     X2_a(iMode),      X2_b(iMode)/)
    CASE(M_ODD)
      BC_type=(/BC_TYPE_ANTISYMM,BC_TYPE_DIRICHLET/)
      BC_val =(/          0.0_wp,      X2_b(iMode)/)
    CASE(M_EVEN)
      BC_type=(/BC_TYPE_SYMMZERO,BC_TYPE_DIRICHLET/)
      BC_val =(/          0.0_wp,      X2_b(iMode)/)
    END SELECT !X1(:,iMode) zero odd even
    CALL X2_base%s%applyBCtoDOF(U0%X2(:,iMode),BC_type,BC_val)
  END DO 
  END ASSOCIATE !X2

  SWRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda ..."
  !initialize Lambda
  LA_gIP(1,:)=0.0_wp !at axis
  DO is=2,LA_base%s%nBase
    spos=LA_base%s%s_IP(is)
    iota_s=EVAL1DPOLY(n_iota_coefs,iota_coefs,spos)
    CALL lambda_Solve(spos,iota_s,U0%X1,U0%X2,LA_gIP(is,:))
  END DO !is
  DO imode=1,LA_base%f%modes
    IF(LA_base%f%zero_odd_even(iMode).EQ.MN_ZERO)THEN
      U0%LA(:,iMode)=0.0_wp !zero mode hsould not be here, but must be zero
    ELSE
      U0%LA(:,iMode)=LA_base%s%initDOF( LA_gIP(:,iMode) )
    END IF!iMode ~ MN_ZERO
  END DO !iMode 
  LA_b = U0%LA(LA_base%s%nbase,:)

  SWRITE(UNIT_stdOut,'(4X,A)') "..DONE."

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
  SDEALLOCATE(X1_b)
  SDEALLOCATE(X2_b)
  SDEALLOCATE(LA_b)
  SDEALLOCATE(X1_a)
  SDEALLOCATE(X2_a)

END SUBROUTINE FinalizeMHD3D

END MODULE MOD_MHD3D
