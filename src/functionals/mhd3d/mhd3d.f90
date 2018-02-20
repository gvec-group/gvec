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
    !-------------------------------------------------------------------------------------------------------------------------------
    LOGICAL :: initialized
    !-------------------------------------------------------------------------------------------------------------------------------
    CONTAINS
      PROCEDURE :: init     => InitMHD3D
      PROCEDURE :: minimize => MinimizeMHD3D
      PROCEDURE :: free     => FinalizeMHD3D
  END TYPE t_functional_mhd3d

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE InitMHD3D(sf) 
  ! MODULES
  USE MOD_MHD3D_Vars
  USE MOD_Globals        , ONLY: TWOPI
  USE MOD_mhdeq_Vars     , ONLY: whichInitEquilibrium
  USE MOD_sgrid          , ONLY: t_sgrid
  USE MOD_fbase          , ONLY: t_fbase,fbase_new
  USE MOD_base           , ONLY: t_base,base_new
  USE MOD_VMEC_Readin    , ONLY: nfp,nFluxVMEC,Phi,xm,xn,lasym
  USE MOD_ReadInTools    , ONLY: GETSTR,GETLOGICAL,GETINT,GETINTARRAY,GETREAL,GETREALALLOCARRAY
  USE MOD_MHD3D_EvalFunc , ONLY: InitializeMHD3D_EvalFunc,EvalEnergy,EvalForce,CheckEvalForce
  USE MOD_Restart_vars   , ONLY: doRestart,RestartFile
  USE MOD_Restart        , ONLY: ReadState
  USE MOD_Analyze        , ONLY: Analyze
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
  INTEGER          :: degGP,mn_nyq(2),mn_nyq_min(2),fac_nyq
  INTEGER          :: nfp_loc,which_hmap 
  INTEGER          :: JacCheck
  REAL(wp)         :: pres_scale
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'
  
  maxIter   = GETINT("maxIter",Proposal=5000)
  outputIter= GETINT("outputIter",Proposal=500)
  logIter   = GETINT("logIter",Proposal=250)
  minimize_tol  =GETREAL("minimize_tol",Proposal=1.0e-12_wp)
  start_dt  =GETREAL("start_dt",Proposal=1.0e-08_wp)

  nElems   = GETINT("sgrid_nElems",Proposal=10)
  grid_type= GETINT("sgrid_grid_type",Proposal=0)
  
  !mandatory global input parameters
  degGP   = GETINT( "degGP",Proposal=4)
  fac_nyq = GETINT( "fac_nyq")
  
  !constants
  
  mu_0    = 2.0e-07_wp*TWOPI
  
  which_init = whichInitEquilibrium ! GETINT("which_init","0")
  init_LA= GETLOGICAL("init_LA",Proposal=.TRUE.)

  PrecondType=GETINT("PrecondType",Proposal=-1)
  
  SELECT CASE(which_init)
  CASE(0)
    init_fromBConly= .TRUE.
    gamm    = GETREAL("GAMMA",Proposal=0.0_wp)
    nfp_loc  = GETINT( "nfp",Proposal=1)
    !hmap
    which_hmap=GETINT("which_hmap",Proposal=1)
    CALL GETREALALLOCARRAY("iota_coefs",iota_coefs,n_iota_coefs,Proposal=(/1.1_wp,0.1_wp/)) !a+b*s+c*s^2...
    CALL GETREALALLOCARRAY("pres_coefs",pres_coefs,n_pres_coefs,Proposal=(/1.0_wp,0.0_wp/)) !a+b*s+c*s^2...
    pres_scale=GETREAL("PRES_SCALE",Proposal=1.0_wp)
    pres_coefs=pres_coefs*pres_scale
    Phi_edge   = GETREAL("PHIEDGE",Proposal=1.0_wp)
    Phi_edge   =-Phi_edge/TWOPI !normalization like in VMEC!!!
  CASE(1) !VMEC init
    init_fromBConly= GETLOGICAL("init_fromBConly",Proposal=.FALSE.)
    gamm = 0.0_wp
    nfp_loc = nfp
    !hmap
    which_hmap=1 !hmap_RZ
    Phi_edge = Phi(nFluxVMEC)
  END SELECT !which_init


  sgammM1=1.0_wp/(gamm-1.0_wp)

  CALL hmap_new(hmap,which_hmap)
  
  
  X1X2_deg     = GETINT(     "X1X2_deg")
  X1X2_cont    = GETINT(     "X1X2_continuity",Proposal=(X1X2_deg-1) )
  X1_mn_max    = GETINTARRAY("X1_mn_max"   ,2 ,Proposal=(/2,0/))
  X2_mn_max    = GETINTARRAY("X2_mn_max"   ,2 ,Proposal=(/2,0/))
  X1_sin_cos   = GETSTR(     "X1_sin_cos"     ,Proposal="_cos_")  !_sin_,_cos_,_sin_cos_
  X2_sin_cos   = GETSTR(     "X2_sin_cos"     ,Proposal="_sin_")
  
  
  LA_deg     = GETINT(     "LA_deg")
  LA_cont    = GETINT(     "LA_continuity",Proposal=-1)
  LA_mn_max  = GETINTARRAY("LA_mn_max", 2 ,Proposal=(/2,0/))
  LA_sin_cos = GETSTR(     "LA_sin_cos"   ,Proposal="_sin_")
  
  IF(fac_nyq.EQ.-1)THEN
    fac_nyq=4
    mn_nyq_min(1)=MAX(1,fac_nyq*MAXVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/)))
    mn_nyq_min(2)=MAX(1,fac_nyq*MAXVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/)))
    mn_nyq  = GETINTARRAY("mn_nyq",2)
    IF(mn_nyq(1).LT.mn_nyq_min(1))THEN
       WRITE(*,*) 'mn_nyq(1) too small, should be >= ',mn_nyq_min(1)
       STOP
    END IF
    IF(mn_nyq(2).LT.mn_nyq_min(2))THEN 
       WRITE(*,*) 'mn_nyq(2) too small, should be >= ',mn_nyq_min(2)
       STOP
    END IF
  ELSE
    mn_nyq(1)=MAX(1,fac_nyq*MAXVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/)))
    mn_nyq(2)=MAX(1,fac_nyq*MAXVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/)))
  END IF
  
  SWRITE(UNIT_stdOut,*)
  SWRITE(UNIT_stdOut,'(A,I4,A,I6," , ",I6,A)')'    fac_nyq = ', fac_nyq,'  ==> interpolation points mn_nyq=( ',mn_nyq(:),' )'
  SWRITE(UNIT_stdOut,*)

  !INITIALIZE GRID  
  CALL sgrid%init(nElems,grid_type)

  !INITIALIZE BASE        !sbase parameter                 !fbase parameter               ...exclude_mn_zero
  CALL base_new(X1_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X1_mn_max,mn_nyq,nfp_loc,X1_sin_cos,.FALSE.)
  CALL base_new(X2_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X2_mn_max,mn_nyq,nfp_loc,X2_sin_cos,.FALSE.)
  CALL base_new(LA_base  ,   LA_deg,  LA_cont,sgrid,degGP , LA_mn_max,mn_nyq,nfp_loc,LA_sin_cos,.TRUE. )

  IF(which_init.EQ.1) THEN !VMEC
    IF(lasym)THEN
      IF((X1_base%f%sin_cos.NE._SINCOS_).OR. &
         (X2_base%f%sin_cos.NE._SINCOS_).OR. &
         (LA_base%f%sin_cos.NE._SINCOS_) ) THEN
        SWRITE(UNIT_stdOut,'(A)')'!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
        SWRITE(UNIT_stdOut,'(A)')'!!!!!!!!   ---->  VMEC was run asymmetric, you should use _sincos_ basis for all variables'
        SWRITE(UNIT_stdOut,'(A)')'!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
      END IF
    END IF
    IF((MAXVAL(INT(xm(:))).GT.MAXVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/))).OR. &
       (MAXVAL(ABS(INT(xn(:))/nfp_loc)).GT.MAXVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/))))THEN
      SWRITE(UNIT_stdOut,'(A)')    '!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
      SWRITE(UNIT_stdOut,'(A,2I6)')'!!!!!!!!   ---->  VMEC was run with a higher mode number (m,n)_max,vmec= ', &
                                    MAXVAL(INT(xm(:))),MAXVAL(ABS(INT(xn(:))/nfp_loc))
      SWRITE(UNIT_stdOut,'(A)')    '!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
    END IF
  END IF

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
  END SELECT !which_init

  !boundary conditions (used in force, in init slightly changed)
  ASSOCIATE(modes        =>X1_base%f%modes, &
            zero_odd_even=>X1_base%f%zero_odd_even)
  ALLOCATE(X1_BC_type(1:2,modes))
  X1_BC_type(BC_EDGE,:)=BC_TYPE_DIRICHLET
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMM     
!      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_NEUMANN
    CASE(M_ODD_FIRST)
      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_ANTISYMM
!      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET
    CASE(M_ODD)
!      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_ANTISYMM
      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET !not too strong for high modes...
    CASE(M_EVEN)
!      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMMZERO
      X1_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET !not too strong for high modes...
    END SELECT !X1(:,iMode) zero odd even
  END DO 
  END ASSOCIATE !X1
  ASSOCIATE(modes        =>X2_base%f%modes, &
            zero_odd_even=>X2_base%f%zero_odd_even)
  ALLOCATE(X2_BC_type(1:2,modes))
  X2_BC_type(BC_EDGE,:)=BC_TYPE_DIRICHLET
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMM     
!      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_NEUMANN
    CASE(M_ODD_FIRST)
      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_ANTISYMM
!      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET
    CASE(M_ODD)
!      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_ANTISYMM
      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET !not too strong for high modes...
    CASE(M_EVEN)
!      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMMZERO
      X2_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET !not too strong for high modes...
    END SELECT !X2(:,iMode) zero odd even
  END DO 
  END ASSOCIATE !X2
  ASSOCIATE(modes        =>LA_base%f%modes, &
            zero_odd_even=>LA_base%f%zero_odd_even)
  ALLOCATE(LA_BC_type(1:2,modes))
  LA_BC_type(BC_EDGE,:)=BC_TYPE_OPEN !no BC for lambda at the edge!
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMMZERO     
!      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_NEUMANN
    CASE(M_ODD_FIRST)
!      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMMZERO     
!      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_ANTISYMM
      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET 
    CASE(M_ODD)
!      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_ANTISYMM
      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET !not too strong for high modes...
    CASE(M_EVEN)
!      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_SYMMZERO
      LA_BC_type(BC_AXIS,iMode)=BC_TYPE_DIRICHLET !not too strong for high modes...
    END SELECT !LA(:,iMode) zero odd even
  END DO 
  END ASSOCIATE !LA
  


  ALLOCATE(U(-1:1))
  CALL U(1)%init((/X1_base%s%nbase,X2_base%s%nbase,LA_base%s%nBase,  &
                   X1_base%f%modes,X2_base%f%modes,LA_base%f%modes/)  )
  ALLOCATE(F(-1:0))
  DO i=-1,0
    CALL U(i)%copy(U(1))
    CALL F(i)%copy(U(1))
  END DO
  ALLOCATE(P(-1:1))
  DO i=-1,1
    CALL P(i)%copy(U(1))
  END DO


  IF(doRestart)THEN
    SWRITE(UNIT_stdOut,'(4X,A)')'... restarting from file ... '
    CALL ReadState(RestartFile,U(0))
    CALL InitSolution(U(0),-1) !only apply BC and recompute lambda
  ELSE 
    CALL InitSolution(U(0),which_init)
  END IF
  CALL U(-1)%set_to(U(0))

 CALL InitializeMHD3D_EvalFunc()
  JacCheck=2
  U(0)%W_MHD3D=EvalEnergy(U(0),.TRUE.,JacCheck)
  IF(JacCheck.EQ.-1)THEN
    CALL Analyze(0)
  END IF
  CALL EvalForce(U(0),.FALSE.,JacCheck, F(0))
  
  CALL CheckEvalForce(U(0),0)
  
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
    get_iMode=GETREAL(TRIM(varstr),Proposal=0.0_wp,quiet_def_in=.TRUE.)
   
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
!> Initialize the solution with the given boundary condition 
!!
!===================================================================================================================================
SUBROUTINE InitSolution(U_init,which_init_in)
! MODULES
  USE MOD_MHD3D_Vars   , ONLY:init_fromBConly
  USE MOD_MHD3D_Vars   , ONLY:X1_base,X1_BC_Type,X1_a,X1_b
  USE MOD_MHD3D_Vars   , ONLY:X2_base,X2_BC_Type,X2_a,X2_b
  USE MOD_MHD3D_Vars   , ONLY:LA_base,init_LA,LA_BC_Type
  USE MOD_sol_var_MHD3D, ONLY:t_sol_var_mhd3d
  USE MOD_lambda_solve,  ONLY:lambda_solve
  USE MOD_VMEC_Vars,     ONLY:Rmnc_spl,Rmns_spl,Zmnc_spl,Zmns_spl
  USE MOD_VMEC_Readin,   ONLY:lasym
  USE MOD_VMEC,          ONLY:VMEC_EvalSplMode
  USE MOD_MHD3D_Profiles,ONLY: Eval_iota
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: which_init_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: U_init
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode,is
  INTEGER  :: BC_type(2)
  REAL(wp) :: BC_val(2)
  REAL(wp) :: spos,iota_s
  REAL(wp) :: X1_gIP(1:X1_base%s%nBase)
  REAL(wp) :: X2_gIP(1:X2_base%s%nBase)
  REAL(wp) :: LA_gIP(1:LA_base%s%nBase,1:LA_base%f%modes)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)') "INTIALIZE SOLUTION..."
  SELECT CASE(which_init_in)
  CASE(-1) !restart
    X1_a(:)=U_init%X1(1,:)
    X2_a(:)=U_init%X2(1,:)
    X1_b(:)=U_init%X1(X1_base%s%nBase,:)
    X2_b(:)=U_init%X2(X2_base%s%nBase,:)
  CASE(0)
    !X1_a,X2_a and X1_b,X2_b already filled from parameter file readin...
  CASE(1) !VMEC
    ASSOCIATE(sin_range    => X1_base%f%sin_range,&
              cos_range    => X1_base%f%cos_range )
    DO imode=cos_range(1)+1,cos_range(2)
      spos=0.0_wp
      X1_a(iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmnc_Spl)
      spos=1.0_wp
      X1_b(iMode:iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmnc_Spl)
    END DO 
    IF(lasym)THEN
      DO imode=sin_range(1)+1,sin_range(2)
        spos=0.0_wp
        X1_a(iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmns_Spl)
        spos=1.0_wp
        X1_b(iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmns_Spl)
      END DO 
    END IF !lasym
    END ASSOCIATE !X1
    ASSOCIATE(sin_range    => X2_base%f%sin_range,&
              cos_range    => X2_base%f%cos_range )
    DO imode=sin_range(1)+1,sin_range(2)
      spos=0.0_wp
      X2_a(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmns_Spl)
      spos=1.0_wp
      X2_b(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmns_Spl)
    END DO 
    IF(lasym)THEN
      DO imode=cos_range(1)+1,cos_range(2)
        spos=0.0_wp
        X2_a(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmnc_Spl)
        spos=1.0_wp
        X2_b(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmnc_Spl)
      END DO 
    END IF !lasym
    END ASSOCIATE !X2
    IF(.NOT.init_fromBConly)THEN !only boundary and axis from VMEC
      ASSOCIATE(s_IP         => X1_base%s%s_IP, &
                nBase        => X1_base%s%nBase, &
                sin_range    => X1_base%f%sin_range,&
                cos_range    => X1_base%f%cos_range )
      DO imode=cos_range(1)+1,cos_range(2)
        DO is=1,nBase
          X1_gIP(is)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,s_IP(is),Rmnc_Spl)
        END DO !is
        U_init%X1(:,iMode)=X1_base%s%initDOF( X1_gIP(:) )
      END DO !imode=cos_range
      IF(lasym)THEN
        DO imode=sin_range(1)+1,sin_range(2)
          DO is=1,nBase
            X1_gIP(is)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,s_IP(is),Rmns_Spl)
          END DO !is
          U_init%X1(:,iMode)=X1_base%s%initDOF( X1_gIP(:) )
        END DO !imode= sin_range
      END IF !lasym
      END ASSOCIATE !X1
      ASSOCIATE(s_IP         => X2_base%s%s_IP, &
                nBase        => X2_base%s%nBase, &
                sin_range    => X2_base%f%sin_range,&
                cos_range    => X2_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        DO is=1,nBase
          X2_gIP(is)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,s_IP(is),Zmns_Spl)
        END DO !is
        U_init%X2(:,iMode)=X2_base%s%initDOF( X2_gIP(:) )
      END DO !imode=sin_range
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          DO is=1,nBase
            X2_gIP(is)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,s_IP(is),Zmnc_Spl)
          END DO !is
          U_init%X2(:,iMode)=X2_base%s%initDOF( X2_gIP(:) )
        END DO !imode= sin_range
      END IF !lasym
      END ASSOCIATE !X2
    END IF !fullIntVmec
  END SELECT !which_init
 
  IF((which_init_in.GE.0).AND.(init_fromBConly))THEN 
    !no restart(=-1) and initialization only 
    !smoothly interpolate between  edge and axis data
    ASSOCIATE(s_IP         =>X1_base%s%s_IP, &
              modes        =>X1_base%f%modes, &
              zero_odd_even=>X1_base%f%zero_odd_even)
    DO imode=1,modes
      SELECT CASE(zero_odd_even(iMode))
      CASE(MN_ZERO,M_ZERO) !X1_a only used here!!
        X1_gIP(:)=(1.0_wp-(s_IP(:)**2))*X1_a(iMode)+(s_IP(:)**2)*X1_b(iMode)  ! meet edge and axis, ~(1-s^2)
      CASE(M_ODD_FIRST)
        X1_gIP(:)=s_IP(:)*X1_b(iMode)      ! first odd mode ~s
      CASE(M_ODD)
        X1_gIP(:)=(s_IP(:)**3)*X1_b(iMode) ! higher odd modes ~s^3
      CASE(M_EVEN)
        X1_gIP(:)=(s_IP(:)**2)*X1_b(iMode)   !even mode ~s^2
      END SELECT !X1(:,iMode) zero odd even
      U_init%X1(:,iMode)=X1_base%s%initDOF( X1_gIP(:) )
    END DO 
    END ASSOCIATE
    
    ASSOCIATE(s_IP         =>X2_base%s%s_IP, &
              modes        =>X2_base%f%modes, &
              zero_odd_even=>X2_base%f%zero_odd_even)
    DO imode=1,modes
      SELECT CASE(zero_odd_even(iMode))
      CASE(MN_ZERO,M_ZERO) !X2_a only used here!!!
        X2_gIP(:)=(1.0_wp-(s_IP(:)**2))*X2_a(iMode)+(s_IP(:)**2)*X2_b(iMode) ! meet edge and axis, ~(1-s^2)
      CASE(M_ODD_FIRST)
        X2_gIP(:)=s_IP(:)*X2_b(iMode)      ! first odd mode ~s
      CASE(M_ODD)
        X2_gIP(:)=(s_IP(:)**3)*X2_b(iMode) ! higher odd modes ~s^3
      CASE(M_EVEN)
        X2_gIP(:)=(s_IP(:)**2)*X2_b(iMode) !even mode ~s^2
      END SELECT !X2(:,iMode) zero odd even
      U_init%X2(:,iMode)=X2_base%s%initDOF( X2_gIP(:))
    END DO 
    END ASSOCIATE
  END IF !init_fromBConly

  !apply strong boundary conditions
  ASSOCIATE(modes        =>X1_base%f%modes, &
            zero_odd_even=>X1_base%f%zero_odd_even)
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      BC_val =(/ X1_a(iMode)    ,      X1_b(iMode)/)
    CASE(M_ODD_FIRST,M_ODD,M_EVEN)
      BC_val =(/          0.0_wp,      X1_b(iMode)/)
    END SELECT !X1(:,iMode) zero odd even
    CALL X1_base%s%applyBCtoDOF(U_init%X1(:,iMode),X1_BC_type(:,iMode),BC_val)
  END DO 
  END ASSOCIATE !X1
  
  ASSOCIATE(modes        =>X2_base%f%modes, &
            zero_odd_even=>X2_base%f%zero_odd_even)
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      BC_val =(/     X2_a(iMode),      X2_b(iMode)/)
    CASE(M_ODD_FIRST,M_ODD,M_EVEN)
      BC_val =(/          0.0_wp,      X2_b(iMode)/)
    END SELECT !X1(:,iMode) zero odd even
    CALL X2_base%s%applyBCtoDOF(U_init%X2(:,iMode),X2_BC_type(:,iMode),BC_val)
  END DO 
  END ASSOCIATE !X2

  IF(init_LA)THEN
    SWRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda from mapping ..."
    !initialize Lambda
    LA_gIP(1,:)=0.0_wp !at axis
    DO is=2,LA_base%s%nBase
      spos=LA_base%s%s_IP(is)
      iota_s=eval_iota(spos)
      CALL lambda_Solve(spos,iota_s,U_init%X1,U_init%X2,LA_gIP(is,:))
    END DO !is
    ASSOCIATE(modes        =>LA_base%f%modes, &
              zero_odd_even=>LA_base%f%zero_odd_even)
    DO imode=1,modes
      IF(zero_odd_even(iMode).EQ.MN_ZERO)THEN
        U_init%LA(:,iMode)=0.0_wp ! (0,0) mode should not be here, but must be zero if its used.
      ELSE
        U_init%LA(:,iMode)=LA_base%s%initDOF( LA_gIP(:,iMode) )
      END IF!iMode ~ MN_ZERO
      BC_val =(/ 0.0_wp, 0.0_wp/)
      CALL LA_base%s%applyBCtoDOF(U_init%LA(:,iMode),LA_BC_type(:,iMode),BC_val)
    END DO !iMode 
    END ASSOCIATE !LA
  ELSE
    !lambda init might not be needed since it has no boundary condition and changes anyway after the update of the mapping...
    SWRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda =0 ..."
    U_init%LA=0.0_wp
  END IF !init_LA

  SWRITE(UNIT_stdOut,'(4X,A)') "... DONE."
  SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE InitSolution


!===================================================================================================================================
!> Compute Equilibrium, iteratively
!!
!===================================================================================================================================
SUBROUTINE MinimizeMHD3D(sf) 
! MODULES
  USE MOD_MHD3D_Vars
  USE MOD_MHD3D_EvalFunc
  USE MOD_Analyze, ONLY:analyze
  USE MOD_Restart, ONLY:WriteState
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iter,nStepDecreased,nSkip_Jac,nSkip_dw
  INTEGER   :: JacCheck,lastoutputIter
  REAL(wp)  :: beta,dt,deltaW,absTol
  REAL(wp)  :: min_dt_out,max_dt_out,min_dw_out,max_dw_out,t_pseudo,Fnorm,Fnorm0,W_MHD3D_0
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)') "MINIMIZE MHD3D FUNCTIONAL..."
  min_dt_out=1.0e+30_wp
  max_dt_out=0.0_wp
  min_dW_out=1.0e+30_wp
  max_dW_out=-1.0e+30_wp
  nSkip_Jac=0
  nSkip_dW =0
  CALL WriteState(U(0),0)

  JacCheck=1 !abort if detJ<0
  CALL EvalAux(           U(0),JacCheck)
  U(0)%W_MHD3D=EvalEnergy(U(0),.FALSE.,JacCheck)
  W_MHD3D_0 = U(0)%W_MHD3D

  abstol=minimize_tol

  CALL EvalForce(         U(0),.FALSE.,JacCheck,F(0))
  Fnorm0=SQRT(SUM(F(0)%norm_2()))
  Fnorm=Fnorm0

  CALL P( -1)%set_to(0.0_wp)

  CALL U( -1)%set_to(U(0))

  beta=0.3_wp  !for damping
  dt=start_dt
  nstepDecreased=0
  t_pseudo=0
  lastOutputIter=0
  iter=1
  SWRITE(UNIT_stdOut,'(A,E11.4,A)')'%%%%%%%%%%  START ITERATION, dt= ',dt, '  %%%%%%%%%%%%%%%%%%%%%%%%%%%'
  DO WHILE(iter.LE.maxIter)
! hirshman method
!    CALL P(0)%AXBY(beta,P(-1),1.0_wp,F(0))
!    CALL P(1)%AXBY(1.0_wp,U(0),dt,P(0)) !overwrites P(1)

!simple gradient
!    CALL P(1)%AXBY(1.0_wp,U(0),dt,F(0)) !overwrites P(1)

!damping (beta >0), U^(n+1)=U^(n)+(1-beta*sqrt(dt))*(U^(n)-U^(n-1))+dt F , beta->(1-beta*sqrt(dt))
               !        =U^(n)*(1+beta) - beta*U^(n-1) =P0 + dt F

!    !for damping (1-beta*sqrt(dt)), beta=20.0_wp  
!    CALL P(0)%AXBY((2.0_wp-beta*sqrt(dt)),U(0),(beta*sqrt(dt)-1.0_wp),U(-1)) !overwrites P(0)

    !for damping   beta=0.6
    CALL P(0)%AXBY((1.0_wp+beta),U(0),-beta,U(-1)) !overwrites P(0)
    CALL P(1)%AXBY(1.0_wp,P(0),dt,F(0)) !overwrites P(1)



    JacCheck=2 !no abort,if detJ<0, JacCheck=-1
    P(1)%W_MHD3D=EvalEnergy(P(1),.TRUE.,JacCheck) 
    IF(JacCheck.EQ.-1)THEN
      dt=0.5_wp*dt
      nstepDecreased=nStepDecreased+1
      nSkip_Jac=nSkip_Jac+1
      SWRITE(UNIT_stdOut,'(8X,I8,A)')iter,'...detJac<0, skip step and decrease stepsize!'
      !do not use P(1), redo the iteration
    ELSE 
      !detJ>0
      deltaW=P(1)%W_MHD3D-U(0)%W_MHD3D!should be <=0, 
      
      IF(deltaW.LE.1.0e-10*W_MHD3D_0)THEN !valid step 
         !LINE SEARCH COULD BE USED HERE !!

!        IF(ABS(deltaW).LE.aborttol)THEN
!          SWRITE(UNIT_stdOut,'(A,A,E11.4)')'Iteration finished, energy stagnates in relative tolerance, ', &
!                                           ' deltaW= ' ,U(0)%W_MHD3D-U(-1)%W_MHD3D
!        IF(Fnorm*dt.LE.reltol*Fnorm0)THEN
        IF(ALL(SQRT(F(0)%norm_2())*dt.LE.abstol))THEN
          SWRITE(UNIT_stdOut,'(74("%")"\n",A,I8,A,2I8,A,E11.4,A,2E11.4,A,E21.14,A,E11.4,A,3E11.4,A)') &
                            '%%%  #ITERATIONS= ',iter,', #skippedIter (Jac/dW)= ',nSkip_Jac,nSkip_dW, &
                    '    %%%\n%%%  t_pseudo= ',t_pseudo,', min/max dt= ',min_dt_out,max_dt_out, &
                   '         %%%\n%%%  W_MHD3D= ',U(0)%W_MHD3D,', deltaW= ' , deltaW , &
          '               %%%\n%%% dU = |Force|*dt= ',SQRT(F(0)%norm_2())*dt, &
   '                        %%%\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
          SWRITE(UNIT_stdOut,'(4x,A)')'==>Iteration finished, |force| in relative tolerance'
          EXIT !DO LOOP
        END IF
        iter=iter+1
        nstepDecreased=0
        min_dt_out=MIN(min_dt_out,dt)
        max_dt_out=MAX(max_dt_out,dt)
        min_dW_out=MIN(min_dW_out,deltaW)
        max_dW_out=MAX(max_dW_out,deltaW)
        IF(MOD(iter,logIter).EQ.0)THEN 
          SWRITE(UNIT_stdOut,'(74("%")"\n",A,I8,A,2I8,A,E11.4,A,2E11.4,A,E21.14,A,E11.4,A,3E11.4,A)') &
                            '%%%  #ITERATIONS= ',iter,', #skippedIter (Jac/dW)= ',nSkip_Jac,nSkip_dW, &
                    '    %%%\n%%%  t_pseudo= ',t_pseudo,', min/max dt= ',min_dt_out,max_dt_out, &
                   '         %%%\n%%%  W_MHD3D= ',U(0)%W_MHD3D,', deltaW= ' , deltaW , &
          '               %%%\n%%% dU = |Force|*dt= ',SQRT(F(0)%norm_2())*dt, &
   '                        %%%\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
          min_dt_out=1.0e+30_wp
          max_dt_out=0.0_wp
          min_dW_out=1.0e+30_wp
          max_dW_out=-1.0e+30_wp
          nSkip_Jac=0
          nSkip_dW =0
        END IF

          
! for hirshman method
!        CALL F(-1)%set_to(F(0))
!        CALL P(-1)%set_to(P(0))
        t_pseudo=t_pseudo+dt
! for simple gradient & damping
        CALL U(-1)%set_to(U(0))
        CALL U(0)%set_to(P(1))
        CALL EvalForce(P(1),.FALSE.,JacCheck,F(0)) !evalAux was already called on P(1)=U(0), so that its set false here.
        Fnorm=SUM(SQRT(F(0)%norm_2()))
! for hirshman method
!        beta=SUM(F(0)%norm_2())/SUM(F(-1)%norm_2())

       !increase time step
        dt=1.001_wp*dt
      ELSE !not a valid step, decrease timestep and skip P(1)
        dt=0.5_wp*dt
        nstepDecreased=nStepDecreased+1
        nSkip_dW=nSkip_dW+1
        SWRITE(UNIT_stdOut,'(8X,I8,A)')iter,'...deltaW>0, skip step and decrease stepsize!'
      END IF
    END IF !JacCheck
   
    IF(nStepDecreased.GT.20) THEN ! 2^20 ~10^6
      SWRITE(UNIT_stdOut,'(A,E21.11)')'Iteration stopped since timestep has been decreased by 2^20: ', dt 
      SWRITE(UNIT_stdOut,fmt_sep)
      RETURN
    END IF
    IF((MOD(iter,outputIter).EQ.0).AND.(lastoutputIter.NE.iter))THEN
      SWRITE(UNIT_stdOut,'(A)')'#######################  VISUALIZATION ##############################'
      CALL Analyze(iter)
      CALL WriteState(U(0),iter)
      CALL CheckEvalForce(U(0),iter)
      lastOutputIter=iter
    END IF
  END DO !iter
  IF(iter.GE.MaxIter)THEN
    SWRITE(UNIT_stdOut,'(A,E21.11)')"maximum iteration count exceeded, not converged" 
  END IF
  SWRITE(UNIT_stdOut,'(A)') "... DONE."
  SWRITE(UNIT_stdOut,fmt_sep)
  CALL Analyze(99999999)
  CALL WriteState(U(0),99999999)
!DEBUG
!  WRITE(FileString,'(A,"_State_",I4.4,"_",I8.8,".dat")')TRIM(ProjectName),OutputLevel,99999999
!  CALL ReadState(FileString,U(-1))
  

END SUBROUTINE MinimizeMHD3D


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeMHD3D(sf) 
! MODULES
  USE MOD_MHD3D_Vars
  USE MOD_MHD3D_EvalFunc,ONLY:FinalizeMHD3D_EvalFunc
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
    CALL P(i)%free()
  END DO
  DO i=-1,0
    CALL F(i)%free()
  END DO
  CALL sgrid%free()

  SDEALLOCATE(U)
  SDEALLOCATE(F)
  SDEALLOCATE(P)
  SDEALLOCATE(X1_BC_type)
  SDEALLOCATE(X2_BC_type)
  SDEALLOCATE(LA_BC_type)
  SDEALLOCATE(X1_b)
  SDEALLOCATE(X2_b)
  SDEALLOCATE(LA_b)
  SDEALLOCATE(X1_a)
  SDEALLOCATE(X2_a)
  SDEALLOCATE(pres_coefs)
  SDEALLOCATE(iota_coefs)

  CALL FinalizeMHD3D_EvalFunc()
END SUBROUTINE FinalizeMHD3D

END MODULE MOD_MHD3D
