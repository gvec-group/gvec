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
MODULE MODgvec_MHD3D
! MODULES
  USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep
  USE MODgvec_c_functional,   ONLY: t_functional
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
      PROCEDURE :: initSolution => InitSolutionMHD3D
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
  USE MODgvec_MHD3D_Vars
  USE MODgvec_Globals        , ONLY: TWOPI
  USE MODgvec_sgrid          , ONLY: t_sgrid
  USE MODgvec_fbase          , ONLY: t_fbase,fbase_new
  USE MODgvec_base           , ONLY: t_base,base_new
  USE MODgvec_hmap           , ONLY: hmap_new
  USE MODgvec_VMEC           , ONLY: InitVMEC
  USE MODgvec_VMEC_vars      , ONLY: switchZeta
  USE MODgvec_VMEC_Readin    , ONLY: nfp,nFluxVMEC,Phi,xm,xn,lasym,mpol,ntor
  USE MODgvec_MHD3D_EvalFunc , ONLY: InitializeMHD3D_EvalFunc
  USE MODgvec_ReadInTools    , ONLY: GETSTR,GETLOGICAL,GETINT,GETINTARRAY,GETREAL,GETREALALLOCARRAY
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
  INTEGER          :: nfp_loc
  INTEGER          :: sign_iota
  INTEGER          :: X1X2_BCtype_axis(0:4),LA_BCtype_axis(0:4)
  INTEGER          :: proposal_mn_max(1:2)=(/2,0/) !!default proposals, changed for VMEC input to automatically match input!
  CHARACTER(LEN=8) :: proposal_X1_sin_cos="_cos_"  !!default proposals, changed for VMEC input to automatically match input!
  CHARACTER(LEN=8) :: proposal_X2_sin_cos="_sin_"  !!default proposals, changed for VMEC input to automatically match input!
  CHARACTER(LEN=8) :: proposal_LA_sin_cos="_sin_"  !!default proposals, changed for VMEC input to automatically match input!
  REAL(wp)         :: pres_scale
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT MHD3D ...'

  which_init = GETINT("whichInitEquilibrium")
  IF(which_init.EQ.1) CALL InitVMEC()

  !-----------MINIMIZER
  MinimizerType= GETINT("MinimizerType",Proposal=0)
  PrecondType  = GETINT("PrecondType",Proposal=-1)

  dW_allowed=GETREAL("dW_allowed",Proposal=1.0e-10) !! for minimizer, accept step if dW<dW_allowed*W_MHD(iter=0) default +10e-10
  IF(MinimizerType.EQ.10)THEN !for hirshman method
    doLineSearch=.FALSE.  !does not work
  ELSE
    doLineSearch=GETLOGICAL("doLineSearch",Proposal=.FALSE.) ! does not improve convergence
  END IF
  maxIter   = GETINT("maxIter",Proposal=5000)
  outputIter= GETINT("outputIter",Proposal=500)
  logIter   = GETINT("logIter",Proposal=250)
  nlogScreen= GETINT("nLogScreen",Proposal=1)
  minimize_tol  =GETREAL("minimize_tol",Proposal=1.0e-12_wp)
  start_dt  =GETREAL("start_dt",Proposal=1.0e-08_wp)
  doCheckDistance=GETLOGICAL("doCheckDistance",Proposal=.FALSE.)
  doCheckAxis=GETLOGICAL("doCheckAxis",Proposal=.TRUE.)
  !-----------


  nElems   = GETINT("sgrid_nElems",Proposal=10)
  grid_type= GETINT("sgrid_grid_type",Proposal=0)

  !mandatory global input parameters
  degGP   = GETINT( "degGP",Proposal=4)
  fac_nyq = GETINT( "fac_nyq")

  !constants
  mu_0    = 2.0e-07_wp*TWOPI


  init_LA= GETLOGICAL("init_LA",Proposal=.TRUE.)

  SELECT CASE(which_init)
  CASE(0)
    init_fromBConly= .TRUE.
    init_BC        = 2
    gamm    = GETREAL("GAMMA",Proposal=0.0_wp)
    nfp_loc  = GETINT( "nfp",Proposal=1)
    !hmap
    which_hmap=GETINT("which_hmap",Proposal=1)
    sign_iota  = GETINT( "sign_iota",Proposal=-1) !if positive in vmec, this should be -1, because of (R,Z,phi) coordinate system
    CALL GETREALALLOCARRAY("iota_coefs",iota_coefs,n_iota_coefs,Proposal=(/1.1_wp,0.1_wp/)) !a+b*s+c*s^2...
    iota_coefs=REAL(sign_iota)*iota_coefs
    CALL GETREALALLOCARRAY("pres_coefs",pres_coefs,n_pres_coefs,Proposal=(/1.0_wp,0.0_wp/)) !a+b*s+c*s^2...
    pres_scale=GETREAL("PRES_SCALE",Proposal=1.0_wp)
    pres_coefs=pres_coefs*pres_scale
    Phi_edge   = GETREAL("PHIEDGE",Proposal=1.0_wp)
    Phi_edge   = Phi_edge/TWOPI !normalization like in VMEC!!!
  CASE(1) !VMEC init
    init_fromBConly= GETLOGICAL("init_fromBConly",Proposal=.FALSE.)
    IF(init_fromBConly)THEN
      !=-1, keep vmec axis and boundary, =0: keep vmec boundary, overwrite axis, =1: keep vmec axis, overwrite boundary, =2: overwrite axis and boundary
      init_BC= GETINT("reinit_BC",Proposal=-1)
    ELSE
      init_BC=-1
    END IF


    proposal_mn_max(:)=(/mpol-1,ntor/)
    IF(lasym)THEN !asymmetric
      proposal_X1_sin_cos="_sincos_"
      proposal_X2_sin_cos="_sincos_"
      proposal_LA_sin_cos="_sincos_"
    END IF
    gamm = 0.0_wp
    nfp_loc = nfp
    !hmap: depends on how vmec data is read:
    IF(switchZeta)THEN
      which_hmap=1 !hmap_RZ
    ELSE
      which_hmap=2 !hmap_RphiZ
    END IF
    Phi_edge = Phi(nFluxVMEC)
  END SELECT !which_init

  init_average_axis= GETLOGICAL("init_average_axis",Proposal=.FALSE.)
  IF(init_average_axis)THEN
    average_axis_move(1) = GETREAL("average_axis_move_X1",Proposal=0.0_wp)
    average_axis_move(2) = GETREAL("average_axis_move_X2",Proposal=0.0_wp)
  END IF

  sgammM1=1.0_wp/(gamm-1.0_wp)

  CALL hmap_new(hmap,which_hmap)


  X1X2_deg     = GETINT(     "X1X2_deg")
  X1X2_cont    = GETINT(     "X1X2_continuity",Proposal=(X1X2_deg-1) )
  X1_mn_max    = GETINTARRAY("X1_mn_max"   ,2 ,Proposal=proposal_mn_max)
  X2_mn_max    = GETINTARRAY("X2_mn_max"   ,2 ,Proposal=proposal_mn_max)
  X1_sin_cos   = GETSTR(     "X1_sin_cos"     ,Proposal=proposal_X1_sin_cos)  !_sin_,_cos_,_sin_cos_
  X2_sin_cos   = GETSTR(     "X2_sin_cos"     ,Proposal=proposal_X2_sin_cos)


  LA_deg     = GETINT(     "LA_deg")
  LA_cont    = GETINT(     "LA_continuity",Proposal=-1)
  LA_mn_max  = GETINTARRAY("LA_mn_max", 2 ,Proposal=proposal_mn_max)
  LA_sin_cos = GETSTR(     "LA_sin_cos"   ,Proposal=proposal_LA_sin_cos)

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
        CALL abort(__STAMP__,&
            '!!!!  VMEC was run asymmetric, you should use _sincos_ basis for all variables')
      END IF
    END IF
    IF((MAXVAL(INT(xm(:))).GT.MINVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/))).OR. &
       (MAXVAL(ABS(INT(xn(:))/nfp_loc)).GT.MINVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/))))THEN
      SWRITE(UNIT_stdOut,'(A)')    '!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
      SWRITE(UNIT_stdOut,'(A,2I6)')'!!!!!!!!   ---->  you use a lower mode number than the VMEC  run  ', &
                                    MAXVAL(INT(xm(:))),MAXVAL(ABS(INT(xn(:))/nfp_loc))
      SWRITE(UNIT_stdOut,'(A)')    '!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
      !  CALL abort(__STAMP__,&
      !'!!!!!  you use a lower mode number than the VMEC  run  (m,n)_max')
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

  IF((init_BC.EQ.0).OR.(init_BC.EQ.2))THEN !READ axis values from input file
    WRITE(UNIT_stdOut,'(4X,A)')'... read axis data for X1:'
    ASSOCIATE(modes=>X1_base%f%modes,sin_range=>X1_base%f%sin_range,cos_range=>X1_base%f%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X1_a(iMode)=get_iMode('X1_a_sin',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1_a(iMode)=get_iMode('X1_a_cos',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
    WRITE(UNIT_stdOut,'(4X,A)')'... read axis data for X2:'
    ASSOCIATE(modes=>X2_base%f%modes,sin_range=>X2_base%f%sin_range,cos_range=>X2_base%f%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X2_a(iMode)=get_iMode('X2_a_sin',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2_a(iMode)=get_iMode('X2_a_cos',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
  END IF
  IF((init_BC.EQ.1).OR.(init_BC.EQ.2))THEN !READ edge values from input file
    WRITE(UNIT_stdOut,'(4X,A)')'... read edge boundary data for X1:'
    ASSOCIATE(modes=>X1_base%f%modes,sin_range=>X1_base%f%sin_range,cos_range=>X1_base%f%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X1_b(iMode)=get_iMode('X1_b_sin',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1_b(iMode)=get_iMode('X1_b_cos',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
    WRITE(UNIT_stdOut,'(4X,A)')'... read edge boundary data for X2:'
    ASSOCIATE(modes=>X2_base%f%modes,sin_range=>X2_base%f%sin_range,cos_range=>X2_base%f%cos_range)
    DO iMode=sin_range(1)+1,sin_range(2)
      X2_b(iMode)=get_iMode('X2_b_sin',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2_b(iMode)=get_iMode('X2_b_cos',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
  END IF !init_BC

  X1X2_BCtype_axis(MN_ZERO    )= GETINT("X1X2_BCtype_axis_mn_zero"    ,Proposal=BC_TYPE_NEUMANN  ) ! stronger: BC_TYPE_SYMM
  X1X2_BCtype_axis(M_ZERO     )= GETINT("X1X2_BCtype_axis_m_zero"     ,Proposal=BC_TYPE_NEUMANN  ) ! stronger: BC_TYPE_SYMM
  X1X2_BCtype_axis(M_ODD_FIRST)= GETINT("X1X2_BCtype_axis_m_odd_first",Proposal=BC_TYPE_DIRICHLET) ! stronger: BC_TYPE_ANTISYMM
  X1X2_BCtype_axis(M_ODD      )= GETINT("X1X2_BCtype_axis_m_odd"      ,Proposal=BC_TYPE_DIRICHLET) ! stronger: BC_TYPE_ANTISYMM
  X1X2_BCtype_axis(M_EVEN     )= GETINT("X1X2_BCtype_axis_m_even"     ,Proposal=BC_TYPE_DIRICHLET) ! stronger: BC_TYPE_SYMMZERO


  !boundary conditions (used in force, in init slightly changed)
  ASSOCIATE(modes        =>X1_base%f%modes, zero_odd_even=>X1_base%f%zero_odd_even)
  ALLOCATE(X1_BC_type(1:2,modes))
  X1_BC_type(BC_EDGE,:)=BC_TYPE_DIRICHLET
  DO imode=1,modes
    X1_BC_type(BC_AXIS,iMode)=X1X2_BCtype_axis(zero_odd_even(iMode))
  END DO
  END ASSOCIATE !X1

  ASSOCIATE(modes        =>X2_base%f%modes, zero_odd_even=>X2_base%f%zero_odd_even)
  ALLOCATE(X2_BC_type(1:2,modes))
  X2_BC_type(BC_EDGE,:)=BC_TYPE_DIRICHLET
  DO imode=1,modes
    X2_BC_type(BC_AXIS,iMode)=X1X2_BCtype_axis(zero_odd_even(iMode))
  END DO
  END ASSOCIATE !X2

  LA_BCtype_axis(MN_ZERO    )= GETINT("LA_BCtype_axis_mn_zero"    ,Proposal=BC_TYPE_DIRICHLET) ! stronger: BC_TYPE_SYMMZERO
  LA_BCtype_axis(M_ZERO     )= GETINT("LA_BCtype_axis_m_zero"     ,Proposal=BC_TYPE_NEUMANN  ) ! stronger: BC_TYPE_SYMM
  LA_BCtype_axis(M_ODD_FIRST)= GETINT("LA_BCtype_axis_m_odd_first",Proposal=BC_TYPE_DIRICHLET) ! stronger: BC_TYPE_ANTISYMM
  LA_BCtype_axis(M_ODD      )= GETINT("LA_BCtype_axis_m_odd"      ,Proposal=BC_TYPE_DIRICHLET) ! stronger: BC_TYPE_ANTISYMM
  LA_BCtype_axis(M_EVEN     )= GETINT("LA_BCtype_axis_m_even"     ,Proposal=BC_TYPE_DIRICHLET) ! stronger:BC_TYPE_SYMMZERO

  ASSOCIATE(modes        =>LA_base%f%modes, zero_odd_even=>LA_base%f%zero_odd_even)
  ALLOCATE(LA_BC_type(1:2,modes))
  LA_BC_type(BC_EDGE,:)=BC_TYPE_OPEN
  DO imode=1,modes
    LA_BC_type(BC_AXIS,iMode)=LA_BCtype_axis(zero_odd_even(iMode))
  END DO
  END ASSOCIATE !LA

  ! ALLOCATE DATA
  ALLOCATE(U(-3:1))
  CALL U(1)%init((/X1_base%s%nbase,X2_base%s%nbase,LA_base%s%nBase,  &
                   X1_base%f%modes,X2_base%f%modes,LA_base%f%modes/)  )
  DO i=-3,0
    CALL U(i)%copy(U(1))
  END DO
  ALLOCATE(F(-1:0))
  DO i=-1,0
    CALL F(i)%copy(U(1))
  END DO
  ALLOCATE(V(-1:1))
  DO i=-1,1
    CALL V(i)%copy(U(1))
  END DO
  ALLOCATE(P(-1:1))
  DO i=-1,1
    CALL P(i)%copy(U(1))
  END DO

  CALL InitializeMHD3D_EvalFunc()

  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)


END SUBROUTINE InitMHD3D



!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE InitSolutionMHD3D(sf)
! MODULES
  USE MODgvec_MHD3D_Vars     , ONLY: which_init,U,F,X1_base,X2_base,LA_base
  USE MODgvec_Restart_vars   , ONLY: doRestart,RestartFile
  USE MODgvec_Restart        , ONLY: RestartFromState
  USE MODgvec_Restart        , ONLY: WriteState
  USE MODgvec_MHD3D_EvalFunc , ONLY: EvalEnergy,EvalForce,CheckEvalForce
  USE MODgvec_Analyze        , ONLY: Analyze
  USE MODgvec_ReadInTools    , ONLY: GETLOGICAL
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: JacCheck,iMode
  LOGICAL              :: boundary_perturb !! false: no boundary perturbation, true: add boundary perturbation X1pert_b,X2pert_b
  REAL(wp),ALLOCATABLE :: X1pert_b(:)      !! fourier modes of the boundary perturbation for X1
  REAL(wp),ALLOCATABLE :: X2pert_b(:)      !! fourier modes of the boundary perturbation for X2
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)') "INTIALIZE SOLUTION..."
  IF(doRestart)THEN
    SWRITE(UNIT_stdOut,'(4X,A)')'... restarting from file ... '
    CALL RestartFromState(RestartFile,U(0))
    !CALL InitSolution(U(0),-1) !would apply BC and recompute lambda
  ELSE
    CALL InitSolution(U(0),which_init)
  END IF

  boundary_perturb=GETLOGICAL('boundary_perturb',Proposal=.FALSE.)
  IF(boundary_perturb)THEN
    ALLOCATE(X1pert_b(1:X1_base%f%modes) )
    ALLOCATE(X2pert_b(1:X2_base%f%modes) )
    X1pert_b=0.0_wp
    X2pert_b=0.0_wp
    !READ boudnary values from input file
    ASSOCIATE(modes=>X1_base%f%modes,sin_range=>X1_base%f%sin_range,cos_range=>X1_base%f%cos_range)
    WRITE(UNIT_stdOut,'(4X,A)')'... read data for X1pert:'
    DO iMode=sin_range(1)+1,sin_range(2)
      X1pert_b(iMode)=get_iMode('X1pert_b_sin',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1pert_b(iMode)=get_iMode('X1pert_b_cos',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
    ASSOCIATE(modes=>X2_base%f%modes,sin_range=>X2_base%f%sin_range,cos_range=>X2_base%f%cos_range)
    WRITE(UNIT_stdOut,'(4X,A)')'... read data for X2pert:'
    DO iMode=sin_range(1)+1,sin_range(2)
      X2pert_b(iMode)=get_iMode('X2pert_b_sin',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2pert_b(iMode)=get_iMode('X2pert_b_cos',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
    CALL AddBoundaryPerturbation(U(0),0.3,X1pert_b,X2pert_b)
    DEALLOCATE(X1pert_b,X2pert_b)
  END IF

  CALL U(-1)%set_to(U(0))

  JacCheck=2
  U(0)%W_MHD3D=EvalEnergy(U(0),.TRUE.,JacCheck)
  CALL WriteState(U(0),0)
  CALL EvalForce(U(0),.FALSE.,JacCheck, F(0))
  SWRITE(UNIT_stdOut,'(8x,A,3E11.4)')'|Force|= ',SQRT(F(0)%norm_2())
  CALL CheckEvalForce(U(0),0)
  CALL Analyze(0)

  SWRITE(UNIT_stdOut,'(4X,A)') "... DONE."
  SWRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE InitSolutionMHD3D


!===================================================================================================================================
!> automatically build the string to be read from parameterfile, varname + m,n mode number, and then read it from parameterfile
!!
!===================================================================================================================================
FUNCTION get_iMode(varname_in,mn_in,nfp_in)
! MODULES
  USE MODgvec_ReadInTools    , ONLY: GETREAL
!$ USE omp_lib
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN) :: varname_in
    INTEGER         ,INTENT(IN) :: mn_in(2),nfp_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL(wp)                    :: get_iMode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
    CHARACTER(LEN=100) :: varstr
!===================================================================================================================================
  WRITE(varstr,'(A,"("I4,";",I4,")")')TRIM(varname_in),mn_in(1),mn_in(2)/nfp_in
  varstr=delete_spaces(varstr)         !quiet on default=0.0
  get_iMode=GETREAL(TRIM(varstr),Proposal=0.0_wp,quiet_def_in=.TRUE.)

  CONTAINS

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
END FUNCTION get_iMode


!===================================================================================================================================
!> Initialize the solution with the given boundary condition
!!
!===================================================================================================================================
SUBROUTINE InitSolution(U_init,which_init_in)
! MODULES
  USE MODgvec_Globals,       ONLY:ProgressBar
  USE MODgvec_MHD3D_Vars   , ONLY:init_fromBConly,init_BC,init_average_axis,average_axis_move
  USE MODgvec_MHD3D_Vars   , ONLY:X1_base,X1_BC_Type,X1_a,X1_b
  USE MODgvec_MHD3D_Vars   , ONLY:X2_base,X2_BC_Type,X2_a,X2_b
  USE MODgvec_MHD3D_Vars   , ONLY:LA_base,init_LA,LA_BC_Type
  USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_mhd3d
  USE MODgvec_lambda_solve,  ONLY:lambda_solve
  USE MODgvec_VMEC_Vars,     ONLY:Rmnc_spl,Rmns_spl,Zmnc_spl,Zmns_spl
  USE MODgvec_VMEC_Vars,     ONLY:lmnc_spl,lmns_spl
  USE MODgvec_VMEC_Readin,   ONLY:lasym
  USE MODgvec_VMEC,          ONLY:VMEC_EvalSplMode
!$ USE omp_lib
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: which_init_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: U_init
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode,is,i_m,i_n
  REAL(wp) :: BC_val(2)
  REAL(wp) :: spos
  REAL(wp) :: StartTime,EndTime
  REAL(wp) :: dl,lint,x1int,x2int
  REAL(wp) :: X1_b_IP(X1_base%f%mn_nyq(1),X1_base%f%mn_nyq(2))
  REAL(wp) :: X2_b_IP(X2_base%f%mn_nyq(1),X2_base%f%mn_nyq(2))
  REAL(wp) :: X1_gIP(1:X1_base%s%nBase)
  REAL(wp) :: X2_gIP(1:X2_base%s%nBase)
  REAL(wp) :: LA_gIP(1:LA_base%s%nBase,1:LA_base%f%modes)
!===================================================================================================================================
  SELECT CASE(which_init_in)
  CASE(-1) !restart
    X1_a(:)=U_init%X1(1,:)
    X2_a(:)=U_init%X2(1,:)
    X1_b(:)=U_init%X1(X1_base%s%nBase,:)
    X2_b(:)=U_init%X2(X2_base%s%nBase,:)
  CASE(0)
    !X1_a,X2_a and X1_b,X2_b already filled from parameter file readin...
  CASE(1) !VMEC
    IF((init_BC.EQ.-1).OR.(init_BC.EQ.1))THEN ! compute  axis from VMEC, else use the one defined in paramterfile
      spos=0.0_wp
      ASSOCIATE(sin_range    => X1_base%f%sin_range, cos_range    => X1_base%f%cos_range )
      DO imode=cos_range(1)+1,cos_range(2)
        X1_a(iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmnc_Spl)
      END DO
      IF(lasym)THEN
        DO imode=sin_range(1)+1,sin_range(2)
          X1_a(iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmns_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X1
      ASSOCIATE(sin_range    => X2_base%f%sin_range, cos_range    => X2_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        X2_a(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmns_Spl)
      END DO
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          X2_a(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmnc_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X2
    END IF
    IF((init_BC.EQ.-1).OR.(init_BC.EQ.0))THEN ! compute edge from VMEC, else use the one defined in paramterfile
      spos=1.0_wp
      ASSOCIATE(sin_range    => X1_base%f%sin_range, cos_range    => X1_base%f%cos_range )
      DO imode=cos_range(1)+1,cos_range(2)
        X1_b(iMode:iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmnc_Spl)
      END DO
      IF(lasym)THEN
        DO imode=sin_range(1)+1,sin_range(2)
          X1_b(iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,spos,Rmns_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X1
      ASSOCIATE(sin_range    => X2_base%f%sin_range, cos_range    => X2_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        X2_b(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmns_Spl)
      END DO
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          X2_b(iMode)  =VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,spos,Zmnc_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X2
    END IF
    IF(init_average_axis)THEN
      ASSOCIATE(m_nyq=>X1_base%f%mn_nyq(1),n_nyq=>X1_base%f%mn_nyq(2))
      X1_b_IP(:,:) = RESHAPE(X1_base%f%evalDOF_IP(0,X1_b),(/m_nyq,n_nyq/))
      X2_b_IP(:,:) = RESHAPE(X2_base%f%evalDOF_IP(0,X2_b),(/m_nyq,n_nyq/))
      DO i_n=1,n_nyq
        !overwrite axis with average axis by center of closed line of the boundary in each poloidal plane:
        !dl=SQRT((X1_b_IP(1,i_n)-X1_b_IP(m_nyq,i_n))**2+(X2_b_IP(1,i_n)-X2_b_IP(m_nyq,i_n))**2)
        !lint=dl
        !x1int=X1_b_IP(1,i_n)*dl
        !x2int=X2_b_IP(1,i_n)*dl
        !DO i_m=2,m_nyq
        !  dl=SQRT((X1_b_IP(i_m,i_n)-X1_b_IP(i_m-1,i_n))**2+(X2_b_IP(i_m,i_n)-X2_b_IP(i_m-1,i_n))**2)
        !  lint=lint+dl
        !  x1int=x1int+X1_b_IP(i_m,i_n)*dl
        !  x2int=x2int+X2_b_IP(i_m,i_n)*dl
        !END DO
        !overwrite axis with centroid  of surface enclosed  by the line of the boundary in each poloidal plane:
        ! c_x= 1/(6A) sum_i (x_i-1+x_i)*(x_i-1*y_i - x_i*y_i-1), A=1/2 sum_i  (x_i-1*y_i - x_i*y_i-1)
        dl=X1_b_IP(m_nyq,i_n)*X2_b_IP(1,i_n)-X1_b_IP(1,i_n)*X2_b_IP(m_nyq,i_n)
        lint=dl
        x1int=(X1_b_IP(m_nyq,i_n)+X1_b_IP(1,i_n))*dl
        x2int=(X2_b_IP(m_nyq,i_n)+X2_b_IP(1,i_n))*dl
        DO i_m=2,m_nyq
          dl=SQRT((X1_b_IP(i_m,i_n)-X1_b_IP(i_m-1,i_n))**2+(X2_b_IP(i_m,i_n)-X2_b_IP(i_m-1,i_n))**2)
          dl=X1_b_IP(i_m-1,i_n)*X2_b_IP(i_m,i_n)-X1_b_IP(i_m,i_n)*X2_b_IP(i_m-1,i_n)
          lint=lint+dl
          x1int=x1int+(X1_b_IP(i_m-1,i_n)+X1_b_IP(i_m,i_n))*dl
          x2int=x2int+(X2_b_IP(i_m-1,i_n)+X2_b_IP(i_m,i_n))*dl
        END DO
        X1_b_IP(:,i_n) = x1int/(3.0_wp*lint) + average_axis_move(1)
        X2_b_IP(:,i_n) = x2int/(3.0_wp*lint) + average_axis_move(2)
      END DO
      X1_a = X1_base%f%initDOF(RESHAPE(X1_b_IP,(/X1_base%f%mn_IP/)))
      X2_a = X2_base%f%initDOF(RESHAPE(X2_b_IP,(/X2_base%f%mn_IP/)))
      END ASSOCIATE
    END IF !init_average_axis
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
      ASSOCIATE(s_IP         => LA_base%s%s_IP, &
                nBase        => LA_base%s%nBase, &
                sin_range    => LA_base%f%sin_range,&
                cos_range    => LA_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        DO is=1,nBase
          LA_gIP(is,iMode)  =VMEC_EvalSplMode(LA_base%f%Xmn(:,iMode),0,s_IP(is),lmns_Spl)
        END DO !is
      END DO !imode= sin_range
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          DO is=1,nBase
            LA_gIP(is,iMode)  =VMEC_EvalSplMode(LA_base%f%Xmn(:,iMode),0,s_IP(is),Lmnc_Spl)
          END DO !is
        END DO !imode=cos_range
      END IF !lasym
      END ASSOCIATE !X1
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
    !CASE(M_ODD_FIRST,M_ODD,M_EVEN)
    CASE DEFAULT
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
    !CASE(M_ODD_FIRST,M_ODD,M_EVEN)
    CASE DEFAULT
      BC_val =(/          0.0_wp,      X2_b(iMode)/)
    END SELECT !X1(:,iMode) zero odd even
    CALL X2_base%s%applyBCtoDOF(U_init%X2(:,iMode),X2_BC_type(:,iMode),BC_val)
  END DO
  END ASSOCIATE !X2

  IF(init_LA)THEN
    CALL CPU_TIME(StartTime)
!$ StartTime=OMP_GET_WTIME()
    SWRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda from mapping ..."
    !initialize Lambda
    CALL ProgressBar(0,LA_base%s%nBase) !init
    DO is=1,LA_base%s%nBase
      spos=LA_base%s%s_IP(is)
      CALL lambda_Solve(spos,U_init%X1,U_init%X2,LA_gIP(is,:))
      CALL ProgressBar(is,LA_base%s%nBase)
    END DO !is
    SWRITE(UNIT_stdOut,'(A)') "... done."
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
    CALL CPU_TIME(EndTime)
!$ EndTime=OMP_GET_WTIME()
    SWRITE(UNIT_stdOut,'(4X,A,F9.2,A)') " init lambda took [ ",EndTime-StartTime," sec]"
  ELSE
    !lambda init might not be needed since it has no boundary condition and changes anyway after the update of the mapping...
    IF(.NOT.init_fromBConly)THEN
      SWRITE(UNIT_stdOut,'(4X,A)') "... lambda initialized with VMEC ..."
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
      SWRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda =0 ..."
      U_init%LA=0.0_wp
    END IF
  END IF !init_LA

END SUBROUTINE InitSolution


!===================================================================================================================================
!> Add boundary perturbation
!!
!===================================================================================================================================
SUBROUTINE AddBoundaryPerturbation(U_init,h,X1pert_b,X2pert_b)
! MODULES
  USE MODgvec_MHD3D_Vars   , ONLY:X1_base,X1_BC_Type,X1_a,X1_b
  USE MODgvec_MHD3D_Vars   , ONLY:X2_base,X2_BC_Type,X2_a,X2_b
  USE MODgvec_MHD3D_Vars   , ONLY:LA_base,LA_BC_Type,init_LA
  USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_mhd3d
  USE MODgvec_lambda_solve,  ONLY:lambda_solve
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp),INTENT(IN) :: h ! depth of perturbation from boundary (0.1..0.3)
  REAL(wp),INTENT(IN) :: X1pert_b(1:X1_base%f%modes)
  REAL(wp),INTENT(IN) :: X2pert_b(1:X2_base%f%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: U_init
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode,is
  REAL(wp) :: BC_val(2),spos
  REAL(wp) :: X1pert_gIP(1:X1_base%s%nBase)
  REAL(wp) :: X2pert_gIP(1:X2_base%s%nBase)
  REAL(wp) :: LA_gIP(1:LA_base%s%nBase,1:LA_base%f%modes)
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)') "ADD BOUNDARY PERTURBATION..."


  ASSOCIATE(s_IP         =>X1_base%s%s_IP, &
            modes        =>X1_base%f%modes )
  DO imode=1,modes
    X1_b(iMode)=X1_b(iMode)+X1pert_b(iMode)
    X1pert_gIP(:)=blend(s_IP)*X1pert_b(iMode)
    U_init%X1(:,iMode)=U_init%X1(:,iMode) + X1_base%s%initDOF( X1pert_gIP(:) )
  END DO
  END ASSOCIATE

  ASSOCIATE(s_IP         =>X2_base%s%s_IP, &
            modes        =>X2_base%f%modes )
  DO imode=1,modes
    X2_b(iMode)=X2_b(iMode)+X2pert_b(iMode)
    X2pert_gIP(:)=blend(s_IP)*X2pert_b(iMode)
    U_init%X2(:,iMode)=U_init%X2(:,iMode) + X2_base%s%initDOF( X2pert_gIP(:))
  END DO
  END ASSOCIATE

  !apply strong boundary conditions
  ASSOCIATE(modes        =>X1_base%f%modes, &
            zero_odd_even=>X1_base%f%zero_odd_even)
  DO imode=1,modes
    SELECT CASE(zero_odd_even(iMode))
    CASE(MN_ZERO,M_ZERO)
      BC_val =(/ X1_a(iMode)    ,      X1_b(iMode)/)
    !CASE(M_ODD_FIRST,M_ODD,M_EVEN)
    CASE DEFAULT
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
    !CASE(M_ODD_FIRST,M_ODD,M_EVEN)
    CASE DEFAULT
      BC_val =(/          0.0_wp,      X2_b(iMode)/)
    END SELECT !X1(:,iMode) zero odd even
    CALL X2_base%s%applyBCtoDOF(U_init%X2(:,iMode),X2_BC_type(:,iMode),BC_val)
  END DO
  END ASSOCIATE !X2

  IF(init_LA)THEN
    SWRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda from mapping ..."
    !initialize Lambda
    DO is=1,LA_base%s%nBase
      spos=LA_base%s%s_IP(is)
      CALL lambda_Solve(spos,U_init%X1,U_init%X2,LA_gIP(is,:))
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
  END IF !init_LA

  SWRITE(UNIT_stdOut,'(4X,A)') "... DONE."
  SWRITE(UNIT_stdOut,fmt_sep)


  CONTAINS

  ELEMENTAL FUNCTION blend(s_in)
    REAL(wp),INTENT(IN) :: s_in !input coordinate [0,1]
    REAL(wp)            :: blend
    !blend= ( MIN(0., (s_in-(1.-h)) ) / h) **4
    !blend= 2. -2./(EXP(8.*(s_in-1.)/h) +1)
    !blend= EXP(-4.*((s_in-1.)/h)**2)
    !blend= s_in
    blend= EXP(-4.*((s_in-1.)/0.6)**2)
  END FUNCTION blend

END SUBROUTINE AddBoundaryPerturbation

!===================================================================================================================================
!> Compute Equilibrium, iteratively
!!
!===================================================================================================================================
SUBROUTINE MinimizeMHD3D(sf)
! MODULES
  USE MODgvec_MHD3D_vars, ONLY: MinimizerType
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  __PERFON('minimizer')
  SELECT CASE(MinimizerType)
  CASE(0,10)
    CALL MinimizeMHD3D_descent(sf)
  CASE DEFAULT
    CALL abort(__STAMP__,&
        "Minimizertype does not exist",MinimizerType,-1.0_wp)
  END SELECT
  __PERFOFF('minimizer')
END SUBROUTINE MinimizeMHD3D

!===================================================================================================================================
!> Compute Equilibrium, iteratively
!!
!===================================================================================================================================
SUBROUTINE MinimizeMHD3D_descent(sf)
! MODULES
  USE MODgvec_MHD3D_Vars
  USE MODgvec_MHD3D_EvalFunc
  USE MODgvec_Analyze, ONLY:analyze
  USE MODgvec_Restart, ONLY:WriteState
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: iter,nStepDecreased,nSkip_Jac,nSkip_dw
  INTEGER   :: JacCheck,lastoutputIter,StartTimeArray(8)
  REAL(wp)  :: dt,deltaW,absTol
  INTEGER,PARAMETER   :: ndamp=10
  REAL(wp)  :: tau(1:ndamp), tau_bar
  REAL(wp)  :: min_dt_out,max_dt_out,min_dw_out,max_dw_out,sum_dW_out,t_pseudo,Fnorm(3),Fnorm0(3),Fnorm_old(3),W_MHD3D_0
  INTEGER   :: logUnit !globally needed for logging
  INTEGER   :: logiter_ramp,logscreen
  LOGICAL   :: restart_iter
  LOGICAL   :: first_iter
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)') "MINIMIZE MHD3D FUNCTIONAL..."


  abstol=minimize_tol


  dt=start_dt
  nstepDecreased=0
  nSkip_Jac=0
  t_pseudo=0
  lastOutputIter=0
  iter=0
  logiter_ramp=1
  logscreen=1

  first_iter=.TRUE.
  restart_iter=.FALSE.

  CALL U(-3)%set_to(U(0)) !initial state, should remain unchanged

  DO WHILE(iter.LE.maxIter)
    IF((first_iter).OR.(restart_iter))THEN
      JacCheck=1 !abort if detJ<0
      CALL EvalAux(           U(0),JacCheck)
      U(0)%W_MHD3D=EvalEnergy(U(0),.FALSE.,JacCheck)
      W_MHD3D_0 = U(0)%W_MHD3D
      CALL EvalForce(         U(0),.FALSE.,JacCheck,F(0))
      Fnorm0=SQRT(F(0)%norm_2())
      Fnorm=Fnorm0
      Fnorm_old=1.1*Fnorm0
      CALL U(-1)%set_to(U(0)) !last state
      CALL U(-2)%set_to(U(0)) !state at last logging interval
      !for hirshman method
      IF(MinimizerType.EQ.10)THEN
        CALL V(-1)%set_to(0.0_wp)
        CALL V( 0)%set_to(0.0_wp)
        tau(1:ndamp)=0.15_wp/dt
      END IF
      min_dt_out=1.0e+30_wp
      max_dt_out=0.0_wp
      min_dW_out=1.0e+30_wp
      max_dW_out=-1.0e+30_wp
      sum_dW_out=0.0_wp
      nSkip_dW =0
      IF(restart_iter) restart_iter=.FALSE.
      IF(first_iter)THEN
        CALL StartLogging()
        first_iter=.FALSE.
      END IF
    END IF !before first iteration or after restart Jac<0

    !COMPUTE NEW SOLUTION P(1) as a prediction

    SELECT CASE(MinimizerType)
    CASE(0) !gradient descent, previously used for minimizerType=0
      CALL P(1)%AXBY(1.0_wp,U(0),dt,F(0)) !overwrites P(1), predicts solution U(1)
    CASE(10) !hirshman method
      !tau is damping parameter
      tau(1:ndamp-1) = tau(2:ndamp) !save old
      tau(ndamp)  = MIN(0.15_wp,ABS(LOG(SUM(Fnorm**2)/SUM(Fnorm_old**2))))/dt  !ln(|F_n|^2/|F_{n-1}|^2), Fnorm=|F_X1|,|F_X2|,|F_LA|
      tau_bar = 0.5*dt*SUM(tau)/REAL(ndamp,wp)   !=1/2 * tauavg
      CALL V(1)%AXBY(((1.0_wp-tau_bar)/(1.0_wp+tau_bar)),V(0),(dt/(1.0_wp+tau_bar)),F(0)) !velocity V(1)
      CALL P(1)%AXBY(1.0_wp,U(0),dt,V(1)) !overwrites P(1), predicst solution U(1)
    END SELECT


    JacCheck=2 !no abort,if detJ<0, JacCheck=-1
    P(1)%W_MHD3D=EvalEnergy(P(1),.TRUE.,JacCheck)
    IF(JacCheck.EQ.-1)THEN
      dt=0.9_wp*dt
      nstepDecreased=nStepDecreased+1
      nSkip_Jac=nSkip_Jac+1
      restart_iter=.TRUE.
      CALL U(0)%set_to(U(-3)) !reset to initial state
      SWRITE(UNIT_stdOut,'(8X,I8,A,E11.4,A)')iter,'...detJac<0, decrease stepsize to dt=',dt,  ' and RESTART simulation!!!!!!!'
    ELSE
      !detJ>0
      deltaW=P(1)%W_MHD3D-U(0)%W_MHD3D!should be <=0,
      IF(deltaW.LE.dW_allowed*W_MHD3D_0)THEN !valid step /hirshman method accept W increase!
        IF(doLineSearch)THEN
          ! LINE SEARCH ... SEEMS THAT IT NEVER REALLY REDUCES THE STEP SIZE... NOT NEEDED?
          CALL U(1)%AXBY(0.5_wp,P(1),0.5_wp,U(0)) !overwrites U(1)
          JacCheck=2 !no abort,if detJ<0, JacCheck=-1, if detJ>0 Jaccheck=1
          U(1)%W_MHD3D=EvalEnergy(U(1),.TRUE.,JacCheck)
          IF((U(1)%W_MHD3D.LT.U(0)%W_MHD3D).AND.(U(1)%W_MHD3D.LT.P(1)%W_MHD3D).AND.(JacCheck.EQ.1))THEN
            CALL P(1)%set_to(U(1)) !accept smaller step
            CALL U(1)%AXBY(0.5_wp,P(1),0.5_wp,U(0)) !overwrites U(1)
            JacCheck=2 !no abort,if detJ<0, JacCheck=-1, if detJ>0 Jaccheck=1
            U(1)%W_MHD3D=EvalEnergy(U(1),.TRUE.,JacCheck)
            IF((U(1)%W_MHD3D.LT.U(0)%W_MHD3D).AND.(U(1)%W_MHD3D.LT.P(1)%W_MHD3D).AND.(JacCheck.EQ.1))THEN
              !SWRITE(UNIT_stdOut,'(8X,I8,A)')iter,' linesearch: 1/4 step!'
              CALL P(1)%set_to(U(1)) !accept smaller step
            ELSE
              !SWRITE(UNIT_stdOut,'(8X,I8,A)')iter,' linesearch: 1/2 step!'
            END IF
          END IF
        END IF !dolinesearch

        IF(ALL(Fnorm.LE.abstol))THEN
          CALL Logging(.FALSE.)
          SWRITE(UNIT_stdOut,'(4x,A)')'==>Iteration finished, |force| in relative tolerance'
          EXIT !DO LOOP
        END IF
        iter=iter+1
        nstepDecreased=0
        min_dt_out=MIN(min_dt_out,dt)
        max_dt_out=MAX(max_dt_out,dt)
        min_dW_out=MIN(min_dW_out,deltaW)
        max_dW_out=MAX(max_dW_out,deltaW)
        sum_dW_out=sum_dW_out+deltaW
        IF(MOD(iter,logIter_ramp).EQ.0)THEN

          CALL Logging(.NOT.((logIter_ramp.GE.logIter).AND.(MOD(logscreen,nLogScreen).EQ.0)))
          IF(.NOT.(logIter_ramp.LT.logIter))THEN !only reset for logIter
            logscreen=logscreen+1
            min_dt_out=1.0e+30_wp
            max_dt_out=0.0_wp
            min_dW_out=1.0e+30_wp
            max_dW_out=-1.0e+30_wp
            sum_dW_out=0.0_wp
            nSkip_dW =0
          END IF
          logIter_ramp=MIN(logIter,logIter_ramp*2)
        END IF

        t_pseudo=t_pseudo+dt
        ! for simple gradient & hirshman
        CALL U(-1)%set_to(U(0))
        CALL U(0)%set_to(P(1))
        ! for hirshman method
        IF(MinimizerType.EQ.10)THEN
          CALL V(-1)%set_to(V(0))
          CALL V(0)%set_to(V(1))
        END IF

        CALL EvalForce(P(1),.FALSE.,JacCheck,F(0)) !evalAux was already called on P(1)=U(0), so that its set false here.
        Fnorm_old=Fnorm
        Fnorm=SQRT(F(0)%norm_2())

      ELSE !not a valid step, decrease timestep and skip P(1)
        dt=0.9_wp*dt
        nstepDecreased=nStepDecreased+1
        nSkip_dW=nSkip_dW+1
        !CALL U(0)%set_to(U(-2))
        restart_iter=.TRUE.
        SWRITE(UNIT_stdOut,'(8X,I8,A,E8.1,A,E11.4)')iter,'...deltaW>',dW_allowed,'*W_MHD3D_0, skip step and decrease stepsize to dt=',dt
      END IF
    END IF !JacCheck

    IF(nStepDecreased.GT.20) THEN ! 2^20 ~10^6
      SWRITE(UNIT_stdOut,'(A,E21.11)')'Iteration stopped since timestep has been decreased by 2^20: ', dt
      SWRITE(UNIT_stdOut,fmt_sep)
      RETURN
    END IF
    IF((MOD(iter,outputIter).EQ.0).AND.(lastoutputIter.NE.iter))THEN
      __PERFON('output')
      SWRITE(UNIT_stdOut,'(A)')'##########################  OUTPUT ##################################'
      CALL Analyze(iter)
      CALL WriteState(U(0),iter)
      CALL CheckEvalForce(U(0),iter)
      SWRITE(UNIT_stdOut,'(A)')'#####################################################################'
      lastOutputIter=iter
      __PERFOFF('output')
    END IF
  END DO !iter
  IF(iter.GE.MaxIter)THEN
    SWRITE(UNIT_stdOut,'(A,E21.11)')"maximum iteration count exceeded, not converged"
  END IF
  SWRITE(UNIT_stdOut,'(A)') "... DONE."
  SWRITE(UNIT_stdOut,fmt_sep)
  CALL Analyze(MIN(iter,MaxIter))
  CALL WriteState(U(0),MIN(iter,MaxIter))
  CALL FinishLogging()
!DEBUG
!  WRITE(FileString,'(A,"_State_",I4.4,"_",I8.8,".dat")')TRIM(ProjectName),OutputLevel,99999999
!  CALL ReadState(FileString,U(-1))


CONTAINS

  !=================================================================================================================================
  !> all screen and logfile tasks, can use all variables from subroutine above
  !!
  !=================================================================================================================================
  SUBROUTINE StartLogging()
  USE MODgvec_Globals,     ONLY: GETFREEUNIT
  USE MODgvec_Output_Vars, ONLY: ProjectName,outputLevel
  USE MODgvec_MHD3D_visu,  ONLY: checkAxis
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  CHARACTER(LEN=255)  :: fileString
  INTEGER             :: TimeArray(8),iLogDat
  REAL(wp)            :: AxisPos(2,2)
  INTEGER,PARAMETER   :: nLogDat=16
  REAL(wp)            :: LogDat(1:nLogDat)
  !=================================================================================================================================
  __PERFON('log_output')
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  SWRITE(UNIT_stdOut,'(A,E11.4,A)')'%%%%%%%%%%  START ITERATION, dt= ',dt, '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  SWRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                 '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
  SWRITE(UNIT_stdOut,'(A,3E21.14)') &
          '%%% dU = |Force|= ',Fnorm(1:3)
  SWRITE(UNIT_stdOut,'(40(" -"))')
  !------------------------------------
  StartTimeArray=TimeArray !save first time stamp

  logUnit=GETFREEUNIT()
  WRITE(FileString,'("logMinimizer_",A,"_",I4.4,".csv")')TRIM(ProjectName),outputLevel
  OPEN(UNIT     = logUnit       ,&
     FILE     = TRIM(FileString) ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' )
  !header
  iLogDat=0
  WRITE(logUnit,'(A)',ADVANCE="NO")'"#iterations","runtime(s)","min_dt","max_dt"'
  WRITE(logUnit,'(A)',ADVANCE="NO")',"W_MHD3D","min_dW","max_dW","sum_dW"'
  WRITE(logUnit,'(A)',ADVANCE="NO")',"normF_X1","normF_X2","normF_LA"'
  LogDat(ilogDat+1:iLogDat+11)=(/0.0_wp,0.0_wp,dt,dt,U(0)%W_MHD3D,0.0_wp,0.0_wp,0.0_wp,Fnorm(1:3)/)
  iLogDat=11
  IF(doCheckDistance) THEN
    WRITE(logUnit,'(A)',ADVANCE="NO")',"max_Dist","avg_Dist"'
    LogDat(iLogDat+1:iLogDat+2)=(/0.0_wp,0.0_wp/)
    iLogDat=iLogDat+2
  END IF!doCheckDistance
  IF(doCheckAxis) THEN
    WRITE(logUnit,'(A)',ADVANCE="NO")',"X1_axis_0","X2_axis_0","X1_axis_1","X2_axis_1"'
    CALL CheckAxis(U(0),2,AxisPos)
    LogDat(iLogDat+1:iLogDat+4)=RESHAPE(AxisPos,(/4/))
    iLogDat=iLogDat+4
  END IF!doCheckAxis
  WRITE(logUnit,'(A)')' '
  !first data line
  WRITE(logUnit,'(*(e23.15,:,","))') logDat(1:iLogDat)
  __PERFOFF('log_output')
  END SUBROUTINE StartLogging

  !=================================================================================================================================
  !> all screen and logfile tasks, can use all variables from subroutine above
  !!
  !=================================================================================================================================
  SUBROUTINE Logging(quiet)
  USE MODgvec_MHD3D_visu, ONLY: checkDistance
  USE MODgvec_MHD3D_visu, ONLY: checkAxis
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: quiet !! True: no screen output
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER             :: TimeArray(8),runtime_ms,iLogDat
  REAL(wp)            :: AxisPos(2,2),maxDist,avgDist
  INTEGER,PARAMETER   :: nLogDat=16
  REAL(wp)            :: LogDat(1:nLogDat)
  !=================================================================================================================================
  __PERFON('log_output')
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  IF(.NOT.quiet)THEN
    SWRITE(UNIT_stdOut,'(80("%"))')
    SWRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                      '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
    SWRITE(UNIT_stdOut,'(A,I8,A,2I8,A,E11.4,A,2E11.4,A,E21.14,A,3E12.4)') &
                      '%%% #ITERATIONS= ',iter,', #skippedIter (Jac/dW)= ',nSkip_Jac,nSkip_dW, &
              '\n%%% t_pseudo= ',t_pseudo,', min/max dt= ',min_dt_out,max_dt_out, &
              '\n%%% W_MHD3D= ',U(0)%W_MHD3D,', min/max/sum deltaW= ' , min_dW_out,max_dW_out,sum_dW_out
    SWRITE(UNIT_stdOut,'(A,3E21.14)') &
                '%%% dU = |Force|= ',Fnorm(1:3)
    !------------------------------------
  END IF!.NOT.quiet
  iLogDat=0
  runtime_ms=MAX(0,SUM((timeArray(5:8)-StartTimearray(5:8))*(/360000,6000,100,1/)))
  LogDat(ilogDat+1:iLogDat+11)=(/REAL(iter,wp),REAL(runtime_ms,wp)/100.0_wp, &
                                min_dt_out,max_dt_out,U(0)%W_MHD3D,min_dW_out,max_dW_out,sum_dW_out, &
                                Fnorm(1:3)/)
  iLogDat=11
  IF(doCheckDistance) THEN
    CALL CheckDistance(U(0),U(-2),maxDist,avgDist)
    CALL U(-2)%set_to(U(0))
    IF(.NOT.quiet)THEN
      SWRITE(UNIT_stdOut,'(A,2E11.4)') &
      '               %%% Dist to last log (max/avg) : ',maxDist,avgDist
    END IF!.NOT.quiet
    LogDat(iLogDat+1:iLogDat+2)=(/maxDist,avgDist/)
    iLogDat=iLogDat+2
  END IF!doCheckDistance
  IF(doCheckAxis) THEN
    CALL CheckAxis(U(0),2,AxisPos)
    IF(.NOT.quiet)THEN
      SWRITE(UNIT_stdOut,'(2(A,2E22.14))') &
        '%%% axis position (X1,X2,zeta=0     ): ',AxisPos(1:2,1), &
      '\n%%% axis position (X1,X2,zeta=pi/nfp): ',AxisPos(1:2,2)
    END IF!.NOT.quiet
    LogDat(iLogDat+1:iLogDat+4)=RESHAPE(AxisPos,(/4/))
    iLogDat=iLogDat+4
  END IF !doCheckAxis
  IF(.NOT.quiet)THEN
    SWRITE(UNIT_stdOut,'(40(" -"))')
  END IF!.NOT.quiet
  WRITE(logUnit,'(*(e23.15,:,","))') logDat(1:iLogDat)
  __PERFOFF('log_output')
  END SUBROUTINE Logging

  !=================================================================================================================================
  !>
  !!
  !=================================================================================================================================
  SUBROUTINE FinishLogging()
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  CLOSE(logUnit)
  END SUBROUTINE FinishLogging

END SUBROUTINE MinimizeMHD3D_descent


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE FinalizeMHD3D(sf)
! MODULES
  USE MODgvec_MHD3D_Vars
  USE MODgvec_MHD3D_EvalFunc,ONLY:FinalizeMHD3D_EvalFunc
  USE MODgvec_VMEC,ONLY: FinalizeVMEC
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
    CALL V(i)%free()
  END DO
  DO i=-1,0
    CALL F(i)%free()
  END DO
  CALL sgrid%free()

  SDEALLOCATE(U)
  SDEALLOCATE(P)
  SDEALLOCATE(V)
  SDEALLOCATE(F)
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
  IF(which_init.EQ.1) CALL FinalizeVMEC()

  CALL hmap%free()
  SDEALLOCATE(hmap)

  SDEALLOCATE(X1_base)
  SDEALLOCATE(X2_base)
  SDEALLOCATE(LA_base)
END SUBROUTINE FinalizeMHD3D

END MODULE MODgvec_MHD3D
