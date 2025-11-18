!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
  USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep,MPIRoot,enter_subregion,exit_subregion
  USE MODgvec_c_functional,   ONLY: t_functional
  IMPLICIT NONE

  PRIVATE
  PUBLIC t_functional_mhd3d

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES
!-----------------------------------------------------------------------------------------------------------------------------------

  TYPE,EXTENDS(t_functional) :: t_functional_mhd3d
    !-------------------------------------------------------------------------------------------------------------------------------
    LOGICAL :: initialized
    !-------------------------------------------------------------------------------------------------------------------------------
    CONTAINS
      PROCEDURE :: init         => InitMHD3D
      PROCEDURE :: initSolution => InitSolutionMHD3D
      PROCEDURE :: minimize     => MinimizeMHD3D
      PROCEDURE :: free         => FinalizeMHD3D
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
  USE MODgvec_base           , ONLY: base_new
  USE MODgvec_boundaryFromFile, ONLY: boundaryFromFile_new
  USE MODgvec_hmap           , ONLY: hmap_new,hmap_new_auxvar
  USE MODgvec_VMEC           , ONLY: InitVMEC
  USE MODgvec_VMEC_vars      , ONLY: vmec_iota_profile,vmec_pres_profile
  USE MODgvec_VMEC_Readin    , ONLY: nfp,nFluxVMEC,Phi,xm,xn,lasym,mpol,ntor !<<< only exists on MPIroot!
  USE MODgvec_MHD3D_EvalFunc , ONLY: InitializeMHD3D_EvalFunc
  USE MODgvec_ReadInTools    , ONLY: GETSTR,GETLOGICAL,GETINT,GETINTARRAY,GETREAL,GETREALALLOCARRAY, GETREALARRAY
  USE MODgvec_MPI            , ONLY: par_BCast,par_barrier
  USE MODgvec_rProfile_poly  , ONLY: t_rProfile_poly
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER          :: i,iMode,nElems,n_sgrid_rho
  INTEGER          :: grid_type
  REAL(wp),ALLOCATABLE :: sgrid_rho(:)
  INTEGER          :: X1X2_deg,X1X2_cont
  INTEGER          :: X1_mn_max(2),X2_mn_max(2)
  INTEGER          :: LA_deg,LA_cont,LA_mn_max(2)
  CHARACTER(LEN=8) :: X1_sin_cos
  CHARACTER(LEN=8) :: X2_sin_cos
  CHARACTER(LEN=8) :: LA_sin_cos
  INTEGER          :: degGP,mn_nyq(2),mn_nyq_min(2),fac_nyq
  INTEGER          :: nfp_loc
  INTEGER          :: X1X2_BCtype_axis(0:4),LA_BCtype_axis(0:4)
  INTEGER          :: proposal_mn_max(1:2)=(/2,0/) !!default proposals, changed for VMEC input to automatically match input!
  CHARACTER(LEN=8) :: proposal_X1_sin_cos="_cos_"  !!default proposals, changed for VMEC input to automatically match input!
  CHARACTER(LEN=8) :: proposal_X2_sin_cos="_sin_"  !!default proposals, changed for VMEC input to automatically match input!
  CHARACTER(LEN=8) :: proposal_LA_sin_cos="_sin_"  !!default proposals, changed for VMEC input to automatically match input!
  REAL(wp)         :: scale_minor_radius
  CHARACTER(LEN=255) ::boundary_filename
  CHARACTER(LEN=8) :: boundary_perturb_type_str !! readin variable for boundary_perturb_type: legacy, cosm
!===================================================================================================================================
  CALL enter_subregion("init-MHD3D")
  CALL par_Barrier(beforeScreenOut='INIT MHD3D ...')

  which_init = GETINT("whichInitEquilibrium",0)
  IF(which_init.EQ.1) CALL InitVMEC()

  !-----------MINIMIZER
  MinimizerType= GETINT("MinimizerType",Proposal=10)
  PrecondType  = GETINT("PrecondType",Proposal=1)

  dW_allowed=GETREAL("dW_allowed",Proposal=1.0e-10_wp) !! for minimizer, accept step if dW<dW_allowed*W_MHD(iter=0) default +10e-10
  maxIter   = GETINT("maxIter",Proposal=5000)
  outputIter= GETINT("outputIter",Proposal=maxIter)
  logIter   = GETINT("logIter",Proposal=maxIter)
  nlogScreen= GETINT("nLogScreen",Proposal=1)
  minimize_tol  =GETREAL("minimize_tol",Proposal=1.0e-12_wp)
  start_dt  =GETREAL("start_dt",Proposal=0.1_wp)
  doCheckDistance=GETLOGICAL("doCheckDistance",Proposal=.FALSE.)
  doCheckAxis=GETLOGICAL("doCheckAxis",Proposal=.TRUE.)
  !-----------

  fac_nyq = GETINT( "fac_nyq",Proposal=4)

  !constants
  mu_0    = 2.0e-07_wp*TWOPI
  gamm    = 0.0_wp  !fixed!
  sgammM1=1.0_wp/(gamm-1.0_wp)

  init_LA= GETLOGICAL("init_LA",Proposal=.TRUE.)

  SELECT CASE(which_init)
  CASE(0)
    init_fromBConly= .TRUE.
    init_BC        = 2

    !hmap
    which_hmap=GETINT("which_hmap",Proposal=1)
    CALL hmap_new(hmap,which_hmap)
    IF(hmap%nfp.NE.-1)THEN
      nfp_loc = hmap%nfp
    ELSE
      nfp_loc = GETINT("nfp")
    END IF

    Phi_edge   = GETREAL("PHIEDGE",Proposal=1.0_wp)
    Phi_edge   = Phi_edge/TWOPI !normalization like in VMEC!!!
    IF(MPIroot)THEN
      CALL InitProfile(sf,"iota",iota_profile)
      CALL InitProfile(sf,"pres",pres_profile)
    END IF

  CASE(1) !VMEC init
    init_fromBConly= GETLOGICAL("init_fromBConly",Proposal=.FALSE.)
    IF(init_fromBConly)THEN
      !=-1, keep vmec axis and boundary, =0: keep vmec boundary, overwrite axis, =1: keep vmec axis, overwrite boundary, =2: overwrite axis and boundary
      init_BC= GETINT("reinit_BC",Proposal=-1)
    ELSE
      init_BC=-1
    END IF

    IF(MPIroot)THEN
      init_with_profile_iota     = GETLOGICAL("init_with_profile_iota", Proposal=.FALSE.)
      IF(init_with_profile_iota)THEN
        CALL InitProfile(sf,"iota",iota_profile)
      ELSE
        iota_profile=vmec_iota_profile
      END IF ! iota from parameterfile
      init_with_profile_pressure = GETLOGICAL("init_with_profile_pressure", Proposal=.FALSE.)
      IF(init_with_profile_pressure)THEN
        CALL InitProfile(sf,"pres",pres_profile)
      ELSE
        pres_profile=vmec_pres_profile
      END IF ! pressure from parameterfile

      proposal_mn_max(:)=(/mpol-1,ntor/)
      IF(lasym)THEN !asymmetric
        proposal_X1_sin_cos="_sincos_"
        proposal_X2_sin_cos="_sincos_"
        proposal_LA_sin_cos="_sincos_"
      END IF
      nfp_loc = nfp
      Phi_edge = Phi(nFluxVMEC)
    END IF !MPIroot

    which_hmap=1 !hmap_RZ
    CALL hmap_new(hmap,which_hmap)
  END SELECT !which_init

  IF(MPIroot)THEN
    Phi_profile=t_rProfile_poly((/0.0_wp,Phi_edge/)) !! choice phi=Phi_edge * s
    !iota= (dChi/ds) / (dPhi/ds) = dchi_ds / Phi_edge  => chi = Phi_edge * int(iota ds)
    chi_profile=iota_profile%antiderivative()
    chi_profile%coefs=chi_profile%coefs*Phi_edge
  END IF !MPIroot

  getBoundaryFromFile=GETINT("getBoundaryFromFile",Proposal=-1)  ! =-1: OFF, get X1b and X2b from parameterfile. 1: get boundary from specific netcdf file
  SELECT CASE(getBoundaryFromFile)
  CASE(-1)
    !do nothing
  CASE(1)
    boundary_filename=GETSTR("boundary_filename")
    scale_minor_radius=GETREAL("scale_minor_radius",Proposal=1.0_wp)
    IF(MPIroot)THEN
      CALL boundaryFromFile_new(BFF,boundary_filename)
      IF(nfp_loc.NE.BFF%nfp) WRITE(UNIT_stdOut,'(6X,A,I4)')'INFO: changed to boundary file NFP= ',BFF%nfp
      nfp_loc=BFF%nfp
      proposal_mn_max(:)=(/BFF%m_max,BFF%n_max/)
      IF(BFF%lasym.EQ.1)THEN !asymmetric
        proposal_X1_sin_cos="_sincos_"
        proposal_X2_sin_cos="_sincos_"
        proposal_LA_sin_cos="_sincos_"
      END IF
    END IF
  END SELECT
  CALL par_BCast(proposal_mn_max,0)
  CALL par_BCast(proposal_X1_sin_cos,0)
  CALL par_BCast(proposal_X2_sin_cos,0)
  CALL par_BCast(proposal_LA_sin_cos,0)
  CALL par_BCast(nfp_loc,0)

  CALL enter_subregion("discretization")
  X1X2_deg     = GETINT(     "X1X2_deg")
  X1X2_cont    = GETINT(     "X1X2_continuity",Proposal=(X1X2_deg-1) )
  X1_mn_max    = GETINTARRAY("X1_mn_max"   ,2 ,Proposal=proposal_mn_max)
  X2_mn_max    = GETINTARRAY("X2_mn_max"   ,2 ,Proposal=proposal_mn_max)
  X1_sin_cos   = GETSTR(     "X1_sin_cos"     ,Proposal=proposal_X1_sin_cos)  !_sin_,_cos_,_sin_cos_
  X2_sin_cos   = GETSTR(     "X2_sin_cos"     ,Proposal=proposal_X2_sin_cos)


  LA_deg     = GETINT(     "LA_deg")
  LA_cont    = GETINT(     "LA_continuity",Proposal=(LA_deg-1))
  LA_mn_max  = GETINTARRAY("LA_mn_max", 2 ,Proposal=proposal_mn_max)
  LA_sin_cos = GETSTR(     "LA_sin_cos"   ,Proposal=proposal_LA_sin_cos)

  IF(fac_nyq.EQ.-1)THEN
    fac_nyq=2
    mn_nyq_min(1)=1+fac_nyq*MAXVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/))
    mn_nyq_min(2)=1+fac_nyq*(MAXVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/))+hmap%n_max)
    mn_nyq  = GETINTARRAY("mn_nyq",2)
    IF(mn_nyq(1).LT.mn_nyq_min(1))THEN
       SWRITE(UNIT_stdOut,'(A,I6)')'WARNING: mn_nyq(1) too small, should be >= ',mn_nyq_min(1)
    END IF
    IF(mn_nyq(2).LT.mn_nyq_min(2))THEN
       SWRITE(UNIT_stdOut,'(A,I6)') 'WARNING: mn_nyq(2) too small, should be >= ',mn_nyq_min(2)
    END IF
  ELSE
    mn_nyq(1)=1+fac_nyq*MAXVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/))
    mn_nyq(2)=1+fac_nyq*(MAXVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/))+hmap%n_max)
  END IF

  SWRITE(UNIT_stdOut,*)
  SWRITE(UNIT_stdOut,'(2(A,I4),A,I6," , ",I6,A)')'    fac_nyq = ', fac_nyq,' hmap%n_max = ',hmap%n_max,'  ==> interpolation points mn_nyq=( ',mn_nyq(:),' )'
  SWRITE(UNIT_stdOut,*)

  grid_type= GETINT("sgrid_grid_type",Proposal=0)
  degGP    = GETINT("degGP",Proposal=MAX(X1X2_deg,LA_deg)+2)

  !INITIALIZE GRID
  IF (grid_type.NE.GRID_TYPE_CUSTOM) THEN
    nElems   = GETINT("sgrid_nElems")
    CALL sgrid%init(nElems,grid_type)
  ELSE
    CALL GETREALALLOCARRAY("sgrid_rho", sgrid_rho, n_sgrid_rho)
    nElems = n_sgrid_rho - 1
    CALL sgrid%init(nElems,grid_type,sgrid_rho)
    SDEALLOCATE(sgrid_rho)
  END IF

  !INITIALIZE BASE        !sbase parameter                 !fbase parameter               ...exclude_mn_zero
  CALL base_new(X1_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X1_mn_max,mn_nyq,nfp_loc,X1_sin_cos,.FALSE.)
  CALL base_new(X2_base  , X1X2_deg,X1X2_cont,sgrid,degGP , X2_mn_max,mn_nyq,nfp_loc,X2_sin_cos,.FALSE.)
  CALL base_new(LA_base  ,   LA_deg,  LA_cont,sgrid,degGP , LA_mn_max,mn_nyq,nfp_loc,LA_sin_cos,.TRUE. )

  CALL hmap_new_auxvar(hmap,X1_base%f%zeta_IP,hmap_auxvar,.FALSE.) !no second derivative needed!

  IF((which_init.EQ.1).AND.MPIroot) THEN !VMEC
    IF(lasym)THEN
      IF((X1_base%f%sin_cos.NE._SINCOS_).OR. &
         (X2_base%f%sin_cos.NE._SINCOS_).OR. &
         (LA_base%f%sin_cos.NE._SINCOS_) ) THEN
        WRITE(UNIT_stdOut,'(A)')'!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
        WRITE(UNIT_stdOut,'(A)')'!!!!!!!!   ---->  VMEC was run asymmetric, you should use _sincos_ basis for all variables'
        WRITE(UNIT_stdOut,'(A)')'!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
        !CALL abort(__STAMP__,&
        !    '!!!!  VMEC was run asymmetric, you should use _sincos_ basis for all variables')
      END IF
    END IF
    IF((MAXVAL(INT(xm(:))).GT.MINVAL((/X1_mn_max(1),X2_mn_max(1),LA_mn_max(1)/))).OR. &
       (MAXVAL(ABS(INT(xn(:))/nfp_loc)).GT.MINVAL((/X1_mn_max(2),X2_mn_max(2),LA_mn_max(2)/))))THEN
      WRITE(UNIT_stdOut,'(A)')    '!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
      WRITE(UNIT_stdOut,'(A,2I6)')'!!!!!!!!   ---->  you use a lower mode number than the VMEC  run  ', &
                                    MAXVAL(INT(xm(:))),MAXVAL(ABS(INT(xn(:))/nfp_loc))
      WRITE(UNIT_stdOut,'(A)')    '!!!!!!!! WARNING: !!!!!!!!!!!!!!!'
      !  CALL abort(__STAMP__,&
      !'!!!!!  you use a lower mode number than the VMEC  run  (m,n)_max')
    END IF
  END IF

  nDOF_X1 = X1_base%s%nBase* X1_base%f%modes
  nDOF_X2 = X2_base%s%nBase* X2_base%f%modes
  nDOF_LA = LA_base%s%nBase* LA_base%f%modes



  !X1X2_BCtype_axis(MN_ZERO    )= GETINT("X1X2_BCtype_axis_mn_zero"    ,Proposal=0 ) !AUTOMATIC,m-dependent
  !X1X2_BCtype_axis(M_ZERO     )= GETINT("X1X2_BCtype_axis_m_zero"     ,Proposal=0 ) !AUTOMATIC,m-dependent
  !X1X2_BCtype_axis(M_ODD_FIRST)= GETINT("X1X2_BCtype_axis_m_odd_first",Proposal=0 ) !AUTOMATIC,m-dependent
  !X1X2_BCtype_axis(M_ODD      )= GETINT("X1X2_BCtype_axis_m_odd"      ,Proposal=0 ) !AUTOMATIC,m-dependent
  !X1X2_BCtype_axis(M_EVEN     )= GETINT("X1X2_BCtype_axis_m_even"     ,Proposal=0 ) !AUTOMATIC,m-dependent
  X1X2_BCtype_axis= 0 !fix to AUTOMATIC, m-dependent

  !boundary conditions (used in force, in init slightly changed)
  ASSOCIATE(modes        =>X1_base%f%modes, zero_odd_even=>X1_base%f%zero_odd_even)
  ALLOCATE(X1_BC_type(1:2,modes))
  X1_BC_type(BC_EDGE,:)=BC_TYPE_DIRICHLET
  DO imode=1,modes
    X1_BC_type(BC_AXIS,iMode)=X1X2_BCtype_axis(zero_odd_even(iMode))
    IF(X1_BC_type(BC_AXIS,iMode).EQ.0)THEN !AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
      X1_BC_type(BC_AXIS,iMode)=-1*MIN(X1_base%s%deg+1,X1_base%f%Xmn(1,iMode))
      IF(MPIroot.AND.(nElems.EQ.1).AND.(X1_base%f%Xmn(1,iMode).GT.X1_base%s%deg)) THEN
        IF(X1_base%f%Xmn(2,iMode).EQ.0) & !warning for all n-modes written once!
           WRITE(UNIT_stdOut,'(4X,A,I4,A)')'WARNING, 1-element spline with BC for m>deg, will ZERO edge coeff. X1_b(m=',&
                                          X1_base%f%Xmn(1,iMode),',n=-n_max,n_max)! (use 2elems instead)'
      END IF
    END IF
  END DO
  END ASSOCIATE !X1

  ASSOCIATE(modes        =>X2_base%f%modes, zero_odd_even=>X2_base%f%zero_odd_even)
  ALLOCATE(X2_BC_type(1:2,modes))
  X2_BC_type(BC_EDGE,:)=BC_TYPE_DIRICHLET
  DO imode=1,modes
    X2_BC_type(BC_AXIS,iMode)=X1X2_BCtype_axis(zero_odd_even(iMode))
    IF(X2_BC_type(BC_AXIS,iMode).EQ.0)THEN ! AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
      X2_BC_type(BC_AXIS,iMode)=-1*MIN(X2_base%s%deg+1,X2_base%f%Xmn(1,iMode))
      IF(MPIroot.AND.(nElems.EQ.1).AND.(X2_base%f%Xmn(1,iMode).GT.X2_base%s%deg)) THEN
        IF(X2_base%f%Xmn(2,iMode).EQ.0) & !warning for all n-modes written once!
          WRITE(UNIT_stdOut,'(4X,A,I4,A)')'WARNING, 1-element spline with BC for m>deg, will ZERO edge coeff. X2_b(m=',&
                                          X2_base%f%Xmn(1,iMode),',n=-n_max,n_max)! (use 2elems instead)'
      END IF
    END IF
  END DO
  END ASSOCIATE !X2

  !LA_BCtype_axis(MN_ZERO    )= GETINT("LA_BCtype_axis_mn_zero"    ,Proposal=BC_TYPE_SYMMZERO )
  !LA_BCtype_axis(M_ZERO     )= GETINT("LA_BCtype_axis_m_zero"     ,Proposal=BC_TYPE_SYMM     )
  !LA_BCtype_axis(M_ODD_FIRST)= GETINT("LA_BCtype_axis_m_odd_first",Proposal=0 ) !AUTOMATIC,m-dependent
  !LA_BCtype_axis(M_ODD      )= GETINT("LA_BCtype_axis_m_odd"      ,Proposal=0 ) !AUTOMATIC,m-dependent
  !LA_BCtype_axis(M_EVEN     )= GETINT("LA_BCtype_axis_m_even"     ,Proposal=0 ) !AUTOMATIC,m-dependent
  LA_BCtype_axis= 0 !fix to AUTOMATIC, m-dependent

  ASSOCIATE(modes        =>LA_base%f%modes, zero_odd_even=>LA_base%f%zero_odd_even)
  ALLOCATE(LA_BC_type(1:2,modes))
  LA_BC_type(BC_EDGE,:)=BC_TYPE_OPEN
  DO imode=1,modes
    LA_BC_type(BC_AXIS,iMode)=LA_BCtype_axis(zero_odd_even(iMode))
    IF(LA_BC_type(BC_AXIS,iMode).EQ.0)THEN ! AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
      LA_BC_type(BC_AXIS,iMode)=-1*MIN(LA_base%s%deg+1,LA_base%f%Xmn(1,iMode))
    END IF
  END DO
  END ASSOCIATE !LA

  CALL exit_subregion("discretization")
  CALL enter_subregion("boundary")
  !INITIALIZATION PARAMETERS (ONLY NECESSARY ON MPIroot)
  IF(MPIroot)THEN
    init_average_axis= GETLOGICAL("init_average_axis",Proposal=.FALSE.)
    IF(init_average_axis)THEN
      average_axis_move(1) = GETREAL("average_axis_move_X1",Proposal=0.0_wp)
      average_axis_move(2) = GETREAL("average_axis_move_X2",Proposal=0.0_wp)
    END IF
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
    IF(getBoundaryFromFile.EQ.-1)THEN
      IF(((init_BC.EQ.1).OR.(init_BC.EQ.2)))THEN !READ edge values from input file
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
    ELSE !getBoundaryFromFile
      CALL BFF%convert_to_modes(X1_base%f,X2_base%f,X1_b,X2_b,scale_minor_radius)
    END IF
  END IF !MPIroot

  boundary_perturb = GETLOGICAL('boundary_perturb', Proposal=.FALSE.)
  boundary_perturb_type_str  = GETSTR("boundary_perturb_type", proposal="legacy")
  IF (boundary_perturb_type_str .EQ. "legacy") THEN
    boundary_perturb_type = BLEND_LEGACY
  ELSE IF (boundary_perturb_type_str .EQ. "cosm") THEN
    boundary_perturb_type = BLEND_COSM
  ELSE
    CALL abort(__STAMP__,&
    'boundary_perturb_type must be "legacy" or "cosm", found '//TRIM(boundary_perturb_type_str),&
    intInfo=boundary_perturb_type,&
    TypeInfo="InvalidParameterError")
  END IF
  boundary_perturb_depth = GETREAL('boundary_perturb_depth', proposal=0.6_wp)
  IF(boundary_perturb)THEN
    ALLOCATE(X1pert_b(1:X1_base%f%modes) )
    ALLOCATE(X2pert_b(1:X2_base%f%modes) )
    X1pert_b=0.0_wp
    X2pert_b=0.0_wp
    !READ boudnary values from input file
    ASSOCIATE(modes=>X1_base%f%modes,sin_range=>X1_base%f%sin_range,cos_range=>X1_base%f%cos_range)
    SWRITE(UNIT_stdOut,'(4X,A)')'... read data for X1pert:'
    DO iMode=sin_range(1)+1,sin_range(2)
      X1pert_b(iMode)=get_iMode('X1pert_b_sin',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X1pert_b(iMode)=get_iMode('X1pert_b_cos',X1_base%f%Xmn(:,iMode),X1_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
    ASSOCIATE(modes=>X2_base%f%modes,sin_range=>X2_base%f%sin_range,cos_range=>X2_base%f%cos_range)
    SWRITE(UNIT_stdOut,'(4X,A)')'... read data for X2pert:'
    DO iMode=sin_range(1)+1,sin_range(2)
      X2pert_b(iMode)=get_iMode('X2pert_b_sin',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      X2pert_b(iMode)=get_iMode('X2pert_b_cos',X2_base%f%Xmn(:,iMode),X2_base%f%nfp)
    END DO !iMode
    END ASSOCIATE
  END IF
  CALL exit_subregion("boundary")
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
  CALL exit_subregion("init-MHD3D")
  CALL par_barrier(afterScreenOut='...DONE')
  SWRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE InitMHD3D

SUBROUTINE InitProfile(sf, var,var_profile)
  ! MODULES
  USE MODgvec_ReadInTools    , ONLY: GETSTR,GETLOGICAL,GETINT,GETINTARRAY,GETREAL,GETREALALLOCARRAY, GETREALARRAY
  USE MODgvec_rProfile_bspl  , ONLY: t_rProfile_bspl
  USE MODgvec_rProfile_poly  , ONLY: t_rProfile_poly
  USE MODgvec_cubic_spline   , ONLY: interpolate_cubic_spline
  USE MODgvec_rProfile_base, ONLY: c_rProfile
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
  CHARACTER(LEN=4), INTENT(IN) :: var
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_rProfile), ALLOCATABLE   :: var_profile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: n_profile_knots
  REAL(wp),ALLOCATABLE :: profile_knots(:)
  INTEGER              :: n_profile_coefs    !! number of polynomial/bspline coeffients for profile profile
  REAL(wp),ALLOCATABLE :: profile_coefs(:)   !! polynomial/bspline coefficients of the profile profile
  CHARACTER(LEN=20)    :: profile_type
  REAL(wp),ALLOCATABLE :: profile_rho2(:)
  REAL(wp),ALLOCATABLE :: profile_vals(:)
  INTEGER              :: n_profile_vals
  INTEGER              :: n_profile_rho2
  REAL(wp)             :: profile_BC_vals(1:2)
  CHARACTER(LEN=10)    :: profile_BC_type(1:2)
  REAL(wp)             :: profile_scale
  CHARACTER(LEN=10)    :: possible_BCs(0:2)
  INTEGER              :: BC(1:2),iBC,jBC
  !===================================================================================================================================
  CALL enter_subregion("get-profile-"//TRIM(var))
  possible_BCs(0)="not_a_knot"
  possible_BCs(1)="1st_deriv"
  possible_BCs(2)="2nd_deriv"

  profile_scale=GETREAL(var//"_scale",Proposal=1.0_wp)
  profile_type  = GETSTR(var//"_type") !make it mandatory
  IF (profile_type.EQ."polynomial") THEN
    CALL GETREALALLOCARRAY(var//"_coefs",profile_coefs,n_profile_coefs) !a+b*s+c*s^2...
    profile_coefs=profile_coefs*profile_scale
    var_profile = t_rProfile_poly(profile_coefs)
  ELSE IF (profile_type.EQ."bspline") THEN
    CALL GETREALALLOCARRAY(var//"_coefs",profile_coefs,n_profile_coefs)
    profile_coefs=profile_coefs*profile_scale
    CALL GETREALALLOCARRAY(var//"_knots",profile_knots,n_profile_knots)
    IF(ABS(profile_knots(1)).GT.1.0d-12) CALL abort(__STAMP__,&
        "First knot position must be =0 for bspline of "//TRIM(var)//" profile!",&
        TypeInfo="InvalidParameterError")
    IF(ABS(profile_knots(n_profile_knots)-1.0_wp).GT.1.0d-12) CALL abort(__STAMP__,&
        "Last knot position must be =1 for bspline of "//TRIM(var)//" profile!",&
        TypeInfo="InvalidParameterError")
    var_profile = t_rProfile_bspl(coefs=profile_coefs,knots=profile_knots)
  ELSE IF (profile_type.EQ."interpolation") THEN
    CALL GETREALALLOCARRAY(var//"_vals",profile_vals, n_profile_vals)
    CALL GETREALALLOCARRAY(var//"_rho2",profile_rho2, n_profile_rho2)
    IF(ABS(profile_rho2(1)).GT.1.0d-12) CALL abort(__STAMP__,&
      "First rho2 position must be =0 for interpolation of "//TRIM(var)//" profile!",&
        TypeInfo="InvalidParameterError")
    IF(ABS(profile_rho2(n_profile_rho2)-1.0_wp).GT.1.0d-12) CALL abort(__STAMP__,&
      "Last rho2 position must be =1 for interpolation of "//TRIM(var)//" profile!",&
        TypeInfo="InvalidParameterError")
    IF (n_profile_vals .NE. n_profile_rho2) THEN
      CALL abort(__STAMP__,&
      'Size of '//var//'_rho2 and '//var//'_vals must be equal!',&
        TypeInfo="InvalidParameterError")
    END IF
    profile_BC_type(1) = GETSTR(var//"_BC_type_axis",Proposal="not_a_knot")
    profile_BC_type(2) = GETSTR(var//"_BC_type_edge",Proposal="not_a_knot")
    BC=-1
    DO iBC=1,2
      DO jBC=0,2
        IF (INDEX(TRIM(profile_BC_type(iBC)),TRIM(possible_BCs(jBC)))>0) THEN
          BC(iBC)=jBC
          EXIT
        END IF
      END DO !jBC
    END DO !iBC
    IF(ANY(BC<0)) THEN
      CALL abort(__STAMP__,&
          "BC_type of profile must be 'not_a_knot', '1st_deriv' or '2nd_deriv' ... got '"//TRIM(profile_BC_type(1))//"' on axis and '"//TRIM(profile_BC_type(2))//"' on edge",&
        TypeInfo="InvalidParameterError")
    END IF

    IF(ANY(BC>0)) THEN
     profile_BC_vals = GETREALARRAY(var//"_BC_vals", 2, Proposal=(/0.0_wp, 0.0_wp/))
     CALL interpolate_cubic_spline(profile_rho2,profile_vals, profile_coefs, profile_knots, BC, profile_BC_vals)
    ELSE
     CALL interpolate_cubic_spline(profile_rho2,profile_vals, profile_coefs, profile_knots, BC)
    END IF

    profile_coefs=profile_coefs*profile_scale
    var_profile = t_rProfile_bspl(profile_knots,profile_coefs)
    SDEALLOCATE(profile_vals)
    SDEALLOCATE(profile_rho2)
  ELSE
    CALL abort(__STAMP__,&
         "Specified "//var//"_type unknown. Expecting 'polynomial', 'bspline' or 'interpolation' ... got '"//TRIM(profile_type)//"'",&
        TypeInfo="InvalidParameterError")
  END IF ! profile type

  SDEALLOCATE(profile_knots)
  SDEALLOCATE(profile_coefs)
  CALL exit_subregion("get-profile-"//TRIM(var))
END SUBROUTINE InitProfile

!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE InitSolutionMHD3D(sf)
! MODULES
  USE MODgvec_MHD3D_Vars     , ONLY: which_init,U,F,init_LA,boundary_perturb,boundary_perturb_depth,boundary_perturb_type
  USE MODgvec_Restart_vars   , ONLY: doRestart,RestartFile
  USE MODgvec_Restart        , ONLY: RestartFromState
  USE MODgvec_Restart        , ONLY: WriteState
  USE MODgvec_MHD3D_EvalFunc , ONLY: InitProfilesGP,EvalEnergy,EvalForce
  USE MODgvec_Analyze        , ONLY: Analyze
  USE MODgvec_ReadInTools    , ONLY: GETLOGICAL
  USE MODgvec_MPI            , ONLY: par_Bcast,par_barrier
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_functional_mhd3d), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: JacCheck
!===================================================================================================================================
  CALL par_barrier(beforeScreenOut="    INITIALIZE SOLUTION...",afterScreenOut="                           ...")
  CALL enter_subregion("init-solution")
  IF(MPIroot) THEN
    IF(doRestart)THEN
      WRITE(UNIT_stdOut,'(4X,A)')'... restarting from file ... '
      CALL RestartFromState(RestartFile,U(0))
      CALL InitSolution(U(0),-1) ! (re-)apply BC and init LA (if init_LA is true)
    ELSE
      CALL InitSolution(U(0),which_init)
    END IF
    IF(boundary_perturb)THEN
      CALL AddBoundaryPerturbation(U(0),boundary_perturb_depth,boundary_perturb_type)
    END IF !boundary_perturb
  END IF !MPIroot
  CALL par_Bcast(U(0)%X1,0)
  CALL par_Bcast(U(0)%X2,0)
  CALL exit_subregion("init-solution")

  IF(init_LA) THEN
    CALL Init_LA_From_Solution(U(0))  !BCast inside
  ELSE
    CALL par_Bcast(U(0)%LA,0)
  END IF

  CALL U(-1)%set_to(U(0))

  JacCheck=2
  CALL InitProfilesGP() !evaluate profiles once at Gauss Points (on MPIroot + BCast)

  CALL enter_subregion("check-solution")
  U(0)%W_MHD3D=EvalEnergy(U(0),.TRUE.,JacCheck)
  IF(JacCheck.EQ.-1)THEN
    CALL Analyze(0)
    CALL abort(__STAMP__,&
        "NEGATIVE JACOBIAN FOUND AFTER INITIALIZATION!",TypeInfo="InitializationError")
  END IF
  CALL WriteState(U(0),0)
  CALL EvalForce(U(0),.FALSE.,JacCheck, F(0))
  SWRITE(UNIT_stdOut,'(8x,A,3E11.4)')'|Force|= ',SQRT(F(0)%norm_2())
  CALL Analyze(0)
  CALL exit_subregion("check-solution")
  CALL par_barrier(afterScreenOut="    ...DONE")
  SWRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE InitSolutionMHD3D


!===================================================================================================================================
!> automatically build the string to be read from parameterfile, varname + m,n mode number, and then read it from parameterfile
!!
!===================================================================================================================================
FUNCTION get_iMode(varname_in,mn_in,nfp_in)
! MODULES
  USE MODgvec_ReadInTools    , ONLY: GETREAL,remove_blanks
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
  SWRITE(varstr,'(A,"("I4,";",I4,")")')TRIM(varname_in),mn_in(1),mn_in(2)/nfp_in
  varstr=remove_blanks(varstr)         !quiet on default=0.0
  get_iMode=GETREAL(TRIM(varstr),Proposal=0.0_wp,quiet_def_in=.TRUE.)
END FUNCTION get_iMode


!===================================================================================================================================
!> Overwrite axis with average axis by center of closed line of the boundary in each poloidal plane
!!
!===================================================================================================================================
SUBROUTINE InitAverageAxis()
! MODULES
  USE MODgvec_MHD3D_Vars   , ONLY:X1_base,X1_a,X1_b
  USE MODgvec_MHD3D_Vars   , ONLY:X2_base,X2_a,X2_b
  USE MODgvec_MHD3D_Vars   , ONLY:average_axis_move
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT Variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL Variables
  REAL(wp) :: X1_b_IP(X1_base%f%mn_nyq(1),X1_base%f%mn_nyq(2))
  REAL(wp) :: X2_b_IP(X2_base%f%mn_nyq(1),X2_base%f%mn_nyq(2))

  INTEGER  :: i_m,i_n
  REAL(wp) :: dl,lint,x1int,x2int
!-----------------------------------------------------------------------------------------------------------------------------------
  ASSOCIATE(m_nyq=>X1_base%f%mn_nyq(1),n_nyq=>X1_base%f%mn_nyq(2))
    X1_b_IP(:,:) = RESHAPE(X1_base%f%evalDOF_IP(0,X1_b),(/m_nyq,n_nyq/))
    X2_b_IP(:,:) = RESHAPE(X2_base%f%evalDOF_IP(0,X2_b),(/m_nyq,n_nyq/))
    DO i_n=1,n_nyq
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
      ! c_x= 1/(6A) sum_i (x_i-1+x_i)*(x_i-1*y_i - x_i*y_i-1), A=1/2 sum_i  (x_i-1*y_i - x_i*y_i-1)
      X1_b_IP(:,i_n) = x1int/(3.0_wp*lint) + average_axis_move(1)
      X2_b_IP(:,i_n) = x2int/(3.0_wp*lint) + average_axis_move(2)
    END DO
    X1_a = X1_base%f%initDOF(RESHAPE(X1_b_IP,(/X1_base%f%mn_IP/)))
    X2_a = X2_base%f%initDOF(RESHAPE(X2_b_IP,(/X2_base%f%mn_IP/)))
  END ASSOCIATE
END SUBROUTINE InitAverageAxis

!===================================================================================================================================
!> Initialize the solution with the given boundary condition
!!
!===================================================================================================================================
SUBROUTINE InitSolution(U_init,which_init_in)
! MODULES
  USE MODgvec_Globals,       ONLY:ProgressBar,getTime
  USE MODgvec_MHD3D_Vars   , ONLY:init_fromBConly,init_BC,init_average_axis,average_axis_move
  USE MODgvec_MHD3D_Vars   , ONLY:X1_base,X1_BC_Type,X1_a,X1_b
  USE MODgvec_MHD3D_Vars   , ONLY:X2_base,X2_BC_Type,X2_a,X2_b
  USE MODgvec_MHD3D_Vars   , ONLY:LA_base,LA_BC_Type,init_LA
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
  INTEGER  :: iMode,i_m
  REAL(wp) :: BC_val(2)
  REAL(wp) :: rhopos
  REAL(wp) :: X1_gIP(1:X1_base%s%nBase,1:X1_base%f%modes)
  REAL(wp) :: X2_gIP(1:X2_base%s%nBase,1:X2_base%f%modes)
  REAL(wp) :: LA_gIP(1:LA_base%s%nBase,1:LA_base%f%modes)
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                       "InitSolution should only be called by MPIroot!")
  X1_gIP=0.0_wp; X2_gIP=0.0_wp; LA_gIP=0.0_wp

  SELECT CASE(which_init_in)
  CASE(-1) !restart
    X1_a(:)=U_init%X1(1,:)
    X2_a(:)=U_init%X2(1,:)
    X1_b(:)=U_init%X1(X1_base%s%nBase,:)
    X2_b(:)=U_init%X2(X2_base%s%nBase,:)
    IF (init_average_axis)THEN
      WRITE(UNIT_stdOut,'(A)') "WARNING: init_average_axis ignored due to restart!"
    END IF
  CASE(0)
    !X1_a,X2_a and X1_b,X2_b already filled from parameter file readin...
    IF(init_average_axis)THEN
      CALL InitAverageAxis()
    END IF !init_average_axis
  CASE(1) !VMEC
    IF((init_BC.EQ.-1).OR.(init_BC.EQ.1))THEN ! compute  axis from VMEC, else use the one defined in paramterfile
      rhopos=0.0_wp
      ASSOCIATE(sin_range    => X1_base%f%sin_range, cos_range    => X1_base%f%cos_range )
      DO imode=cos_range(1)+1,cos_range(2)
        X1_a(iMode:iMode)= VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,(/rhopos/),Rmnc_Spl)
      END DO
      IF(lasym)THEN
        DO imode=sin_range(1)+1,sin_range(2)
          X1_a(iMode:iMode)=X1_a(iMode:imode)  +VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,(/rhopos/),Rmns_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X1
      ASSOCIATE(sin_range    => X2_base%f%sin_range, cos_range    => X2_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        X2_a(iMode:iMode)=VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,(/rhopos/),Zmns_Spl)
      END DO
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          X2_a(iMode:iMode)=X2_a(iMode:iMode)  +VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,(/rhopos/),Zmnc_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X2
    END IF
    IF((init_BC.EQ.-1).OR.(init_BC.EQ.0))THEN ! compute edge from VMEC, else use the one defined in paramterfile
      rhopos=1.0_wp
      ASSOCIATE(sin_range    => X1_base%f%sin_range, cos_range    => X1_base%f%cos_range )
      DO imode=cos_range(1)+1,cos_range(2)
        X1_b(iMode:iMode)=VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,(/rhopos/),Rmnc_Spl)
      END DO
      IF(lasym)THEN
        DO imode=sin_range(1)+1,sin_range(2)
          X1_b(iMode:iMode)=X1_b(iMode:iMode)  +VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,(/rhopos/),Rmns_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X1
      ASSOCIATE(sin_range    => X2_base%f%sin_range, cos_range    => X2_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        X2_b(iMode:iMode)=VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,(/rhopos/),Zmns_Spl)
      END DO
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          X2_b(iMode:iMode)=X2_b(iMode:iMode)  +VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,(/rhopos/),Zmnc_Spl)
        END DO
      END IF !lasym
      END ASSOCIATE !X2
    END IF
    IF(init_average_axis)THEN
      CALL InitAverageAxis()
    END IF !init_average_axis
    IF(.NOT.init_fromBConly)THEN !only boundary and axis from VMEC
      ASSOCIATE(s_IP         => X1_base%s%s_IP, &
                nBase        => X1_base%s%nBase, &
                sin_range    => X1_base%f%sin_range,&
                cos_range    => X1_base%f%cos_range )
      DO imode=cos_range(1)+1,cos_range(2)
        X1_gIP(:,iMode)  =VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,s_IP,Rmnc_Spl)
      END DO !imode=cos_range
      IF(lasym)THEN
        DO imode=sin_range(1)+1,sin_range(2)
          X1_gIP(:,iMode)=X1_gIP(:,iMode)  +VMEC_EvalSplMode(X1_base%f%Xmn(:,iMode),0,s_IP,Rmns_Spl)
        END DO !imode= sin_range
      END IF !lasym
      END ASSOCIATE !X1
      ASSOCIATE(s_IP         => X2_base%s%s_IP, &
                nBase        => X2_base%s%nBase, &
                sin_range    => X2_base%f%sin_range,&
                cos_range    => X2_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        X2_gIP(:,iMode)=VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,s_IP,Zmns_Spl)
      END DO !imode=sin_range
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          X2_gIP(:,iMode)=X2_gIP(:,iMode)  +VMEC_EvalSplMode(X2_base%f%Xmn(:,iMode),0,s_IP,Zmnc_Spl)
        END DO !imode= sin_range
      END IF !lasym
      END ASSOCIATE !X2
      ASSOCIATE(s_IP         => LA_base%s%s_IP, &
                nBase        => LA_base%s%nBase, &
                sin_range    => LA_base%f%sin_range,&
                cos_range    => LA_base%f%cos_range )
      DO imode=sin_range(1)+1,sin_range(2)
        LA_gIP(:,iMode)=VMEC_EvalSplMode(LA_base%f%Xmn(:,iMode),0,s_IP,lmns_Spl)
      END DO !imode= sin_range
      IF(lasym)THEN
        DO imode=cos_range(1)+1,cos_range(2)
          LA_gIP(:,iMode)=LA_gIP(:,iMode)  +VMEC_EvalSplMode(LA_base%f%Xmn(:,iMode),0,s_IP,lmnc_Spl)
        END DO !imode=cos_range
      END IF !lasym
      END ASSOCIATE !LA
     END IF !fullIntVmec
  END SELECT !which_init


  IF((which_init_in.NE.-1).AND.(init_fromBConly))THEN
    !no restart(=-1) and initialization only
    !smoothly interpolate between  edge and axis data
    ASSOCIATE(s_IP         =>X1_base%s%s_IP, &
              modes        =>X1_base%f%modes, &
              zero_odd_even=>X1_base%f%zero_odd_even)
    DO imode=1,modes
      SELECT CASE(zero_odd_even(iMode))
      CASE(MN_ZERO,M_ZERO) !X1_a only used here!!
        X1_gIP(:,iMode)=(1.0_wp-(s_IP(:)**2))*X1_a(iMode)+(s_IP(:)**2)*X1_b(iMode)  ! meet edge and axis, ~(1-s^2)
      CASE(M_ODD_FIRST)
        X1_gIP(:,iMode)=s_IP(:)*X1_b(iMode)      ! first odd mode ~s
      CASE DEFAULT ! ~s^m
        X1_gIP(:,iMode)=s_IP(:)*X1_b(iMode)
        DO i_m=2,X1_base%f%Xmn(1,iMode)
          X1_gIP(:,iMode)=X1_gIP(:,iMode)*s_IP(:)
        END DO
      END SELECT !X1(:,iMode) zero odd even
    END DO !iMode
    END ASSOCIATE

    ASSOCIATE(s_IP         =>X2_base%s%s_IP, &
              modes        =>X2_base%f%modes, &
              zero_odd_even=>X2_base%f%zero_odd_even)
    DO imode=1,modes
      SELECT CASE(zero_odd_even(iMode))
      CASE(MN_ZERO,M_ZERO) !X2_a only used here!!!
        X2_gIP(:,iMode)=(1.0_wp-(s_IP(:)**2))*X2_a(iMode)+(s_IP(:)**2)*X2_b(iMode) ! meet edge and axis, ~(1-s^2)
      CASE(M_ODD_FIRST)
        X2_gIP(:,iMode)=s_IP(:)*X2_b(iMode)      ! first odd mode ~s
      CASE DEFAULT ! ~s^m
        X2_gIP(:,iMode)=s_IP(:)*X2_b(iMode)
        DO i_m=2,X2_base%f%Xmn(1,iMode)
          X2_gIP(:,iMode)=X2_gIP(:,iMode)*s_IP(:)
        END DO
      END SELECT !X2(:,iMode) zero odd even
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
    IF(which_init_in.NE.-1) U_init%X1(:,iMode)=X1_base%s%initDOF( X1_gIP(:,iMode) )
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
    IF(which_init_in.NE.-1) U_init%X2(:,iMode)=X2_base%s%initDOF( X2_gIP(:,iMode) )
    CALL X2_base%s%applyBCtoDOF(U_init%X2(:,iMode),X2_BC_type(:,iMode),BC_val)
  END DO
  END ASSOCIATE !X2

  ASSOCIATE(modes        =>LA_base%f%modes, &
            zero_odd_even=>LA_base%f%zero_odd_even)
  IF((which_init_in.NE.-1).AND.(.NOT.init_LA)) THEN
    !lambda init might not be needed since it has no boundary condition and changes anyway after the update of the mapping...
    IF(.NOT.init_fromBConly)THEN
      WRITE(UNIT_stdOut,'(4X,A)') "... lambda initialized with VMEC ..."
    ELSE
      WRITE(UNIT_stdOut,'(4X,A)') "... initialize lambda =0 ..."
      LA_gIP=0.0_wp
    END IF
  END IF !!which_init_in>-1 and not init_LA
  !always apply strong BC
  DO imode=1,modes
    IF(zero_odd_even(iMode).EQ.MN_ZERO)THEN
      U_init%LA(:,iMode)=0.0_wp ! (0,0) mode should not be here, but must be zero if its used.
    ELSE
      BC_val =(/ 0.0_wp, 0.0_wp/)
      IF(which_init_in.NE.-1) U_init%LA(:,iMode)=LA_base%s%initDOF( LA_gIP(:,iMode) )
      CALL LA_base%s%applyBCtoDOF(U_init%LA(:,iMode),LA_BC_type(:,iMode),BC_val)
    END IF!iMode ~ MN_ZERO
  END DO !iMode

  END ASSOCIATE !LA
END SUBROUTINE InitSolution


!===================================================================================================================================
!> Initialize LAMBDA FROM U_init%X1,%X2 and iota profile, this computation is distributed over MPIranks
!!
!===================================================================================================================================
SUBROUTINE Init_LA_from_Solution(U_init)
! MODULES
  USE MODgvec_Globals,       ONLY:ProgressBar,getTime,myRank,nRanks
  USE MODgvec_MHD3D_Vars   , ONLY:X1_base,X2_base,LA_base,LA_BC_Type,hmap, Phi_profile, chi_profile
  USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_mhd3d
  USE MODgvec_lambda_solve,  ONLY:lambda_solve
  USE MODgvec_MPI           ,ONLY:par_reduce,par_BCast
  USE MODgvec_hmap          ,ONLY:hmap_new_auxvar,PP_T_HMAP_AUXVAR
!$ USE omp_lib
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: U_init
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode,is,ns_str,ns_end,iRank,nBase
  REAL(wp) :: BC_val(2),rhopos
  REAL(wp) :: StartTime,EndTime
  REAL(wp),DIMENSION(1:LA_base%s%nBase):: PhiPrime,chiPrime
  REAL(wp) :: LA_gIP(1:LA_base%s%nBase,1:LA_base%f%modes)
#ifdef PP_WHICH_HMAP
  TYPE(PP_T_HMAP_AUXVAR),ALLOCATABLE  :: hmap_xv(:) !! auxiliary variables for hmap
#else
  CLASS(PP_T_HMAP_AUXVAR),ALLOCATABLE  :: hmap_xv(:) !! auxiliary variables for hmap
#endif
!===================================================================================================================================
  StartTime=GetTime()
  SWRITE(UNIT_stdOut,'(4X,A)') "... Initialize lambda from mapping ..."
  CALL enter_subregion("reinit-lambda")
  nBase        = LA_base%s%nBase
  ASSOCIATE(modes        => LA_base%f%modes, &
            s_IP         => LA_base%s%s_IP, &
            zero_odd_even=> LA_base%f%zero_odd_even, &
            modes_str    => LA_base%f%modes_str, &
            modes_end    => LA_base%f%modes_end, &
            offset_modes => LA_Base%f%offset_modes )
  !evaluate profiles only in MPIroot!
  IF(MPIroot)THEN
    DO is=1,nBase
      rhopos=MIN(1.0_wp-1.0e-12_wp,MAX(1.0e-4_wp,s_IP(is))) !exclude axis
      phiPrime(is)=Phi_profile%eval_at_rho(rhopos,deriv=1)
      chiPrime(is)=chi_profile%eval_at_rho(rhopos,deriv=1)
    END DO
  END IF !MPIroot
  CALL par_BCast(phiPrime,0)
  CALL par_BCast(chiPrime,0)
  !initialize Lambda, radially parallel
  ns_str = (nBase*(myRank  ))/nRanks+1
  ns_end = (nBase*(myRank+1))/nRanks
  LA_gIP=0.0_wp
  CALL ProgressBar(0,ns_end) !init

  CALL hmap_new_auxvar(hmap,X1_base%f%x_IP(2,:),hmap_xv,.FALSE.) !no 2nd derivative needed
  DO is=ns_str,ns_end
    rhopos=MIN(1.0_wp-1.0e-12_wp,MAX(1.0e-4_wp,s_IP(is))) !exclude axis
    CALL lambda_Solve(rhopos,hmap,hmap_xv,X1_base,X2_base,LA_base%f,U_init%X1,U_init%X2,LA_gIP(is,:),phiPrime(is),chiPrime(is))
    CALL ProgressBar(is,ns_end)
  END DO !is
  DEALLOCATE(hmap_xv)
!!!  CALL par_reduce(LA_gIP,'SUM',0)
!!!  IF(MPIroot)THEN
!!!    DO iMode=1,modes
!!!      IF(zero_odd_even(iMode).EQ.MN_ZERO)THEN
!!!        U_init%LA(:,iMode)=0.0_wp ! (0,0) mode should not be here, but must be zero if its used.
!!!      ELSE
!!!        U_init%LA(:,iMode)=LA_base%s%initDOF( LA_gIP(:,iMode) )
!!!      END IF!iMode ~ MN_ZERO
!!!      BC_val =(/ 0.0_wp, 0.0_wp/)
!!!      CALL LA_base%s%applyBCtoDOF(U_init%LA(:,iMode),LA_BC_type(:,iMode),BC_val)
!!!    END DO !iMode=1,modes
!!!  END IF
!!!  CALL par_BCast(U_init%LA,0)
  !reduce radially, different mode sets to different MPIranks (should be a gatherv)
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      CALL par_Reduce(LA_gIP(1:nbase,offset_modes(iRank)+1:offset_modes(iRank+1)),'SUM',iRank)
  END DO
  DO iMode=modes_str,modes_end
    IF(zero_odd_even(iMode).EQ.MN_ZERO)THEN
      U_init%LA(1:nBase,iMode)=0.0_wp ! (0,0) mode should not be here, but must be zero if its used.
    ELSE
      U_init%LA(1:nBase,iMode)=LA_base%s%initDOF( LA_gIP(1:nBase,iMode) )
    END IF!iMode ~ MN_ZERO
    BC_val =(/ 0.0_wp, 0.0_wp/)
    CALL LA_base%s%applyBCtoDOF(U_init%LA(:,iMode),LA_BC_type(:,iMode),BC_val)
  END DO !iMode=modes_str, modes_end
  ! broadcast result: different mode ranges to different MPIranks
  DO iRank=0,nRanks-1
    IF(offset_modes(iRank+1)-offset_modes(iRank).GT.0) &
      CALL par_Bcast(U_init%LA(1:nBase,offset_modes(iRank)+1:offset_modes(iRank+1)),iRank)
  END DO
  END ASSOCIATE !LA
  EndTime=GetTime()
  CALL exit_subregion("reinit-lambda")
  SWRITE(UNIT_stdOut,'(4X,A,F9.2,A)') " init lambda took [ ",EndTime-StartTime," sec]"
END SUBROUTINE Init_LA_from_solution

!===================================================================================================================================
!> Add boundary perturbation
!!
!===================================================================================================================================
SUBROUTINE AddBoundaryPerturbation(U_init, depth, blend_type)
! MODULES
  USE MODgvec_MHD3D_Vars   , ONLY:X1_base,X1_BC_Type,X1_a,X1_b,X1pert_b
  USE MODgvec_MHD3D_Vars   , ONLY:X2_base,X2_BC_Type,X2_a,X2_b,X2pert_b
  USE MODgvec_sol_var_MHD3D, ONLY:t_sol_var_mhd3d
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp),INTENT(IN) :: depth ! depth of perturbation from boundary (0.1..0.3)
  INTEGER, INTENT(IN) :: blend_type ! 0/BLEND_LEGACY: legacy Gaussian, 1/BLEND_COSM: new cosine
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: U_init
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode
  REAL(wp) :: BC_val(2)
  REAL(wp) :: X1pert_gIP(1:X1_base%s%nBase)
  REAL(wp) :: X2pert_gIP(1:X2_base%s%nBase)
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                       "AddBoundaryPerturbation should only be called by MPIroot!")
  WRITE(UNIT_stdOut,'(4X,A)') "ADD BOUNDARY PERTURBATION..."

  ASSOCIATE(s_IP         =>X1_base%s%s_IP, &
            modes        =>X1_base%f%modes )
  DO imode=1,modes
    X1_b(iMode)=X1_b(iMode)+X1pert_b(iMode)
    X1pert_gIP(:)=blend(s_IP, depth, X1_base%f%Xmn(1, iMode))*X1pert_b(iMode)
    U_init%X1(:,iMode)=U_init%X1(:,iMode) + X1_base%s%initDOF( X1pert_gIP(:) )
  END DO
  END ASSOCIATE

  ASSOCIATE(s_IP         =>X2_base%s%s_IP, &
            modes        =>X2_base%f%modes )
  DO imode=1,modes
    X2_b(iMode)=X2_b(iMode)+X2pert_b(iMode)
    X2pert_gIP(:)=blend(s_IP, depth, X2_base%f%Xmn(1, iMode))*X2pert_b(iMode)
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

  WRITE(UNIT_stdOut,'(4X,A)') "... DONE."
  WRITE(UNIT_stdOut,fmt_sep)

  CONTAINS

  ELEMENTAL FUNCTION blend(s_in, depth, m)
    USE MODgvec_Globals, ONLY: wp, PI
    USE MODgvec_MHD3D_Vars, ONLY: BLEND_LEGACY
    REAL(wp),INTENT(IN) :: s_in !input coordinate [0,1]
    REAL(wp)            :: blend
    REAL(wp),INTENT(IN) :: depth
    INTEGER,INTENT(IN)  :: m     ! exponent for cosine blending (poloidal mode number)
    ASSOCIATE(shift => 1.0_wp - depth)
      IF (blend_type == BLEND_LEGACY) THEN
        blend = EXP(-4.0_wp * ((s_in - 1.0_wp) / depth)**2)
      ELSE IF (s_in .GE. shift) THEN
        blend = (COS(((s_in - shift) / (1.0_wp - shift) - 1.0_wp) * PI) / 2.0_wp + 0.5_wp)**m
      ELSE IF (m .EQ. 0) THEN
        blend = 1.0_wp
      ELSE
        blend = 0.0_wp
      END IF
    END ASSOCIATE
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
  CALL enter_subregion("minimize")
  SELECT CASE(MinimizerType)
  CASE(0,10)
    CALL MinimizeMHD3D_descent(sf)
  CASE DEFAULT
    CALL abort(__STAMP__,&
        "requested MinimizeType does not exist, expecting 0 or 10",intinfo=MinimizerType,&
        TypeInfo="InvalidParameterError")
  END SELECT
  CALL exit_subregion("minimize")
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
  USE MODgvec_MHD3D_visu, ONLY:WriteSFLoutfile
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
  REAL(wp)  :: min_dt_out,max_dt_out,min_dw_out,max_dw_out,sum_dW_out,t_pseudo,Fnorm(3),Vnorm(3),Fnorm0(3),Fnorm_old(3),W_MHD3D_0
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
  Vnorm=0.0_wp
  logiter_ramp=1
  logscreen=1

  first_iter=.TRUE.
  restart_iter=.FALSE.

  CALL U(-3)%set_to(U(0)) !initial state, should remain unchanged

  DO WHILE(iter.LT.maxIter)
    IF((first_iter).OR.(restart_iter))THEN
      JacCheck=1 !abort if detJ<0
      CALL EvalAux(           U(0),JacCheck)
      U(0)%W_MHD3D=EvalEnergy(U(0),.FALSE.,JacCheck)
      W_MHD3D_0 = U(0)%W_MHD3D
      CALL EvalForce(         U(0),.FALSE.,JacCheck,F(0))
      Fnorm0=SQRT(F(0)%norm_2())
      Fnorm=Fnorm0
      Fnorm_old=1.1_wp*Fnorm0
      CALL U(-1)%set_to(U(0)) !last state
      CALL U(-2)%set_to(U(0)) !state at last logging interval
      !for hirshman method
      IF(MinimizerType.EQ.10)THEN
        CALL V(-1)%set_to(0.0_wp)
        CALL V( 0)%set_to(0.0_wp)
        tau(1:ndamp)=0.15_wp/dt
        tau_bar = 0.075_wp
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
      tau_bar = 0.5_wp*dt*SUM(tau)/REAL(ndamp,wp)   !=1/2 * tauavg
      CALL V(1)%AXBY(((1.0_wp-tau_bar)/(1.0_wp+tau_bar)),V(0),(dt/(1.0_wp+tau_bar)),F(0)) !velocity V(1)
      CALL P(1)%AXBY(1.0_wp,U(0),dt,V(1)) !overwrites P(1), predicst solution U(1)
      Vnorm=SQRT(V(1)%norm_2())
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

        IF(ALL(Fnorm.LE.abstol))THEN
          CALL Logging(.FALSE.)
          SWRITE(UNIT_stdOut,'(4x,A)')'==>Iteration finished, |force| in relative tolerance'
          EXIT !DO LOOP
        END IF
        iter=iter+1
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

      ELSE !not a valid step, decrease timestep and skip P(1)
        dt=0.9_wp*dt
        nstepDecreased=nStepDecreased+1
        nSkip_dW=nSkip_dW+1
        !CALL U(0)%set_to(U(-2))
        restart_iter=.TRUE.
        SWRITE(UNIT_stdOut,'(8X,I8,A,E8.1,A,E11.4)')iter,'...deltaW>',dW_allowed,'*W_MHD3D_0, skip step and decrease stepsize to dt=',dt
      END IF
    END IF !JacCheck

    IF(nStepDecreased.GT.130) THEN ! 0.9^130 ~10^-6
      SWRITE(UNIT_stdOut,'(A,E21.11)')'Iteration stopped since timestep has been decreased by 0.9^130: ', dt
      SWRITE(UNIT_stdOut,fmt_sep)
      RETURN
    END IF
    IF((MOD(iter,outputIter).EQ.0).AND.(lastoutputIter.NE.iter))THEN
      __PERFON('output')
      SWRITE(UNIT_stdOut,'(A)')'##########################  OUTPUT ##################################'
      CALL Analyze(iter)
      CALL WriteState(U(0),iter)
      SWRITE(UNIT_stdOut,'(A)')'#####################################################################'
      lastOutputIter=iter
      __PERFOFF('output')
    END IF
  END DO !iter
  IF(iter.GE.MaxIter)THEN
    SWRITE(UNIT_stdOut,'(A,E21.11)')"maximum iteration count exceeded"
  END IF
  SWRITE(UNIT_stdOut,'(A)') "... DONE."
  SWRITE(UNIT_stdOut,fmt_sep)
  IF(lastoutputIter.NE.iter)THEN
    CALL Analyze(MIN(iter,MaxIter))
    CALL WriteState(U(0),MIN(iter,MaxIter))
  END IF
  CALL FinishLogging()
  CALL writeSFLoutfile(U(0),MIN(iter,MaxIter))


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
  REAL(wp)            :: AxisPos(2,2),W_MHD3D
  INTEGER,PARAMETER   :: nLogDat=20
  REAL(wp)            :: LogDat(1:nLogDat)
  !=================================================================================================================================
  IF(.NOT.MPIroot) RETURN
  __PERFON('log_output')
  W_MHD3D=U(0)%W_MHD3D
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  WRITE(UNIT_stdOut,'(A,E11.4,A)')'%%%%%%%%%%  START ITERATION, dt= ',dt, '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                 '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
  WRITE(UNIT_stdOut,'(A,3E21.14)') &
          '%%% dU = |Force|= ',Fnorm(1:3)
  IF(MinimizerType.EQ.10) THEN
    WRITE(UNIT_stdOut,'(A,E11.4,A,3E11.4)') &
          '%%% accel.GD: tau= ',tau_bar,' |vel|= ',Vnorm(1:3)
  END IF

  WRITE(UNIT_stdOut,'(40(" -"))')
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
  LogDat(ilogDat+1:iLogDat+11)=(/0.0_wp,0.0_wp,dt,dt,W_MHD3D,0.0_wp,0.0_wp,0.0_wp,Fnorm(1:3)/)
  iLogDat=11
  IF(MinimizerType.EQ.10) THEN
    WRITE(logUnit,'(A)',ADVANCE="NO")',"tau","normV_X1","normV_X2","normV_LA"'
    LogDat(ilogDat+1:iLogDat+4)=(/tau_bar,Vnorm(1:3)/)
    iLogDat=iLogDat+4
  END IF
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
  REAL(wp)            :: AxisPos(2,2),maxDist,avgDist,W_MHD3D
  INTEGER,PARAMETER   :: nLogDat=20
  REAL(wp)            :: LogDat(1:nLogDat)
  !=================================================================================================================================
  IF(.NOT.MPIroot) RETURN
  __PERFON('log_output')
  CALL DATE_AND_TIME(values=TimeArray) ! get System time
  W_MHD3D=U(0)%W_MHD3D
  IF(.NOT.quiet)THEN
    WRITE(UNIT_stdOut,'(80("%"))')
    WRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                      '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
    WRITE(UNIT_stdOut,'(A,I8,A,2I8,A,E11.4,A,2E11.4,A,E21.14,A,3E12.4)') &
                      '%%% #ITERATIONS= ',iter,', #skippedIter (Jac/dW)= ',nSkip_Jac,nSkip_dW, &
              '\n%%% t_pseudo= ',t_pseudo,', min/max dt= ',min_dt_out,max_dt_out, &
              '\n%%% W_MHD3D= ',W_MHD3D,', min/max/sum deltaW= ' , min_dW_out,max_dW_out,sum_dW_out
    WRITE(UNIT_stdOut,'(A,3E21.14)') &
                '%%% dU = |Force|= ',Fnorm(1:3)
    !------------------------------------
  END IF!.NOT.quiet
  iLogDat=0
  runtime_ms=MAX(0,SUM((timeArray(5:8)-StartTimearray(5:8))*(/360000,6000,100,1/)))
  LogDat(ilogDat+1:iLogDat+11)=(/REAL(iter,wp),REAL(runtime_ms,wp)/100.0_wp, &
                                min_dt_out,max_dt_out,W_MHD3D,min_dW_out,max_dW_out,sum_dW_out, &
                                Fnorm(1:3)/)
  iLogDat=11
  IF(MinimizerType.EQ.10) THEN
    IF(.NOT.quiet)THEN
      WRITE(UNIT_stdOut,'(A,E11.4,A,3E11.4)') &
            '%%% accel.GD: tau= ',tau_bar,' |vel|= ',Vnorm(1:3)
    END IF!.NOT.quiet
    LogDat(ilogDat+1:iLogDat+4)=(/tau_bar,Vnorm(1:3)/)
    iLogDat=iLogDat+4
  END IF
  IF(doCheckDistance) THEN
    CALL CheckDistance(U(0),U(-2),maxDist,avgDist)
    CALL U(-2)%set_to(U(0))
    IF(.NOT.quiet)THEN
      WRITE(UNIT_stdOut,'(A,2E11.4)') &
      '               %%% Dist to last log (max/avg) : ',maxDist,avgDist
    END IF!.NOT.quiet
    LogDat(iLogDat+1:iLogDat+2)=(/maxDist,avgDist/)
    iLogDat=iLogDat+2
  END IF!doCheckDistance
  IF(doCheckAxis) THEN
    CALL CheckAxis(U(0),2,AxisPos)
    IF(.NOT.quiet)THEN
      WRITE(UNIT_stdOut,'(2(A,2E22.14))') &
        '%%% axis position (X1,X2,zeta=0     ): ',AxisPos(1:2,1), &
      '\n%%% axis position (X1,X2,zeta=pi/nfp): ',AxisPos(1:2,2)
    END IF!.NOT.quiet
    LogDat(iLogDat+1:iLogDat+4)=RESHAPE(AxisPos,(/4/))
    iLogDat=iLogDat+4
  END IF !doCheckAxis
  IF(.NOT.quiet)THEN
    WRITE(UNIT_stdOut,'(40(" -"))')
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
  CALL enter_subregion("finalize-MHD3D")
  IF(ALLOCATED(X1_base)) CALL X1_base%free()
  IF(ALLOCATED(X2_base)) CALL X2_base%free()
  IF(ALLOCATED(LA_base)) CALL LA_base%free()

  DO i=-1,1
    IF(ALLOCATED(U)) CALL U(i)%free()
    IF(ALLOCATED(P)) CALL P(i)%free()
    IF(ALLOCATED(V)) CALL V(i)%free()
  END DO
  DO i=-1,0
    IF(ALLOCATED(F)) CALL F(i)%free()
  END DO
  CALL sgrid%free()
  IF(ALLOCATED(BFF)) THEN
    CALL BFF%free()
    DEALLOCATE(BFF)
  END IF

  SDEALLOCATE(U)
  SDEALLOCATE(P)
  SDEALLOCATE(V)
  SDEALLOCATE(F)
  SDEALLOCATE(X1_BC_type)
  SDEALLOCATE(X2_BC_type)
  SDEALLOCATE(LA_BC_type)
  SDEALLOCATE(X1_b)
  SDEALLOCATE(X2_b)
  SDEALLOCATE(X1pert_b)
  SDEALLOCATE(X2pert_b)
  SDEALLOCATE(LA_b)
  SDEALLOCATE(X1_a)
  SDEALLOCATE(X2_a)

  SDEALLOCATE(iota_profile)
  SDEALLOCATE(pres_profile)
  SDEALLOCATE(Phi_profile)
  SDEALLOCATE(chi_profile)

  CALL FinalizeMHD3D_EvalFunc()
  IF(which_init.EQ.1) CALL FinalizeVMEC()

  SDEALLOCATE(hmap)
  SDEALLOCATE(hmap_auxvar)

  SDEALLOCATE(X1_base)
  SDEALLOCATE(X2_base)
  SDEALLOCATE(LA_base)
  CALL exit_subregion("finalize-MHD3D")
END SUBROUTINE FinalizeMHD3D

END MODULE MODgvec_MHD3D
