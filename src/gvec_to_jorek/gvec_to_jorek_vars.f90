!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================

!===================================================================================================================================
!>
!!# Module ** gvec_to_jorek Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_jorek_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
USE MODgvec_base   ,ONLY: t_base
USE MODgvec_fbase  ,ONLY: t_fbase
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
!INPUT VARIABLES
CHARACTER(LEN=255) :: gvecfileName           !< name of GVEC file
CHARACTER(LEN=255) :: fileNameOut            !< name of output file
INTEGER            :: Ns_out                 !< number of equidistant points in radial s-direction (includes axis and edge!)
INTEGER            :: npfactor               !< factor theta,zeta resolution Ntheta=Factor*m_max, Nzeta=MAX(1,Factor*n_max)
INTEGER            :: SFLcoord               !< which angular coordinates to choose: =0: GVEC coord. (no SFL), =1: PEST SFL, =2: BOOZER SFL
!INTEGER            :: factorSFL              !< factor for SFL coordinates, mn_max_sfl=mn_max*factorSFL, default=3
REAL(wp)            :: factorField            !< factor for output field representation, mn_max_out=mn_max*factorField, default=1
REAL(wp)            :: s_max                 !< radial range goes from [0,1]*smax, thats THE RADIAL LIKE COORDINATE. 0 < smax<=1.0 , default=1
CHARACTER(LEN=700) :: cmdline                !< full command line stored
LOGICAL            :: generate_test_data     !< Determine whether to generate fourier representation, or test data for JOREK
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTEGER               :: nfp_out            !< number of field periods
INTEGER               :: asym_out           !< =0: symmetric configuration (R~cos,Z~sin,lambda~sin), =1 asymmetric
INTEGER               :: mn_max_out(2)      !< maximum number of modes in m,n
INTEGER               :: fac_nyq_fields     !< nyquist factor for field fourier representations constructed for JOREK
INTEGER               :: Nthet_out          !< total number of points in theta direction theta[0,2pi (
INTEGER               :: Nzeta_out          !< total number of points in zeta direction zeta[0,-2pi/NFP( opposite sign compared to GVEC!! SET INTERNALLY
REAL(wp),ALLOCATABLE  :: s_pos(:)           !< positions in s for evaluation s=sqrt(phi/phiEdge), size (Ns_out)
REAL(wp),ALLOCATABLE  :: thet_pos(:)        !< positions in theta for evaluation, size (Nthet_out)
REAL(wp),ALLOCATABLE  :: zeta_pos(:)        !< positions in zeta for evaluation , size (Nzeta_out)
INTEGER               :: n_modes            !< total number of toroidal modes in output (from fbase_zeta)
INTEGER               :: sin_range(2)       !< start/end position of sin of toroidal modes in output (from fbase_zeta)
INTEGER               :: cos_range(2)       !< start/end position of cos of toroidal modes in output (from fbase_zeta)

!base needed for evaluation on increased mode numbers / integration points
CLASS(t_Base),ALLOCATABLE           :: out_base        !< same full base for  all output fields
TYPE(t_fbase),ALLOCATABLE          :: X1_fbase_nyq    !< same as X1_base_r%f exept integration points
TYPE(t_fbase),ALLOCATABLE          :: X2_fbase_nyq    !< same as X2_base_r%f exept integration points
TYPE(t_fbase),ALLOCATABLE          :: LA_fbase_nyq    !< same as LA_base_r%f exept integration points
TYPE(t_fbase),ALLOCATABLE          :: fbase_zeta      !< base for doing a fourier transform in zeta only, for 1D toroidal representation
!1D data - Unnecessary for initial JOREK import
!INTEGER,PARAMETER     :: nVar1D = 0          !< number of variables in 1d profiles
!INTEGER,PARAMETER     :: SPOS__    = 1
!INTEGER,PARAMETER     :: PHI__     = 2
!INTEGER,PARAMETER     :: DPHIDS__  = 3
!INTEGER,PARAMETER     :: CHI__     = 4
!INTEGER,PARAMETER     :: DCHIDS__  = 5
!INTEGER,PARAMETER     :: IOTA__    = 6
!INTEGER,PARAMETER     :: PRESSURE__= 7
!INTEGER,PARAMETER     :: ITOR__    = 8
!INTEGER,PARAMETER     :: IPOL__    = 9
!INTEGER,PARAMETER     :: FAVG__    =10
!INTEGER,PARAMETER     :: FMIN__    =11
!INTEGER,PARAMETER     :: FMAX__    =12
!CHARACTER(LEN=50),DIMENSION(nVar1D),PARAMETER :: StrVarNames1D(nVar1D)=(/ CHARACTER(LEN=50) :: &
!                           's'            & ! 1 : position s =sqrt(phi/phiEdge) [0,1]
!                          ,'Phi'          & ! 2 : toroidal flux
!                          ,'dPhi_ds'      & ! 3 : derivative of toroidal flux to s coordinate
!                          ,'Chi'          & ! 4 : poloidal flux
!                          ,'dChi_ds'      & ! 5 : derivative of poloidal flux to s coordinate
!                          ,'iota'         & ! 6 : iota profile
!                          ,'Pressure'     & ! 7 : pressure
!                          ,'Itor'         & ! 8 : Toroidal current
!                          ,'Ipol'         & ! 9 : Poloidal current
!                          ,'Favg'         & !10 : Only tokamaks(n=0!), toroidal magnetic field strength is F/R (averaged over theta)
!                          ,'Fmin'         & !11 : F(s) is averaged over theta, Fmin(s) = min(F(s,theta))
!                          ,'Fmax'         & !12 : F(s) is averaged over theta, Fmax(s) = max(F(s,theta))
!                                    /)
!REAL(wp),ALLOCATABLE  :: data_1D(:,:)        !< 1D profiles size (nVar1D,Ns_out)

! 2D scalar data
INTEGER,PARAMETER     :: nVarScalar2D   = 48         !< number of variables in 2D data
INTEGER, PARAMETER    :: R__            = 1
INTEGER, PARAMETER    :: R_S__          = 2
INTEGER, PARAMETER    :: R_T__          = 3
INTEGER, PARAMETER    :: R_ST__         = 4
INTEGER, PARAMETER    :: Z__            = 5
INTEGER, PARAMETER    :: Z_S__          = 6
INTEGER, PARAMETER    :: Z_T__          = 7
INTEGER, PARAMETER    :: Z_ST__         = 8
INTEGER, PARAMETER    :: P2D__          = 9
INTEGER, PARAMETER    :: P2D_S__        = 10  !<
INTEGER, PARAMETER    :: P2D_T__        = 11
INTEGER, PARAMETER    :: P2D_ST__       = 12
INTEGER, PARAMETER    :: A_R2D__        = 13
INTEGER, PARAMETER    :: A_R2D_S__      = 14
INTEGER, PARAMETER    :: A_R2D_T__      = 15
INTEGER, PARAMETER    :: A_R2D_ST__     = 16
INTEGER, PARAMETER    :: A_Z2D__        = 17
INTEGER, PARAMETER    :: A_Z2D_S__      = 18
INTEGER, PARAMETER    :: A_Z2D_T__      = 19
INTEGER, PARAMETER    :: A_Z2D_ST__     = 20
INTEGER, PARAMETER    :: A_phi2D__      = 21
INTEGER, PARAMETER    :: A_phi2D_S__    = 22
INTEGER, PARAMETER    :: A_phi2D_T__    = 23
INTEGER, PARAMETER    :: A_phi2D_ST__   = 24
INTEGER, PARAMETER    :: B_R2D__        = 25
INTEGER, PARAMETER    :: B_R2D_S__      = 26 !<
INTEGER, PARAMETER    :: B_R2D_T__      = 27
INTEGER, PARAMETER    :: B_R2D_ST__     = 28 !<
INTEGER, PARAMETER    :: B_Z2D__        = 29
INTEGER, PARAMETER    :: B_Z2D_S__      = 30
INTEGER, PARAMETER    :: B_Z2D_T__      = 31
INTEGER, PARAMETER    :: B_Z2D_ST__     = 32
INTEGER, PARAMETER    :: B_phi2D__      = 33  !<
INTEGER, PARAMETER    :: B_phi2D_S__    = 34  !<
INTEGER, PARAMETER    :: B_phi2D_T__    = 35  !<
INTEGER, PARAMETER    :: B_phi2D_ST__   = 36  !<
INTEGER, PARAMETER    :: J_R2D__        = 37  !<<<
INTEGER, PARAMETER    :: J_R2D_S__      = 38  !<<<
INTEGER, PARAMETER    :: J_R2D_T__      = 39  !<<<
INTEGER, PARAMETER    :: J_R2D_ST__     = 40  !<<<
INTEGER, PARAMETER    :: J_Z2D__        = 41  !<<<
INTEGER, PARAMETER    :: J_Z2D_S__      = 42  !<<<
INTEGER, PARAMETER    :: J_Z2D_T__      = 43  !<<<
INTEGER, PARAMETER    :: J_Z2D_ST__     = 44  !<<<
INTEGER, PARAMETER    :: J_phi2D__      = 45  !<<<
INTEGER, PARAMETER    :: J_phi2D_S__    = 46  !<<<
INTEGER, PARAMETER    :: J_phi2D_T__    = 47  !<<<
INTEGER, PARAMETER    :: J_phi2D_ST__   = 48  !<<<

CHARACTER(LEN=50),DIMENSION(nVarScalar2D),PARAMETER :: StrVarNamesScalar2D(nVarScalar2D)=(/ CHARACTER(LEN=50) :: &
                        'R'             & ! 1  : Major radius
                        ,'R_s'          & ! 2  : radial derivative of major radius
                        ,'R_t'          & ! 3  : poloidal derivative of major radius
                        ,'R_st'         & ! 4  : cross derivative of major radius
                        ,'Z'            & ! 5  : Vertical position
                        ,'Z_s'          & ! 6  : radial derivative of vertical position
                        ,'Z_t'          & ! 7  : poloidal derivative of vertical position
                        ,'Z_st'         & ! 8  : cross derivative of vertical position
                        ,'P'            & ! 9  : Pressure
                        ,'P_s'          & ! 10 : radial derivative of pressure
                        ,'P_t'          & ! 11 : poloidal derivative of pressure
                        ,'P_st'         & ! 12 : cross derivative of pressure
                        ,'A_R'          & ! 13 : X Component of  vector potential
                        ,'A_R_s'        & ! 14 : radial derivative of radial vector potential
                        ,'A_R_t'        & ! 15 : poloidal derivative of radial vector potential
                        ,'A_R_st'       & ! 16 : cross derivative of radial vector potential
                        ,'A_Z'          & ! 17 : Y Component of  vector potential
                        ,'A_Z_s'        & ! 18 : radial derivative of radial vector potential
                        ,'A_Z_t'        & ! 19 : poloidal derivative of radial vector potential
                        ,'A_Z_st'       & ! 20 : cross derivative of radial vector potential
                        ,'A_phi'        & ! 21 : Vertical vector potential
                        ,'A_phi_s'      & ! 22 : radial derivative of vertical vector potential
                        ,'A_phi_t'      & ! 23 : poloidal derivative of vertical vector potential
                        ,'A_phi_st'     & ! 24 : cross derivative of vertical vector potential
                        ,'B_R'          & ! 25 : X Component of  magnetic field
                        ,'B_R_s'        & ! 26 : radial derivative of radial magnetic field
                        ,'B_R_t'        & ! 27 : poloidal derivative of radial magnetic field
                        ,'B_R_st'       & ! 28 : cross derivative of radial magnetic field
                        ,'B_Z'          & ! 29 : Y Component of  magnetic field
                        ,'B_Z_s'        & ! 30 : radial derivative of radial magnetic field
                        ,'B_Z_t'        & ! 31 : poloidal derivative of radial magnetic field
                        ,'B_Z_st'       & ! 32 : cross derivative of radial magnetic field
                        ,'B_phi'        & ! 33 : Vertical magnetic field
                        ,'B_phi_s'      & ! 34 : radial derivative of vertical magnetic field
                        ,'B_phi_t'      & ! 35 : poloidal derivative of vertical magnetic field
                        ,'B_phi_st'     & ! 36 : cross derivative of vertical magnetic field
                        ,'J_R'          & ! 37 : X Component of  current density
                        ,'J_R_s'        & ! 38 : radial derivative of radial current density
                        ,'J_R_t'        & ! 39 : poloidal derivative of radial current density
                        ,'J_R_st'       & ! 40 : cross derivative of radial current density
                        ,'J_Z'          & ! 41 : Y Component of  current density
                        ,'J_Z_s'        & ! 42 : radial derivative of radial current density
                        ,'J_Z_t'        & ! 43 : poloidal derivative of radial current density
                        ,'J_Z_st'       & ! 44 : cross derivative of radial current density
                        ,'J_phi'        & ! 45 : Vertical current density
                        ,'J_phi_s'      & ! 46 : radial derivative of vertical current density
                        ,'J_phi_t'      & ! 47 : poloidal derivative of vertical current density
                        ,'J_phi_st'     & ! 48 : cross derivative of vertical current density
                            /)
REAL(wp),ALLOCATABLE  :: data_scalar2D(:,:,:,:)    !< Size (Nthet_out,Ns_out,1:n_modes,nVarScalar2D)

INTEGER,DIMENSION(nVarScalar2D)  :: map_vars_3D_2D  !< map variables from 3D data container to 2D+fourier data container

!3D scalar data
INTEGER,PARAMETER     :: nVarScalar3D = 49           !< number of variables in 3D data
INTEGER,PARAMETER     :: S__          = 1
INTEGER,PARAMETER     :: THET__       = 2
INTEGER,PARAMETER     :: ZETA__       = 3
INTEGER,PARAMETER     :: X1__         = 4
INTEGER,PARAMETER     :: X1_S__       = 5
INTEGER,PARAMETER     :: X1_T__       = 6
INTEGER,PARAMETER     :: X1_ST__      = 7
INTEGER,PARAMETER     :: X2__         = 8
INTEGER,PARAMETER     :: X2_S__       = 9
INTEGER,PARAMETER     :: X2_T__       = 10
INTEGER,PARAMETER     :: X2_ST__      = 11
INTEGER,PARAMETER     :: P__          = 12
INTEGER,PARAMETER     :: P_S__        = 13
INTEGER,PARAMETER     :: A_R__        = 14
INTEGER,PARAMETER     :: A_R_S__      = 15
INTEGER,PARAMETER     :: A_R_T__      = 16
INTEGER,PARAMETER     :: A_R_ST__     = 17
INTEGER,PARAMETER     :: A_Z__        = 18
INTEGER,PARAMETER     :: A_Z_S__      = 19
INTEGER,PARAMETER     :: A_Z_T__      = 20
INTEGER,PARAMETER     :: A_Z_ST__     = 21
INTEGER,PARAMETER     :: A_phi__      = 22
INTEGER,PARAMETER     :: A_phi_S__    = 23
INTEGER,PARAMETER     :: A_phi_T__    = 24
INTEGER,PARAMETER     :: A_phi_ST__   = 25
INTEGER,PARAMETER     :: B_R__        = 26
INTEGER,PARAMETER     :: B_R_S__      = 27
INTEGER,PARAMETER     :: B_R_T__      = 28
INTEGER,PARAMETER     :: B_R_ST__     = 29
INTEGER,PARAMETER     :: B_Z__        = 30
INTEGER,PARAMETER     :: B_Z_S__      = 31
INTEGER,PARAMETER     :: B_Z_T__      = 32
INTEGER,PARAMETER     :: B_Z_ST__     = 33
INTEGER,PARAMETER     :: B_phi__      = 34
INTEGER,PARAMETER     :: B_phi_S__    = 35
INTEGER,PARAMETER     :: B_phi_T__    = 36
INTEGER,PARAMETER     :: B_phi_ST__   = 37
INTEGER,PARAMETER     :: J_R__        = 38
INTEGER,PARAMETER     :: J_R_S__      = 39
INTEGER,PARAMETER     :: J_R_T__      = 40
INTEGER,PARAMETER     :: J_R_ST__     = 41
INTEGER,PARAMETER     :: J_Z__        = 42
INTEGER,PARAMETER     :: J_Z_S__      = 43
INTEGER,PARAMETER     :: J_Z_T__      = 44
INTEGER,PARAMETER     :: J_Z_ST__     = 45
INTEGER,PARAMETER     :: J_phi__      = 46
INTEGER,PARAMETER     :: J_phi_S__    = 47
INTEGER,PARAMETER     :: J_phi_T__    = 48
INTEGER,PARAMETER     :: J_phi_ST__   = 49
CHARACTER(LEN=50),DIMENSION(nVarScalar3D),PARAMETER :: StrVarNamesScalar3D(nVarScalar3D)=(/ CHARACTER(LEN=50) :: &
                          'S'               & ! 1  : Radial coordinate
                          ,'THET'           & ! 2  : Poloidal coordinate
                          ,'ZETA'           & ! 3  : Toroidal coordinate
                          ,'X1(R)'          & ! 4  : for Torus map (hmap=1), R=X1
                          ,'X1_s(R)'        & ! 5  : radial derivative of R
                          ,'X1_t(R)'        & ! 6  : poloidal derivative of R
                          ,'X1_st(R)'       & ! 7  : cross derivative of R
                          ,'X2(Z)'          & ! 8  : for Torus map (hmap=1), Z=X2
                          ,'X2_s(Z)'        & ! 9  : radial derivative of Z
                          ,'X2_t(Z)'        & ! 10 : poloidal derivative of Z
                          ,'X2_st(Z)'       & ! 11 : cross derivative of Z
                          ,'P'              & ! 12 : Pressure
                          ,'P_s'            & ! 13 : radial derivative of Pressure
                          ,'A_R'            & ! 14 : X component vector potential component
                          ,'A_R_s'          & ! 15 : X component vector potential component s derivative
                          ,'A_R_t'          & ! 16 : X component vector potential component t derivative
                          ,'A_R_st'         & ! 17 : X component vector potential component st derivative
                          ,'A_Z'            & ! 18 : Y component vector potential component
                          ,'A_Z_s'          & ! 19 : Y component vector potential component s derivative
                          ,'A_Z_t'          & ! 20 : Y component vector potential component t derivative
                          ,'A_Z_st'         & ! 21 : Y component vector potential component st derivative
                          ,'A_phi'          & ! 22 : vertical vector potential component
                          ,'A_phi_s'        & ! 23 : vertical vector potential component s derivative
                          ,'A_phi_t'        & ! 24 : vertical vector potential component t derivative
                          ,'A_phi_st'       & ! 25 : vertical vector potential component st derivative
                          ,'B_R'            & ! 26 : X component magnetic field component
                          ,'B_R_s'          & ! 27 : X component magnetic field component s derivative
                          ,'B_R_t'          & ! 28 : X component magnetic field component t derivative
                          ,'B_R_st'         & ! 29 : X component magnetic field component st derivative
                          ,'B_Z'            & ! 30 : Y component magnetic field component
                          ,'B_Z_s'          & ! 31 : Y component magnetic field component s derivative
                          ,'B_Z_t'          & ! 32 : Y component magnetic field component t derivative
                          ,'B_Z_st'         & ! 33 : Y component magnetic field component st derivative
                          ,'B_phi'          & ! 34 : vertical magnetic field component
                          ,'B_phi_s'        & ! 35 : vertical magnetic field component s derivative
                          ,'B_phi_t'        & ! 36 : vertical magnetic field component t derivative
                          ,'B_phi_st'       & ! 37 : vertical magnetic field component st derivative
                          ,'J_R'            & ! 38 : X component current density component
                          ,'J_R_s'          & ! 39 : X component current density component s derivative
                          ,'J_R_t'          & ! 40 : X component current density component t derivative
                          ,'J_R_st'         & ! 41 : X component current density component st derivative
                          ,'J_Z'            & ! 42 : Y component current density component
                          ,'J_Z_s'          & ! 43 : Y component current density component s derivative
                          ,'J_Z_t'          & ! 44 : Y component current density component t derivative
                          ,'J_Z_st'         & ! 45 : Y component current density component st derivative
                          ,'J_phi'          & ! 46 : vertical current density component
                          ,'J_phi_s'        & ! 47 : vertical current density component s derivative
                          ,'J_phi_t'        & ! 48 : vertical current density component t derivative
                          ,'J_phi_st'       & ! 49 : vertical current density component st derivative
                             /)
REAL(wp),ALLOCATABLE  :: data_scalar3D(:,:,:,:)    !< Size (Nthet_out,Nzeta_out,Ns_out,nVar3D)

!!3D vector data - Unnecessary for initial JOREK import
!INTEGER,PARAMETER     :: nVarVector3D=0           !< number of variables in 3D data
!INTEGER,PARAMETER     :: BFIELD__     = 1
!INTEGER,PARAMETER     :: ECOV_S__     = 2
!INTEGER,PARAMETER     :: ECOV_THETA__ = 3
!INTEGER,PARAMETER     :: ECOV_ZETA__  = 4
!CHARACTER(LEN=50),DIMENSION(nVarVector3D),PARAMETER :: StrVarNamesVector3D(nVarVector3D)=(/ CHARACTER(LEN=50) :: &
!                           'Bfield'           & ! 1 : magnetic field vector    , (x,y,z) cartesian components
!                          ,'ecov_s'           & ! 2 : covariant vector in s    , (x,y,z) cartesian components
!                          ,'ecov_theta'       & ! 3 : covariant vector in theta, (x,y,z) cartesian components
!                          ,'ecov_zeta'        & ! 4 : covariant vector in zeta , (x,y,z) cartesian components
!                                    /)
!REAL(wp),ALLOCATABLE  :: data_vector3D(:,:,:,:,:)    !< Size (3,Nthet_out,Nzeta_out,Ns_out,nVarVector3D)

!===================================================================================================================================

CONTAINS





END MODULE MODgvec_gvec_to_jorek_Vars
