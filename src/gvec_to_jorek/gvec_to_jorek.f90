!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **gvec_to_jorek**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_jorek
! MODULES
USE MODgvec_Globals, ONLY:wp,UNIT_stdOut,fmt_sep
IMPLICIT NONE
PRIVATE

INTERFACE get_cla_gvec_to_jorek
  MODULE PROCEDURE get_cla_gvec_to_jorek
END INTERFACE

INTERFACE init_gvec_to_jorek
  MODULE PROCEDURE init_gvec_to_jorek
END INTERFACE

INTERFACE gvec_to_jorek_prepare
  MODULE PROCEDURE gvec_to_jorek_prepare
END INTERFACE

INTERFACE gvec_to_jorek_writeToFile
  MODULE PROCEDURE gvec_to_jorek_writeToFile_ASCII
END INTERFACE

INTERFACE finalize_gvec_to_jorek
  MODULE PROCEDURE finalize_gvec_to_jorek
END INTERFACE

PUBLIC::get_cla_gvec_to_jorek
PUBLIC::init_gvec_to_jorek
PUBLIC::gvec_to_jorek_prepare
PUBLIC::gvec_to_jorek_writeToFile
PUBLIC::finalize_gvec_to_jorek

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get command line arguments
!!
!===================================================================================================================================
SUBROUTINE get_CLA_gvec_to_jorek()
! MODULES
USE MODgvec_cla
USE MODgvec_gvec_to_jorek_Vars, ONLY: gvecfileName,FileNameOut
USE MODgvec_gvec_to_jorek_Vars, ONLY: Ns_out,npfactor,s_max,factorField,Nthet_out,SFLcoord,cmdline, generate_test_data
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=STRLEN)   :: f_str
CHARACTER(LEN=24)       :: execname="convert_gvec_to_jorek"
LOGICAL                 :: commandFailed
CHARACTER(LEN=6),DIMENSION(0:2),PARAMETER :: SFLcoordName=(/" GVEC "," PEST ","BOOZER"/)
!===================================================================================================================================

  !USING CLAF90 module to get command line arguments!
  CALL cla_init()

  CALL cla_register('-r',  '--rpoints', &
       'Number of radial points in s=[0,1] for output [MANDATORY, >1]', cla_int,'2') !must be provided
  CALL cla_register('-n',  '--npfactor', &
       'Number of angular points, computed from max. mode numbers (=npfactor*mn_max_in) [DEFAULT = 4]',cla_int,'4')
  CALL cla_register('-p',  '--polpoints', &
       'Number of poloidal points, if specified overwrites factor*m_max  [OPTIONAL]',cla_int,'-1')
!  CALL cla_register('-t',  '--torpoints', &
!       'Number of toroidal points, if specified overwrites factor*n_max  [OPTIONAL]',cla_int,'-1')
  CALL cla_register('-s',  '--sflcoord', &
       'which angular coordinates to choose: =0: GVEC coord. (no SFL), =1: PEST SFL, =2: BOOZER SFL [DEFAULT = 0]', &
       cla_int,'0')
  CALL cla_register('-f',  '--factorfield', &
       'factor(real number) on max. mode numbers (mn_max_out=factorfield*mn_max_in) in output representation. [DEFAULT = 1.0]',&
             cla_float,'1.0')
  CALL cla_register('-x',  '--smax', &
       'radial range goes from [0,1]*smax, thats THE RADIAL LIKE COORDINATE.0 < smax<=1.0 [DEFAULT = 1.0]',&
             cla_float,'1.0')
  CALL cla_register('-g',  '--generate_test_data', &
       'determine whether test data is generated  [DEFAULT = FALSE]', cla_int,'.false.')
  !positional argument
  CALL cla_posarg_register('gvecfile.dat', &
       'Input filename of GVEC restart file [MANDATORY]',  cla_char,'xxx') !
  CALL cla_posarg_register('outfile.dat', &
       'Output filename  [OPTIONAL, DEFAULT: gvec2jorek_nameofgvecfile]',  cla_char,'yyy') !

  CALL cla_validate(execname)
  CALL cla_get('-r',NS_out)
  CALL cla_get('-n',npfactor)
  CALL cla_get('-p',Nthet_out)
!!!!!  CALL cla_get('-t',Nzeta_out)
  CALL cla_get('-s',SFLcoord)
  CALL cla_get('-f',factorField)
  CALL cla_get('-x',s_max)
  CALL cla_get('-g',generate_test_data)
  CALL cla_get('gvecfile.dat',f_str)
  gvecfilename=TRIM(f_str)
  CALL cla_get('outfile.dat',f_str)
  FileNameOut=TRIM(f_str)

  commandFailed=.FALSE.
  IF(.NOT.((cla_key_present('-r')).AND.(Ns_out.GE.2))) THEN
    IF (.not. generate_test_data) THEN
      IF(.NOT.commandFailed) CALL cla_help(execname)
      commandFailed=.TRUE.
      SWRITE(UNIT_StdOut,*) " ==> [-r,--rpoints] argument is MANDATORY and must be >1 !!!"
    END IF
  END IF
  IF((SFLcoord.LT.0).OR.(SFLcoord.GT.2)) THEN
    IF(.NOT.commandFailed) CALL cla_help(execname)
    commandFailed=.TRUE.
    SWRITE(UNIT_StdOut,*) " ==> [-s,--sflcoord] argument  must be 0,1,2 !!!"
  END IF
  IF((INDEX(gvecfilename,'xxx').NE.0))THEN
    IF(.NOT.commandFailed) CALL cla_help(execname)
    commandFailed=.TRUE.
    SWRITE(UNIT_StdOut,*) " ==> input gvec filename is MANDATORY must be specified as first positional argument!!!"
  END IF
  IF((INDEX(FileNameOut,'yyy').NE.0))THEN
    FileNameOut="gvec2jorek_"//TRIM(gvecfilename)
  END IF
  IF((s_max.GT.1.0_wp).OR.(s_max.LE.0.0_wp)) THEN
    commandFailed=.TRUE.
    SWRITE(UNIT_StdOut,*) " ==> input parameter smax must be 0.0<smax <=1.0!!!"
  END IF
  IF(commandFailed) STOP

  SWRITE(UNIT_stdOut,'(A)')     ' INPUT PARAMETERS:'
  SWRITE(UNIT_stdOut,'(A,I6)')  '  * Number of radial points        : ',Ns_out
  SWRITE(UNIT_stdOut,'(A,I4)')  '  * npfactor points from modes     : ',npfactor
  IF(Nthet_out.NE.-1) THEN
    SWRITE(UNIT_stdOut,'(A,I4)')'  * number of points in theta      : ',Nthet_out
  END IF
!  IF(Nzeta_out.NE.-1) THEN
!    SWRITE(UNIT_stdOut,'(A,I4)')'  * number of points in zeta       : ',Nzeta_out
!  END IF
  SWRITE(UNIT_stdOut,'(A,E11.4)')'  * factor for output modes        : ',factorField
  SWRITE(UNIT_stdOut,'(A,I4,A)')'  * SFL coordinates flag           : ',SFLcoord,' ( '//SFLcoordName(SFLcoord)//' )'
  SWRITE(UNIT_stdOut,'(A,A)')   '  * generate test data for JOREK   : ',MERGE(".false.",".true. ",generate_test_data)
  SWRITE(UNIT_stdOut,'(A,A)')   '  * GVEC input file                : ',TRIM(gvecfileName)
  SWRITE(UNIT_stdOut,'(A,A)')   '  * output file name               : ',TRIM(FileNameOut)
  SWRITE(UNIT_stdOut,fmt_sep)
  IF(SFLcoord.NE.0)THEN
    SWRITE(UNIT_StdOut,*) "SFLcoord=",SFLcoord,", but only GVEC coordinates (=0) is currently implemented!"
    STOP
  END IF

  CALL GET_COMMAND(cmdline)

END SUBROUTINE get_cla_gvec_to_jorek

!===================================================================================================================================
!> Initialize Module
!!
!===================================================================================================================================
SUBROUTINE init_gvec_to_jorek()
! MODULES
USE MODgvec_Globals,ONLY: TWOPI
USE MODgvec_ReadState         ,ONLY: ReadState
USE MODgvec_ReadState_vars    ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_ReadState_vars    ,ONLY: LA_r,X1_r,X2_r
!USE MODgvec_transform_sfl_vars,ONLY: X1sfl_base,X1sfl,X2sfl_base,X2sfl ,GZsfl_base,GZsfl
!USE MODgvec_transform_sfl     ,ONLY: BuildTransform_SFL
USE MODgvec_gvec_to_jorek_vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: i
REAL(wp) :: r
INTEGER  :: seed=1139           ! Random seed for test data
REAL(wp) :: phi_direction=1     ! direction of phi in JOREK and GVEC is clockwise, so direction does not need to be flipped
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT GVEC-TO-CASTOR3D ...'

  ! Initialize grid variables from GVEC restart file
  CALL ReadState(TRIM(gvecfileName))

  mn_max_out(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
  mn_max_out(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))
  nfp_out          = X1_base_r%f%nfp
  ! increase modal respresentation of the computed fields
  mn_max_out(1) = NINT(mn_max_out(1)*factorField)
  mn_max_out(2) = NINT(mn_max_out(2)*factorField)
  fac_nyq_fields=4 !hard coded for now

  IF((X1_base_r%f%sin_cos.EQ._COS_).AND.(X2_base_r%f%sin_cos.EQ._SIN_).AND.(LA_base_r%f%sin_cos.EQ._SIN_))THEN
    asym_out = 0 !R~cos,Z~sin,lambda~sin
  ELSE
    asym_out = 1 !full fourier
  END IF

  ! Initialize sample points in s, theta, zeta. s and theta can be randomly sampled for testing purposes.
  ALLOCATE(s_pos(Ns_out))
  !ALLOCATE(data_1D(nVar1D,Ns_out))
  IF (generate_test_data) THEN
    call RANDOM_SEED(seed)
    DO i=1,Ns_out
      call RANDOM_NUMBER(r)
      WRITE(*,*) "Random number: ", r
      s_pos(i) = (1.0 - 1.0e-12_wp - 1.0e-08_wp) * r + 1.0e-08_wp
      !s_pos(i) = r
    END DO
    IF (Ns_out .eq. 1) s_pos(1) = 1.0
  ELSE
    s_pos(1)=1.0e-08_wp !avoid axis
    DO i=2,Ns_out-1
        s_pos(i) = REAL(i-1,wp)/REAL(Ns_out-1,wp)
    END DO !i
    s_pos(Ns_out)=1. - 1.0e-12_wp !avoid edge
  END IF

  IF(Nthet_out.EQ.-1) THEN !overwrite with default from factorFourier
    Nthet_out = npfactor*mn_max_out(1)
  END IF
  IF(Nthet_out .LT. 4*mn_max_out(1)) WRITE(UNIT_StdOut,'(A)')'WARNING: number of poloidal points for output should be >=4*m_max!'
  Nzeta_out = MAX(1,fac_nyq_fields*mn_max_out(2)) !if n=0, output 1 point

  ALLOCATE(thet_pos(Nthet_out))
  ALLOCATE(zeta_pos(Nzeta_out))
  IF (generate_test_data) THEN
   DO i=1,Nthet_out
     call RANDOM_NUMBER(r)
     thet_pos(i)=r
   END DO
   IF (Nthet_out .eq. 1) thet_pos(1) =0.0
  ELSE
    DO i=1,Nthet_out
      thet_pos(i)=(REAL((i-1),wp))/REAL(Nthet_out,wp)
    END DO
  END IF
  DO i=1,Nzeta_out
    zeta_pos(i)=phi_direction * (TWOPI*REAL((i-0.5),wp))/REAL((Nzeta_out*nfp_out),wp)
  END DO

  CALL Init_Base(mn_max_out,fac_nyq_fields)

  n_modes      = fbase_zeta%modes
  sin_range(:) = fbase_zeta%sin_range(:)
  cos_range(:) = fbase_zeta%cos_range(:)
  ALLOCATE(data_scalar2D(Nthet_out,  Ns_out,n_modes,nVarScalar2D))
  ALLOCATE(data_scalar3D(Nthet_out,Nzeta_out,Ns_out,nVarscalar3D))
  !ALLOCATE(data_vector3D(3,Nthet_out,Nzeta_out,Ns_out,nVarvector3D))

  SWRITE(UNIT_stdOut,'(A,3I6)')'  Number OF N_s,N_theta,N_zeta evaluation points:',Ns_out,Nthet_out,Nzeta_out
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)


  SELECT CASE(SFLcoord)
  CASE(0) ! GVEC coordinates - toroidal coordinate is the cylindrical toroidal direction
    CALL gvec_to_jorek_prepare(X1_base_r,X1_r,X2_base_r,X2_r,LA_base_r,LA_r)
  CASE DEFAULT
    SWRITE(UNIT_StdOut,*)'This SFLcoord is not yet implemented',SFLcoord
    STOP
  END SELECT
END SUBROUTINE init_gvec_to_jorek

!===================================================================================================================================
!> initialize base classes declared in _vars module, needed for computation of output fields
!!
!===================================================================================================================================
SUBROUTINE Init_Base(mn_max,fac_nyq)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut
USE MODgvec_base   ,ONLY: base_new
USE MODgvec_fbase  ,ONLY: fbase_new,sin_cos_map
USE MODgvec_ReadState_vars  ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_gvec_to_jorek_vars, ONLY: X1_fbase_nyq,X2_fbase_nyq,LA_fbase_nyq,out_base,fbase_zeta,Nzeta_out
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER      ,INTENT(IN) :: mn_max(2)                                     !< maximum number for new variables in SFL coordinates
  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: mn_nyq(2)
!===================================================================================================================================
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max)
  SWRITE(UNIT_StdOut,'(2(A,2I6))')'INITIALIZE OUTPUT BASE, mn_max_out=',mn_max,', mn_int=',mn_nyq

  ! Initialize basis for field_representation based on existing grid representation
  CALL base_new(out_base,  X1_base_r%s%deg,        &
                           X1_base_r%s%continuity, &
                           X1_base_r%s%grid,       &
                           X1_base_r%s%degGP,      &
                           mn_max,mn_nyq,          &
                           X1_base_r%f%nfp,        &
                           '_sincos_   ', &  !full basis
                           .False.) !do not exclude m=n=0

  CALL fbase_new(fbase_zeta, (/0, mn_max(2)/), (/1, Nzeta_out/), X1_base_r%f%nfp, "_sincos_", .false.)
  ! Initialize bases for existing grid at higher number of integration points, based on nyquist condition
  CALL fbase_new( X1_fbase_nyq, X1_base_r%f%mn_max,  mn_nyq, &
                                X1_base_r%f%nfp, &
                    sin_cos_map(X1_base_r%f%sin_cos), &
                                X1_base_r%f%exclude_mn_zero)

  CALL fbase_new( X2_fbase_nyq, X2_base_r%f%mn_max,  mn_nyq, &
                                X2_base_r%f%nfp, &
                    sin_cos_map(X2_base_r%f%sin_cos), &
                                X2_base_r%f%exclude_mn_zero)

  CALL fbase_new(LA_fbase_nyq,  LA_base_r%f%mn_max,  mn_nyq, &
                                LA_base_r%f%nfp, &
                    sin_cos_map(LA_base_r%f%sin_cos), &
                                LA_base_r%f%exclude_mn_zero)
END SUBROUTINE Init_Base


!===================================================================================================================================
!> prepare all data to be written
!!
!===================================================================================================================================
SUBROUTINE gvec_to_jorek_prepare(X1_base_in,X1_in,X2_base_in,X2_in,LG_base_in,LG_in)
! MODULES
USE MODgvec_gvec_to_jorek_Vars
USE MODgvec_Globals,        ONLY: CROSS,TWOPI,ProgressBar
USE MODgvec_ReadState_Vars, ONLY: profiles_1d,hmap_r,sbase_prof !for profiles
USE MODgvec_Base,           ONLY: t_base
USE MODgvec_fBase,          ONLY: t_fbase, fbase_new


IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CLASS(t_base) ,INTENT(IN) :: X1_base_in,X2_base_in,LG_base_in
REAL(wp)      ,INTENT(IN) :: X1_in(1:X1_base_in%s%nBase,1:X1_base_in%f%modes)
REAL(wp)      ,INTENT(IN) :: X2_in(1:X2_base_in%s%nBase,1:X2_base_in%f%modes)
REAL(wp)      ,INTENT(IN) :: LG_in(1:LG_base_in%s%nBase,1:LG_base_in%f%modes) ! is either LA if SFLcoord=0/1 or G if SFLcoord=2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i_s,ithet,izeta,iVar,jVar                                        ! Indices for enumeration
REAL(wp)                                 :: spos,xp(2),sqrtG                                                 ! Current s, theta, zeta, position and metric tensor
REAL(wp)                                 :: dX1ds,dX1dthet,dX1dzeta,d2X1dsdthet                              ! X coordinate and derivatives
REAL(wp)                                 :: dX2ds,dX2dthet,dX2dzeta,d2X2dsdthet                              ! Z coordinate and derivatives
REAL(wp)                                 :: dLAdthet,dLAdzeta                                                ! SFL transformation
REAL(wp)                                 :: Phi_int,dPhids_int,Chi_int,dChids_int,iota_int                   ! Flux variables and derivatives
REAL(wp)                                 :: P_int,dPds_int                                                   ! Pressure and derivatives
REAL(wp)                                 :: A_R_int, dA_Rds_int, dA_Rdthet_int, d2A_Rdsdthet_int             ! Vector potential components and derivatives
REAL(wp)                                 :: A_Z_int, dA_Zds_int, dA_Zdthet_int, d2A_Zdsdthet_int
REAL(wp)                                 :: A_phi_int, dA_phids_int, dA_phidthet_int, d2A_phidsdthet_int
REAL(wp)                                 :: B_R_int, dB_Rds_int, dB_Rdthet_int, d2B_Rdsdthet_int             ! Magnetic field components and derivatives
REAL(wp)                                 :: B_Z_int, dB_Zds_int, dB_Zdthet_int, d2B_Zdsdthet_int
REAL(wp)                                 :: B_phi_int, dB_phids_int, dB_phidthet_int, d2B_phidsdthet_int
REAL(wp)                                 :: J_R_int, dJ_Rds_int, dJ_Rdthet_int, d2J_Rdsdthet_int             ! Current density components and derivatives
REAL(wp)                                 :: J_Z_int, dJ_Zds_int, dJ_Zdthet_int, d2J_Zdsdthet_int
REAL(wp)                                 :: J_phi_int, dJ_phids_int, dJ_phidthet_int, d2J_phidsdthet_int
REAL(wp)                                 :: X1_int,X2_int,G_int,dGds,dGdthet,dGdzeta
REAL(wp)                                 :: AR_diff,AZ_diff,Aphi_diff,BR_diff,BZ_diff,Bphi_diff              ! Diagnostics for the convergence of the generated field representation
REAL(wp)                                 :: AR_diff_max,AZ_diff_max,Aphi_diff_max,BR_diff_max,BZ_diff_max,Bphi_diff_max
REAL(wp),DIMENSION(3)                    :: qvec,e_s,e_thet,e_zeta                                           ! Vectors for local covariant coordinate system
REAL(wp),DIMENSION(3)                    :: Acart,A_orig,Bcart,B_orig                                        ! Test variables for magnetic field and position
REAL(wp)                                 :: Bthet, Bzeta
REAL(wp),DIMENSION(3)                    :: grad_s,grad_thet,grad_zeta,grad_R,grad_Z                         ! Contravariant coordinates

! 2D theta/zeta fourier representation of variables in GVEC
REAL(wp),DIMENSION(1:X1_base_in%f%modes) :: X1_s,dX1ds_s
REAL(wp),DIMENSION(1:X2_base_in%f%modes) :: X2_s,dX2ds_s
REAL(wp),DIMENSION(1:out_base%f%modes)   :: A_R_s, A_Rds_int, A_Z_s, A_Zds_int, A_phi_s, A_phids_int
REAL(wp),DIMENSION(1:out_base%f%modes)   :: B_R_s, B_Rds_int, B_Z_s, B_Zds_int, B_phi_s, B_phids_int
REAL(wp),DIMENSION(1:out_base%f%modes)   :: J_R_s, J_Rds_int, J_Z_s, J_Zds_int, J_phi_s, J_phids_int
REAL(wp),DIMENSION(1:LG_base_in%f%modes) :: LG_s


! Variables for new GVEC fourier representations needed for JOREK inputs
REAL(wp),DIMENSION(out_base%s%nBase,out_base%f%modes):: A_R, A_Z, A_phi !< data (1:nBase,1:modes) of A_* in GVEC coords
REAL(wp),DIMENSION(out_base%s%nBase,out_base%f%modes):: B_R, B_Z, B_phi !< data (1:nBase,1:modes) of B_* in GVEC coords
REAL(wp),DIMENSION(out_base%s%nBase,out_base%f%modes):: J_R, J_Z, J_phi !< data (1:nBase,1:modes) of J_* in GVEC coords

!===================================================================================================================================
! ------------------------------------------------------------------------------------------------------------
! ------- CALCULATE 2D FOURIER REPRESENTATION OF Vector Potential, Magnetic Field, and Current Density -------
! ------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)')'PREPARE FIELD BASES ...'
CALL get_field(2,4,A_R)
CALL get_field(2,3,A_Z)
CALL get_field(2,5,A_phi)
CALL get_field(1,4,B_R)
CALL get_field(1,3,B_Z)
CALL get_field(1,5,B_phi)
CALL get_field(3,4,J_R)
CALL get_field(3,3,J_Z)
CALL get_field(3,5,J_phi)

! -----------------------------------------------------------------------------
! ---------- GENERATE 3D POINTS FROM GVEC REPRESENTATION ----------------------
! -----------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)')'PREPARE 3D DATA FOR GVEC-TO-JOREK ...'
CALL ProgressBar(0,Ns_out) !init
DO i_s=1,Ns_out
  !!spos  = s_pos(i_s)
  !! SCALE DOMAIN TO s=s_logical*s_max
  spos = s_pos(i_s)*s_max

  Phi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,1))
  dPhids_int  = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  Chi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,2)) !Chi representation inconsistent - see ReadState routine
  dChids_int  = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2)) !Chi representation inconsistent - see ReadState routine
  iota_int    = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,3))
  P_int       = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,4))
  dPds_int    = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,4))

  dLAdthet= 0.0_wp !only changed for SFLcoord=0
  dLAdzeta= 0.0_wp !only changed for SFLcoord=0
  G_int   = 0.0_wp !only changed for SFLcoords=2
  dGds    = 0.0_wp !only changed for SFLcoords=2
  dGdthet = 0.0_wp !only changed for SFLcoords=2
  dGdzeta = 0.0_wp !only changed for SFLcoords=2

  !interpolate radially
  X1_s(   :) = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes,       0,X1_in(:,:))
  dX1ds_s(:) = X1_base_in%s%evalDOF2D_s(spos,X1_base_in%f%modes, DERIV_S,X1_in(:,:))

  X2_s(   :) = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes,       0,X2_in(:,:))
  dX2ds_s(:) = X2_base_in%s%evalDOF2D_s(spos,X2_base_in%f%modes, DERIV_S,X2_in(:,:))
  IF(SFLcoord.EQ.0)THEN !GVEC coordinates
    LG_s(  :) = LG_base_in%s%evalDOF2D_s(spos,LG_base_in%f%modes,       0,LG_in(:,:))
  ELSE
    SWRITE(UNIT_StdOut,*)'This SFLcoord is not valid',SFLcoord
    STOP
  END IF

  ! Generate representations of the fields in (R,Z,phi)
  A_R_s(:)       =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,         0, A_R(:,:))
  A_Rds_int(:)   =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,   DERIV_S, A_R(:,:))
  A_Z_s(:)       =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,         0, A_Z(:,:))
  A_Zds_int(:)   =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,   DERIV_S, A_Z(:,:))
  A_phi_s(:)     =   out_base%s%evalDOF2D_s(spos, out_base%f%modes,         0, A_phi(:,:))
  A_phids_int(:) =   out_base%s%evalDOF2D_s(spos, out_base%f%modes,   DERIV_S, A_phi(:,:))
  B_R_s(:)       =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,         0, B_R(:,:))
  B_Rds_int(:)   =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,   DERIV_S, B_R(:,:))
  B_Z_s(:)       =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,         0, B_Z(:,:))
  B_Zds_int(:)   =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,   DERIV_S, B_Z(:,:))
  B_phi_s(:)     =   out_base%s%evalDOF2D_s(spos, out_base%f%modes,         0, B_phi(:,:))
  B_phids_int(:) =   out_base%s%evalDOF2D_s(spos, out_base%f%modes,   DERIV_S, B_phi(:,:))
  J_R_s(:)       =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,         0, J_R(:,:))
  J_Rds_int(:)   =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,   DERIV_S, J_R(:,:))
  J_Z_s(:)       =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,         0, J_Z(:,:))
  J_Zds_int(:)   =   out_base%s%evalDOF2D_s(spos,   out_base%f%modes,   DERIV_S, J_Z(:,:))
  J_phi_s(:)     =   out_base%s%evalDOF2D_s(spos, out_base%f%modes,         0, J_phi(:,:))
  J_phids_int(:) =   out_base%s%evalDOF2D_s(spos, out_base%f%modes,   DERIV_S, J_phi(:,:))

  BR_diff=0.0_wp; BZ_diff=0.0_wp; Bphi_diff=0.0_wp; AR_diff=0.0_wp; AZ_diff=0.0_wp; Aphi_diff=0.0_wp
  BR_diff_max=0.0_wp; BZ_diff_max=0.0_wp; Bphi_diff_max=0.0_wp; AR_diff_max=0.0_wp; AZ_diff_max=0.0_wp; Aphi_diff_max=0.0_wp
!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE) COLLAPSE(2)                                                                  &
!$OMP   REDUCTION(+:BR_diff, BZ_diff, Bphi_diff, AR_diff, AZ_diff, Aphi_diff)                                                  &
!$OMP   REDUCTION(max: BR_diff_max, BZ_diff_max, Bphi_diff_max, AR_diff_max,AZ_diff_max,Aphi_diff_max)                         &
!$OMP   PRIVATE(izeta,ithet,X1_int,dX1ds,dX1dthet,dX1dzeta,d2X1dsdthet, X2_int,dX2ds,dX2dthet,dX2dzeta, d2X2dsdthet,           &
!$OMP           xp,qvec,e_s,e_thet,e_zeta,sqrtG,                                                                               &
!$OMP           A_R_int,dA_Rds_int,dA_Rdthet_int,d2A_Rdsdthet_int,                                                             &
!$OMP           A_Z_int,dA_Zds_int,dA_Zdthet_int,d2A_Zdsdthet_int,                                                             &
!$OMP           A_phi_int,dA_phids_int,dA_phidthet_int,d2A_phidsdthet_int,                                                     &
!$OMP           B_R_int,dB_Rds_int,dB_Rdthet_int,d2B_Rdsdthet_int,                                                             &
!$OMP           B_Z_int,dB_Zds_int,dB_Zdthet_int,d2B_Zdsdthet_int,                                                             &
!$OMP           B_phi_int,dB_phids_int,dB_phidthet_int,d2B_phidsdthet_int,                                                     &
!$OMP           J_R_int,dJ_Rds_int,dJ_Rdthet_int,d2J_Rdsdthet_int,                                                             &
!$OMP           J_Z_int,dJ_Zds_int,dJ_Zdthet_int,d2J_Zdsdthet_int,                                                             &
!$OMP           J_phi_int,dJ_phids_int,dJ_phidthet_int,d2J_phidsdthet_int,                                                     &
!$OMP           Acart,A_orig,grad_s,grad_thet,grad_zeta,grad_R,grad_Z,                                                   &
!$OMP           Bthet,Bzeta,Bcart,B_orig)                                                                                      &
!$OMP   FIRSTPRIVATE(dLAdthet,dLAdzeta,G_int,dGds,dGdthet,dGdzeta)                                                             &
!$OMP   SHARED(i_s,Nzeta_out,Nthet_out,spos,s_pos,s_max, thet_pos,zeta_pos,X1_base_in,X2_base_in,LG_base_in,out_base,               &
!$OMP          hmap_r,X1_s,dX1ds_s,X2_s,dX2ds_s,LG_s,SFLcoord, Phi_int, dPhids_int, iota_int, Chi_int, dChids_int, P_int, dPds_int, &
!$OMP          A_R_s, A_Rds_int,A_Z_s, A_Zds_int,A_phi_s, A_phids_int, B_R_s, B_Rds_int,B_Z_s, B_Zds_int,B_phi_s, B_phids_int,      &
!$OMP          J_R_s, J_Rds_int,J_Z_s, J_Zds_int,J_phi_s, J_phids_int,                                                              &
!$OMP          data_scalar3D)
  !interpolate in the angles
  DO izeta=1,Nzeta_out; DO ithet=1,Nthet_out
    xp=(/TWOPI * thet_pos(ithet),zeta_pos(izeta)/)

    X1_int      = X1_base_in%f%evalDOF_x(xp,          0, X1_s  )
    dX1ds       = X1_base_in%f%evalDOF_x(xp,          0,dX1ds_s)
    dX1dthet    = X1_base_in%f%evalDOF_x(xp, DERIV_THET, X1_s  )
    d2X1dsdthet = X1_base_in%f%evalDOF_x(xp, DERIV_THET,dX1ds_s)
    dX1dzeta    = X1_base_in%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

    X2_int      = X2_base_in%f%evalDOF_x(xp,          0, X2_s  )
    dX2ds       = X2_base_in%f%evalDOF_x(xp,          0,dX2ds_s)
    dX2dthet    = X2_base_in%f%evalDOF_x(xp, DERIV_THET, X2_s  )
    d2X2dsdthet = X2_base_in%f%evalDOF_x(xp, DERIV_THET,dX2ds_s)
    dX2dzeta    = X2_base_in%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

    ! Get A components
    A_R_int            = out_base%f%evalDOF_x(xp,          0,   A_R_s)
    dA_Rds_int         = out_base%f%evalDOF_x(xp,          0, A_Rds_int)
    dA_Rdthet_int      = out_base%f%evalDOF_x(xp, DERIV_THET,   A_R_s)
    d2A_Rdsdthet_int   = out_base%f%evalDOF_x(xp, DERIV_THET, A_Rds_int)
    A_Z_int            = out_base%f%evalDOF_x(xp,          0,   A_Z_s)
    dA_Zds_int         = out_base%f%evalDOF_x(xp,          0, A_Zds_int)
    dA_Zdthet_int      = out_base%f%evalDOF_x(xp, DERIV_THET,   A_Z_s)
    d2A_Zdsdthet_int   = out_base%f%evalDOF_x(xp, DERIV_THET, A_Zds_int)
    A_phi_int          = out_base%f%evalDOF_x(xp,          0,   A_phi_s)
    dA_phids_int       = out_base%f%evalDOF_x(xp,          0, A_phids_int)
    dA_phidthet_int    = out_base%f%evalDOF_x(xp, DERIV_THET,   A_phi_s)
    d2A_phidsdthet_int = out_base%f%evalDOF_x(xp, DERIV_THET, A_phids_int)

    ! Get B components
    B_R_int            = out_base%f%evalDOF_x(xp,          0,   B_R_s)
    dB_Rds_int         = out_base%f%evalDOF_x(xp,          0, B_Rds_int)
    dB_Rdthet_int      = out_base%f%evalDOF_x(xp, DERIV_THET,   B_R_s)
    d2B_Rdsdthet_int   = out_base%f%evalDOF_x(xp, DERIV_THET, B_Rds_int)
    B_Z_int            = out_base%f%evalDOF_x(xp,          0,   B_Z_s)
    dB_Zds_int         = out_base%f%evalDOF_x(xp,          0, B_Zds_int)
    dB_Zdthet_int      = out_base%f%evalDOF_x(xp, DERIV_THET,   B_Z_s)
    d2B_Zdsdthet_int   = out_base%f%evalDOF_x(xp, DERIV_THET, B_Zds_int)
    B_phi_int          = out_base%f%evalDOF_x(xp,          0,   B_phi_s)
    dB_phids_int       = out_base%f%evalDOF_x(xp,          0, B_phids_int)
    dB_phidthet_int    = out_base%f%evalDOF_x(xp, DERIV_THET,   B_phi_s)
    d2B_phidsdthet_int = out_base%f%evalDOF_x(xp, DERIV_THET, B_phids_int)

    ! Get J components
    J_R_int            = out_base%f%evalDOF_x(xp,          0,   J_R_s)
    dJ_Rds_int         = out_base%f%evalDOF_x(xp,          0, J_Rds_int)
    dJ_Rdthet_int      = out_base%f%evalDOF_x(xp, DERIV_THET,   J_R_s)
    d2J_Rdsdthet_int   = out_base%f%evalDOF_x(xp, DERIV_THET, J_Rds_int)
    J_Z_int            = out_base%f%evalDOF_x(xp,          0,   J_Z_s)
    dJ_Zds_int         = out_base%f%evalDOF_x(xp,          0, J_Zds_int)
    dJ_Zdthet_int      = out_base%f%evalDOF_x(xp, DERIV_THET,   J_Z_s)
    d2J_Zdsdthet_int   = out_base%f%evalDOF_x(xp, DERIV_THET, J_Zds_int)
    J_phi_int          = out_base%f%evalDOF_x(xp,          0,   J_phi_s)
    dJ_phids_int       = out_base%f%evalDOF_x(xp,          0, J_phids_int)
    dJ_phidthet_int    = out_base%f%evalDOF_x(xp, DERIV_THET,   J_phi_s)
    d2J_phidsdthet_int = out_base%f%evalDOF_x(xp, DERIV_THET, J_phids_int)

    ! Get straight field line transformation
    IF(SFLcoord.EQ.0)THEN !GVEC coordinates (else=0)
      G_int    = LG_base_in%f%evalDOF_x(xp, 0, LG_s)
      dLAdthet = LG_base_in%f%evalDOF_x(xp, DERIV_THET, LG_s)
      dLAdzeta = LG_base_in%f%evalDOF_x(xp, DERIV_ZETA, LG_s)
    END IF

    ! --- Compare generated field representations of A and B with calculations from
    ! --- the original representation to ensure convergence
    ! Get the covariant basis vectors
    qvec     = (/ X1_int, X2_int, xp(2) /) !(X1,X2,zeta)
    e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   , 0.0    /)) !dxvec/ds
    e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet, 0.0    /)) !dxvec/dthet
    e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta, 1.0_wp /)) !dxvec/dzeta
    sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

    ! Get contravarian basis vectors
    grad_s    = CROSS(e_thet,e_zeta) /sqrtG
    grad_thet = CROSS(e_zeta,e_s   ) /sqrtG
    grad_zeta = CROSS(e_s   ,e_thet) /sqrtG

    ! Get Grad R and Grad Z pol - WARNING: this implementation only works for PEST coordinates
    grad_R = dX1ds * grad_s + dX1dthet * grad_thet + dX1dzeta * grad_zeta
    grad_Z = dX2ds * grad_s + dX2dthet * grad_thet + dX2dzeta * grad_zeta

    ! Calculate ave/max error in (R,Z, phi) magnetic field
    Bthet = ((iota_int-dLAdzeta )*dPhids_int)   !/sqrtG
    Bzeta = ((1.0_wp  +dLAdthet )*dPhids_int)   !/sqrtG
    Bcart(:) =  ( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG
    B_orig(1) =  Bcart(1) * grad_R(1)    + Bcart(2) * grad_R(2)    + Bcart(3) * grad_R(3)
    B_orig(2) =  Bcart(1) * grad_Z(1)    + Bcart(2) * grad_Z(2)    + Bcart(3) * grad_Z(3)
    B_orig(3) =  X1_int * (Bcart(1) * grad_zeta(1) + Bcart(2) * grad_zeta(2) + Bcart(3) * grad_zeta(3))
    BR_diff   = BR_diff   + ABS((B_orig(1) - B_R_int));   BR_diff_max   = MAX(BR_diff_max,   ABS(B_orig(1) - B_R_int))
    BZ_diff   = BZ_diff   + ABS((B_orig(2) - B_Z_int));   BZ_diff_max   = MAX(BZ_diff_max,   ABS(B_orig(2) - B_Z_int))
    Bphi_diff = Bphi_diff + ABS((B_orig(3) - B_phi_int)); Bphi_diff_max = MAX(Bphi_diff_max, ABS(B_orig(3) - B_phi_int))

    ! Calculate ave/max error in (R,Z,phi) vector potential
    Acart(:)  = (Phi_int * grad_thet(:) - (G_int * dPhids_int) * grad_s(:) - chi_int * grad_zeta(:))
    A_orig(1) =  Acart(1) * grad_R(1)    + Acart(2) * grad_R(2)    + Acart(3) * grad_R(3)
    A_orig(2) =  Acart(1) * grad_Z(1)    + Acart(2) * grad_Z(2)    + Acart(3) * grad_Z(3)
    A_orig(3) =  X1_int * (Acart(1) * grad_zeta(1) + Acart(2) * grad_zeta(2) + Acart(3) * grad_zeta(3))
    AR_diff   = AR_diff   + ABS((A_orig(1) - A_R_int));   AR_diff_max   = MAX(AR_diff_max,   ABS(A_orig(1) - A_R_int))
    AZ_diff   = AZ_diff   + ABS((A_orig(2) - A_Z_int));   AZ_diff_max   = MAX(AZ_diff_max,   ABS(A_orig(2) - A_Z_int))
    Aphi_diff = Aphi_diff + ABS((A_orig(3) - A_phi_int)); Aphi_diff_max = MAX(Aphi_diff_max, ABS(A_orig(3) - A_phi_int))

    !==========
    ! save data

    !!! SCALED DOMAIN: NEW LOGICAL DOMAIN REMAINS s_logical=[0,1],
    !!!                -->  position to evaluated was s = s_logical*s_max ,  d/ds_logical = ds/ds_logical * d/ds = s_max *d/ds

    !!! data_scalar3D(ithet,izeta,i_s, S__)          = spos
    data_scalar3D(ithet,izeta,i_s, S__)          = s_pos(i_s) !! =s_new,
    data_scalar3D(ithet,izeta,i_s, THET__)       = thet_pos(ithet)
    data_scalar3D(ithet,izeta,i_s, ZETA__)       = zeta_pos(izeta)
    data_scalar3D(ithet,izeta,i_s, X1__)         = X1_int
    data_scalar3D(ithet,izeta,i_s, X1_S__)       = dX1ds *s_max !! scale s-derivative with new domain size d/ds_new=d/ds*ds/ds_new
    data_scalar3D(ithet,izeta,i_s, X1_T__)       = dX1dthet
    data_scalar3D(ithet,izeta,i_s, X1_ST__)      = d2X1dsdthet *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, X2__)         = X2_int
    data_scalar3D(ithet,izeta,i_s, X2_S__)       = dX2ds *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, X2_T__)       = dX2dthet
    data_scalar3D(ithet,izeta,i_s, X2_ST__)      = d2X2dsdthet *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, P__)          = P_int
    data_scalar3D(ithet,izeta,i_s, P_S__)        = dPds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, A_R__)        = A_R_int
    data_scalar3D(ithet,izeta,i_s, A_R_S__)      = dA_Rds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, A_R_T__)      = dA_Rdthet_int
    data_scalar3D(ithet,izeta,i_s, A_R_ST__)     = d2A_Rdsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, A_Z__)        = A_Z_int
    data_scalar3D(ithet,izeta,i_s, A_Z_S__)      = dA_Zds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, A_Z_T__)      = dA_Zdthet_int
    data_scalar3D(ithet,izeta,i_s, A_Z_ST__)     = d2A_Zdsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, A_phi__)      = A_phi_int
    data_scalar3D(ithet,izeta,i_s, A_phi_S__)    = dA_phids_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, A_phi_T__)    = dA_phidthet_int
    data_scalar3D(ithet,izeta,i_s, A_phi_ST__)   = d2A_phidsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, B_R__)        = B_R_int
    data_scalar3D(ithet,izeta,i_s, B_R_S__)      = dB_Rds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, B_R_T__)      = dB_Rdthet_int
    data_scalar3D(ithet,izeta,i_s, B_R_ST__)     = d2B_Rdsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, B_Z__)        = B_Z_int
    data_scalar3D(ithet,izeta,i_s, B_Z_S__)      = dB_Zds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, B_Z_T__)      = dB_Zdthet_int
    data_scalar3D(ithet,izeta,i_s, B_Z_ST__)     = d2B_Zdsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, B_phi__)      = B_phi_int
    data_scalar3D(ithet,izeta,i_s, B_phi_S__)    = dB_phids_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, B_phi_T__)    = dB_phidthet_int
    data_scalar3D(ithet,izeta,i_s, B_phi_ST__)   = d2B_phidsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, J_R__)        = J_R_int
    data_scalar3D(ithet,izeta,i_s, J_R_S__)      = dJ_Rds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, J_R_T__)      = dJ_Rdthet_int
    data_scalar3D(ithet,izeta,i_s, J_R_ST__)     = d2J_Rdsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, J_Z__)        = J_Z_int
    data_scalar3D(ithet,izeta,i_s, J_Z_S__)      = dJ_Zds_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, J_Z_T__)      = dJ_Zdthet_int
    data_scalar3D(ithet,izeta,i_s, J_Z_ST__)     = d2J_Zdsdthet_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, J_phi__)      = J_phi_int
    data_scalar3D(ithet,izeta,i_s, J_phi_S__)    = dJ_phids_int *s_max !! scale s-derivative with new domain size
    data_scalar3D(ithet,izeta,i_s, J_phi_T__)    = dJ_phidthet_int
    data_scalar3D(ithet,izeta,i_s, J_phi_ST__)   = d2J_phidsdthet_int *s_max !! scale s-derivative with new domain size
    !==========

  END DO ; END DO !izeta,ithet
!$OMP END PARALLEL DO
  CALL ProgressBar(i_s,Ns_out)
END DO !i_s=1,Ns_out
BR_diff = BR_diff / REAL(Ns_out*Nthet_out*Nzeta_out,wp);BZ_diff = BZ_diff / REAL(Ns_out*Nthet_out*Nzeta_out,wp);Bphi_diff = Bphi_diff / REAL(Ns_out*Nthet_out*Nzeta_out,wp)
AR_diff = AR_diff / REAL(Ns_out*Nthet_out*Nzeta_out,wp);AZ_diff = AZ_diff / REAL(Ns_out*Nthet_out*Nzeta_out,wp);Aphi_diff = Aphi_diff / REAL(Ns_out*Nthet_out*Nzeta_out,wp)
SWRITE(UNIT_stdOut,'(A,6E16.7)')'AVE/MAX diff in B(R, Z, phi) :', BR_diff, BR_diff_max, BZ_diff, BZ_diff_max, Bphi_diff, Bphi_diff_max
SWRITE(UNIT_stdOut,'(A,6E16.7)')'MIN/MAX B(R, Z, phi)         :',MINVAL(data_scalar3D(:,:,:,B_R__)),MAXVAL(data_scalar3D(:,:,:,B_R__)),&
                                                                 MINVAL(data_scalar3D(:,:,:,B_Z__)),MAXVAL(data_scalar3D(:,:,:,B_Z__)),&
                                                                 MINVAL(data_scalar3D(:,:,:,B_phi__)),MAXVAL(data_scalar3D(:,:,:,B_phi__))
SWRITE(UNIT_stdOut,'(A,6E16.7)')'AVE/MAX diff in A(R, Z, phi) :', AR_diff,AR_diff_max,AZ_diff,AZ_diff_max,Aphi_diff,Aphi_diff_max
SWRITE(UNIT_stdOut,'(A,6E16.7)')'MIN/MAX A(R, Z, phi)         :',MINVAL(data_scalar3D(:,:,:,A_R__)),MAXVAL(data_scalar3D(:,:,:,A_R__)),&
                                                                 MINVAL(data_scalar3D(:,:,:,A_Z__)),MAXVAL(data_scalar3D(:,:,:,A_Z__)),&
                                                                 MINVAL(data_scalar3D(:,:,:,A_phi__)),MAXVAL(data_scalar3D(:,:,:,A_phi__))

! -----------------------------------------------------------------------------
! ---------------- CONVERT TO 1D TOROIDAL REPRESENTATION ----------------------
! -----------------------------------------------------------------------------
SWRITE(Unit_stdOut, *)  "Writing 3D variables to toroidal fourier representation..."
map_vars_3D_2D(R__     )      = X1__
map_vars_3D_2D(R_S__   )      = X1_S__
map_vars_3D_2D(R_T__   )      = X1_T__
map_vars_3D_2D(R_ST__  )      = X1_ST__
map_vars_3D_2D(Z__     )      = X2__
map_vars_3D_2D(Z_S__   )      = X2_S__
map_vars_3D_2D(Z_T__   )      = X2_T__
map_vars_3D_2D(Z_ST__  )      = X2_ST__
map_vars_3D_2D(P2D__   )      = P__
map_vars_3D_2D(P2D_S__ )      = P_S__
map_vars_3D_2D(P2D_T__ )      = -1 !leave zero
map_vars_3D_2D(P2D_ST__)      = -1 !leave zero
map_vars_3D_2D(A_R2D__ )      =A_R__
map_vars_3D_2D(A_R2D_S__)     =A_R_S__
map_vars_3D_2D(A_R2D_T__)     =A_R_T__
map_vars_3D_2D(A_R2D_ST__)    =A_R_ST__
map_vars_3D_2D(A_Z2D__)       =A_Z__
map_vars_3D_2D(A_Z2D_S__)     =A_Z_S__
map_vars_3D_2D(A_Z2D_T__)     =A_Z_T__
map_vars_3D_2D(A_Z2D_ST__)    =A_Z_ST__
map_vars_3D_2D(A_phi2D__)     =A_phi__
map_vars_3D_2D(A_phi2D_S__)   =A_phi_S__
map_vars_3D_2D(A_phi2D_T__)   =A_phi_T__
map_vars_3D_2D(A_phi2D_ST__)  =A_phi_ST__
map_vars_3D_2D(B_R2D__)       =B_R__
map_vars_3D_2D(B_R2D_S__)     =B_R_S__
map_vars_3D_2D(B_R2D_T__)     =B_R_T__
map_vars_3D_2D(B_R2D_ST__)    =B_R_ST__
map_vars_3D_2D(B_Z2D__)       =B_Z__
map_vars_3D_2D(B_Z2D_S__)     =B_Z_S__
map_vars_3D_2D(B_Z2D_T__)     =B_Z_T__
map_vars_3D_2D(B_Z2D_ST__)    =B_Z_ST__
map_vars_3D_2D(B_phi2D__)     =B_phi__
map_vars_3D_2D(B_phi2D_S__)   =B_phi_S__
map_vars_3D_2D(B_phi2D_T__)   =B_phi_T__
map_vars_3D_2D(B_phi2D_ST__)  =B_phi_ST__
map_vars_3D_2D(J_R2D__)       =J_R__
map_vars_3D_2D(J_R2D_S__)     =J_R_S__
map_vars_3D_2D(J_R2D_T__)     =J_R_T__
map_vars_3D_2D(J_R2D_ST__)    =J_R_ST__
map_vars_3D_2D(J_Z2D__)       =J_Z__
map_vars_3D_2D(J_Z2D_S__)     =J_Z_S__
map_vars_3D_2D(J_Z2D_T__)     =J_Z_T__
map_vars_3D_2D(J_Z2D_ST__)    =J_Z_ST__
map_vars_3D_2D(J_phi2D__)     =J_phi__
map_vars_3D_2D(J_phi2D_S__)   =J_phi_S__
map_vars_3D_2D(J_phi2D_T__)   =J_phi_T__
map_vars_3D_2D(J_phi2D_ST__)  =J_phi_ST__

data_scalar2D=0.0_wp
DO iVar=1,nVarScalar2D
  jVar=map_vars_3D_2D(iVar)
  IF(jVar.LE.0)CYCLE !do not set these
  DO ithet=1, Nthet_out
    DO i_s=1, Ns_out
      ! Store DOFs for writing to file
      data_scalar2D(ithet, i_s, 1:n_modes,iVar) = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,jVar))
    END DO ! ithet=1, Nthet_out
  END DO ! i_s=1, Ns_out
END DO !iVar=1,nVarScalar2D


SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE gvec_to_jorek_prepare

!===================================================================================================================================
!> compute different fields depending on the input parameters field_type and vector_component,
!!
!===================================================================================================================================
!SUBROUTINE Get_Field_base(mn_max,fac_nyq,field_type,vector_component,sgrid_in,out_base, field_out)
SUBROUTINE Get_Field(field_type,vector_component, field_out)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,CROSS,TWOPI,PI,ProgressBar
!USE MODgvec_LinAlg
!USE MODgvec_base   ,ONLY: t_base,base_new
!USE MODgvec_sGrid  ,ONLY: t_sgrid
!USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map
USE MODgvec_ReadState_vars  ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_ReadState_vars  ,ONLY: LA_r,X1_r,X2_r
USE MODgvec_ReadState_Vars  ,ONLY: profiles_1d,hmap_r,sbase_prof !for profiles
USE MODgvec_gvec_to_jorek_vars, ONLY: X1_fbase_nyq,X2_fbase_nyq,LA_fbase_nyq,out_base

IMPLICIT NONE

!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!  INTEGER      ,INTENT(IN) :: mn_max(2)                                     !< maximum number for new variables in SFL coordinates
!  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
                                                                            !< n_IP=fac_nyq*max(mn_max_in)
  INTEGER      ,INTENT(IN) :: field_type                                    !< field to be transformed (1=B, 2=A)
  INTEGER      ,INTENT(IN) :: vector_component                              !< vector component to be transformed (1=R, 2=Z, 3=phi)
!  CLASS(t_sgrid), INTENT(IN   ),TARGET :: sgrid_in                          !< change grid for G_base_out
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!  CLASS(t_Base),ALLOCATABLE,INTENT(INOUT) :: out_base                 !< new fourier basis of function Gthet,Gzeta
  REAL(wp),INTENT(INOUT) :: field_out(out_base%s%nBase,out_base%f%modes)  !< coefficients of toroidal vector potential
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,modes,i_mn,mn_IP                  !< enumerators
  INTEGER               :: BCtype_axis(0:4),BCaxis
  REAL(wp)              :: spos
  REAL(wp)              :: Phi_int,dPhids_int,iota_int,Chi_int
  REAL(wp)              :: dPhids_int_eps,iota_int_eps
  REAL(wp)              :: theta, zeta, sqrtG
  REAL(wp)              :: xp(2),qvec(3),e_s(3), e_thet(3),e_zeta(3)                                      !< position and covariant coordinate
  REAL(wp)              :: grad_s(3), grad_thet(3),grad_zeta(3),grad_R(3),grad_Z(3)                       !< contravariant coordinates and grad R/Z vectors
  REAL(wp)              :: Acart(3), Bcart(3), Bthet, Bzeta                                               !< cartesian vector potential, magnetic field and toroidal/poloidal components
  REAL(wp)              :: Jcart(3), B_ds(3), B_dthet(3), B_dzeta(3), grad_Bcart(3, 3)                    !< cartesion current density and gradient of magnetic field components

  REAL(wp)                            :: X1_s(  1:X1_base_r%f%modes),  X1_s_eps(  1:X1_base_r%f%modes)
  REAL(wp)                            :: dX1ds_s(1:X1_base_r%f%modes), dX1ds_s_eps(1:X1_base_r%f%modes)
  REAL(wp)                            :: X2_s(  1:X2_base_r%f%modes),  X2_s_eps(  1:X2_base_r%f%modes)
  REAL(wp)                            :: dX2ds_s(1:X2_base_r%f%modes), dX2ds_s_eps(1:X2_base_r%f%modes)
  REAL(wp)                            :: LA_s(   1:LA_base_r%f%modes), LA_s_eps(   1:LA_base_r%f%modes)

  REAL(wp),DIMENSION(1:out_base%f%mn_IP) :: dLAdthet_IP,dLAdzeta_IP, dLAdthet_IP_eps,dLAdzeta_IP_eps
  REAL(wp),DIMENSION(1:out_base%f%mn_IP) :: X1_IP,dX1ds_IP,dX1dthet_IP,dX1dzeta_IP,X1_IP_eps,dX1ds_IP_eps,dX1dthet_IP_eps,dX1dzeta_IP_eps
  REAL(wp),DIMENSION(1:out_base%f%mn_IP) :: X2_IP,dX2ds_IP,dX2dthet_IP,dX2dzeta_IP,X2_IP_eps,dX2ds_IP_eps,dX2dthet_IP_eps,dX2dzeta_IP_eps
  REAL(wp),DIMENSION(1:out_base%f%mn_IP) :: LA_IP, LA_IP_eps,field_out_IP
!  TYPE(t_fbase),ALLOCATABLE          :: X1_fbase_nyq
!  TYPE(t_fbase),ALLOCATABLE          :: X2_fbase_nyq
!  TYPE(t_fbase),ALLOCATABLE          :: LA_fbase_nyq

  ! Variables used in local interpolations for finite difference calculation of current density
  REAL(wp) :: X1_int,dX1ds,dX1dthet,dX1dzeta
  REAL(wp) :: X2_int,dX2ds,dX2dthet,dX2dzeta
  REAL(wp) :: dLAdthet,dLAdzeta

  REAL(wp) :: eps=1.0e-08                                   !< Small displacement for finite difference operations,
                                                            !  local variables appended with _eps are used in finite different operations
  REAL(wp) :: sgn
!===================================================================================================================================
  ! use maximum number of integration points from maximum mode number in both directions
  SWRITE(UNIT_StdOut,'(2(A,I4))')'GET FIELD,  field_type=',field_type,', vector_component=',vector_component
  mn_IP        = out_base%f%mn_IP  !total number of integration points
  modes        = out_base%f%modes  !number of modes in output
  nBase        = out_base%s%nBase  !number of radial points in output

  ! Loop over radial coordinate and evaluate modes of field_out
  DO is=1,nBase
    ! Avoid magnetic axis and plasma boundary
    spos=MIN(MAX(1.0e-08_wp,out_base%s%s_IP(is)),1.0_wp-1.0e-12_wp) !interpolation points for q_in

    ! Evaluate grid position, derivatives and field variables at integration points and finite difference eps points
    Phi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,1))
    dPhids_int  = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
    iota_int    = sbase_prof%evalDOF_s(spos,        0,profiles_1d(:,3))
    Chi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,2))

    !interpolate radially
    X1_s(:)        = X1_base_r%s%evalDOF2D_s(spos    ,X1_base_r%f%modes,      0,X1_r(:,:))
    dX1ds_s(:)     = X1_base_r%s%evalDOF2D_s(spos    ,X1_base_r%f%modes,DERIV_S,X1_r(:,:))
    X2_s(:)        = X2_base_r%s%evalDOF2D_s(spos    ,X2_base_r%f%modes,      0,X2_r(:,:))
    dX2ds_s(:)     = X2_base_r%s%evalDOF2D_s(spos    ,X2_base_r%f%modes,DERIV_S,X2_r(:,:))

    ! Interpolate finite difference point in radial direction - direction of finite step is changed for last element to stay inside the domain
    sgn = 1.0_wp
    IF (is .eq. nBase) sgn=-1.0_wp
    dPhids_int_eps = sbase_prof%evalDOF_s(spos+sgn*eps, DERIV_S ,profiles_1d(:,1))
    iota_int_eps   = sbase_prof%evalDOF_s(spos+sgn*eps,        0,profiles_1d(:,3))
    X1_s_eps(:)    = X1_base_r%s%evalDOF2D_s(spos+sgn*eps,X1_base_r%f%modes,      0,X1_r(:,:))
    dX1ds_s_eps(:) = X1_base_r%s%evalDOF2D_s(spos+sgn*eps,X1_base_r%f%modes,DERIV_S,X1_r(:,:))
    X2_s_eps(:)    = X2_base_r%s%evalDOF2D_s(spos+sgn*eps,X2_base_r%f%modes,      0,X2_r(:,:))
    dX2ds_s_eps(:) = X2_base_r%s%evalDOF2D_s(spos+sgn*eps,X2_base_r%f%modes,DERIV_S,X2_r(:,:))

    LA_s(:)        = LA_base_r%s%evalDOF2D_s(spos     ,LA_base_r%f%modes,         0,LA_r(:,:))
    LA_s_eps(:)    = LA_base_r%s%evalDOF2D_s(spos+sgn*eps,LA_base_r%f%modes,      0,LA_r(:,:))

    ! evaluate at integration points and finite difference eps from points
    X1_IP       = X1_fbase_nyq%evalDOF_IP(         0, X1_s(  :)); X1_IP_eps       = X1_fbase_nyq%evalDOF_IP(         0, X1_s_eps(  :))
    dX1ds_IP    = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s(:)); dX1ds_IP_eps    = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s_eps(:))
    dX1dthet_IP = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s(  :)); dX1dthet_IP_eps = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s_eps(  :))
    dX1dzeta_IP = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s(  :)); dX1dzeta_IP_eps = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s_eps(  :))

    X2_IP       = X2_fbase_nyq%evalDOF_IP(         0, X2_s(  :)); X2_IP_eps       = X2_fbase_nyq%evalDOF_IP(         0, X2_s_eps(  :))
    dX2ds_IP    = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s(:)); dX2ds_IP_eps    = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s_eps(:))
    dX2dthet_IP = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s(  :)); dX2dthet_IP_eps = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s_eps(  :))
    dX2dzeta_IP = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s(  :)); dX2dzeta_IP_eps = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s_eps(  :))

    LA_IP(:)       = LA_fbase_nyq%evalDOF_IP(         0,LA_s(:)); LA_IP_eps(:)       = LA_fbase_nyq%evalDOF_IP(         0,LA_s_eps(:))
    dLAdthet_IP(:) = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:)); dLAdthet_IP_eps(:) = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s_eps(:))
    dLAdzeta_IP(:) = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s(:)); dLAdzeta_IP_eps(:) = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s_eps(:))

    ! Loop over surface points and evaluate field_out
    DO i_mn=1,mn_IP
      theta = X1_fbase_nyq%x_IP(1, i_mn)
      zeta  = X1_fbase_nyq%x_IP(2, i_mn)

      ! Get the covariant basis vectors
      qvec     = (/ X1_IP(i_mn), X2_IP(i_mn), zeta /) !(X1,X2,zeta)
      e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds_IP(i_mn)   ,dX2ds_IP(i_mn)   , 0.0    /)) !dxvec/ds
      e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet_IP(i_mn),dX2dthet_IP(i_mn), 0.0    /)) !dxvec/dthet
      e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta_IP(i_mn),dX2dzeta_IP(i_mn), 1.0_wp /)) !dxvec/dzeta
      sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

      ! Get contravarian basis vectors
      grad_s    = CROSS(e_thet,e_zeta) /sqrtG
      grad_thet = CROSS(e_zeta,e_s   ) /sqrtG
      grad_zeta = CROSS(e_s   ,e_thet) /sqrtG

      ! Get Grad R and Grad Z pol - WARNING: this implementation only works for PEST coordinates
      grad_R = dX1ds_IP(i_mn) * grad_s + dX1dthet_IP(i_mn) * grad_thet + dX1dzeta_IP(i_mn) * grad_zeta
      grad_Z = dX2ds_IP(i_mn) * grad_s + dX2dthet_IP(i_mn) * grad_thet + dX2dzeta_IP(i_mn) * grad_zeta

      ! Get A and X in cartesian coordinates
      Acart(:)  = (Phi_int * grad_thet(:) - (LA_IP(i_mn) * dPhids_int) * grad_s(:) - Chi_int * grad_zeta)

      ! Calculate cartesian magnetic field
      Bthet = (iota_int-dLAdzeta_IP(i_mn) )*dPhids_int   !/sqrtG
      Bzeta = (1.0_wp  +dLAdthet_IP(i_mn) )*dPhids_int   !/sqrtG
      Bcart(:) =  ( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG

      SELECT CASE(field_type)
      CASE(1)  ! Magnetic Field
        SELECT CASE(vector_component)
        CASE(1)
          ! Get Vertical Magnetic Field - B_X
          field_out_IP(i_mn) = Bcart(1)
        CASE(2)
          ! Get Vertical Magnetic Field - B_Y
          field_out_IP(i_mn) = Bcart(2)
        CASE(3)
          ! Get Vertical Magnetic Field - B_Z
          field_out_IP(i_mn) = Bcart(3)
        CASE(4)
          ! Get Radial Magnetic Field - B_R
          field_out_IP(i_mn) = SUM(Bcart(1:3) * grad_R(1:3))
        CASE(5)
          ! Get Toroidal Magnetic Field - B_phi = R * (B.grad(zeta))
          field_out_IP(i_mn) = X1_IP(i_mn)*SUM(Bcart(1:3) * grad_zeta(1:3))
        CASE DEFAULT
          SWRITE(UNIT_StdOut,*) "Invalid vector component selected: ", vector_component
        END SELECT
      CASE(2)  ! Vector Potential
        SELECT CASE(vector_component)
        CASE(1)
          ! Get Vertical Flux - A_X
          field_out_IP(i_mn) = Acart(1)
        CASE(2)
          ! Get Vertical Flux - A_Y
          field_out_IP(i_mn) = Acart(2)
        CASE(3)
          ! Get Vertical Flux - A_Z
          field_out_IP(i_mn) = Acart(3)
        CASE(4)
          ! Get Radial Flux - A_R
          field_out_IP(i_mn) = SUM(Acart(1:3) * grad_R(1:3))
        CASE(5)
          ! Get Poloidal Flux - A_phi = R * (A.grad(zeta))
          field_out_IP(i_mn) = X1_IP(i_mn)*SUM(Acart(1:3) * grad_zeta(1:3))
        CASE DEFAULT
          SWRITE(UNIT_StdOut,*) "Invalid vector component selected: ", vector_component
        END SELECT
      CASE(3)  ! Current density
        ! Calculate ds derivative of B
        qvec     = (/ X1_IP_eps(i_mn), X2_IP_eps(i_mn), zeta /) !(X1,X2,zeta)
        e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds_IP_eps(i_mn)   ,dX2ds_IP_eps(i_mn)   , 0.0_wp /)) !dxvec/ds
        e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet_IP_eps(i_mn),dX2dthet_IP_eps(i_mn), 0.0_wp /)) !dxvec/dthet
        e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta_IP_eps(i_mn),dX2dzeta_IP_eps(i_mn), 1.0_wp /)) !dxvec/dzeta
        sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

        Bthet = (iota_int_eps-dLAdzeta_IP_eps(i_mn) )*dPhids_int_eps   !/sqrtG
        Bzeta = (1.0_wp  +dLAdthet_IP_eps(i_mn) )*dPhids_int_eps   !/sqrtG
        B_ds(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / (sgn*eps)      ! calculating dBx_ds, dBy_ds, dBz_ds

        ! Calculate dtheta derivative of B
        xp = (/theta+eps, zeta/)
        X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
        dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
        dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
        dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

        X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
        dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
        dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
        dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

        dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
        dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

        qvec   = (/X1_int,X2_int,zeta/)
        e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   , 0.0_wp /)) !dxvec/ds
        e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet, 0.0_wp /)) !dxvec/dthet
        e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta, 1.0_wp /)) !dxvec/dzeta
        sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

        Bthet = (iota_int-dLAdzeta )*dPhids_int   !/sqrtG
        Bzeta = (1.0_wp  +dLAdthet )*dPhids_int   !/sqrtG
        B_dthet(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / eps      ! calculating dBx_dtheta, dBy_dtheta, dBz_dtheta

        ! Calculate dzeta derivative of B
        xp = (/theta, zeta+eps/)
        X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
        dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
        dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
        dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

        X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
        dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
        dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
        dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

        dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
        dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

        qvec   = (/X1_int,X2_int,xp(2)/) !xp(2)=zeta +eps  here!
        e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   , dX2ds   , 0.0_wp /)) !dxvec/ds
        e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet, dX2dthet, 0.0_wp /)) !dxvec/dthet
        e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta, dX2dzeta, 1.0_wp /)) !dxvec/dzeta
        sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

        Bthet = (iota_int-dLAdzeta )*dPhids_int   !/sqrtG
        Bzeta = (1.0_wp  +dLAdthet )*dPhids_int   !/sqrtG
        B_dzeta(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / eps      ! calculating dBx_dzeta, dBy_dzeta, dBz_dzeta

        ! Calculate B derivatives by finite difference
        grad_Bcart(1, :) = B_ds(1) * grad_s(:) + B_dthet(1) * grad_thet(:) + B_dzeta(1) * grad_zeta(:)   ! grad_BX
        grad_Bcart(2, :) = B_ds(2) * grad_s(:) + B_dthet(2) * grad_thet(:) + B_dzeta(2) * grad_zeta(:)   ! grad_BY
        grad_Bcart(3, :) = B_ds(3) * grad_s(:) + B_dthet(3) * grad_thet(:) + B_dzeta(3) * grad_zeta(:)   ! grad_BZ

        ! Calculate current cartesian components
        Jcart(1) = grad_Bcart(3, 2) - grad_Bcart(2, 3)   ! dBZ_dY - dBY_dZ
        Jcart(2) = grad_Bcart(1, 3) - grad_Bcart(3, 1)   ! dBX_dZ - dBZ_dX
        Jcart(3) = grad_Bcart(2, 1) - grad_Bcart(1, 2)   ! dBY_dX - dBX_dY

        SELECT CASE(vector_component)
        CASE(1)
          ! Get Vertical Current Density - J_X
          field_out_IP(i_mn) = Jcart(1)
        CASE(2)
          ! Get Vertical Current Density - J_Y
          field_out_IP(i_mn) = Jcart(2)
        CASE(3)
          ! Get Vertical Current Density - J_Z
          field_out_IP(i_mn) = Jcart(3)
        CASE(4)
          ! Get Radial Current Density - J_R
          field_out_IP(i_mn) = SUM(Jcart(1:3) * grad_R(1:3))
        CASE(5)
          ! Get Toroidal Current Density - J_phi = R * (J.grad(zeta))
          field_out_IP(i_mn) = X1_IP(i_mn)*SUM(Jcart(1:3) * grad_zeta(1:3))
        CASE DEFAULT
          SWRITE(UNIT_StdOut,*) "Invalid vector component selected: ", vector_component
        END SELECT

      CASE DEFAULT
        SWRITE(UNIT_StdOut,*) "Invalid field type selected: ", field_type
      END SELECT

    END DO ! i_mn

    ! Convert interation points into fourier mode representation
    field_out(is,:) = out_base%f%initDOF(field_out_IP(:))
  END DO ! is

  ! Convert radial fourier representation into radial spline

  ! SETTING boundary conditions at the axis (standard low order BC).
  !    Note: B_R and B_Z on the axis should be zero for mode m=n=0, but this should already be in the data and is not imposed here
  BCtype_axis(MN_ZERO    )= BC_TYPE_NEUMANN   ! derivative zero at axis
  BCtype_axis(M_ZERO     )= BC_TYPE_NEUMANN   ! derivative zero at axis
  BCtype_axis(M_ODD_FIRST)= BC_TYPE_DIRICHLET !=0 at axis m>0 modes should not contribute
  BCtype_axis(M_ODD      )= BC_TYPE_DIRICHLET !=0 at axis
  BCtype_axis(M_EVEN     )= BC_TYPE_DIRICHLET !=0 at axis
  ! for smooth axis BC use instead:
  ! BCtype_axis=0
  DO iMode=1, modes
    field_out(:, iMode) = out_base%s%initDOF(field_out(:, iMode))
    BCaxis=BCtype_axis(out_base%f%zero_odd_even(iMode))
    IF(BCaxis.EQ.0)THEN !AUTOMATIC, m-dependent BC, for m>deg, switch off all DOF up to deg+1
      BCaxis=-1*MIN(out_base%s%deg+1,out_base%f%Xmn(1,iMode))
    END IF
    CALL out_base%s%applyBCtoDOF(field_out(:,iMode), &
                                 (/BCaxis,BC_TYPE_OPEN/),(/0.0_wp,0.0_wp/))
  END DO

END SUBROUTINE Get_Field

!===================================================================================================================================
!> write data to file
!!
!===================================================================================================================================
SUBROUTINE gvec_to_jorek_writeToFile_ASCII()
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MODgvec_gvec_to_jorek_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: ioUnit,iVar
INTEGER       :: date_time_values(8)
CHARACTER(LEN=30) :: curr_date_time
!===================================================================================================================================
  CALL DATE_AND_TIME(VALUES=date_time_values)
  WRITE(curr_date_time,'(1X,I4.4,"-",I2.2,"-",I2.2,2X,I2.2,":",I2.2,":",I2.2)') date_time_values(1:3),date_time_values(5:7)

  WRITE(UNIT_stdOut,'(A)')'WRITING NEW JOREK FILE    "'//TRIM(FileNameOut)//'" , date: '//TRIM(curr_date_time)//' ... '
  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit       ,&
     FILE     = TRIM(FileNameOut) ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' )

!HEADER
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## GVEC-TO-JOREK file, VERSION: 1.0                                                                 '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## data is written on equidistant points in s,theta,zeta coordinates,                               '
  WRITE(ioUnit,'(A100)')'## * radially outward coordinate s=sqrt(phi_tor/phi_tor_edge) in [0,1]                              '
  WRITE(ioUnit,'(A100)')'##   s(1:Ns) , with  s(1)=0, s(Ns)=1                                                                '
  WRITE(ioUnit,'(A100)')'## * poloidal angle theta in [0,2pi] , sign: theta ~ atan(z/sqrt(x^2+y^2))                          '
  WRITE(ioUnit,'(A100)')'##   theta(1:Ntheta)  with theta(1)=0, theta(Ntheta)=2pi*(Ntheta-1)*/Ntheta                         '
  WRITE(ioUnit,'(A100)')'## * toroidal angle zeta in [0,2pi/nfp], sign: zeta ~ atan(y/x)  (opposite to GVEC definition!)     '
  WRITE(ioUnit,'(A100)')'##   zeta(1:Nzeta)  with zeta(1)=0, zeta(Nzeta)=2pi/nfp*(Nzeta-1)*/Nzeta                            '
  WRITE(ioUnit,'(A100)')'## * Angular coordinates can represent GVEC coordinates, which are not SFL (straight-field line)    '
  WRITE(ioUnit,'(A100)')'##   coordinates or can be SFL coordinates, either PEST or BOOZER. See global parameter "SFLcoord"  '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## 2D arrays containing toroidal Fourier coefficients for GVEC fields are used for the import       '
  WRITE(ioUnit,'(A100)')'## 3D test data can be generated instead if the -g option has been used                                                                                                 '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'##   WARNING: Note that the change in the coordinate system is important!                           '
  WRITE(ioUnit,'(A100)')'##            GVEC and JOREK both use left handed coodinate systems so:                             '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'##                          (x,y,z)=(Rcos(zeta),-Rsin(zeta),Z)                                      '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## Global variables:                                                                                '
  WRITE(ioUnit,'(A100)')'## * SFLcoord    : =0: GVEC coords (not SFL), =1: PEST SFL coords. , =2: BOOZER SFL coords.         '
  WRITE(ioUnit,'(A100)')'## * nfp         : number of toroidal field periods (toroidal angle [0,2pi/nfp])                    '
  WRITE(ioUnit,'(A100)')'## * asym        :  =0: symmetric cofiguration (R~cos,Z~sin), 1: asymmetric                         '
  WRITE(ioUnit,'(A100)')'## * m_max       : maximum number of poloidal modes in R,Z,lambda variables                         '
  WRITE(ioUnit,'(A100)')'## * n_max       : maximum number of toroidal modes in R,Z,lambda variables                         '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## 2D arrays of scalar fourier coefficients (1:Ntheta,1:Ns)                                         '
  WRITE(ioUnit,'(A100)')'## * R              : major radius                                                                  '
  WRITE(ioUnit,'(A100)')'## * R_s            : radial derivative of major radius                                             '
  WRITE(ioUnit,'(A100)')'## * R_t            : poloidal derivative of major radius                                           '
  WRITE(ioUnit,'(A100)')'## * R_st           : cross derivative of major radius                                              '
  WRITE(ioUnit,'(A100)')'## * Z              : vertical position                                                             '
  WRITE(ioUnit,'(A100)')'## * Z_s            : radial derivative of vertical position                                        '
  WRITE(ioUnit,'(A100)')'## * Z_t            : poloidal derivative of vertical position                                      '
  WRITE(ioUnit,'(A100)')'## * Z_st           : cross derivative of vertical position                                         '
  WRITE(ioUnit,'(A100)')'## * P              : pressure                                                                      '
  WRITE(ioUnit,'(A100)')'## * P_s            : radial derivative of pressure                                                 '
  WRITE(ioUnit,'(A100)')'## * P_t            : poloidal derivative of pressure                                               '
  WRITE(ioUnit,'(A100)')'## * P_st           : cross derivative of pressure                                                  '
  WRITE(ioUnit,'(A100)')'## * A_R            : X component of  vector potential                                              '
  WRITE(ioUnit,'(A100)')'## * A_R_s          : radial derivative of X component of  vector potential                         '
  WRITE(ioUnit,'(A100)')'## * A_R_t          : poloidal derivative of X component of  vector potential                       '
  WRITE(ioUnit,'(A100)')'## * A_R_st         : cross derivative of X component of  vector potential                          '
  WRITE(ioUnit,'(A100)')'## * A_Z            : Y Component vector potential                                                  '
  WRITE(ioUnit,'(A100)')'## * A_Z_s          : radial derivative of Y component of  vector potential                         '
  WRITE(ioUnit,'(A100)')'## * A_Z_t          : poloidal derivative of Y component of  vector potential                       '
  WRITE(ioUnit,'(A100)')'## * A_Z_st         : cross derivative of Y component of  vector potential                          '
  WRITE(ioUnit,'(A100)')'## * A_phi          : Vertical vector potential                                                     '
  WRITE(ioUnit,'(A100)')'## * A_phi_s        : radial derivative of Vertical vector potential                                '
  WRITE(ioUnit,'(A100)')'## * A_phi_t        : poloidal derivative of Vertical vector potential                              '
  WRITE(ioUnit,'(A100)')'## * A_phi_st       : cross derivative of Vertical vector potential                                 '
  WRITE(ioUnit,'(A100)')'## * B_R            : X component of  magnetic field                                                '
  WRITE(ioUnit,'(A100)')'## * B_R_s          : radial derivative of X component of  magnetic field                           '
  WRITE(ioUnit,'(A100)')'## * B_R_t          : poloidal derivative of X component of  magnetic field                         '
  WRITE(ioUnit,'(A100)')'## * B_R_st         : cross derivative of X component of  magnetic field                            '
  WRITE(ioUnit,'(A100)')'## * B_Z            : Y Component magnetic field                                                    '
  WRITE(ioUnit,'(A100)')'## * B_Z_s          : radial derivative of Y component of  magnetic field                           '
  WRITE(ioUnit,'(A100)')'## * B_Z_t          : poloidal derivative of Y component of  magnetic field                         '
  WRITE(ioUnit,'(A100)')'## * B_Z_st         : cross derivative of Y component of  magnetic field                            '
  WRITE(ioUnit,'(A100)')'## * B_phi          : Vertical magnetic field                                                       '
  WRITE(ioUnit,'(A100)')'## * B_phi_s        : radial derivative of Vertical magnetic field                                  '
  WRITE(ioUnit,'(A100)')'## * B_phi_t        : poloidal derivative of Vertical magnetic field                                '
  WRITE(ioUnit,'(A100)')'## * B_phi_st       : cross derivative of Vertical magnetic field                                   '
  WRITE(ioUnit,'(A100)')'## * J_R            : X component of  current density                                               '
  WRITE(ioUnit,'(A100)')'## * J_R_s          : radial derivative of X component of  current density                          '
  WRITE(ioUnit,'(A100)')'## * J_R_t          : poloidal derivative of X component of  current density                        '
  WRITE(ioUnit,'(A100)')'## * J_R_st         : cross derivative of X component of  current density                           '
  WRITE(ioUnit,'(A100)')'## * J_Z            : Y Component current density                                                   '
  WRITE(ioUnit,'(A100)')'## * J_Z_s          : radial derivative of Y component of  current density                          '
  WRITE(ioUnit,'(A100)')'## * J_Z_t          : poloidal derivative of Y component of  current density                        '
  WRITE(ioUnit,'(A100)')'## * J_Z_st         : cross derivative of Y component of  current density                           '
  WRITE(ioUnit,'(A100)')'## * J_phi          : Vertical current density                                                      '
  WRITE(ioUnit,'(A100)')'## * J_phi_s        : radial derivative of Vertical current density                                 '
  WRITE(ioUnit,'(A100)')'## * J_phi_t        : poloidal derivative of Vertical current density                               '
  WRITE(ioUnit,'(A100)')'## * J_phi_st       : cross derivative of Vertical current density                                  '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## 3D arrays of scalars (1:Ntheta,1:Nzeta,1:Ns)                                                     '
  WRITE(ioUnit,'(A100)')'## * s              : radial coordinate                                                             '
  WRITE(ioUnit,'(A100)')'## * t              : poloidal coordinate                                                           '
  WRITE(ioUnit,'(A100)')'## * p              : toroidal coordinate                                                           '
  WRITE(ioUnit,'(A100)')'## * X1 (R)         : coordinate R=sqrt(x^2+y^2) ( called X1 in GVEC, only=R for hmap=1)            '
  WRITE(ioUnit,'(A100)')'## * X1_s (R)       : radial derivative of X1                                                       '
  WRITE(ioUnit,'(A100)')'## * X1_t (R)       : poloidal derivative of X1                                                     '
  WRITE(ioUnit,'(A100)')'## * X1_st (R)      : cross derivative of X1                                                        '
  WRITE(ioUnit,'(A100)')'## * X2 (Z)         : coordinate Z=z ( called X2 in GVEC, only=Z for hmap=1)                        '
  WRITE(ioUnit,'(A100)')'## * X2_s (Z)       : radial derivative of X2                                                       '
  WRITE(ioUnit,'(A100)')'## * X2_t (Z)       : poloidal derivative of X2                                                     '
  WRITE(ioUnit,'(A100)')'## * X2_st (Z)      : cross derivative of X2                                                        '
  WRITE(ioUnit,'(A100)')'## * P              : pressure                                                                      '
  WRITE(ioUnit,'(A100)')'## * P_s            : radial derivative of pressure                                                 '
  WRITE(ioUnit,'(A100)')'## * A_R            : X component of  vector potential                                              '
  WRITE(ioUnit,'(A100)')'## * A_R_s          : radial derivative of X component of  vector potential                         '
  WRITE(ioUnit,'(A100)')'## * A_R_t          : poloidal derivative of X component of  vector potential                       '
  WRITE(ioUnit,'(A100)')'## * A_R_st         : cross derivative of X component of  vector potential                          '
  WRITE(ioUnit,'(A100)')'## * A_Z            : Y Component vector potential                                                  '
  WRITE(ioUnit,'(A100)')'## * A_Z_s          : radial derivative of Y component of  vector potential                         '
  WRITE(ioUnit,'(A100)')'## * A_Z_t          : poloidal derivative of Y component of  vector potential                       '
  WRITE(ioUnit,'(A100)')'## * A_Z_st         : cross derivative of Y component of  vector potential                          '
  WRITE(ioUnit,'(A100)')'## * A_phi          : Vertical vector potential                                                     '
  WRITE(ioUnit,'(A100)')'## * A_phi_s        : radial derivative of Vertical vector potential                                '
  WRITE(ioUnit,'(A100)')'## * A_phi_t        : poloidal derivative of Vertical vector potential                              '
  WRITE(ioUnit,'(A100)')'## * A_phi_st       : cross derivative of Vertical vector potential                                 '
  WRITE(ioUnit,'(A100)')'## * B_R            : X component of  magnetic field                                                '
  WRITE(ioUnit,'(A100)')'## * B_R_s          : radial derivative of X component of  magnetic field                           '
  WRITE(ioUnit,'(A100)')'## * B_R_t          : poloidal derivative of X component of  magnetic field                         '
  WRITE(ioUnit,'(A100)')'## * B_R_st         : cross derivative of X component of  magnetic field                            '
  WRITE(ioUnit,'(A100)')'## * B_Z            : Y Component magnetic field                                                    '
  WRITE(ioUnit,'(A100)')'## * B_Z_s          : radial derivative of Y component of  magnetic field                           '
  WRITE(ioUnit,'(A100)')'## * B_Z_t          : poloidal derivative of Y component of  magnetic field                         '
  WRITE(ioUnit,'(A100)')'## * B_Z_st         : cross derivative of Y component of  magnetic field                            '
  WRITE(ioUnit,'(A100)')'## * B_phi          : Vertical magnetic field                                                       '
  WRITE(ioUnit,'(A100)')'## * B_phi_s        : radial derivative of Vertical magnetic field                                  '
  WRITE(ioUnit,'(A100)')'## * B_phi_t        : poloidal derivative of Vertical magnetic field                                '
  WRITE(ioUnit,'(A100)')'## * B_phi_st       : cross derivative of Vertical magnetic field                                   '
  WRITE(ioUnit,'(A100)')'## * J_R            : X component of  current density                                               '
  WRITE(ioUnit,'(A100)')'## * J_R_s          : radial derivative of X component of  current density                          '
  WRITE(ioUnit,'(A100)')'## * J_R_t          : poloidal derivative of X component of  current density                        '
  WRITE(ioUnit,'(A100)')'## * J_R_st         : cross derivative of X component of  current density                           '
  WRITE(ioUnit,'(A100)')'## * J_Z            : Y Component current density                                                   '
  WRITE(ioUnit,'(A100)')'## * J_Z_s          : radial derivative of Y component of  current density                          '
  WRITE(ioUnit,'(A100)')'## * J_Z_t          : poloidal derivative of Y component of  current density                        '
  WRITE(ioUnit,'(A100)')'## * J_Z_st         : cross derivative of Y component of  current density                           '
  WRITE(ioUnit,'(A100)')'## * J_phi          : Vertical current density                                                      '
  WRITE(ioUnit,'(A100)')'## * J_phi_s        : radial derivative of Vertical current density                                 '
  WRITE(ioUnit,'(A100)')'## * J_phi_t        : poloidal derivative of Vertical current density                               '
  WRITE(ioUnit,'(A100)')'## * J_phi_st       : cross derivative of Vertical current density                                  '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(2A)')  '## CALLED AS: ',TRIM(cmdline)
  WRITE(ioUnit,'(2A)')  '## CALLED ON: ',TRIM(curr_date_time)
  WRITE(ioUnit,'(A100)')'####################################################################################################'
  WRITE(ioUnit,'(A)')'##<< number of grid points: 1:Ns (radial), 1:Ntheta (poloidal),1:Nzeta (toroidal) '
  WRITE(ioUnit,'(*(I8,:,1X))')Ns_out,Nthet_out,Nzeta_out
  WRITE(ioUnit,'(A)')'##<< global: SFLcoord,nfp,  asym, m_max, n_max, n_modes, sin_min, sin_max, cos_min, cos_max'
  WRITE(ioUnit,'(12X,*(I6,:,1X))')SFLcoord,nfp_out,asym_out,mn_max_out(1:2), n_modes,sin_range(1:2),cos_range(1:2)
  IF (generate_test_data) THEN
    ! Write 3D data only
    DO iVar=1,nVarScalar3D
      WRITE(ioUnit,'(A)',ADVANCE='NO')'##<< 3D scalar variable (1:Ntheta,1:Nzeta,1:Ns), Variable name: '
      WRITE(ioUNIT,'(A)')' "'//TRIM(StrVarNamesScalar3D(iVar))//'"'
      WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar3D(1:Nthet_out,1:Nzeta_out,1:Ns_out,iVar)
    END DO !iVar=1,nVarScalar3D
  ELSE
    ! Write 2D data only
    DO iVar=1,nVarScalar2D
      WRITE(ioUnit,'(A)',ADVANCE='NO')'##<< 2D scalar variable fourier modes (1:Ntheta,1:Ns), Variable name: '
      WRITE(ioUNIT,'(A)')' "'//TRIM(StrVarNamesScalar2D(iVar))//'"'
      WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar2D(1:Nthet_out,1:Ns_out,1:n_modes,iVar)
    END DO !iVar=1,nVarScalar2D
  END IF

  CLOSE(ioUnit)

  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE gvec_to_jorek_writeToFile_ASCII


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE finalize_gvec_to_jorek
! MODULES
USE MODgvec_gvec_to_jorek_Vars
USE MODgvec_readState, ONLY: finalize_readState
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  CALL Finalize_ReadState()
  SDEALLOCATE(s_pos)
  SDEALLOCATE(thet_pos)
  SDEALLOCATE(zeta_pos)
  !SDEALLOCATE(data_1D)
  SDEALLOCATE(data_scalar3D)
  SDEALLOCATE(data_scalar2D)
  !SDEALLOCATE(data_vector3D)

  CALL out_base%free()
  DEALLOCATE(out_base)

  CALL X1_fbase_nyq%free()
  CALL X2_fbase_nyq%free()
  CALL LA_fbase_nyq%free()
  CALL fbase_zeta%free()
  DEALLOCATE(X1_fbase_nyq)
  DEALLOCATE(X2_fbase_nyq)
  DEALLOCATE(LA_fbase_nyq)
  DEALLOCATE(fbase_zeta)


END SUBROUTINE finalize_gvec_to_jorek

END MODULE MODgvec_gvec_to_jorek
