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
USE MODgvec_gvec_to_jorek_Vars, ONLY: gvecfileName,FileNameOut, Ns_out,npfactor,factorSFL,Nthet_out,Nzeta_out,SFLcoord,cmdline, generate_test_data
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
       'Number of angular points, computed from max. mode numbers (=npfactor*mn_max) [DEFAULT = 4]',cla_int,'4') 
  CALL cla_register('-p',  '--polpoints', &
       'Number of poloidal points, if specified overwrites factor*m_max  [OPTIONAL]',cla_int,'-1') 
  CALL cla_register('-t',  '--torpoints', &
       'Number of toroidal points, if specified overwrites factor*n_max  [OPTIONAL]',cla_int,'-1') 
  CALL cla_register('-s',  '--sflcoord', &
       'which angular coordinates to choose: =0: GVEC coord. (no SFL), =1: PEST SFL, =2: BOOZER SFL [DEFAULT = 0]', cla_int,'0')
  CALL cla_register('-f',  '--factorsfl', &
       'factor on maximum mode numbers (mn_max=factorsfl*mn_max_in) for SFL coords.  [DEFAULT = 4]', cla_int,'4')
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
  CALL cla_get('-t',Nzeta_out)
  CALL cla_get('-s',SFLcoord)
  CALL cla_get('-f',factorSFL)
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
  IF(commandFailed) STOP

  SWRITE(UNIT_stdOut,'(A)')     ' INPUT PARAMETERS:' 
  SWRITE(UNIT_stdOut,'(A,I6)')  '  * Number of radial points        : ',Ns_out
  SWRITE(UNIT_stdOut,'(A,I4)')  '  * npfactor points from modes     : ',npfactor
  IF(Nthet_out.NE.-1) THEN
    SWRITE(UNIT_stdOut,'(A,I4)')'  * number of points in theta      : ',Nthet_out
  END IF
  IF(Nzeta_out.NE.-1) THEN
    SWRITE(UNIT_stdOut,'(A,I4)')'  * number of points in zeta       : ',Nzeta_out
  END IF
  SWRITE(UNIT_stdOut,'(A,I4,A)')'  * SFL coordinates flag           : ',SFLcoord,' ( '//SFLcoordName(SFLcoord)//' )'
  SWRITE(UNIT_stdOut,'(A,I4)')  '  * factor for modes of SFL coords : ',factorSFL
  SWRITE(UNIT_stdOut,'(A,A)')   '  * generate test data for JOREK   : ',MERGE(".false.",".true. ",generate_test_data)
  SWRITE(UNIT_stdOut,'(A,A)')   '  * GVEC input file                : ',TRIM(gvecfileName)
  SWRITE(UNIT_stdOut,'(A,A)')   '  * output file name               : ',TRIM(FileNameOut)
  SWRITE(UNIT_stdOut,fmt_sep)

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
USE MODgvec_transform_sfl     ,ONLY: BuildTransform_SFL
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

  ! Initialise grid variables from GVEC restart file
  CALL ReadState(TRIM(gvecfileName))

  mn_max_out(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
  mn_max_out(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))
  nfp_out          = X1_base_r%f%nfp

  IF((X1_base_r%f%sin_cos.EQ._COS_).AND.(X2_base_r%f%sin_cos.EQ._SIN_).AND.(LA_base_r%f%sin_cos.EQ._SIN_))THEN
    asym_out = 0 !R~cos,Z~sin,lambda~sin
  ELSE
    asym_out = 1 !full fourier
  END IF

  IF(SFLcoord.NE.0)THEN
   SWRITE(UNIT_StdOut,*) "Only GVEC coordinates are currently handled!"
  END IF
  
  ! Initialise sample points in s, theta, zeta. s and theta can be randomly sampled for testing purposes.
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
  IF(Nzeta_out.EQ.-1) THEN !overwrite with default from factorFourier
    Nzeta_out = MAX(1,npfactor*mn_max_out(2)) !if n=0, output 1 point
  END IF
  !check
  IF(Nthet_out .LT. 4*mn_max_out(1)) WRITE(UNIT_StdOut,'(A)')'WARNING: number of poloidal points for output should be >=4*m_max!'
  IF(Nzeta_out .LT. 4*mn_max_out(2)) WRITE(UNIT_StdOut,'(A)')'WARNING: number of toroidal points for output should be >=4*n_max!'

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

  n_modes = X1_base_r%f%modes
  ALLOCATE(data_scalar2D(Nthet_out, Ns_out, n_modes*nVarScalar2D))
  ALLOCATE(data_scalar3D(  Nthet_out,Nzeta_out,Ns_out,nVarscalar3D))
  !ALLOCATE(data_vector3D(3,Nthet_out,Nzeta_out,Ns_out,nVarvector3D))

  SWRITE(UNIT_stdOut,'(A,3I6)')'  Number OF N_s,N_theta,N_zeta evaluation points:',Ns_out,Nthet_out,Nzeta_out
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)

  SELECT CASE(SFLcoord)
  CASE(0) ! PEST - toroidal coordinate is the cylindrical toroidal direction
    CALL gvec_to_jorek_prepare(X1_base_r,X1_r,X2_base_r,X2_r,LA_base_r,LA_r)
  CASE DEFAULT
    SWRITE(UNIT_StdOut,*)'This SFLcoord is not valid',SFLcoord
    STOP
  END SELECT
END SUBROUTINE init_gvec_to_jorek


!===================================================================================================================================
!> prepare all data to be written
!!
!===================================================================================================================================
SUBROUTINE gvec_to_jorek_prepare(X1_base_in,X1_in,X2_base_in,X2_in,LG_base_in,LG_in)
! MODULES
USE MODgvec_gvec_to_jorek_Vars 
USE MODgvec_Globals,        ONLY: CROSS,TWOPI,ProgressBar
USE MODgvec_ReadState_Vars, ONLY: profiles_1d,hmap_r,sbase_prof !for profiles
USE MODgvec_Base,           ONLY: t_base, t_fbase, fbase_new
USE MODgvec_get_field
USE MODgvec_sGrid  ,        ONLY: t_sgrid


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
INTEGER                                  :: i_s,ithet,izeta                                                  ! Indices for enumeration
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
REAL(wp),DIMENSION(3)                    :: Acart,A_orig,Bcart,B_orig,Xcart                                  ! Test variables for magnetic field and position
REAL(wp)                                 :: Bthet, Bzeta
REAL(wp),DIMENSION(3)                    :: grad_s,grad_thet,grad_zeta,grad_R,grad_Z                         ! Contravariant coordinates

! 2D theta/zeta fourier representation of variables in GVEC
REAL(wp),DIMENSION(1:X1_base_in%f%modes) :: X1_s,dX1ds_s                                                     
REAL(wp),DIMENSION(1:X2_base_in%f%modes) :: X2_s,dX2ds_s 
REAL(wp),DIMENSION(:), ALLOCATABLE       :: A_R_s, A_Rds_int, A_Z_s, A_Zds_int, A_phi_s, A_phids_int
REAL(wp),DIMENSION(:), ALLOCATABLE       :: B_R_s, B_Rds_int, B_Z_s, B_Zds_int, B_phi_s, B_phids_int
REAL(wp),DIMENSION(:), ALLOCATABLE       :: J_R_s, J_Rds_int, J_Z_s, J_Zds_int, J_phi_s, J_phids_int
REAL(wp),DIMENSION(1:LG_base_in%f%modes) :: LG_s,dGds_s

! 1D zeta fourier representation of variables in JOREK
REAL(wp), DIMENSION(:), ALLOCATABLE     :: X1_DOFs, X1_S_DOFs, X1_T_DOFs, X1_ST_DOFs, X2_DOFs, X2_S_DOFs, X2_T_DOFs, X2_ST_DOFs 
REAL(wp), DIMENSION(:), ALLOCATABLE     :: P_DOFs, P_S_DOFs, P_T_DOFs, P_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: A_R_DOFs, A_R_S_DOFs, A_R_T_DOFs, A_R_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: A_Z_DOFs, A_Z_S_DOFs, A_Z_T_DOFs, A_Z_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: A_phi_DOFs, A_phi_S_DOFs, A_phi_T_DOFs, A_phi_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: B_R_DOFs, B_R_S_DOFs, B_R_T_DOFs, B_R_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: B_Z_DOFs, B_Z_S_DOFs, B_Z_T_DOFs, B_Z_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: B_phi_DOFs, B_phi_S_DOFs, B_phi_T_DOFs, B_phi_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: J_R_DOFs, J_R_S_DOFs, J_R_T_DOFs, J_R_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: J_Z_DOFs, J_Z_S_DOFs, J_Z_T_DOFs, J_Z_ST_DOFs
REAL(wp), DIMENSION(:), ALLOCATABLE     :: J_phi_DOFs, J_phi_S_DOFs, J_phi_T_DOFs, J_phi_ST_DOFs
CLASS(t_fBase),ALLOCATABLE              :: fbase_zeta                                                                           ! Basis for 1D toroidal representation

! Variables for new GVEC fourier representations needed for JOREK inputs
INTEGER                     :: fac_nyq_fields                            !< nyquist factor for field fourier representations constructed for JOREK
INTEGER                     :: mn_max_fields(2)                          !< maximum number for new variables in SFL coordinates
TYPE(t_sgrid)               :: sgrid_fields                              !< uniform grid for generating new fourier representations
CLASS(t_base),  ALLOCATABLE :: A_R_base, A_Z_base, A_phi_base            !< container for base of variable A_*
REAL(wp),       ALLOCATABLE :: A_R(:,:), A_Z(:,:), A_phi(:,:)            !< data (1:nBase,1:modes) of A_* in GVEC coords
CLASS(t_base),  ALLOCATABLE :: B_R_base, B_Z_base, B_phi_base            !< container for base of variables B_*
REAL(wp),       ALLOCATABLE :: B_R(:,:), B_Z(:,:), B_phi(:,:)            !< data (1:nBase,1:modes) of B_* in GVEC coords
CLASS(t_base),  ALLOCATABLE :: J_R_base, J_Z_base, J_phi_base            !< container for base of variables J_*
REAL(wp),       ALLOCATABLE :: J_R(:,:), J_Z(:,:), J_phi(:,:)            !< data (1:nBase,1:modes) of J_* in GVEC coords

!===================================================================================================================================
! ------------------------------------------------------------------------------------------------------------
! ------- CALCULATE 2D FOURIER REPRESENTATION OF Vector Potential, Magnetic Field, and Current Density -------
! ------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)')'PREPARE FIELD BASES ...'
mn_max_fields = X1_base_in%f%mn_max
CALL sgrid_fields%copy(X1_base_in%s%grid)
fac_nyq_fields=4
CALL get_field_base(mn_max_fields,fac_nyq_fields,2,4,sgrid_fields,A_R_base,A_R)
ALLOCATE(A_R_s(1:A_R_base%f%modes), A_Rds_int(1:A_R_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,2,3,sgrid_fields,A_Z_base,A_Z)
ALLOCATE(A_Z_s(1:A_Z_base%f%modes), A_Zds_int(1:A_Z_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,2,5,sgrid_fields,A_phi_base,A_phi)
ALLOCATE(A_phi_s(1:A_phi_base%f%modes), A_phids_int(1:A_phi_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,1,4,sgrid_fields,B_R_base,B_R)
ALLOCATE(B_R_s(1:B_R_base%f%modes), B_Rds_int(1:B_R_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,1,3,sgrid_fields,B_Z_base,B_Z)
ALLOCATE(B_Z_s(1:B_Z_base%f%modes), B_Zds_int(1:B_Z_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,1,5,sgrid_fields,B_phi_base,B_phi)
ALLOCATE(B_phi_s(1:B_phi_base%f%modes), B_phids_int(1:B_phi_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,3,4,sgrid_fields,J_R_base,J_R)
ALLOCATE(J_R_s(1:J_R_base%f%modes), J_Rds_int(1:J_R_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,3,3,sgrid_fields,J_Z_base,J_Z)
ALLOCATE(J_Z_s(1:J_Z_base%f%modes), J_Zds_int(1:J_Z_base%f%modes))
CALL get_field_base(mn_max_fields,fac_nyq_fields,3,5,sgrid_fields,J_phi_base,J_phi)
ALLOCATE(J_phi_s(1:J_phi_base%f%modes), J_phids_int(1:J_phi_base%f%modes))

! -----------------------------------------------------------------------------
! ---------- GENERATE 3D POINTS FROM GVEC REPRESENTATION ----------------------
! -----------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A)')'PREPARE 3D DATA FOR GVEC-TO-JOREK ...'
CALL ProgressBar(0,Ns_out) !init
DO i_s=1,Ns_out
  spos          = s_pos(i_s)

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
  A_R_s(:)       =   A_R_base%s%evalDOF2D_s(spos,   A_R_base%f%modes,         0, A_R(:,:))
  A_Rds_int(:)   =   A_R_base%s%evalDOF2D_s(spos,   A_R_base%f%modes,   DERIV_S, A_R(:,:))
  A_Z_s(:)       =   A_Z_base%s%evalDOF2D_s(spos,   A_Z_base%f%modes,         0, A_Z(:,:))
  A_Zds_int(:)   =   A_Z_base%s%evalDOF2D_s(spos,   A_Z_base%f%modes,   DERIV_S, A_Z(:,:))
  A_phi_s(:)     =   A_phi_base%s%evalDOF2D_s(spos, A_phi_base%f%modes,         0, A_phi(:,:))
  A_phids_int(:) =   A_phi_base%s%evalDOF2D_s(spos, A_phi_base%f%modes,   DERIV_S, A_phi(:,:))
  B_R_s(:)       =   B_R_base%s%evalDOF2D_s(spos,   B_R_base%f%modes,         0, B_R(:,:))
  B_Rds_int(:)   =   B_R_base%s%evalDOF2D_s(spos,   B_R_base%f%modes,   DERIV_S, B_R(:,:))
  B_Z_s(:)       =   B_Z_base%s%evalDOF2D_s(spos,   B_Z_base%f%modes,         0, B_Z(:,:))
  B_Zds_int(:)   =   B_Z_base%s%evalDOF2D_s(spos,   B_Z_base%f%modes,   DERIV_S, B_Z(:,:))
  B_phi_s(:)     =   B_phi_base%s%evalDOF2D_s(spos, B_phi_base%f%modes,         0, B_phi(:,:))
  B_phids_int(:) =   B_phi_base%s%evalDOF2D_s(spos, B_phi_base%f%modes,   DERIV_S, B_phi(:,:))
  J_R_s(:)       =   J_R_base%s%evalDOF2D_s(spos,   J_R_base%f%modes,         0, J_R(:,:))
  J_Rds_int(:)   =   J_R_base%s%evalDOF2D_s(spos,   J_R_base%f%modes,   DERIV_S, J_R(:,:))
  J_Z_s(:)       =   J_Z_base%s%evalDOF2D_s(spos,   J_Z_base%f%modes,         0, J_Z(:,:))
  J_Zds_int(:)   =   J_Z_base%s%evalDOF2D_s(spos,   J_Z_base%f%modes,   DERIV_S, J_Z(:,:))
  J_phi_s(:)     =   J_phi_base%s%evalDOF2D_s(spos, J_phi_base%f%modes,         0, J_phi(:,:))
  J_phids_int(:) =   J_phi_base%s%evalDOF2D_s(spos, J_phi_base%f%modes,   DERIV_S, J_phi(:,:))
  
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
!$OMP           Acart,A_orig,Xcart,grad_s,grad_thet,grad_zeta,grad_R,grad_Z,                                                   &
!$OMP           Bthet,Bzeta,Bcart,B_orig)                                                                                      &
!$OMP   FIRSTPRIVATE(dLAdthet,dLAdzeta,G_int,dGds,dGdthet,dGdzeta)                                                             &
!$OMP   SHARED(i_s,Nzeta_out,Nthet_out,spos, thet_pos,zeta_pos,X1_base_in,X2_base_in,LG_base_in,                               &
!$OMP          A_R_base, A_Z_base, A_phi_base,B_R_base, B_Z_base, B_phi_base,J_R_base, J_Z_base,J_phi_base,                    &
!$OMP          hmap_r,X1_s,dX1ds_s,X2_s,dX2ds_s,LG_s,dGds_s,SFLcoord, Phi_int, dPhids_int, iota_int, Chi_int, dChids_int, P_int, dPds_int,      &
!$OMP          A_R_s, A_Rds_int,A_Z_s, A_Zds_int,A_phi_s, A_phids_int, B_R_s, B_Rds_int,B_Z_s, B_Zds_int,B_phi_s, B_phids_int,                  &
!$OMP          J_R_s, J_Rds_int,J_Z_s, J_Zds_int,J_phi_s, J_phids_int,                                                                          &
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
    A_R_int            = A_R_base%f%evalDOF_x(xp,          0,   A_R_s)
    dA_Rds_int         = A_R_base%f%evalDOF_x(xp,          0, A_Rds_int)
    dA_Rdthet_int      = A_R_base%f%evalDOF_x(xp, DERIV_THET,   A_R_s)
    d2A_Rdsdthet_int   = A_R_base%f%evalDOF_x(xp, DERIV_THET, A_Rds_int)
    A_Z_int            = A_Z_base%f%evalDOF_x(xp,          0,   A_Z_s)
    dA_Zds_int         = A_Z_base%f%evalDOF_x(xp,          0, A_Zds_int)
    dA_Zdthet_int      = A_Z_base%f%evalDOF_x(xp, DERIV_THET,   A_Z_s)
    d2A_Zdsdthet_int   = A_Z_base%f%evalDOF_x(xp, DERIV_THET, A_Zds_int)
    A_phi_int          = A_phi_base%f%evalDOF_x(xp,          0,   A_phi_s)
    dA_phids_int       = A_phi_base%f%evalDOF_x(xp,          0, A_phids_int)
    dA_phidthet_int    = A_phi_base%f%evalDOF_x(xp, DERIV_THET,   A_phi_s)
    d2A_phidsdthet_int = A_phi_base%f%evalDOF_x(xp, DERIV_THET, A_phids_int)
    
    ! Get B components
    B_R_int            = B_R_base%f%evalDOF_x(xp,          0,   B_R_s)
    dB_Rds_int         = B_R_base%f%evalDOF_x(xp,          0, B_Rds_int)
    dB_Rdthet_int      = B_R_base%f%evalDOF_x(xp, DERIV_THET,   B_R_s)
    d2B_Rdsdthet_int   = B_R_base%f%evalDOF_x(xp, DERIV_THET, B_Rds_int)
    B_Z_int            = B_Z_base%f%evalDOF_x(xp,          0,   B_Z_s)
    dB_Zds_int         = B_Z_base%f%evalDOF_x(xp,          0, B_Zds_int)
    dB_Zdthet_int      = B_Z_base%f%evalDOF_x(xp, DERIV_THET,   B_Z_s)
    d2B_Zdsdthet_int   = B_Z_base%f%evalDOF_x(xp, DERIV_THET, B_Zds_int)
    B_phi_int          = B_phi_base%f%evalDOF_x(xp,          0,   B_phi_s)
    dB_phids_int       = B_phi_base%f%evalDOF_x(xp,          0, B_phids_int)
    dB_phidthet_int    = B_phi_base%f%evalDOF_x(xp, DERIV_THET,   B_phi_s)
    d2B_phidsdthet_int = B_phi_base%f%evalDOF_x(xp, DERIV_THET, B_phids_int)

    ! Get J components
    J_R_int            = J_R_base%f%evalDOF_x(xp,          0,   J_R_s)
    dJ_Rds_int         = J_R_base%f%evalDOF_x(xp,          0, J_Rds_int)
    dJ_Rdthet_int      = J_R_base%f%evalDOF_x(xp, DERIV_THET,   J_R_s)
    d2J_Rdsdthet_int   = J_R_base%f%evalDOF_x(xp, DERIV_THET, J_Rds_int)
    J_Z_int            = J_Z_base%f%evalDOF_x(xp,          0,   J_Z_s)
    dJ_Zds_int         = J_Z_base%f%evalDOF_x(xp,          0, J_Zds_int)
    dJ_Zdthet_int      = J_Z_base%f%evalDOF_x(xp, DERIV_THET,   J_Z_s)
    d2J_Zdsdthet_int   = J_Z_base%f%evalDOF_x(xp, DERIV_THET, J_Zds_int)
    J_phi_int          = J_phi_base%f%evalDOF_x(xp,          0,   J_phi_s)
    dJ_phids_int       = J_phi_base%f%evalDOF_x(xp,          0, J_phids_int)
    dJ_phidthet_int    = J_phi_base%f%evalDOF_x(xp, DERIV_THET,   J_phi_s)
    d2J_phidsdthet_int = J_phi_base%f%evalDOF_x(xp, DERIV_THET, J_phids_int)

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
    data_scalar3D(ithet,izeta,i_s, S__)          = spos 
    data_scalar3D(ithet,izeta,i_s, THET__)       = thet_pos(ithet) 
    data_scalar3D(ithet,izeta,i_s, ZETA__)       = zeta_pos(izeta)
    data_scalar3D(ithet,izeta,i_s, X1__)         = X1_int 
    data_scalar3D(ithet,izeta,i_s, X1_S__)       = dX1ds 
    data_scalar3D(ithet,izeta,i_s, X1_T__)       = dX1dthet 
    data_scalar3D(ithet,izeta,i_s, X1_ST__)      = d2X1dsdthet 
    data_scalar3D(ithet,izeta,i_s, X2__)         = X2_int
    data_scalar3D(ithet,izeta,i_s, X2_S__)       = dX2ds
    data_scalar3D(ithet,izeta,i_s, X2_T__)       = dX2dthet
    data_scalar3D(ithet,izeta,i_s, X2_ST__)      = d2X2dsdthet
    data_scalar3D(ithet,izeta,i_s, P__)          = P_int
    data_scalar3D(ithet,izeta,i_s, P_S__)        = dPds_int
    data_scalar3D(ithet,izeta,i_s, A_R__)        = A_R_int
    data_scalar3D(ithet,izeta,i_s, A_R_S__)      = dA_Rds_int
    data_scalar3D(ithet,izeta,i_s, A_R_T__)      = dA_Rdthet_int
    data_scalar3D(ithet,izeta,i_s, A_R_ST__)     = d2A_Rdsdthet_int
    data_scalar3D(ithet,izeta,i_s, A_Z__)        = A_Z_int
    data_scalar3D(ithet,izeta,i_s, A_Z_S__)      = dA_Zds_int
    data_scalar3D(ithet,izeta,i_s, A_Z_T__)      = dA_Zdthet_int
    data_scalar3D(ithet,izeta,i_s, A_Z_ST__)     = d2A_Zdsdthet_int
    data_scalar3D(ithet,izeta,i_s, A_phi__)      = A_phi_int
    data_scalar3D(ithet,izeta,i_s, A_phi_S__)    = dA_phids_int
    data_scalar3D(ithet,izeta,i_s, A_phi_T__)    = dA_phidthet_int
    data_scalar3D(ithet,izeta,i_s, A_phi_ST__)   = d2A_phidsdthet_int
    data_scalar3D(ithet,izeta,i_s, B_R__)        = B_R_int
    data_scalar3D(ithet,izeta,i_s, B_R_S__)      = dB_Rds_int
    data_scalar3D(ithet,izeta,i_s, B_R_T__)      = dB_Rdthet_int
    data_scalar3D(ithet,izeta,i_s, B_R_ST__)     = d2B_Rdsdthet_int
    data_scalar3D(ithet,izeta,i_s, B_Z__)        = B_Z_int
    data_scalar3D(ithet,izeta,i_s, B_Z_S__)      = dB_Zds_int
    data_scalar3D(ithet,izeta,i_s, B_Z_T__)      = dB_Zdthet_int
    data_scalar3D(ithet,izeta,i_s, B_Z_ST__)     = d2B_Zdsdthet_int
    data_scalar3D(ithet,izeta,i_s, B_phi__)      = B_phi_int
    data_scalar3D(ithet,izeta,i_s, B_phi_S__)    = dB_phids_int
    data_scalar3D(ithet,izeta,i_s, B_phi_T__)    = dB_phidthet_int
    data_scalar3D(ithet,izeta,i_s, B_phi_ST__)   = d2B_phidsdthet_int
    data_scalar3D(ithet,izeta,i_s, J_R__)        = J_R_int
    data_scalar3D(ithet,izeta,i_s, J_R_S__)      = dJ_Rds_int
    data_scalar3D(ithet,izeta,i_s, J_R_T__)      = dJ_Rdthet_int
    data_scalar3D(ithet,izeta,i_s, J_R_ST__)     = d2J_Rdsdthet_int
    data_scalar3D(ithet,izeta,i_s, J_Z__)        = J_Z_int
    data_scalar3D(ithet,izeta,i_s, J_Z_S__)      = dJ_Zds_int
    data_scalar3D(ithet,izeta,i_s, J_Z_T__)      = dJ_Zdthet_int
    data_scalar3D(ithet,izeta,i_s, J_Z_ST__)     = d2J_Zdsdthet_int
    data_scalar3D(ithet,izeta,i_s, J_phi__)      = J_phi_int
    data_scalar3D(ithet,izeta,i_s, J_phi_S__)    = dJ_phids_int
    data_scalar3D(ithet,izeta,i_s, J_phi_T__)    = dJ_phidthet_int
    data_scalar3D(ithet,izeta,i_s, J_phi_ST__)   = d2J_phidsdthet_int
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
CALL fbase_new(fbase_zeta, (/0, X1_base_in%f%mn_max(2)/), (/1, Nzeta_out/), X1_base_in%f%nfp, "_sincos_", .false.)
ALLOCATE(X1_DOFs(fbase_zeta%modes),       X1_S_DOFs(fbase_zeta%modes),   X1_T_DOFs(fbase_zeta%modes), X1_ST_DOFs(fbase_zeta%modes))
ALLOCATE(X2_DOFs(fbase_zeta%modes),       X2_S_DOFs(fbase_zeta%modes),   X2_T_DOFs(fbase_zeta%modes), X2_ST_DOFs(fbase_zeta%modes))
ALLOCATE(P_DOFs(fbase_zeta%modes),        P_S_DOFs(fbase_zeta%modes),   P_T_DOFs(fbase_zeta%modes), P_ST_DOFs(fbase_zeta%modes))
ALLOCATE(A_R_DOFs(fbase_zeta%modes),     A_R_S_DOFs(fbase_zeta%modes),   A_R_T_DOFs(fbase_zeta%modes), A_R_ST_DOFs(fbase_zeta%modes))
ALLOCATE(A_Z_DOFs(fbase_zeta%modes),     A_Z_S_DOFs(fbase_zeta%modes),   A_Z_T_DOFs(fbase_zeta%modes), A_Z_ST_DOFs(fbase_zeta%modes))
ALLOCATE(A_phi_DOFs(fbase_zeta%modes),     A_phi_S_DOFs(fbase_zeta%modes),   A_phi_T_DOFs(fbase_zeta%modes), A_phi_ST_DOFs(fbase_zeta%modes))
ALLOCATE(B_R_DOFs(fbase_zeta%modes),     B_R_S_DOFs(fbase_zeta%modes),   B_R_T_DOFs(fbase_zeta%modes), B_R_ST_DOFs(fbase_zeta%modes))
ALLOCATE(B_Z_DOFs(fbase_zeta%modes),     B_Z_S_DOFs(fbase_zeta%modes),   B_Z_T_DOFs(fbase_zeta%modes), B_Z_ST_DOFs(fbase_zeta%modes))
ALLOCATE(B_phi_DOFs(fbase_zeta%modes),     B_phi_S_DOFs(fbase_zeta%modes),   B_phi_T_DOFs(fbase_zeta%modes), B_phi_ST_DOFs(fbase_zeta%modes))
ALLOCATE(J_R_DOFs(fbase_zeta%modes),     J_R_S_DOFs(fbase_zeta%modes),   J_R_T_DOFs(fbase_zeta%modes), J_R_ST_DOFs(fbase_zeta%modes))
ALLOCATE(J_Z_DOFs(fbase_zeta%modes),     J_Z_S_DOFs(fbase_zeta%modes),   J_Z_T_DOFs(fbase_zeta%modes), J_Z_ST_DOFs(fbase_zeta%modes))
ALLOCATE(J_phi_DOFs(fbase_zeta%modes),     J_phi_S_DOFs(fbase_zeta%modes),   J_phi_T_DOFs(fbase_zeta%modes), J_phi_ST_DOFs(fbase_zeta%modes))
n_modes = fbase_zeta%modes
sin_range(:) = fbase_zeta%sin_range(:)
cos_range(:) = fbase_zeta%cos_range(:)
DO ithet=1, Nthet_out
  DO i_s=1, Ns_out
    ! Zero all DOFs
    X1_DOFS(:)       = 0    
    X1_S_DOFS(:)     = 0    
    X1_T_DOFS(:)     = 0    
    X1_ST_DOFS(:)    = 0    
    X2_DOFS(:)       = 0    
    X2_S_DOFS(:)     = 0    
    X2_T_DOFS(:)     = 0    
    X2_ST_DOFS(:)    = 0
    P_DOFS(:)        = 0    
    P_S_DOFS(:)      = 0    
    P_T_DOFS(:)      = 0    
    P_ST_DOFS(:)     = 0
    A_R_DOFS(:)      = 0    
    A_R_S_DOFS(:)    = 0    
    A_R_T_DOFS(:)    = 0    
    A_R_ST_DOFS(:)   = 0
    A_Z_DOFS(:)      = 0    
    A_Z_S_DOFS(:)    = 0    
    A_Z_T_DOFS(:)    = 0    
    A_Z_ST_DOFS(:)   = 0
    A_phi_DOFS(:)    = 0    
    A_phi_S_DOFS(:)  = 0    
    A_phi_T_DOFS(:)  = 0    
    A_phi_ST_DOFS(:) = 0
    B_R_DOFS(:)      = 0    
    B_R_S_DOFS(:)    = 0    
    B_R_T_DOFS(:)    = 0    
    B_R_ST_DOFS(:)   = 0
    B_Z_DOFS(:)      = 0    
    B_Z_S_DOFS(:)    = 0    
    B_Z_T_DOFS(:)    = 0    
    B_Z_ST_DOFS(:)   = 0
    B_phi_DOFS(:)    = 0    
    B_phi_S_DOFS(:)  = 0    
    B_phi_T_DOFS(:)  = 0    
    B_phi_ST_DOFS(:) = 0
    J_R_DOFS(:)      = 0    
    J_R_S_DOFS(:)    = 0    
    J_R_T_DOFS(:)    = 0    
    J_R_ST_DOFS(:)   = 0
    J_Z_DOFS(:)      = 0    
    J_Z_S_DOFS(:)    = 0    
    J_Z_T_DOFS(:)    = 0    
    J_Z_ST_DOFS(:)   = 0
    J_phi_DOFS(:)    = 0    
    J_phi_S_DOFS(:)  = 0    
    J_phi_T_DOFS(:)  = 0    
    J_phi_ST_DOFS(:) = 0
    
    ! Calculate toroidal fourier representation of DOFs on node
    X1_DOFS(:)       = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X1__))
    X1_S_DOFS(:)     = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X1_S__))
    X1_T_DOFS(:)     = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X1_T__))
    X1_ST_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X1_ST__))
    X2_DOFS(:)       = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X2__))
    X2_S_DOFS(:)     = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X2_S__))
    X2_T_DOFS(:)     = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X2_T__))
    X2_ST_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,X2_ST__))
    P_DOFS(:)        = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,P__))         ! Leave poloidal derivative terms at 0 for radial profiles
    P_S_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,P_S__))
    A_R_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_R__))
    A_R_S_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_R_S__))
    A_R_T_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_R_T__))
    A_R_ST_DOFS(:)   = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_R_ST__))
    A_Z_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_Z__))
    A_Z_S_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_Z_S__))
    A_Z_T_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_Z_T__))
    A_Z_ST_DOFS(:)   = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_Z_ST__))
    A_phi_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_phi__))
    A_phi_S_DOFS(:)  = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_phi_S__))
    A_phi_T_DOFS(:)  = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_phi_T__))
    A_phi_ST_DOFS(:) = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,A_phi_ST__))
    B_R_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_R__))
    B_R_S_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_R_S__))
    B_R_T_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_R_T__))
    B_R_ST_DOFS(:)   = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_R_ST__))
    B_Z_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_Z__))
    B_Z_S_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_Z_S__))
    B_Z_T_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_Z_T__))
    B_Z_ST_DOFS(:)   = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_Z_ST__))
    B_phi_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_phi__))
    B_phi_S_DOFS(:)  = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_phi_S__))
    B_phi_T_DOFS(:)  = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_phi_T__))
    B_phi_ST_DOFS(:) = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,B_phi_ST__))
    J_R_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_R__))
    J_R_S_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_R_S__))
    J_R_T_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_R_T__))
    J_R_ST_DOFS(:)   = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_R_ST__))
    J_Z_DOFS(:)      = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_Z__))
    J_Z_S_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_Z_S__))
    J_Z_T_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_Z_T__))
    J_Z_ST_DOFS(:)   = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_Z_ST__))
    J_phi_DOFS(:)    = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_phi__))
    J_phi_S_DOFS(:)  = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_phi_S__))
    J_phi_T_DOFS(:)  = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_phi_T__))
    J_phi_ST_DOFS(:) = fbase_zeta%initDOF(data_scalar3D(ithet,:,i_s,J_phi_ST__))

    ! Store DOFs for writing to file
    data_scalar2D(ithet, i_s, 1:R__*n_modes)                                    = X1_DOFS(:)
    data_scalar2D(ithet, i_s, (R_S__-1)*n_modes+1:R_S__*n_modes)                = X1_S_DOFS(:)
    data_scalar2D(ithet, i_s, (R_T__-1)*n_modes+1:R_T__*n_modes)                = X1_T_DOFS(:)
    data_scalar2D(ithet, i_s, (R_ST__-1)*n_modes+1:R_ST__*n_modes)              = X1_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (Z__-1)*n_modes+1:Z__*n_modes)                    = X2_DOFS(:)
    data_scalar2D(ithet, i_s, (Z_S__-1)*n_modes+1:Z_S__*n_modes)                = X2_S_DOFS(:)
    data_scalar2D(ithet, i_s, (Z_T__-1)*n_modes+1:Z_T__*n_modes)                = X2_T_DOFS(:)
    data_scalar2D(ithet, i_s, (Z_ST__-1)*n_modes+1:Z_ST__*n_modes)              = X2_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (P2D__-1)*n_modes+1:P2D__*n_modes)                = P_DOFS(:)
    data_scalar2D(ithet, i_s, (P2D_S__-1)*n_modes+1:P2D_S__*n_modes)            = P_S_DOFS(:)
    data_scalar2D(ithet, i_s, (P2D_T__-1)*n_modes+1:P2D_T__*n_modes)            = P_T_DOFS(:)
    data_scalar2D(ithet, i_s, (P2D_ST__-1)*n_modes+1:P2D_ST__*n_modes)          = P_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (A_R2D__-1)*n_modes+1:A_R2D__*n_modes)            = A_R_DOFS(:)
    data_scalar2D(ithet, i_s, (A_R2D_S__-1)*n_modes+1:A_R2D_S__*n_modes)        = A_R_S_DOFS(:)
    data_scalar2D(ithet, i_s, (A_R2D_T__-1)*n_modes+1:A_R2D_T__*n_modes)        = A_R_T_DOFS(:)
    data_scalar2D(ithet, i_s, (A_R2D_ST__-1)*n_modes+1:A_R2D_ST__*n_modes)      = A_R_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (A_Z2D__-1)*n_modes+1:A_Z2D__*n_modes)            = A_Z_DOFS(:)
    data_scalar2D(ithet, i_s, (A_Z2D_S__-1)*n_modes+1:A_Z2D_S__*n_modes)        = A_Z_S_DOFS(:)
    data_scalar2D(ithet, i_s, (A_Z2D_T__-1)*n_modes+1:A_Z2D_T__*n_modes)        = A_Z_T_DOFS(:)
    data_scalar2D(ithet, i_s, (A_Z2D_ST__-1)*n_modes+1:A_Z2D_ST__*n_modes)      = A_Z_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (A_phi2D__-1)*n_modes+1:A_phi2D__*n_modes)        = A_phi_DOFS(:)
    data_scalar2D(ithet, i_s, (A_phi2D_S__-1)*n_modes+1:A_phi2D_S__*n_modes)    = A_phi_S_DOFS(:)
    data_scalar2D(ithet, i_s, (A_phi2D_T__-1)*n_modes+1:A_phi2D_T__*n_modes)    = A_phi_T_DOFS(:)
    data_scalar2D(ithet, i_s, (A_phi2D_ST__-1)*n_modes+1:A_phi2D_ST__*n_modes)  = A_phi_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (B_R2D__-1)*n_modes+1:B_R2D__*n_modes)            = B_R_DOFS(:)
    data_scalar2D(ithet, i_s, (B_R2D_S__-1)*n_modes+1:B_R2D_S__*n_modes)        = B_R_S_DOFS(:)
    data_scalar2D(ithet, i_s, (B_R2D_T__-1)*n_modes+1:B_R2D_T__*n_modes)        = B_R_T_DOFS(:)
    data_scalar2D(ithet, i_s, (B_R2D_ST__-1)*n_modes+1:B_R2D_ST__*n_modes)      = B_R_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (B_Z2D__-1)*n_modes+1:B_Z2D__*n_modes)            = B_Z_DOFS(:)
    data_scalar2D(ithet, i_s, (B_Z2D_S__-1)*n_modes+1:B_Z2D_S__*n_modes)        = B_Z_S_DOFS(:)
    data_scalar2D(ithet, i_s, (B_Z2D_T__-1)*n_modes+1:B_Z2D_T__*n_modes)        = B_Z_T_DOFS(:)
    data_scalar2D(ithet, i_s, (B_Z2D_ST__-1)*n_modes+1:B_Z2D_ST__*n_modes)      = B_Z_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (B_phi2D__-1)*n_modes+1:B_phi2D__*n_modes)        = B_phi_DOFS(:)
    data_scalar2D(ithet, i_s, (B_phi2D_S__-1)*n_modes+1:B_phi2D_S__*n_modes)    = B_phi_S_DOFS(:)
    data_scalar2D(ithet, i_s, (B_phi2D_T__-1)*n_modes+1:B_phi2D_T__*n_modes)    = B_phi_T_DOFS(:)
    data_scalar2D(ithet, i_s, (B_phi2D_ST__-1)*n_modes+1:B_phi2D_ST__*n_modes)  = B_phi_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (J_R2D__-1)*n_modes+1:J_R2D__*n_modes)            = J_R_DOFS(:)
    data_scalar2D(ithet, i_s, (J_R2D_S__-1)*n_modes+1:J_R2D_S__*n_modes)        = J_R_S_DOFS(:)
    data_scalar2D(ithet, i_s, (J_R2D_T__-1)*n_modes+1:J_R2D_T__*n_modes)        = J_R_T_DOFS(:)
    data_scalar2D(ithet, i_s, (J_R2D_ST__-1)*n_modes+1:J_R2D_ST__*n_modes)      = J_R_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (J_Z2D__-1)*n_modes+1:J_Z2D__*n_modes)            = J_Z_DOFS(:)
    data_scalar2D(ithet, i_s, (J_Z2D_S__-1)*n_modes+1:J_Z2D_S__*n_modes)        = J_Z_S_DOFS(:)
    data_scalar2D(ithet, i_s, (J_Z2D_T__-1)*n_modes+1:J_Z2D_T__*n_modes)        = J_Z_T_DOFS(:)
    data_scalar2D(ithet, i_s, (J_Z2D_ST__-1)*n_modes+1:J_Z2D_ST__*n_modes)      = J_Z_ST_DOFS(:)
    data_scalar2D(ithet, i_s, (J_phi2D__-1)*n_modes+1:J_phi2D__*n_modes)        = J_phi_DOFS(:)
    data_scalar2D(ithet, i_s, (J_phi2D_S__-1)*n_modes+1:J_phi2D_S__*n_modes)    = J_phi_S_DOFS(:)
    data_scalar2D(ithet, i_s, (J_phi2D_T__-1)*n_modes+1:J_phi2D_T__*n_modes)    = J_phi_T_DOFS(:)
    data_scalar2D(ithet, i_s, (J_phi2D_ST__-1)*n_modes+1:J_phi2D_ST__*n_modes)  = J_phi_ST_DOFS(:)
  END DO ! ithet=1, Nthet_out
END DO ! i_s=1, Ns_out


SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)

DEALLOCATE(X1_DOFs,  X1_S_DOFs,   X1_T_DOFs, X1_ST_DOFs)
DEALLOCATE(X2_DOFs,  X2_S_DOFs,   X2_T_DOFs, X2_ST_DOFs)
DEALLOCATE(P_DOFs,   P_S_DOFs,   P_T_DOFs,   P_ST_DOFs)
DEALLOCATE(A_R_DOFs, A_R_S_DOFs,  A_R_T_DOFs, A_R_ST_DOFs)
DEALLOCATE(A_Z_DOFs, A_Z_S_DOFs,  A_Z_T_DOFs, A_Z_ST_DOFs)
DEALLOCATE(A_phi_DOFs, A_phi_S_DOFs,  A_phi_T_DOFs, A_phi_ST_DOFs)
DEALLOCATE(B_R_DOFs, B_R_S_DOFs,  B_R_T_DOFs, B_R_ST_DOFs)
DEALLOCATE(B_Z_DOFs, B_Z_S_DOFs,  B_Z_T_DOFs, B_Z_ST_DOFs)
DEALLOCATE(B_phi_DOFs, B_phi_S_DOFs,  B_phi_T_DOFs, B_phi_ST_DOFs)
DEALLOCATE(J_R_DOFs, J_R_S_DOFs,  J_R_T_DOFs, J_R_ST_DOFs)
DEALLOCATE(J_Z_DOFs, J_Z_S_DOFs,  J_Z_T_DOFs, J_Z_ST_DOFs)
DEALLOCATE(J_phi_DOFs, J_phi_S_DOFs,  J_phi_T_DOFs, J_phi_ST_DOFs)
DEALLOCATE(A_R_s, A_Rds_int)
DEALLOCATE(A_Z_s, A_Zds_int)
DEALLOCATE(A_phi_s, A_phids_int)
DEALLOCATE(B_R_s, B_Rds_int)
DEALLOCATE(B_Z_s, B_Zds_int)
DEALLOCATE(B_phi_s, B_phids_int)
DEALLOCATE(J_R_s, J_Rds_int)
DEALLOCATE(J_Z_s, J_Zds_int)
DEALLOCATE(J_phi_s, J_phids_int)
END SUBROUTINE gvec_to_jorek_prepare

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
  if (generate_test_data) then
    ! Write 3D data only
    DO iVar=1,nVarScalar3D
      WRITE(ioUnit,'(A)',ADVANCE='NO')'##<< 3D scalar variable (1:Ntheta,1:Nzeta,1:Ns), Variable name: '
      WRITE(ioUNIT,'(A)')' "'//TRIM(StrVarNamesScalar3D(iVar))//'"'
      WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar3D(1:Nthet_out,1:Nzeta_out,1:Ns_out,iVar) 
    END DO !iVar=1,nVarScalar3D
  else
    ! Write 2D data only
    DO iVar=1,nVarScalar2D
      WRITE(ioUnit,'(A)',ADVANCE='NO')'##<< 2D scalar variable fourier modes (1:Ntheta,1:Ns), Variable name: '
      WRITE(ioUNIT,'(A)')' "'//TRIM(StrVarNamesScalar2D(iVar))//'"'
      if (iVar .eq. 1) then
        WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar2D(1:Nthet_out,1:Ns_out,1:n_modes) 
      else
        WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar2D(1:Nthet_out,1:Ns_out,(iVar-1)*n_modes+1:iVar*n_modes) 
      endif
    END DO !iVar=1,nVarScalar2D
  endif

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

END SUBROUTINE finalize_gvec_to_jorek

END MODULE MODgvec_gvec_to_jorek
