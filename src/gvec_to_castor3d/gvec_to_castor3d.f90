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
!!# Module **gvec_to_castor3d**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_castor3d
! MODULES
USE MODgvec_Globals, ONLY:wp,UNIT_stdOut,fmt_sep
IMPLICIT NONE
PRIVATE

INTERFACE get_cla_gvec_to_castor3d
  MODULE PROCEDURE get_cla_gvec_to_castor3d
END INTERFACE

INTERFACE init_gvec_to_castor3d
  MODULE PROCEDURE init_gvec_to_castor3d
END INTERFACE

INTERFACE gvec_to_castor3d_prepare
  MODULE PROCEDURE gvec_to_castor3d_prepare
END INTERFACE

INTERFACE gvec_to_castor3d_writeToFile
  MODULE PROCEDURE gvec_to_castor3d_writeToFile_ASCII
END INTERFACE

INTERFACE finalize_gvec_to_castor3d
  MODULE PROCEDURE finalize_gvec_to_castor3d
END INTERFACE

PUBLIC::get_cla_gvec_to_castor3d
PUBLIC::init_gvec_to_castor3d
PUBLIC::gvec_to_castor3d_prepare
PUBLIC::gvec_to_castor3d_writeToFile
PUBLIC::finalize_gvec_to_castor3d

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get command line arguments 
!!
!===================================================================================================================================
SUBROUTINE get_CLA_gvec_to_castor3d() 
! MODULES
USE MODgvec_cla
USE MODgvec_gvec_to_castor3d_Vars, ONLY: fileName, Ns_out,npfactor,factorSFL,Nthet_out,Nzeta_out,SFLcoord
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=STRLEN)   :: f_str
CHARACTER(LEN=24)       :: execname="convert_gvec_to_castor3d"
LOGICAL                 :: commandFailed
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
  !positional argument
  CALL cla_posarg_register('gvecfile.dat', &
       'Filename of GVEC restart file [MANDATORY]',  cla_char,'xxx') !    

  CALL cla_validate(execname)
  CALL cla_get('-r',NS_out)
  CALL cla_get('-n',npfactor)
  CALL cla_get('-p',Nthet_out)
  CALL cla_get('-t',Nzeta_out)
  CALL cla_get('-s',SFLcoord)
  CALL cla_get('-f',factorSFL)
  CALL cla_get('gvecfile.dat',f_str)
  filename=TRIM(f_str)

  commandFailed=.FALSE.
  IF(.NOT.((cla_key_present('-r')).AND.(Ns_out.GE.2))) THEN
    IF(.NOT.commandFailed) CALL cla_help(execname)
    commandFailed=.TRUE.
    SWRITE(UNIT_StdOut,*) " ==> [-r,--rpoints] argument is MANDATORY and must be >1 !!!"
  END IF
  IF((SFLcoord.LT.0).OR.(SFLcoord.GT.2)) THEN
    IF(.NOT.commandFailed) CALL cla_help(execname)
    commandFailed=.TRUE.
    SWRITE(UNIT_StdOut,*) " ==> [-s,--sflcoord] argument  must be 0,1,2 !!!"
  END IF 
  IF((INDEX(filename,'xxx').NE.0))THEN
    IF(.NOT.commandFailed) CALL cla_help(execname)
    commandFailed=.TRUE.
    SWRITE(UNIT_StdOut,*) " ==> filename gvecfile.dat is MANDATORY must be specified !!!"
  END IF
  IF(commandFailed) STOP

  SWRITE(UNIT_stdOut,'(A)')   ' INPUT PARAMETERS:' 
  SWRITE(UNIT_stdOut,'(A,I6)')'  * Number of radial points        : ',Ns_out
  SWRITE(UNIT_stdOut,'(A,I4)')'  * npfactor points from modes     : ',npfactor
  SWRITE(UNIT_stdOut,'(A,I4)')'  * SFL coordinates flag           : ',SFLcoord
  SWRITE(UNIT_stdOut,'(A,I4)')'  * factor for modes of SFL coords : ',factorSFL
  SWRITE(UNIT_stdOut,'(A,A)') '  * GVEC input file                : ',TRIM(fileName)
  SWRITE(UNIT_stdOut,fmt_sep)

!  nArgs=COMMAND_ARGUMENT_COUNT()
!  
!  commandFailed=.FALSE.
!  IF(nArgs.EQ.3)THEN
!    CALL GET_COMMAND_ARGUMENT(1,tmp)
!    READ(tmp,*,IOSTAT=err)Ns_out
!    IF(err.NE.0) commandFailed=.TRUE.
!    CALL GET_COMMAND_ARGUMENT(2,tmp)
!    READ(tmp,*,IOSTAT=err)npfactor
!    IF(err.NE.0) commandFailed=.TRUE.
!    CALL GET_COMMAND_ARGUMENT(3,FileName)
!  ELSE
!    commandFailed=.TRUE.
!  END IF
!  IF(commandfailed)  STOP ' GVEC TO CASTOR3D: command not correct: "./executable Ns FourierFactor gvec_file.dat"'
!  IF(Ns_out.LT.2) STOP ' GVEC TO CASTOR3D: choose number in radial dierction Ns>=2' 
!  IF(npfactor.GT.16) STOP ' GVEC TO CASTOR3D: choose fourierFactor  <=16' 

!  ppos = INDEX(FileName,'.dat',back=.TRUE.)
!  IF ( ppos  .GT. 0 )THEN
!     FileNameOut = FileName(1:ppos-1)
!  ELSE
!    STOP ' did not find .dat extension in input file'
!  END IF

END SUBROUTINE get_cla_gvec_to_castor3d

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE init_gvec_to_castor3d() 
! MODULES
USE MODgvec_Globals,ONLY: TWOPI
USE MODgvec_ReadState         ,ONLY: ReadState
USE MODgvec_ReadState_vars    ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_ReadState_vars    ,ONLY: LA_r,X1_r,X2_r 
USE MODgvec_transform_sfl_vars,ONLY: X1sfl_base,X1sfl,X2sfl_base,X2sfl
USE MODgvec_transform_sfl     ,ONLY: BuildTransform_SFL
USE MODgvec_gvec_to_castor3d_vars
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
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT GVEC-TO-CASTOR3D ...'

  CALL ReadState(TRIM(fileName))


  mn_max_out(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
  mn_max_out(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))
  nfp_out          = X1_base_r%f%nfp

  IF((X1_base_r%f%sin_cos.EQ._COS_).AND.(X2_base_r%f%sin_cos.EQ._SIN_).AND.(LA_base_r%f%sin_cos.EQ._SIN_))THEN
    asym_out = 0 !R~cos,Z~sin,lambda~sin
  ELSE
    asym_out = 1 !full fourier
  END IF

  IF(SFLcoord.NE.0)THEN
   mn_max_out=mn_max_out*factorSFL !*SFLfactor on modes
   CALL buildTransform_SFL(X1_base_r,X1_r,X2_base_r,X2_r,LA_base_r,LA_r,Ns_out,mn_max_out,SFLcoord)
  END IF
  

  ALLOCATE(s_pos(Ns_out))
  ALLOCATE(data_1D(nVar1D,Ns_out))

  s_pos(1)=1.0e-06_wp !avoid axis
  DO i=2,Ns_out-1
      s_pos(i) = REAL(i-1,wp)/REAL(Ns_out-1,wp)
  END DO !i
  s_pos(Ns_out)=1. - 1.0e-12_wp !avoid edge
  IF(Nthet_out.EQ.-1) THEN !overwrite with default from factorFouier
    Nthet_out = npfactor*mn_max_out(1)
  END IF
  IF(Nzeta_out.EQ.-1) THEN !overwrite with default from factorFourier
    Nzeta_out = MAX(1,npfactor*mn_max_out(2)) !if n=0, output 1 point
  END IF
  !check
  IF(Nthet_out .LE. mn_max_out(1)) STOP 'number of poloidal points for output should be >m_max!'
  IF(Nzeta_out .LE. mn_max_out(2)) STOP 'number of torodial points for output must be >n_max!'

  ALLOCATE(thet_pos(Nthet_out))
  ALLOCATE(zeta_pos(Nzeta_out))
  DO i=1,Nthet_out
    thet_pos(i)=(TWOPI*REAL((i-1),wp))/REAL(Nthet_out) 
  END DO
  !zeta goes in opposite direction!!! -> iota and phi,phi' have opposite sign
  DO i=1,Nzeta_out
    zeta_pos(i)=-(TWOPI*REAL((i-1),wp))/REAL((Nzeta_out*nfp_out),wp)
  END DO

  ALLOCATE(data_scalar3D(  Nthet_out,Nzeta_out,Ns_out,nVarscalar3D))
  ALLOCATE(data_vector3D(3,Nthet_out,Nzeta_out,Ns_out,nVarvector3D))

  SWRITE(UNIT_stdOut,'(A,3I6)')'  Number OF N_s,N_theta,N_zeta evaluation points:',Ns_out,Nthet_out,Nzeta_out
  SWRITE(UNIT_stdOut,'(A)')'... DONE'
  SWRITE(UNIT_stdOut,fmt_sep)

  SELECT CASE(SFLcoord)
  CASE(0)
    CALL gvec_to_castor3D_prepare(X1_base_r,X1_r,X2_base_r,X2_r,LA_base_r,LA_r)
  CASE(1)
    CALL gvec_to_castor3D_prepare(X1sfl_base,X1sfl,X2sfl_base,X2sfl,LA_base_r,LA_r) !LA not needed, used as placeholder
  CASE(2)
    !CALL gvec_to_castor3D_prepare(X1_base_r%s,X1sfl_base,X1sfl,X2sfl_base,X2sfl,Gsfl_base,Gsfl)
    STOP 'BOOZER not yet implemented'
  CASE DEFAULT
    SWRITE(UNIT_StdOut,*)'This SFLcoord is not valid',SFLcoord
    STOP
  END SELECT
END SUBROUTINE init_gvec_to_castor3d


!===================================================================================================================================
!> prepare all data to be written
!!
!===================================================================================================================================
SUBROUTINE gvec_to_castor3d_prepare(X1_base_in,X1_in,X2_base_in,X2_in,LG_base_in,LG_in)
! MODULES
USE MODgvec_gvec_to_castor3d_Vars 
USE MODgvec_Globals,        ONLY: CROSS,TWOPI,PI
USE MODgvec_ReadState_Vars, ONLY: profiles_1d,hmap_r,sbase_prof !for profiles
USE MODgvec_Base,           ONLY: t_base
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
INTEGER                                 :: i_s,iMode,ithet,izeta
REAL(wp)                                :: spos,xp(2),sqrtG
REAL(wp)                                :: dX1ds,dX1dthet,dX1dzeta
REAL(wp)                                :: dX2ds,dX2dthet,dX2dzeta
REAL(wp)                                :: dLAdthet,dLAdzeta
REAL(wp)                                :: Phi_int,dPhids_int,Chi_int,dChids_int,iota_int
REAL(wp)                                :: pressure_int,F_loc,Favg_int,Fmin_int,Fmax_int
REAL(wp)                                :: X1_int,X2_int,G_int,dGds,dGdthet,dGdzeta
REAL(wp)                                :: Ipol_int,Itor_int
REAL(wp),DIMENSION(3)                   :: qvec,e_s,e_thet,e_zeta,Bfield,grad_zeta
REAL(wp),DIMENSION(1:X1_base_in%f%modes) :: X1_s,dX1ds_s
REAL(wp),DIMENSION(1:X2_base_in%f%modes) :: X2_s,dX2ds_s
REAL(wp),DIMENSION(1:LG_base_in%f%modes) :: LG_s,dGds_s
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'PREPARE DATA FOR GVEC-TO-CASTOR3D ...'
DO i_s=1,Ns_out
  IF(MOD(i_s,MAX(1,Ns_out/100)).EQ.0) THEN
    SWRITE(UNIT_stdOut,'(4X,I4,A4,I4,A13,A1)',ADVANCE='NO')i_s, ' of ',NS_out,' evaluated...',ACHAR(13)
  END IF
  spos          = s_pos(i_s)
  data_1D(SPOS__,i_s)=s_pos(i_s)

  Phi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,1))
  dPhids_int  = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  Chi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,2)) !Chi not yet working
  !dChids_int  = sbase_profs%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2)) !Chi not yet working
  iota_int    = sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,3))
  dChids_int  = dPhids_int*iota_int 
  pressure_int= sbase_prof%evalDOF_s(spos, 0,profiles_1d(:,4))

  
  Fmin_int=+1.0e12
  Fmax_int=-1.0e12
  Favg_int = 0.
  Itor_int = 0.
  Ipol_int = 0.

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
  ELSEIF(SFLcoord.EQ.2)THEN !BOOZER
    LG_s(  :) = LG_base_in%s%evalDOF2D_s(spos,LG_base_in%f%modes,       0,LG_in(:,:))
    dGds_s(:) = LG_base_in%s%evalDOF2D_s(spos,LG_base_in%f%modes, DERIV_S,LG_in(:,:))
  END IF
  !interpolate in the angles
  DO izeta=1,Nzeta_out; DO ithet=1,Nthet_out
    xp=(/thet_pos(ithet),zeta_pos(izeta)/)


    X1_int   = X1_base_in%f%evalDOF_x(xp,          0, X1_s  )
    dX1ds    = X1_base_in%f%evalDOF_x(xp,          0,dX1ds_s)
    dX1dthet = X1_base_in%f%evalDOF_x(xp, DERIV_THET, X1_s  )
    dX1dzeta = X1_base_in%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )
    
    X2_int   = X2_base_in%f%evalDOF_x(xp,          0, X2_s  )
    dX2ds    = X2_base_in%f%evalDOF_x(xp,          0,dX2ds_s)
    dX2dthet = X2_base_in%f%evalDOF_x(xp, DERIV_THET, X2_s  )
    dX2dzeta = X2_base_in%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )
    
    IF(SFLcoord.EQ.0)THEN !GVEC coordinates (else=0)
      dLAdthet = LG_base_in%f%evalDOF_x(xp, DERIV_THET, LG_s)
      dLAdzeta = LG_base_in%f%evalDOF_x(xp, DERIV_ZETA, LG_s)
    END IF
    
    IF(SFLcoord.EQ.2)THEN !BOOZER coordinates (else=0)
      G_int   = LG_base_in%f%evalDOF_x(xp,         0, LG_s)
      dGds    = LG_base_in%f%evalDOF_x(xp,         0, dGds_s)
      dGdthet = LG_base_in%f%evalDOF_x(xp, DERIV_THET, LG_s)
      dGdzeta = LG_base_in%f%evalDOF_x(xp, DERIV_ZETA, LG_s)
    END IF
    
    qvec=(/X1_int,X2_int,xp(2)/) !(X1,X2,zeta)
    e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,        dGds   /)) !dxvec/ds
    e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,        dGdthet/)) !dxvec/dthet
    e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp +dGdzeta/)) !dxvec/dzeta
    sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))
   
    !check: J = e_s*(e_thet x e_zeta) 
    !sqrtG_check = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 
    !WRITE(*,*)'CHECK sqrtG',sqrtG,sqrtG_check,sqrtG-sqrtG_check

    Bfield(:) = (  e_thet(:)*(iota_int-dLAdzeta )  & 
                 + e_zeta(:)*(1.0_wp  +dLAdthet ) )*(dPhids_int/sqrtG)

    !grad_s(:)   = CROSS(e_thet,e_zeta)/sqrtG
    !grad_thet(:)= CROSS(e_zeta,e_s   )/sqrtG
    grad_zeta(:)= CROSS(e_s   ,e_thet)/sqrtG


    !F-profile, only makes sense for tokamak configurations (n=0): F=-phi' R*(1+dlambda_dtheta)/(Jac*|\nabla\zeta|)
    F_loc = -dPhids_int*X1_int*(1.0_wp+dLAdthet)/(sqrtG*SQRT(SUM(grad_zeta(:)**2)))

    Favg_int = Favg_int + F_loc
    Fmin_int = MIN(Fmin_int,F_loc)
    Fmax_int = MAX(Fmax_int,F_loc)

    !poloidal and toroidal current profiles, line integral: integration over one angle /average over other...
    !Itor= int_0^2pi B_theta dtheta = (nfp/2pi) int_0^2pi int_0^(2pi/nfp) B_theta dtheta dzeta
    !Ipol= int_0^2pi B_zeta  dzeta  = nfp* int_0^(2pi/nfp) B_zeta dzeta = (nfp/2pi) int_0^2pi int_0^(2pi/nfp) B_zeta  dtheta dzeta
    Itor_int = Itor_int+ SUM(Bfield(:)*e_thet(:))   !B_theta=B.e_thet 
    Ipol_int = Ipol_int+ SUM(Bfield(:)*e_zeta(:))   !B_zeta =B.e_zeta

    !========== 
    ! save data
    data_scalar3D(ithet,izeta,i_s, X1__)   = X1_int 
    data_scalar3D(ithet,izeta,i_s, X2__)   = X2_int
    data_scalar3D(ithet,izeta,i_s, GZETA__) =(-1.0_wp)*G_int                             !sign change of zeta coordinate! 
    data_scalar3D(ithet,izeta,i_s, BSUPT__) = (-iota_int+dLAdzeta )*(-dPhids_int/sqrtG)  !iota,dzeta,dphids all change sign  
    data_scalar3D(ithet,izeta,i_s, BSUPZ__) =  (1.0_wp+dLAdthet )*(-dPhids_int/sqrtG)    !sign change of zeta coordinate! 

    data_vector3D(:,ithet,izeta,i_s,BFIELD__    )  = Bfield
    data_vector3D(:,ithet,izeta,i_s,ECOV_S__    )  = e_s   
    data_vector3D(:,ithet,izeta,i_s,ECOV_THETA__)  = e_thet
    data_vector3D(:,ithet,izeta,i_s,ECOV_ZETA__ )  = -e_zeta  !sign change of zeta coordinate! 
    !========== 

  END DO ; END DO !izeta,ithet
  Favg_int = Favg_int/REAL((Nthet_out*Nzeta_out),wp)

  Itor_int = Itor_int*REAL(nfp_out)/(TWOPI*REAL((Nthet_out*Nzeta_out),wp))
  Ipol_int = Ipol_int*REAL(nfp_out)/(TWOPI*REAL((Nthet_out*Nzeta_out),wp))


  !========== 
  !save data
  data_1D( PHI__     ,i_s) = -TWOPI*Phi_int      !sign change of zeta coordinate! 
  data_1D( DPHIDS__  ,i_s) = -TWOPI*dPhids_int   !sign change of zeta coordinate!
  data_1D( CHI__     ,i_s) = TWOPI*Chi_int     
  data_1D( DCHIDS__  ,i_s) = TWOPI*dChids_int  
  data_1D( IOTA__    ,i_s) = -iota_int           !sign change of zeta coordinate!
  data_1D( PRESSURE__,i_s) = pressure_int
  data_1D( ITOR__    ,i_s) = -Itor_int*PI*1.0e+7_wp  !*(2pi)^2/mu_0 = pi*10^7, sign change due to LHS output coordinate system**
  data_1D( IPOL__    ,i_s) =  Ipol_int*PI*1.0e+7_wp  !*(2pi)^2/mu_0 , two sign changes, one from zeta, one from LHS coords.**
  data_1D( FAVG__    ,i_s) = Favg_int    
  data_1D( FMIN__    ,i_s) = Fmin_int    
  data_1D( FMAX__    ,i_s) = Fmax_int    
  !========== 
  ! ** LHS coordinate system u=theta, v=-zeta. covariant field representation with the currents:
  ! B=I_tor(s,u,v) grad u + I_pol(s,u,v) grad v 
  !  = I_tor(s,u,v) (-drvec/dv x drvec/ds)*sqrtg + Ipol(s,u,v) (-drvec/ds x drvec/du)*sqrtg
  !=> I_tor(s):  int_0^2pi B. drvec/du du = int_0^2pi Itor(s,u,0) (-1) du = - Itor(s) = - int_0^2pi B.drvec/dtheta dtheta
  !=> I_pol(s):  int_0^2pi B. drvec/dv dv = int_0^2pi Ipol(s,u,0) (-1) dv = - Ipol(s) = + int_0^2pi B.drvec/dzeta dzeta
  !   
  !  magnetic field in contra-variant components: (Phi=toroidal flux, chi'=iota*phi')
  ! B = Phi'/sqrtg drvec/dv + chi'/sqrtg drvec/du 
  ! toroidal flux:  int_0^s int_0^pi B. grad v ds du  = int_0^s int_0^pi B. (-drvec/ds x drvec/du)*sqrtg ds du    
  !                = -int_0^s Phi' ds  = -Phi(s)

END DO !i_s=1,Ns_out 
PhiEdge=data_1D(PHI__,Ns_out)
ChiEdge=data_1D(CHI__,Ns_out)

SWRITE(UNIT_stdOut,'(A)')'... DONE'
SWRITE(UNIT_stdOut,fmt_sep)
END SUBROUTINE gvec_to_castor3d_prepare


!===================================================================================================================================
!> write data to file
!!
!===================================================================================================================================
SUBROUTINE gvec_to_castor3d_writeToFile_ASCII()
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MODgvec_gvec_to_castor3d_Vars 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255) :: fileNameOut !< name of output file
INTEGER            :: ioUnit,iVar,i_s
!===================================================================================================================================
  FileNameOut='gvec2castor3d_'//TRIM(FileName)
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITING NEW CASTOR3D FILE    "'//TRIM(FileNameOut)//'" ...'
  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit       ,&
     FILE     = TRIM(FileNameOut) ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' ) 

!HEADER
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## GVEC-TO-CASTOR3D file, VERSION: 2.0                                                              '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## data is written on equidistant points in s,theta,zeta coordinates,                               '
  WRITE(ioUnit,'(A100)')'## * radially outward coordinate s=sqrt(phi_tor/phi_tor_edge) in [0,1]                              '
  WRITE(ioUnit,'(A100)')'##   s(1:Ns) , with  s(1)=0, s(Ns)=1                                                                '
  WRITE(ioUnit,'(A100)')'## * poloidal angle theta in [0,2pi] , sign: theta ~ atan(z/sqrt(x^2+y^2))                          '
  WRITE(ioUnit,'(A100)')'##   theta(1:Ntheta)  with theta(1)=0, theta(Ntheta)=2pi*(Ntheta-1)*/Ntheta                         '
  WRITE(ioUnit,'(A100)')'## * toroidal angle zeta in [0,2pi/nfp], sign: zeta ~ atan(y/x)  (opposite to GVEC definition!)     '
  WRITE(ioUnit,'(A100)')'##   zeta(1:Nzeta)  with zeta(1)=0, zeta(Nzeta)=2pi/nfp*(Nzeta-1)*/Nzeta                            '
  WRITE(ioUnit,'(A100)')'## * Angular coordinates can represent GVEC coordinates, with are not SFL (straight-field line)     '
  WRITE(ioUnit,'(A100)')'##   coordinates or can be SFL coordinates, either PEST or BOOZER. See global parameter "SFLcoord"  '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## Global variables:                                                                                '
  WRITE(ioUnit,'(A100)')'## * SFLcoord    : =0: GVEC coords (not SFL), =1: PEST SFL coords. , =2: BOOZER SFL coords.         '
  WRITE(ioUnit,'(A100)')'## * nfp         : number of toroidal field periods (toroidal angle [0,2pi/nfp])                    '
  WRITE(ioUnit,'(A100)')'## * asym        :  =0: symmetric cofiguration (R~cos,Z~sin), 1: asymmetric                         '
  WRITE(ioUnit,'(A100)')'## * m_max       : maximum number of poloidal modes in R,Z,lambda variables                         '
  WRITE(ioUnit,'(A100)')'## * n_max       : maximum number of toroidal modes in R,Z,lambda variables                         '
  WRITE(ioUnit,'(A100)')'## * PhiEdge     : Toroidal Flux at the last flux surface [T*m^2]                                   '
  WRITE(ioUnit,'(A100)')'## * ChiEdge     : Poloidal Flux at the last flux surface [T*m^2]                                   '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## Variables of arrays of radial variables (1:Ns)                                                   '
  WRITE(ioUnit,'(A100)')'## * Phi(s)      : toroidal magnetic flux [T*m^2]                                                   '
  WRITE(ioUnit,'(A100)')'## * dPhi_ds(s)  : derivative of toroidal magnetic flux versus s coordinate                         '
  WRITE(ioUnit,'(A100)')'## * Chi(s)      : poloidal magnetic flux [T*m^2]                                                   '
  WRITE(ioUnit,'(A100)')'## * dChi_ds(s)  : derivative of poloidal magnetic flux versus s coordinate                         '
  WRITE(ioUnit,'(A100)')'## * iota(s)     : rotational transform, dChi_ds/dPhi_ds [-]                                        '
  WRITE(ioUnit,'(A100)')'## * pressure(s) : pressure profile [N/(m^2)]                                                       '
  WRITE(ioUnit,'(A100)')'## * Itor(s)     : toroidal current profile [A]   = int_0^2pi B_theta dtheta                       ' 
  WRITE(ioUnit,'(A100)')'## * Ipol(s)     : poloidal current profile [A]   = int_0^2pi B_zeta dzeta                         ' 
  WRITE(ioUnit,'(A100)')'## * Favg(s)     : For tokamaks, poloidal current profile (input GS solver)                         '
  WRITE(ioUnit,'(A100)')'## * Fmin/Fmax(s): minimum /maximum of F(s)                                                         ' 
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## 3D arrays of scalars (1:Ntheta,1:Nzeta,1:Ns)                                                     '
  WRITE(ioUnit,'(A100)')'## * X1 (R)      : coordinate R=sqrt(x^2+y^2) ( called X1 in GVEC, only=R for hmap=1)               '
  WRITE(ioUnit,'(A100)')'## * X2 (Z)      : coordinate Z=z ( called X2 in GVEC, only=Z for hmap=1)                           '
  WRITE(ioUnit,'(A100)')'## * Gzeta       : map to geometric toroidal angle phi = zeta + Gzeta ( spans RHS coords {R,phi,Z} )'
  WRITE(ioUnit,'(A100)')'## * B^theta     : theta component of magnetic field, B^theta = B. grad(theta)                      '
  WRITE(ioUnit,'(A100)')'## * B^zeta      :  zeta component of magnetic field, B^zeta  = B. grad(zeta)                       '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## 3D arrays of vectors (1:3,1:Ntheta,1:Nzeta,1:Ns)                                                 '
  WRITE(ioUnit,'(A100)')'## * Bfield      : vector of magnetic field, components in cartesian coordinates (x,y,z)            '
  WRITE(ioUnit,'(A100)')'## * ecov_s      : covariant vector in s, components in cartesian coordinates (x,y,z)               '
  WRITE(ioUnit,'(A100)')'## * ecov_theta  : covariant vector in theta, components in cartesian coordinates (x,y,z)           '
  WRITE(ioUnit,'(A100)')'## * ecov_zeta   : covariant vector in zeta, components in cartesian coordinates (x,y,z)            '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'####################################################################################################'
  WRITE(ioUnit,*)
  WRITE(ioUnit,'(A)')'##<< number of grid points: 1:Ns (radial), 1:Ntheta (poloidal),1:Nzeta (toroidal) '
  WRITE(ioUnit,'(*(I8,:,1X))')Ns_out,Nthet_out,Nzeta_out
  WRITE(ioUnit,'(A)')'##<< global: SFLcoord,nfp,  asym, m_max, n_max'
  WRITE(ioUnit,'(12X,*(I6,:,1X))')SFLcoord,nfp_out,asym_out,mn_max_out(1:2)
  WRITE(ioUnit,'(A)')'##<< global: PhiEdge, ChiEdge'
  WRITE(ioUnit,'(*(E23.15,:,1X))')PhiEdge,ChiEdge
  WRITE(ioUnit,'(A,I4,A)',ADVANCE='NO')'##<< 1D profiles (',nVar1D,',1:Ns), variable names : '
  DO iVar=1,nVar1D-1
    WRITE(ioUnit,'(A,1X,(1X))',ADVANCE='NO' )  '"'//TRIM(StrVarNames1D(iVar))//'", '
  END DO
  WRITE(ioUnit,'(A)')  '"'//TRIM(StrVarNames1D(nVar1D))//'"'
  DO i_s=1,Ns_out
    WRITE(ioUnit,'(*(e23.15,:,1X))') data_1d(1:nVar1D,i_s) 
  END DO
  DO iVar=1,nVarScalar3D
    WRITE(ioUnit,'(A)',ADVANCE='NO')'##<< 3D scalar variable (1:Ntheta,1:Nzeta,1:Ns), Variable name: '
    WRITE(ioUNIT,'(A)')' "'//TRIM(StrVarNamesScalar3D(iVar))//'"'
    WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar3D(1:Nthet_out,1:Nzeta_out,1:Ns_out,iVar) 
  END DO !iVar=1,nVarScalar3D
  DO iVar=1,nVarVector3D
    WRITE(ioUnit,'(A)',ADVANCE='NO')'##<< 3D vector variable, cartesian components (1:3,1:Ntheta,1:Nzeta,1:Ns),Variable name: '
    WRITE(ioUNIT,'(A)')' "'//TRIM(StrVarNamesVector3D(iVar))//'"'
    WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_vector3d(1:3,1:Nthet_out,1:Nzeta_out,1:Ns_out,iVar) 
  END DO

  

  CLOSE(ioUnit)
  !write profiles again in csv
  OPEN(UNIT     = ioUnit       ,&
     FILE     = 'profiles.csv' ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' ) 
  DO iVar=1,nVar1D-1
    WRITE(ioUnit,'(A,1X,(1X))',ADVANCE='NO' )  '"'//TRIM(StrVarNames1D(iVar))//'", '
  END DO
  WRITE(ioUnit,'(A)')  '"'//TRIM(StrVarNames1D(nVar1D))//'"'
  DO i_s=1,Ns_out
    WRITE(ioUnit,'(*(e23.15,:,","))') data_1d(1:nVar1D,i_s) 
  END DO
  CLOSE(ioUnit)
  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE gvec_to_castor3d_writeToFile_ASCII


!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE finalize_gvec_to_castor3d 
! MODULES
USE MODgvec_gvec_to_castor3d_Vars 
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
  SDEALLOCATE(data_1D) 
  SDEALLOCATE(data_scalar3D) 
  SDEALLOCATE(data_vector3D) 

END SUBROUTINE finalize_gvec_to_castor3d

END MODULE MODgvec_gvec_to_castor3d
