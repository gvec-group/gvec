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

PUBLIC::init_gvec_to_castor3d
PUBLIC::gvec_to_castor3d_prepare
PUBLIC::gvec_to_castor3d_writeToFile
PUBLIC::finalize_gvec_to_castor3d

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize Module 
!!
!===================================================================================================================================
SUBROUTINE init_gvec_to_castor3d(fileName,Ns_in,factorFourier) 
! MODULES
USE MODgvec_Globals,ONLY: TWOPI
USE MODgvec_ReadState,ONLY: ReadState
USE MODgvec_ReadState_vars,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_gvec_to_castor3d_vars
USE MODgvec_ReadState_vars,ONLY: LA_r,X1_r 
USE MODgvec_transform_sfl,ONLY:test_sfl
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileName !< name of GVEC file
INTEGER         , INTENT(IN) :: Ns_in !< number of equidistant points in radial s-direction (includes axis and edge!)
INTEGER         , INTENT(IN) :: factorFourier !< factor theta,zeta resolution Ntheta=Factor*m_max, Nzeta=MAX(1,Factor*n_max)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: i
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT GVEC-TO-CASTOR3D ...'

  CALL ReadState(fileName)

  CALL test_sfl(LA_base_r,LA_r,X1_base_r,X1_r)

  mn_max_out(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
  mn_max_out(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))
  nfp_out          = X1_base_r%f%nfp

  Ns_out    = Ns_in
  IF((X1_base_r%f%sin_cos.EQ._COS_).AND.(X2_base_r%f%sin_cos.EQ._SIN_).AND.(LA_base_r%f%sin_cos.EQ._SIN_))THEN
    asym_out = 0 !R~cos,Z~sin,lambda~sin
  ELSE
    asym_out = 1 !full fourier
  END IF
  ALLOCATE(s_pos(Ns_out))
  ALLOCATE(data_1D(nVar1D,Ns_out))

  s_pos(1)=1.0e-08_wp !avoid axis
  DO i=2,Ns_out-1
      s_pos(i) = REAL(i-1,wp)/REAL(Ns_out-1,wp)
  END DO !i
  s_pos(Ns_out)=1. - 1.0e-12_wp !aviod edge
  Nthet_out = factorFourier*mn_max_out(1)
  Nzeta_out = MAX(1,factorFourier*mn_max_out(2)) !if n=0, output 1 point

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

END SUBROUTINE init_gvec_to_castor3d


!===================================================================================================================================
!> prepare all data to be written
!!
!===================================================================================================================================
SUBROUTINE gvec_to_castor3d_prepare()
! MODULES
USE MODgvec_gvec_to_castor3d_Vars 
USE MODgvec_Globals,        ONLY: CROSS,TWOPI
USE MODgvec_ReadState_Vars, ONLY: profiles_1d,hmap_r,X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
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
REAL(wp)                                :: X1_int,X2_int,LA_int
REAL(wp)                                :: Ipol_int,Itor_int
REAL(wp),DIMENSION(3)                   :: qvec,e_s,e_thet,e_zeta,Bfield,grad_zeta
REAL(wp),DIMENSION(1:X1_base_r%f%modes) :: X1_s,dX1ds_s
REAL(wp),DIMENSION(1:X2_base_r%f%modes) :: X2_s,dX2ds_s
REAL(wp),DIMENSION(1:LA_base_r%f%modes) :: LA_s
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'PREPARE DATA FOR GVEC-TO-CASTOR3D ...'
DO i_s=1,Ns_out
  IF(MOD(i_s,MAX(1,Ns_out/100)).EQ.0) THEN
    SWRITE(UNIT_stdOut,'(4X,I4,A4,I4,A13,A1)',ADVANCE='NO')i_s, ' of ',NS_out,' evaluated...',ACHAR(13)
  END IF
  spos          = s_pos(i_s)
  data_1D(SPOS__,i_s)=s_pos(i_s)

  Phi_int     = X1_base_r%s%evalDOF_s(spos,       0 ,profiles_1d(:,1))
  dPhids_int  = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  Chi_int     = X1_base_r%s%evalDOF_s(spos,       0 ,profiles_1d(:,2)) !Chi not yet working
  dChids_int  = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2)) !Chi not yet working
  iota_int    = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))
  pressure_int= X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,4))

  
  Fmin_int=+1.0e12
  Fmax_int=-1.0e12
  Favg_int = 0.
  Itor_int = 0.
  Ipol_int = 0.

  !interpolate radially
  DO iMode=1,X1_base_r%f%modes
    X1_s( iMode)  = X1_base_r%s%evalDOF_s(spos,      0,X1_r(:,iMode))
    dX1ds_s(iMode)= X1_base_r%s%evalDOF_s(spos,DERIV_S,X1_r(:,iMode))
  END DO !iMode
  DO iMode=1,X2_base_r%f%modes
    X2_s(   iMode)= X1_base_r%s%evalDOF_s(spos,      0,X2_r(:,iMode))
    dX2ds_s(iMode)= X1_base_r%s%evalDOF_s(spos,DERIV_S,X2_r(:,iMode))
  END DO !iMode
  DO iMode=1,LA_base_r%f%modes
    LA_s(iMode)   = LA_base_r%s%evalDOF_s(spos,      0,LA_r(:,iMode))
  END DO !iMode
  !interpolate in the angles
  DO izeta=1,Nzeta_out; DO ithet=1,Nthet_out
    xp=(/thet_pos(ithet),zeta_pos(izeta)/)


    X1_int   = X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
    dX1ds    = X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
    dX1dthet = X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
    dX1dzeta = X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )
    
    X2_int   = X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
    dX2ds    = X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
    dX2dthet = X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
    dX2dzeta = X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )
    
    LA_int   = LA_base_r%f%evalDOF_x(xp,          0, LA_s)
    dLAdthet = LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
    dLAdzeta = LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)
    
    qvec=(/X1_int,X2_int,xp(2)/)
    e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,0.0_wp/)) !dxvec/ds
    e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/)) !dxvec/dthet
    e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/)) !dxvec/dzeta
    sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))
   
    !check: J = e_s*(e_thet x e_zeta) 
    !sqrtG_check = hmap_r%eval_Jh(qvec)*(dX1ds*dX2dthet -dX2ds*dX1dthet) 
    !WRITE(*,*)'CHECK sqrtG',sqrtG,sqrtG_check,sqrtG-sqrtG_check

    Bfield(:) = (  e_thet(:)*(iota_int-dLAdzeta )  & 
                 + e_zeta(:)*(1.0_wp+dLAdthet   ) )*(dPhids_int/sqrtG)

    !grad_s(:)   = CROSS(e_thet,e_zeta)/sqrtG
    !grad_thet(:)= CROSS(e_zeta,e_s   )/sqrtG
    grad_zeta(:)= CROSS(e_s   ,e_thet)/sqrtG


    !F-profile, only makes sense for tokamak configurations (n=0): F=-phi' R*(1+dlambda_dtheta)/(Jac*|\nabla\zeta|)
    F_loc = -dPhids_int*X1_int*(1+dLAdthet)/(sqrtG*SQRT(SUM(grad_zeta(:)**2)))

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
    data_scalar3D(ithet,izeta,i_s, LA__)   = LA_int

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
  data_1D( FAVG__    ,i_s) = Favg_int    
  data_1D( FMIN__    ,i_s) = Fmin_int    
  data_1D( FMAX__    ,i_s) = Fmax_int    
  data_1D( ITOR__    ,i_s) = Itor_int !/mu_0?   
  data_1D( IPOL__    ,i_s) =-Ipol_int !/mu_0?    !sign change of zeta coordinate!  
  !========== 
  
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
SUBROUTINE gvec_to_castor3d_writeToFile_ASCII(fileNameOut)
! MODULES
USE MODgvec_Globals,ONLY:Unit_stdOut,GETFREEUNIT
USE MODgvec_gvec_to_castor3d_Vars 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileNameOut !< name of output file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: ioUnit,iVar,i_s
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'   WRITING NEW CASTOR3D FILE    "'//TRIM(FileNameOut)//'" ...'
  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit       ,&
     FILE     = TRIM(FileNameOut) ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' ) 

!HEADER
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## GVEC-TO-CASTOR3D file, VERSION: 1.1                                                              '
  WRITE(ioUnit,'(A100)')'## -------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A100)')'## data is written on equidistant points in s,theta,zeta coordinates,                               '
  WRITE(ioUnit,'(A100)')'## * radially outward coordinate s=sqrt(phi_tor/phi_tor_edge) in [0,1]                              '
  WRITE(ioUnit,'(A100)')'##   s(1:Ns) , with  s(1)=0, s(Ns)=1                                                                '
  WRITE(ioUnit,'(A100)')'## * poloidal angle theta in [0,2pi] , sign: theta ~ atan(z/sqrt(x^2+y^2))                          '
  WRITE(ioUnit,'(A100)')'##   theta(1:Ntheta)  with theta(1)=0, theta(Ntheta)=2pi*(Ntheta-1)*/Ntheta                         '
  WRITE(ioUnit,'(A100)')'## * toroidal angle zeta in [0,2pi/nfp], sign: zeta ~ atan(y/x)  (opposite to GVEC definition!)     '
  WRITE(ioUnit,'(A100)')'##   zeta(1:Nzeta)  with zeta(1)=0, zeta(Nzeta)=2pi/nfp*(Nzeta-1)*/Nzeta                            '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## Global variables:                                                                                '
  WRITE(ioUnit,'(A100)')'## * nfp         : number of toroidal field periods (toroidal angle [0,2pi/nfp])                    '
  WRITE(ioUnit,'(A100)')'## * asym        :  =0: symmetric cofiguration (R~cos,Z~sin,lambda~sin), 1: asymmetric              '
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
  WRITE(ioUnit,'(A100)')'## * Favg(s)     : For tokamaks, poloidal current profile (input GS solver)                         '
  WRITE(ioUnit,'(A100)')'## * Fmin/Fmax(s): minimum /maximum of F(s)                                                         ' 
  WRITE(ioUnit,'(A100)')'## * Itor(s)     : toroidal current profile [??]   = int_0^2pi B_theta dtheta                       ' 
  WRITE(ioUnit,'(A100)')'## * Ipol(s)     : poloidal current profile [??]   = int_0^2pi B_zeta dzeta                         ' 
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## 3D arrays of scalars (1:Ntheta,1:Nzeta,1:Ns)                                                     '
  WRITE(ioUnit,'(A100)')'## * X1 (R)      : coordinate R=sqrt(x^2+y^2) ( called X1 in GVEC, only=R for hmap=1)               '
  WRITE(ioUnit,'(A100)')'## * X2 (Z)      : coordinate Z=z ( called X2 in GVEC, only=Z for hmap=1)                           '
  WRITE(ioUnit,'(A100)')'## * LA          : lambda variable, straight-field line angle theta*=theta+lambda (PEST coordinates)'
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'## 3D arrays of vectors (1:3,1:Ntheta,1:Nzeta,1:Ns)                                                 '
  WRITE(ioUnit,'(A100)')'## * Bfield      : vector of magneitc field, components in cartesian coordinates (x,y,z)            '
  WRITE(ioUnit,'(A100)')'## * ecov_s      : covariant vector in s, components in cartesian coordinates (x,y,z)               '
  WRITE(ioUnit,'(A100)')'## * ecov_theta  : covariant vector in theta, components in cartesian coordinates (x,y,z)           '
  WRITE(ioUnit,'(A100)')'## * ecov_zeta   : covariant vector in zeta, components in cartesian coordinates (x,y,z)            '
  WRITE(ioUnit,'(A100)')'##                                                                                                  '
  WRITE(ioUnit,'(A100)')'####################################################################################################'
  WRITE(ioUnit,*)
  WRITE(ioUnit,'(A)')'##<< number of grid points: 1:Ns (radial), 1:Ntheta (poloidal),1:Nzeta (toroidal) '
  WRITE(ioUnit,'(*(I8,:,1X))')Ns_out,Nthet_out,Nzeta_out
  WRITE(ioUnit,'(A)')'##<< global: nfp, asym, m_max, n_max'
  WRITE(ioUnit,'(*(I8,:,1X))')nfp_out,asym_out,mn_max_out(1:2)
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
