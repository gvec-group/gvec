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
USE MODgvec_Globals, ONLY:wp
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
SUBROUTINE init_gvec_to_castor3d(fileName,Nfactor) 
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,fmt_sep
USE MODgvec_Globals,ONLY: TWOPI
USE MODgvec_ReadState,ONLY: ReadState
USE MODgvec_ReadState_vars,ONLY: sgrid_r,X1_base_r,X2_base_r,LA_base_r
USE MODgvec_gvec_to_castor3d_vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*), INTENT(IN) :: fileName !< name of GVEC file
INTEGER         , INTENT(IN) :: Nfactor(3) !< factor for s,theta,zeta resolution 
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iElem,i
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(A)')'INIT GVEC-TO-CASTOR3D ...'

  CALL ReadState(fileName)

  mn_max_out(1)    = MAXVAL((/X1_base_r%f%mn_max(1),X2_base_r%f%mn_max(1),LA_base_r%f%mn_max(1)/))
  mn_max_out(2)    = MAXVAL((/X1_base_r%f%mn_max(2),X2_base_r%f%mn_max(2),LA_base_r%f%mn_max(2)/))
  nfp_out          = X1_base_r%f%nfp

  Ns_out    = Nfactor(1)*sgrid_r%nElems +1
  ALLOCATE(s_pos(Ns_out))
  ALLOCATE(data_1D(nVar1D,Ns_out))

  s_pos(1)=1.0e-08_wp !avoid axis
  DO iElem=1,sgrid_r%nElems
    DO i=1,Nfactor(1)
      s_pos((iElem-1)*Nfactor(1)+1+i) = sgrid_r%sp(iElem-1)+REAL(i,wp)/REAL(Nfactor(1),wp)*sgrid_r%ds(iElem)
    END DO !i=1,Nfactor(1)
  END DO !iElem=1,sgrid_r%nElems
  s_pos(Ns_out)=1. - 1.0e-12_wp !aviod edge
  Nthet_out = Nfactor(2)*mn_max_out(1)
  Nzeta_out = MAX(1,Nfactor(3)*mn_max_out(2)) !if n=0, output 1 point

  ALLOCATE(thet_pos(Nthet_out))
  ALLOCATE(zeta_pos(Nzeta_out))
  DO i=1,Nthet_out
    thet_pos(i)=(TWOPI*REAL((i-1),wp))/REAL(Nthet_out) 
  END DO
  DO i=1,Nzeta_out
    zeta_pos(i)=(TWOPI*REAL((i-1),wp))/REAL((Nzeta_out*nfp_out),wp)
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
USE MODgvec_Globals,        ONLY: CROSS
USE MODgvec_ReadState_Vars, ONLY: profiles_1d,hmap_r,X1_base_r,X2_base_r,LA_base_r,X1_r,X2_r,LA_r
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: i_s,iMode,ithet,izeta
REAL(wp)                                :: spos,xp(2),F_loc,qvec(3),sqrtG
REAL(wp)                                :: grad_zeta(3) 
REAL(wp)                                :: dX1ds,dX1dthet,dX1dzeta
REAL(wp)                                :: dX2ds,dX2dthet,dX2dzeta
REAL(wp)                                :: dLAdthet,dLAdzeta
REAL(wp),DIMENSION(1:X1_base_r%f%modes) :: X1_s,dX1ds_s
REAL(wp),DIMENSION(1:X2_base_r%f%modes) :: X2_s,dX2ds_s
REAL(wp),DIMENSION(1:LA_base_r%f%modes) :: LA_s
!===================================================================================================================================

DO i_s=1,Ns_out
  spos          = s_pos(i_s)
  data_1D(SPOS__,i_s)=s_pos(i_s)
  ASSOCIATE( Phi_int     => data_1D( PHI__     ,i_s) &
            ,dPhids_int  => data_1D( DPHIDS__  ,i_s) &
            ,Chi_int     => data_1D( CHI__     ,i_s) &
            ,dChids_int  => data_1D( DCHIDS__  ,i_s) &
            ,iota_int    => data_1D( IOTA__    ,i_s) &
            ,pressure_int=> data_1D( PRESSURE__,i_s) &
            ,Favg_int    => data_1D( FAVG__    ,i_s) &
            ,Fmin_int    => data_1D( FMIN__    ,i_s) &
            ,Fmax_int    => data_1D( FMAX__    ,i_s) &
           )

  Phi_int     = X1_base_r%s%evalDOF_s(spos,       0 ,profiles_1d(:,1))
  dPhids_int  = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
  Chi_int     = X1_base_r%s%evalDOF_s(spos,       0 ,profiles_1d(:,2)) !Chi not yet working
  dChids_int  = X1_base_r%s%evalDOF_s(spos, DERIV_S ,profiles_1d(:,2)) !Chi not yet working
  iota_int    = X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,3))
  pressure_int= X1_base_r%s%evalDOF_s(spos, 0,profiles_1d(:,4))

  
  Fmin_int=+1.0e12
  Fmax_int=-1.0e12
  Favg_int = 0.

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

    ASSOCIATE( X1_int => data_scalar3D(ithet,izeta,i_s, X1__) & 
              ,X2_int => data_scalar3D(ithet,izeta,i_s, X2__) & 
              ,LA_int => data_scalar3D(ithet,izeta,i_s, LA__) &
              ,e_s    => data_vector3D(:,ithet,izeta,i_s,ECOV_S__     )  &
              ,e_thet => data_vector3D(:,ithet,izeta,i_s,ECOV_THETA__ )  &
              ,e_zeta => data_vector3D(:,ithet,izeta,i_s,ECOV_ZETA__  )  &
              ,Bfield => data_vector3D(:,ithet,izeta,i_s,BFIELD__     )  &
             )

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


    END ASSOCIATE !data_3D,ecov*
  END DO ; END DO !izeta,ithet
  Favg_int = Favg_int/REAL((Nthet_out*Nzeta_out),wp)
  END ASSOCIATE ! data_1D
END DO !i_s=1,Ns_out 

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

  WRITE(ioUnit,'(A)')'## GVEC-TO-CASTOR3D file, VERSION: 1.0'
  WRITE(ioUnit,'(A)')'##<< number of grid points: 1:Ns (radial), 1:Ntheta (poloidal),1:Nzeta (toroidal) ###########################'
  WRITE(ioUnit,'(*(I8,:,1X))')Ns_out,Nthet_out,Nzeta_out
  WRITE(ioUnit,'(A)')'##<< global: nfp,m_max,n_max       ##########################################################################'
  WRITE(ioUnit,'(*(I8,:,1X))')nfp_out,mn_max_out(1:2)
  WRITE(ioUnit,'(A,I4,A)')'##<< 1D profiles (',nVar1D,',1:Ns), variable names :  ##############################################'
  WRITE(ioUNIT,'(A)',ADVANCE='NO') '##  '
  DO iVar=1,nVar1D-1
    WRITE(ioUnit,'(A,1X,(1X))',ADVANCE='NO' )  '"'//TRIM(StrVarNames1D(iVar))//'"'
  END DO
  WRITE(ioUnit,'(A)')  '"'//TRIM(StrVarNames1D(nVar1D))//'"'
  DO i_s=1,Ns_out
    WRITE(ioUnit,'(*(e23.15,:,1X))') data_1d(1:nVar1D,i_s) 
  END DO
  DO iVar=1,nVarScalar3D
    WRITE(ioUnit,'(A)')'##<< 3D scalar variable (1:Ntheta,1:Nzeta,1:Ns), Variable name:  ##########################################'
    WRITE(ioUNIT,'(A)')'##  "'//TRIM(StrVarNamesScalar3D(iVar))//'"'
    WRITE(ioUnit,'(*(6(e23.15,:,1X),/))') data_scalar3D(1:Nthet_out,1:Nzeta_out,1:Ns_out,iVar) 
  END DO !iVar=1,nVarScalar3D
  DO iVar=1,nVarVector3D
    WRITE(ioUnit,'(A)')'##<< 3D vector variable, cartesian components (1:3,1:Ntheta,1:Nzeta,1:Ns),Variable name: ##################'
    WRITE(ioUNIT,'(A)')'##  "'//TRIM(StrVarNamesVector3D(iVar))//'"'
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
