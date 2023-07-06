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
!!# Module **Read in of State file**
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_ReadState
! MODULES
USE MODgvec_Globals, ONLY:wp,GETFREEUNIT,abort
IMPLICIT NONE
PRIVATE

INTERFACE ReadState
  MODULE PROCEDURE ReadStateFileFromASCII
END INTERFACE

INTERFACE Finalize_ReadState
  MODULE PROCEDURE Finalize_ReadState
END INTERFACE

PUBLIC::ReadState
PUBLIC::Finalize_ReadState
!===================================================================================================================================

CONTAINS



!===================================================================================================================================
!> read an input solution and initialize U(0) (X1,X2,LA) of size X1/X2/LA_base , from an ascii .dat file 
!! if size of grid/X1/X2/LA  not equal X1/X2/X3_base
!! interpolate readin solution to the current base of Uin
!!
!===================================================================================================================================
SUBROUTINE ReadStateFileFromASCII(fileString)
! MODULES
USE MODgvec_ReadState_Vars
USE MODgvec_Globals,ONLY: Unit_stdOut,GETFREEUNIT,abort
USE MODgvec_sgrid,  ONLY: t_sgrid
USE MODgvec_base,   ONLY: t_base, base_new
USE MODgvec_sbase,  ONLY: sbase_new
USE MODgvec_fbase,  ONLY: sin_cos_map 
USE MODgvec_hmap,  ONLY: hmap_new
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=*)    , INTENT(IN   ) :: fileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL              :: file_exists
  INTEGER              :: ioUnit,iMode,is,nElems_r,grid_type_r,nfp_r,degGP_r,mn_nyq_r(2),which_hmap_r
  INTEGER              :: X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  INTEGER              :: X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  INTEGER              :: LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r
  INTEGER,ALLOCATABLE  :: X1_mn_r(:,:),X2_mn_r(:,:),LA_mn_r(:,:)
  REAL(wp),ALLOCATABLE :: sp_r(:),profiles_IP(:,:)
  INTEGER              :: X1_mn_max_r(2),X2_mn_max_r(2),LA_mn_max_r(2)

  INTEGER              :: iGP
  REAL                 :: chi_int, iota_int, phi_edge,  ds
  REAL(wp),ALLOCATABLE :: w_GP(:),s_GP(:), chi_IP(:)
  
!===================================================================================================================================
  WRITE(UNIT_stdOut,'(A)')'   READ STATEFILE    "'//TRIM(FileString)//'" ...'

  INQUIRE(FILE=TRIM(FileString), EXIST=file_exists)
  
  IF(.NOT.file_exists) CALL abort(__STAMP__, &
        TRIM("STATEFILE: "//TRIM(FileString)//" DOES NOT EXIST!!"))
  ioUnit=GETFREEUNIT()
  OPEN(UNIT     = ioUnit         ,&
     FILE     = TRIM(FileString) ,&
     STATUS   = 'OLD'            ,&
     ACTION   = 'READ'           ,&
     ACCESS   = 'SEQUENTIAL' ) 

  READ(ioUnit,*) !## MHD3D Solution file
  READ(ioUnit,*) outputLevel_r,fileID_r
  READ(ioUnit,*) !## grid: nElems, grid_type
  READ(ioUnit,*) nElems_r,grid_type_r
  ALLOCATE(sp_r(0:nElems_r))

  READ(ioUnit,*) !## grid: sp(0:nElems)
  READ(ioUnit,*)sp_r(:)
  READ(ioUnit,*) !## global: nfp, degGP, mn_nyq
  READ(ioUnit,*) nfp_r, degGP_r,mn_nyq_r,which_hmap_r
  READ(ioUnit,*) !## X1_base: 
  READ(ioUnit,*) X1_nBase_r,X1_deg_r,X1_cont_r,X1_modes_r,X1_sin_cos_r,X1_excl_mn_zero_r
  READ(ioUnit,*) !## X2_base:                 
  READ(ioUnit,*) X2_nBase_r,X2_deg_r,X2_cont_r,X2_modes_r,X2_sin_cos_r,X2_excl_mn_zero_r
  READ(ioUnit,*) !## LA_base:                 
  READ(ioUnit,*) LA_nBase_r,LA_deg_r,LA_cont_r,LA_modes_r,LA_sin_cos_r,LA_excl_mn_zero_r

  ALLOCATE(X1_r(1:X1_nbase_r,1:X1_modes_r))
  ALLOCATE(X2_r(1:X2_nbase_r,1:X2_modes_r))
  ALLOCATE(LA_r(1:LA_nbase_r,1:LA_modes_r))

  ALLOCATE(profiles_IP(1:X1_nbase_r,5))
  ALLOCATE(X1_mn_r(2,1:X1_modes_r))
  ALLOCATE(X2_mn_r(2,1:X2_modes_r))
  ALLOCATE(LA_mn_r(2,1:LA_modes_r))
  READ(ioUnit,*) !## X1: 
  DO iMode=1,X1_modes_r
    READ(ioUnit,*)X1_mn_r(:,iMode),X1_r(:,iMode)
  END DO
  READ(ioUnit,*) !## X2: 
  DO iMode=1,X2_modes_r
    READ(ioUnit,*)X2_mn_r(:,iMode),X2_r(:,iMode)
  END DO
  READ(ioUnit,*) !## LA: 
  DO iMode=1,LA_modes_r
    READ(ioUnit,*)LA_mn_r(:,iMode),LA_r(:,iMode)
  END DO
  READ(ioUnit,*) !## profiles at X1_base IP points : spos,phi,chi,iota,pressure 
  DO is=1,X1_nbase_r
    READ(ioUnit,*)profiles_IP(is,:)
  END DO
  READ(ioUnit,*) !## a_minor,r_major,volume 
  READ(ioUnit,*)a_minor,r_major,volume

  CLOSE(ioUnit)

  CALL hmap_new(hmap_r,which_hmap_r)
  IF((hmap_r%nfp.NE.-1).AND.(hmap_r%nfp.NE.nfp_r)) CALL abort(__STAMP__,&
                        "nfp from restart file does not match to nfp used in hmap.")
  ! check if input has changed:

  CALL sgrid_r%init(nElems_r,grid_type_r)

  !needed to build base of restart file
  X1_mn_max_r = (/MAXVAL(X1_mn_r(1,:)),MAXVAL(X1_mn_r(2,:))/nfp_r/)
  X2_mn_max_r = (/MAXVAL(X2_mn_r(1,:)),MAXVAL(X2_mn_r(2,:))/nfp_r/)
  LA_mn_max_r = (/MAXVAL(LA_mn_r(1,:)),MAXVAL(LA_mn_r(2,:))/nfp_r/)

  CALL sbase_new(sbase_prof,X1_deg_r,X1_cont_r,sgrid_r,degGP_r) !up to now, same base as X1 is used for profiles!
  CALL base_new(X1_base_r,X1_deg_r,X1_cont_r,sgrid_r,degGP_r,X1_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(X1_sin_cos_r),(X1_excl_mn_zero_r.EQ.1))
  CALL base_new(X2_base_r,X2_deg_r,X2_cont_r,sgrid_r,degGP_r,X2_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(X2_sin_cos_r),(X2_excl_mn_zero_r.EQ.1))
  CALL base_new(LA_base_r,LA_deg_r,LA_cont_r,sgrid_r,degGP_r,LA_mn_max_r,mn_nyq_r,nfp_r, &
                sin_cos_map(LA_sin_cos_r),(LA_excl_mn_zero_r.EQ.1))

  ALLOCATE(profiles_1d(1:X1_nbase_r,4))
  !convert to spline DOF
  profiles_1d(:,1) =X1_base_r%s%initDOF( profiles_IP(:,1+1) ) !phi
  !profiles_1d(:,2) =X1_base_r%s%initDOF( profiles_IP(:,1+2) ) !chi, NOT SAVE TO BE USED
  profiles_1d(:,3) =X1_base_r%s%initDOF( profiles_IP(:,1+3) ) !iota
  profiles_1d(:,4) =X1_base_r%s%initDOF( profiles_IP(:,1+4) ) !pressure

  ! Calculate Chi profile by quadrature Chi(s) = \int{x=0}^s Phi' iota dx , Phi'(s)=Phi_edge*s^2
  ALLOCATE(s_GP(1:degGP_r+1),w_GP(1:degGP_r+1))
  ALLOCATE(chi_IP(1:X1_nbase_r))
  phi_edge = sbase_prof%evalDOF_s(1.0_wp, 0 ,profiles_1d(:,1))
  chi_IP(1)=0.0_wp
  chi_int = 0.0_wp
  DO is=2, X1_nbase_r
    ds = profiles_IP(is, 1) - profiles_IP(is-1,1)

    w_GP(1:(degGP_r+1))=                       0.5_wp * sbase_prof%w_GPloc(:)        * ds
    s_GP(1:(degGP_r+1))= profiles_IP(is-1,1) + 0.5_wp * (sbase_prof%xi_GP(:)+1.0_wp) * ds

    DO iGP=1, degGP_r+1
      iota_int  = sbase_prof%evalDOF_s(s_GP(iGP), 0 ,profiles_1d(:,3))
      chi_int = chi_int + iota_int * 2.0_wp * s_GP(iGP) * phi_edge * w_GP(iGP)
    END DO
    chi_IP(is) = chi_int
  END DO
  profiles_1d(:, 2) = X1_base_r%s%initDOF(chi_IP)

  DEALLOCATE(s_GP, w_GP, chi_IP)
  
  DEALLOCATE(sp_r,profiles_IP)
  DEALLOCATE(X1_mn_r)
  DEALLOCATE(X2_mn_r)
  DEALLOCATE(LA_mn_r)


  WRITE(UNIT_stdOut,'(A)')'...DONE.'
END SUBROUTINE ReadStateFileFromASCII 

!===================================================================================================================================
!> Finalize Module
!!
!===================================================================================================================================
SUBROUTINE Finalize_ReadState 
! MODULES
USE MODgvec_ReadState_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  DEALLOCATE(X1_r)
  DEALLOCATE(X2_r)
  DEALLOCATE(LA_r)
  DEALLOCATE(profiles_1d)
  CALL sgrid_r%free()
  CALL hmap_r%free()
  CALL sbase_prof%free()
  CALL X1_base_r%free()
  CALL X2_base_r%free()
  CALL LA_base_r%free()
  DEALLOCATE(hmap_r)
  DEALLOCATE(sbase_prof)
  DEALLOCATE(X1_base_r)
  DEALLOCATE(X2_base_r)
  DEALLOCATE(LA_base_r)
END SUBROUTINE Finalize_ReadState

END MODULE MODgvec_ReadState
