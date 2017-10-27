!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>

! This file is a modified version from HOPR (github.com/fhindenlang/hopr).
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

!===================================================================================================================================
!>
!!# Module **VMEC readin**
!!
!!
!===================================================================================================================================
MODULE MOD_VMEC_Readin

  USE MOD_Globals,ONLY:wp,twoPi
  IMPLICIT NONE
  PUBLIC

!> mu_0 (Permeability)
  REAL(wp), PARAMETER :: mu0   = 2.0_wp * twoPi * 1.0E-7_wp         ! V s / (A m)
!> version number
  REAL(wp), PARAMETER :: version = 3.75_wp

  INTEGER :: nFluxVMEC       !< number of flux surfaces in VMEC output
  INTEGER :: mn_mode         !< number of modes (theta,zeta)
! INTEGER :: mn_mode_nyq     !< number of modes (Nyquist, 2D)
  INTEGER :: nfp             !< number of field periods
  LOGICAL :: lasym=.FALSE.   !< if lasym=0, solution is (rmnc,zmns,lmns),
                             !< if lasym=1 solution is (rmnc,rmns,zmnc,zmns,lmnc,lmns)
  INTEGER :: mPol            !< poloidal mode number
  INTEGER :: nTor            !< toroidal mode number
  INTEGER :: signgs          !< signum of sqrtG
  INTEGER :: mnmax           !< maximum mode number over m,n
! INTEGER :: mnmax_nyq       !< maximum mode number over m_nyq,n_nyq
  REAL(wp) :: b0             !< magnetic field on axis B_0
  REAL(wp) :: rMajor         !< major radius
  REAL(wp) :: volume         !< volume

!! vector parameters
  REAL(wp), ALLOCATABLE :: xm(:)      !< poloidal mode numbers 
  REAL(wp), ALLOCATABLE :: xn(:)      !< toroidal mode numbers
! REAL(wp), ALLOCATABLE :: xm_nyq(:)  !< poloidal mode numbers (Nyquist)
! REAL(wp), ALLOCATABLE :: xn_nyq(:)  !< toroidal mode numbers (Nyquist) 
  REAL(wp), ALLOCATABLE :: iotaf(:)   !< iota on full mesh
  REAL(wp), ALLOCATABLE :: presf(:)   !< pressure on full mesh
  REAL(wp), ALLOCATABLE :: phi(:)     !< toroidal flux on full mesh
  REAL(wp), ALLOCATABLE :: chi(:)     !< poloidal flux on full mesh
  REAL(wp), ALLOCATABLE :: rmnc(:, :) !< R(iMode,iFlux) (cosine components on full mesh)
  REAL(wp), ALLOCATABLE :: rmns(:, :) 
  REAL(wp), ALLOCATABLE :: zmnc(:, :)
  REAL(wp), ALLOCATABLE :: zmns(:, :) !< Z(iMode,iFlux) (sine components on full mesh)
  REAL(wp), ALLOCATABLE :: lmnc(:, :)
  REAL(wp), ALLOCATABLE :: lmns(:, :) !< lambda(iMode,iFlux) (sine components (read on half mesh, needs to be interpolated on full mesh))

! REAL(wp), ALLOCATABLE :: gmnc(:, :) !< jacobian (cosine components (read on half mesh, interpolated on full
                                      !< mesh, mnMode_nyqist ))

INTERFACE ReadVMEC
  MODULE PROCEDURE ReadVMEC
END INTERFACE ReadVMEC

CONTAINS


!===================================================================================================================================
!> READ VMEC "wout" datafile, needs netcdf library!
!! 
!===================================================================================================================================
SUBROUTINE ReadVMEC(fileName)
!===================================================================================================================================
  USE MOD_VMEC_Vars
  IMPLICIT NONE
  INCLUDE "netcdf.inc"
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
  CHARACTER(LEN = *), INTENT(IN) :: fileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: aStat, aError, ioError, ncid, id
!===================================================================================================================================

    WRITE(*,'(4X,A)')'VMEC READ WOUT FILE...'
    !! open NetCDF input file
    ioError = NF_OPEN(TRIM(fileName), NF_NOWRITE, ncid)
    IF (ioError /= 0) THEN
      WRITE(*,*) " Cannot open ", TRIM(fileName), " in ReadVmecOutput!"
      CALL EXIT(2)
    END IF

    !! get array dimensions
    !! radial dimension
!    ioError = NF_INQ_DIMID(ncid, "radius", id)
!    ioError = ioError + NF_INQ_DIMLEN(ncid, id, nFluxVMEC)
    !! number of fourier components of r, z, lambda
    ioError = ioError + NF_INQ_DIMID(ncid, "mn_mode", id)
    ioError = ioError + NF_INQ_DIMLEN(ncid, id, mn_mode)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading mn_mode'
    !! number of fourier components of b_u, b_v, b_s
!    ioError = ioError + NF_INQ_DIMID(ncid, "mn_mode_nyq", id)
!    ioError = ioError + NF_INQ_DIMLEN(ncid, id, mn_mode_nyq)

    !! get number of field periods
    ioError = ioError + NF_INQ_VARID(ncid, "nfp", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, nfp)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading n-fp'
    !! get dimension of s
    ioError = ioError + NF_INQ_VARID(ncid, "ns", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, nFluxVMEC)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading n-s'
    !! get poloidal mode number
    ioError = ioError + NF_INQ_VARID(ncid, "mpol", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, mPol)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading mpol'
    !! get toroidal mode number
    ioError = ioError + NF_INQ_VARID(ncid, "ntor", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, nTor)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading ntor'
    !! get mnmax
    ioError = ioError + NF_INQ_VARID(ncid, "mnmax", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, mnmax)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading mnmax'
    !! get mnmax_nyq
!    ioError = ioError + NF_INQ_VARID(ncid, "mnmax_nyq", id)
!    ioError = ioError + NF_GET_VAR_INT(ncid, id, mnmax_nyq)
    !! get iasym
    ioError = ioError + NF_INQ_VARID(ncid, "lasym__logical__", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, lasym)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading lasym'
    !! get B_0
    ioError = ioError + NF_INQ_VARID(ncid, "b0", id)
    ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, b0)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading b0'
    !! get mnmax_nyq
    !! check the sign of b0
    IF (b0 < 0) THEN
      WRITE(*,'(4X,A)') "  VMEC run with b0 < 0 !!!"
    END IF
    !! get signgs
    ioError = ioError + NF_INQ_VARID(ncid, "signgs", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, signgs)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading signgs'
    !! get Rmajor_p
    ioError = ioError + NF_INQ_VARID(ncid, "Rmajor_p", id)
    ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, rMajor)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading Rmajor_p'
    !! get volume_p
    ioError = ioError + NF_INQ_VARID(ncid, "volume_p", id)
    ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, volume)
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading volume_p'

    IF (ioError /= 0) THEN
      WRITE(*,*) " Cannot read ", TRIM(fileName), "!"
      WRITE(*,*) " Possible wrong file format!"
      CALL EXIT(3)
    END IF

    !! allocate memory for fourier arrays
    aError = 0
    ALLOCATE(xm(mn_mode), xn(mn_mode), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(xm_nyq(mn_mode_nyq), xn_nyq(mn_mode_nyq), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(iotaf(nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(presf(nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(phi(nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(chi(nFluxVMEC), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(phipf(nFluxVMEC), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(rmnc(mn_mode, nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(zmns(mn_mode, nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(lmns(mn_mode, nFluxVMEC), stat = aStat)
    aError = aError + aStat
    IF(lasym)THEN
      ALLOCATE(rmns(mn_mode, nFluxVMEC), stat = aStat)
      aError = aError + aStat
      ALLOCATE(zmnc(mn_mode, nFluxVMEC), stat = aStat)
      aError = aError + aStat
      ALLOCATE(lmnc(mn_mode, nFluxVMEC), stat = aStat)
      aError = aError + aStat
    END IF

    IF (aError /= 0) THEN
      WRITE(*,*) "Allocation failure in subroutine ReadVmecOutput!"
      CALL EXIT(4)
    END IF

    !! read x_m
    ioError = NF_INQ_VARID(ncid, "xm", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /), (/ mn_mode /),&
         xm(:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading xm'
    !! read x_n
    ioError = NF_INQ_VARID(ncid, "xn", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /), (/ mn_mode /),&
         xn(:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading xn'
    !! read x_m^nyq
!    ioError = NF_INQ_VARID(ncid, "xm_nyq", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!         (/ mn_mode_nyq /), xm_nyq(:))
    !! read x_n^nyq
!    ioError = NF_INQ_VARID(ncid, "xn_nyq", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!         (/ mn_mode_nyq /), xn_nyq(:))
    !! read iotaf
    ioError = NF_INQ_VARID(ncid, "iotaf", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), iotaf(:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading iotaf'
    !! read presf
    ioError = NF_INQ_VARID(ncid, "presf", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), presf(:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading presf'
    !! read phi
    ioError = NF_INQ_VARID(ncid, "phi", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), phi(:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading phi'

    !! scale toroidal flux to get internal VMEC Phi
    phi(:nFluxVMEC) = REAL(signgs, wp) * phi(:nFluxVMEC) / TwoPi

    !! read chi
    ioError = NF_INQ_VARID(ncid, "chi", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), chi(:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading chi'

    !! scale poloidal flux to get internal VMEC chi, WITHOUT SIGNGS (so that iota=chi'/phi' >0 )
    chi(:nFluxVMEC) = chi(:nFluxVMEC) / TwoPi

    !! read phipf
!    ioError = NF_INQ_VARID(ncid, "phipf", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!         (/ nFluxVMEC /), phipf(:))
!    !! scale toroidal flux to get internal VMEC Phi
!    phipf(:nFluxVMEC) = REAL(signgs, wp) * phipf(:nFluxVMEC) / TwoPi
    !! read R_mn
    ioError = NF_INQ_VARID(ncid, "rmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), rmnc(:, 1:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading Rmnc'
    !! read Z_mn
    ioError = NF_INQ_VARID(ncid, "zmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), zmns(:, 1:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading Zmns'
    !! read lambda_mn on HALF MESH
    ioError = NF_INQ_VARID(ncid, "lmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), lmns(:, 1:))
    IF (ioError /= 0)  STOP 'VMEC READIN: problem reading lmns'
    IF(lasym)THEN
      ioError = NF_INQ_VARID(ncid, "rmns", id)
      ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
           nFluxVMEC /), rmns(:, 1:))
      IF (ioError /= 0)  STOP 'VMEC READIN: problem reading Rmns'
      !! read Z_mn
      ioError = NF_INQ_VARID(ncid, "zmnc", id)
      ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
           nFluxVMEC /), zmnc(:, 1:))
      IF (ioError /= 0)  STOP 'VMEC READIN: problem reading Zmnc'
      !! read lambda_mn on HALF MESH
      ioError = NF_INQ_VARID(ncid, "lmnc", id)
      ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
           nFluxVMEC /), lmnc(:, 1:))
      IF (ioError /= 0)  STOP 'VMEC READIN: problem reading lmnc'
    END IF
    !! read jacobian_mn on HALF MESH!!
!    ioError = NF_INQ_VARID(ncid, "gmnc", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
!         mn_mode_nyq, nFluxVMEC /), gmnc(:, 1:))

    IF (ioError /= 0) THEN
      WRITE(*,*) " Cannot read variables from ", TRIM(fileName), "!"
      WRITE(*,*) " Possible wrong data!"
      CALL EXIT(5)
    END IF

    ioError = NF_CLOSE(ncid)
    WRITE(*,'(4X,A)')'...DONE.'

  END SUBROUTINE ReadVMEC

END MODULE MOD_VMEC_Readin
