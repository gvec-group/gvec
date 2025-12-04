!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **VMEC readin**
!!
!!
!===================================================================================================================================
MODULE MODgvec_VMEC_Readin

  USE MODgvec_Globals,ONLY:wp,twoPi,UNIT_stdOut,abort
  IMPLICIT NONE

  PUBLIC

!> mu_0 (Permeability)
  REAL(wp), PARAMETER :: mu0   = 2.0_wp * twoPi * 1.0E-7_wp         ! V s / (A m)
!> version number
  REAL(wp), PARAMETER :: version = 3.75_wp

  INTEGER :: nFluxVMEC       !< number of flux surfaces in VMEC output
  INTEGER :: mn_mode         !< number of modes (theta,zeta)
  INTEGER :: mn_mode_nyq     !< number of modes (Nyquist, 2D)
  INTEGER :: nfp             !< number of field periods
  LOGICAL :: lasym=.FALSE.   !< if lasym=0, solution is (rmnc,zmns,lmns),
                             !< if lasym=1 solution is (rmnc,rmns,zmnc,zmns,lmnc,lmns)
  LOGICAL :: lrfp=.FALSE.    !< if lrfp=0, profiles use normalized toroidal flux, iota is used to compute pol. flux
                             !< if lrfp=1, profiles use normalized poloidal flux, q is used to compute tor. flux
  INTEGER :: mPol            !< poloidal mode number
  INTEGER :: nTor            !< toroidal mode number
  INTEGER :: signgs=-1       !< signum of sqrtG
!  INTEGER :: mnmax           !< maximum mode number over m,n
!  INTEGER :: mnmax_nyq       !< maximum mode number over m_nyq,n_nyq
  REAL(wp) :: b0             !< magnetic field on axis B_0
  REAL(wp) :: aMinor         !< aMinor
  REAL(wp) :: rMajor         !< major radius
  REAL(wp) :: volume         !< volume

!! vector parameters
  REAL(wp), ALLOCATABLE :: xm(:)      !< poloidal mode numbers
  REAL(wp), ALLOCATABLE :: xn(:)      !< toroidal mode numbers
  REAL(wp), ALLOCATABLE :: xm_nyq(:)  !< poloidal mode numbers (Nyquist)
  REAL(wp), ALLOCATABLE :: xn_nyq(:)  !< toroidal mode numbers (Nyquist)
  REAL(wp), ALLOCATABLE :: iotaf(:)   !< iota profile on full grid
  REAL(wp), ALLOCATABLE :: presf(:)   !< pressure on full grid
  REAL(wp), ALLOCATABLE :: phi(:)     !< toroidal flux on full mesh
  REAL(wp), ALLOCATABLE :: chi(:)     !< poloidal flux on full mesh
  REAL(wp), ALLOCATABLE :: rmnc(:, :) !< R(iMode,iFlux) (cosine components on full mesh)
  REAL(wp), ALLOCATABLE :: rmns(:, :) !< R(iMode,iFlux) (sine components on full mesh), if lasym=T
  REAL(wp), ALLOCATABLE :: zmns(:, :) !< Z(iMode,iFlux) (sine components on full mesh)
  REAL(wp), ALLOCATABLE :: zmnc(:, :) !< Z(iMode,iFlux) (cosine components on full mesh), if lasym=T
  CHARACTER(LEN=4)      :: lambda_grid!< if lambda read from "full" or "half"
  REAL(wp), ALLOCATABLE :: lmns(:, :) !< lambda(iMode,iFlux) (sine components )
  REAL(wp), ALLOCATABLE :: lmnc(:, :) !< lambda(iMode,iFlux) (cosine components), if lasym=T

  REAL(wp), ALLOCATABLE :: gmnc(:, :) !< jacobian (cosine components (read on half mesh, interpolated on full
                                      !< mesh, mnMode_nyqist ))
  REAL(wp), ALLOCATABLE :: bmnc(:, :) !< |B| (cosine components (read on half mesh, interpolated on full
                                      !< mesh, mnMode_nyqist ))
  REAL(wp), ALLOCATABLE :: gmns(:, :) !< jacobian (sine components (read on half mesh, interpolated on full
                                      !< mesh, mnMode_nyqist ))
  REAL(wp), ALLOCATABLE :: bmns(:, :) !< |B| (sine components (read on half mesh, interpolated on full
                                      !< mesh, mnMode_nyqist ))
#if NETCDF
INTERFACE ReadVMEC
  MODULE PROCEDURE ReadVMEC
END INTERFACE ReadVMEC
#endif /*NETCDF*/

INTERFACE FinalizeReadVMEC
  MODULE PROCEDURE FinalizeReadVMEC
END INTERFACE FinalizeReadVMEC

CONTAINS

!===================================================================================================================================
!> READ VMEC "wout" datafile, needs netcdf library!
!!
!===================================================================================================================================
SUBROUTINE ReadVMEC(fileName,file_Format)
  USE MODgvec_Globals,ONLY:MPIroot
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
  CHARACTER(LEN = *), INTENT(IN) :: fileName
  INTEGER           , INTENT(IN) :: file_Format
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: ok
!===================================================================================================================================
  IF(.NOT. MPIroot)RETURN
  SELECT CASE(file_Format)
  CASE(0)
#if NETCDF
    ! read VMEC 2000 output (netcdf)
    WRITE(UNIT_stdOut,'(4X,A)')'VMEC READ WOUT FILE "'//TRIM(fileName)//'" in NETCDF format ...'
    CALL ReadVmec_NETCDF(fileName)
#else
    CALL abort(__STAMP__,&
        "cannot read VMEC file, since code is compiled with BUILD_NETCDF=OFF")
#endif
  CASE(1)
    WRITE(UNIT_stdOut,'(4X,A)')'VMEC READ WOUT FILE "'//TRIM(fileName)//'" in NEMEC ASCII format ...'
    CALL ReadNEMEC(fileName,1,ok)
    IF(ok.NE.0) CALL abort(__STAMP__,&
        "Problems with VMEC readin from NEMEC format (VMECwoutfile_format=1), maybe file is binary")
  CASE(2)
    WRITE(UNIT_stdOut,'(4X,A)')'VMEC READ WOUT FILE "'//TRIM(fileName)//'" in NEMEC BINARY format ...'
    CALL ReadNEMEC(fileName,0,ok)
    IF(ok.NE.0) CALL abort(__STAMP__,&
        "Problems with VMEC readin from NEMEC format (VMECwoutfile_format=2), maybe file is ascii")
  END SELECT !file_format
END SUBROUTINE ReadVMEC

!===================================================================================================================================
!> READ VMEC "wout" datafile generated by NEMEC, routine provided by Erika Strumberger, IPP Garching
!! can be either ascii or binary
!===================================================================================================================================
SUBROUTINE ReadNEMEC(fileName,itype,ok)
  USE MODgvec_globals,ONLY:GETFREEUNIT
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
  CHARACTER(LEN = *), INTENT(IN) :: fileName
  INTEGER, INTENT(IN)           :: itype  !! =0: binary, =1: ascii
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(OUT)           :: ok
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: nrho,mpnt,nsin !nfp
  INTEGER              :: mpol1,iasym
  INTEGER              :: inUnit,iMode
  INTEGER              :: ierr,j,m,n,nmin0
  REAL(wp)             :: enfp,enrho,empol,entor,empnt,eiasym
  REAL(wp)             :: ds
  REAL(wp)             :: gam,phiedge

  REAL(wp),ALLOCATABLE :: hiota(:),hpres(:),hbuco(:),hbvco(:)
  REAL(wp),ALLOCATABLE :: hmass(:),hphip(:),fphi(:),hvp(:)
  REAL(wp),ALLOCATABLE :: hoverr(:),fjcuru(:),fjcurv(:),hspecw(:)

  REAL(wp),ALLOCATABLE :: frmnc(:,:,:),frmns(:,:,:)
  REAL(wp),ALLOCATABLE :: fzmnc(:,:,:),fzmns(:,:,:)
  REAL(wp),ALLOCATABLE :: hbsmnc_dw(:,:,:),hbsmns_dw(:,:,:)
  REAL(wp),ALLOCATABLE :: hbumnc_dw(:,:,:),hbumns_dw(:,:,:)
  REAL(wp),ALLOCATABLE :: hbvmnc_dw(:,:,:),hbvmns_dw(:,:,:)
  REAL(wp),ALLOCATABLE :: hbumnc_up(:,:,:),hbumns_up(:,:,:)
  REAL(wp),ALLOCATABLE :: hbvmnc_up(:,:,:),hbvmns_up(:,:,:)
  REAL(wp),ALLOCATABLE :: flmnc(:,:,:),flmns(:,:,:)

  REAL(wp),ALLOCATABLE :: fsve(:),hsve(:)    ! radial mesh
!===================================================================================================================================
  inUnit=GETFREEUNIT()
  ok = 0

  ! test input
  IF(itype.LT.0)THEN
    WRITE(UNIT_stdOut,*) '********** USER error **********',itype
    WRITE(UNIT_stdOut,*)'problem reading file:',TRIM(filename)
    ok = -1 ; RETURN
  END IF

  ! open equilibrium file
  IF(itype.EQ.0) THEN
    OPEN(inUnit,file=trim(filename),form='unformatted', status='old',iostat=ierr)
  ELSE
    OPEN(inUnit,file=trim(filename),form='formatted', status='old',iostat=ierr)
  END IF
  IF(ierr.NE.0) THEN
    WRITE(UNIT_stdOut,*) '********** USER error **********'
    WRITE(UNIT_stdOut,*)'ierr = ',ierr
    WRITE(UNIT_stdOut,*)'could not open file:',TRIM(filename)
    ok = -1 ; RETURN
  END IF

! --- read dimensions
  IF(itype.EQ.0) THEN
    READ(inUnit) gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge
  ELSE
    READ(inUnit,*) gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge
  END IF
  nfp    = NINT(enfp)
  nrho   = NINT(enrho)
  mpol   = NINT(empol)
  ntor   = NINT(entor)
  mpnt   = NINT(empnt)
  iasym  = NINT(eiasym)

  WRITE(UNIT_stdOut,'(6X,A)')        '-----------------------------------'
  WRITE(UNIT_stdOut,'(6X,A)')        'NEMEC file parameters:'
  WRITE(UNIT_stdOut,'(6X,A,I6)')     '  nfp     : ',nfp
  WRITE(UNIT_stdOut,'(6X,A,I6)')     '  nrho    : ',nrho
  WRITE(UNIT_stdOut,'(6X,A,I6)')     '  mpol    : ',mpol
  WRITE(UNIT_stdOut,'(6X,A,I6)')     '  ntor    : ',ntor
  WRITE(UNIT_stdOut,'(6X,A,E23.15)') '  gamma   : ',gam
  WRITE(UNIT_stdOut,'(6X,A,E23.15)') '  phiEdge : ',phiEdge
  IF(iasym.EQ.0)THEN
    WRITE(UNIT_stdOut,'(6X,A,I6)')   '  iasym=0... SYMMETRIC (half fourier modes)'
  ELSE
    WRITE(UNIT_stdOut,'(6X,A,I6)')   '  iasym=1... ASYMMETRIC'
  END IF
  WRITE(UNIT_stdOut,'(6X,A)')        '-----------------------------------'

  nsin   = nrho-1
  mpol1  = mpol-1
  ds     = 1./REAL(nsin,wp)

  ALLOCATE(frmnc(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(frmns(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(fzmnc(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(fzmns(0:mpol1,-ntor:ntor,0:nsin))

  ALLOCATE(hbumnc_up(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbumns_up(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbvmnc_up(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbvmns_up(0:mpol1,-ntor:ntor,0:nsin))

  ALLOCATE(hbsmnc_dw(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbsmns_dw(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbumnc_dw(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbumns_dw(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbvmnc_dw(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(hbvmns_dw(0:mpol1,-ntor:ntor,0:nsin))

  ALLOCATE(flmnc(0:mpol1,-ntor:ntor,0:nsin))
  ALLOCATE(flmns(0:mpol1,-ntor:ntor,0:nsin))

  ALLOCATE(hiota(1:nsin),hpres(1:nsin),hbuco(1:nsin),hbvco(1:nsin))
  ALLOCATE(hmass(1:nsin),hphip(1:nsin),fphi(1:nsin),hvp(1:nsin))
  ALLOCATE(hoverr(1:nsin),fjcuru(1:nsin),fjcurv(1:nsin))
  ALLOCATE(hspecw(1:nsin))

  ALLOCATE(fsve(0:nsin),hsve(1:nsin))

! --- read NEMEC output
  IF(itype.EQ.0) THEN

    DO j=0,nsin
      DO m=0,mpol1
        nmin0=-ntor
        IF(m.EQ.0) nmin0=0
        DO n=nmin0,ntor
! --- full mesh
          READ(inUnit) frmnc(m,n,j),fzmns(m,n,j),         &
                      frmns(m,n,j),fzmnc(m,n,j),         &
! --- half mesh
                      hbumnc_up(m,n,j),hbvmnc_up(m,n,j), &
                      hbumns_up(m,n,j),hbvmns_up(m,n,j), &
! --- full mesh
                      flmns(m,n,j),flmnc(m,n,j),         &
! --- half mesh
                      hbumnc_dw(m,n,j),hbvmnc_dw(m,n,j), &
                      hbsmns_dw(m,n,j),                  &
                      hbumns_dw(m,n,j),hbvmns_dw(m,n,j), &
                      hbsmnc_dw(m,n,j)
        END DO
      END DO
    END DO

! --- half mesh
    READ(inUnit) (hiota(j),hmass(j),hpres(j),hphip(j),hbuco(j), &
                hbvco(j),fphi(j),hvp(j),hoverr(j),fjcuru(j),   &
                fjcurv(j),hspecw(j),j=1,nsin)
  ELSE !itype = 1

    DO j=0,nsin
      DO m=0,mpol1
        nmin0=-ntor
        IF(m.EQ.0) nmin0=0
        DO n=nmin0,ntor
! --- full mesh
          READ(inUnit,*) frmnc(m,n,j),fzmns(m,n,j),       &
                      frmns(m,n,j),fzmnc(m,n,j),         &
! --- half mesh
                      hbumnc_up(m,n,j),hbvmnc_up(m,n,j), &
                      hbumns_up(m,n,j),hbvmns_up(m,n,j), &
! --- full mesh
                      flmns(m,n,j),flmnc(m,n,j),         &
! --- half mesh
                      hbumnc_dw(m,n,j),hbvmnc_dw(m,n,j), &
                      hbsmns_dw(m,n,j),                  &
                      hbumns_dw(m,n,j),hbvmns_dw(m,n,j), &
                      hbsmnc_dw(m,n,j)
        END DO
      END DO
    END DO

! --- half mesh
    READ(inUnit,*) (hiota(j),hmass(j),hpres(j),hphip(j),hbuco(j),   &
                  hbvco(j),fphi(j),hvp(j),hoverr(j),fjcuru(j),     &
                  fjcurv(j),hspecw(j),j=1,nsin)
  END IF

  CLOSE(inUnit)

! --- grid in s
  fsve(0)     = 0.
  fsve(nsin)  = 1.
  DO j=1,nsin-1
! --- half mesh
    hsve(j)  = (j-0.5)*ds
! --- full mesh
    fsve(j)   = j*ds
  END DO
  hsve(nsin) = (nsin-0.5)*ds

!!!  !test output
!!!!  WRITE(UNIT_StdOut,"(6x,'hsve',10x'hiota',9x,'hmass',9x,'hpres',9x,'hphip', &
!!!!                      9x,'hbuco',9x,'hbvco',10x,'fphi',11x,'hvp',9x,'hoverr', &
!!!!                      8x,'hspecw')")
!!!  DO j=1,nsin
!!!    WRITE(UNIT_StdOut,"(13(2x,e12.4))") hsve(j),hiota(j),hmass(j),hpres(j),hphip(j), &
!!!                                        hbuco(j),hbvco(j),fphi(j),hvp(j),hoverr(j),hspecw(j)
!!!  END DO
!!!  WRITE(UNIT_StdOut,"(6x,'fsve',10x,'fjcuru',8x,'fjcurv')")
!!!  DO j=1,nsin
!!!    WRITE(UNIT_StdOut,"(3(2x,e12.4))") fsve(j),fjcuru(j),fjcurv(j)
!!!  END DO

! --- translate to GVEC datastructure

  IF(gam.GT.1.0e-04) &
    CALL abort(__STAMP__, &
                      "readNEMEC: currently, only gamma=0 is supported")

  IF(SUM(ABS(fsve(1:nsin)*phiEdge-fphi(1:nsin)))/REAL(nsin,wp).GT.1.0e-08) &
    CALL abort(__STAMP__, &
                      "readNEMEC: toroidal flux does not seem equidistant...")

  !not needed further
  DEALLOCATE(hbumnc_up,hbumns_up,hbvmnc_up,hbvmns_up)
  DEALLOCATE(hbsmnc_dw,hbsmns_dw,hbumnc_dw,hbumns_dw,hbvmnc_dw,hbvmns_dw)
  DEALLOCATE(hbuco,hbvco,hmass,hphip,hvp,hoverr,fphi,fjcuru,fjcurv,hspecw)


  nFluxVMEC=nrho

  mn_mode = ntor+1 + (mpol-1)*(2*ntor+1)

  mn_mode_nyq=mn_mode !nyquist do not exist here

  lasym = (iasym.EQ.1)

  CALL alloc_all() !needs nFluxVMEC,mn_mode,mn_mode_nyq and lasym

  ! parameters that are not read?! set to zero here
  b0=0.
  aMinor=0.
  rMajor=0.
  volume=0.
  gmnc=0.
  bmnc=0.

  !use same loop as for read in to set modes m,n
  iMode=0
  DO m=0,mpol1
    nmin0=-ntor
    IF(m.EQ.0) nmin0=0
    DO n=nmin0,ntor
      iMode=iMode+1
      xm(iMode)=REAL(m,wp) ; xn(iMode)=REAL(n,wp)
    END DO !n
  END DO !m
  IF(iMode.NE.mn_mode) STOP 'mn_mode not correct'

  lambda_grid="full"

  DO iMode=1,mn_mode
    rmnc(iMode,1:nFluxVMEC)=frmnc(NINT(xm(iMode)),NINT(xn(iMode)),0:nsin)
    zmns(iMode,1:nFluxVMEC)=fzmns(NINT(xm(iMode)),NINT(xn(iMode)),0:nsin)
    lmns(iMode,1:nFluxVMEC)=flmns(NINT(xm(iMode)),NINT(xn(iMode)),0:nsin)
  END DO
  IF(lasym)THEN
    DO iMode=1,mn_mode
      rmns(iMode,1:nFluxVMEC)=frmns(NINT(xm(iMode)),NINT(xn(iMode)),0:nsin)
      zmnc(iMode,1:nFluxVMEC)=fzmnc(NINT(xm(iMode)),NINT(xn(iMode)),0:nsin)
      lmnc(iMode,1:nFluxVMEC)=flmnc(NINT(xm(iMode)),NINT(xn(iMode)),0:nsin)
    END DO !iMode
    gmns=0.
    bmns=0.
  END IF

  !dont forget nfp on toroidal modes
  xn=REAL(nfp)*xn

  xm_nyq=xm; xn_nyq=xn

  phi= fsve*(-phiedge) /TWOPI
  chi = 0.

  CALL halfToFull(nFluxVMEC,hiota,iotaf)

  CALL halfToFull(nFluxVMEC,hpres,presf)
  presf=presf/(2.0e-07_wp*TWOPI) !/mu_0

  DEALLOCATE(fsve,hsve,hiota,hpres)
  DEALLOCATE(frmnc,frmns,fzmnc,fzmns)
  DEALLOCATE(flmnc,flmns)

END SUBROUTINE ReadNEMEC


!===================================================================================================================================
!> Fit disrete data along flux surfaces as spline and then interpolate to full data
!!
!===================================================================================================================================
SUBROUTINE HalfToFull(nFlux,y_half,y_full)
! MODULES
  USE MODgvec_cubic_spline, ONLY:t_cubspl
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: nFlux          !! number of flux surfaces
  REAL(wp),INTENT(IN) :: y_half(2:nFlux)  !! value number 2 at axis+1/2 grid point,last at 1-1/2dx
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
! LOCAL VARIABLES
  REAL(wp),INTENT(OUT):: y_full(1:nFlux) !! values at full grid (first on axis)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: iFlux
  REAL(wp)          :: y_half_Spl(nFlux+1)  ! point values fitted on half grid
  REAL(wp)          :: rho(1:nFlux)
  REAL(wp)          :: rho_half(1:nFlux+1)
  INTEGER           :: iFlag
  CHARACTER(len=100):: message
  TYPE(t_cubspl),ALLOCATABLE    :: spl_half ! spline on half grid
!===================================================================================================================================
  rho(1)=0.
  DO iFlux=2,nFlux
    rho(iFlux)=REAL(iFlux-1,wp)/REAL(nFlux-1,wp) !equidistant grid [0,1]
    rho_half(iFlux)=SQRT(0.5_wp*(rho(iFlux)+rho(iFlux-1))) !halfgrid scaled with sqrt
  END DO
  rho=SQRT(rho) !full grid scaled with sqrt
  !add end points
  rho_half(1)=0.0_wp
  rho_half(nFlux+1)=1.0_wp

  y_half_Spl=0.
  DO iFlux=2,nFlux
      y_half_Spl(iFlux)=y_half(iFlux)
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
  y_Half_Spl(1)=(y_Half_Spl(2)*rho_half(3)**2-y_Half_Spl(3)*rho_half(2)**2) /(rho_half(3)**2-rho_half(2)**2)
  !Extrapolate to Edge
  y_Half_Spl(nFlux+1)= ( y_half_Spl(nFlux  )*(rho_half(nFlux+1)-rho_half(nFlux-1))     &
                        -y_half_Spl(nFlux-1)*(rho_half(nFlux+1)-rho_half(nFlux  )) )   &
                      /(  rho_half(nFlux)   -rho_half(nFlux-1) )
  spl_half=t_cubspl(rho_half,y_half_spl(:), BC=(/1,0/))
  iflag=0
  message=''
  y_full(:)=spl_half%eval(rho,0)
  !redo extrapolation with rho
  y_full(1)  =(y_full(2)*rho(3)**2-y_full(3)*rho(2)**2) /(rho(3)**2-rho(2)**2)

END SUBROUTINE HalfToFull

#if NETCDF
!===================================================================================================================================
!> READ VMEC "wout" datafile in netcdf format, needs netcdf library!
!! lrfp
!===================================================================================================================================
SUBROUTINE ReadVMEC_NETCDF(fileName)
  !USE netcdf
  IMPLICIT NONE
  INCLUDE "netcdf.inc"
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
  CHARACTER(LEN = *), INTENT(IN) :: fileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: ioError, ncid, id,int_logical
!===================================================================================================================================

  !! open NetCDF input file
  ioError = NF_OPEN(TRIM(fileName), NF_NOWRITE, ncid)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
        " Cannot open "//TRIM(fileName)//" in ReadVmecOutput!")

  !! get array dimensions
  !! radial dimension
!  ioError = NF_INQ_DIMID(ncid, "radius", id)
!  ioError = ioError + NF_INQ_DIMLEN(ncid, id, nFluxVMEC)
  !! number of fourier components of r, z, lambda
  ioError =  NF_INQ_DIMID(ncid, "mn_mode", id)
  ioError = ioError + NF_INQ_DIMLEN(ncid, id, mn_mode)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                          'VMEC READIN: problem reading mn_mode' )
  !! number of fourier components of b_u, b_v, b_s
  ioError = NF_INQ_DIMID(ncid, "mn_mode_nyq", id)
  ioError = ioError + NF_INQ_DIMLEN(ncid, id, mn_mode_nyq)
  IF (ioError .NE. 0) WRITE(UNIT_stdOut,*) 'INFO VMEC READIN: problem reading mn_mode_nyq'

  !! get number of field periods
  ioError = NF_INQ_VARID(ncid, "nfp", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, nfp)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                          'VMEC READIN: problem reading n-fp' )
  !! get dimension of s
  ioError = NF_INQ_VARID(ncid, "ns", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, nFluxVMEC)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                          'VMEC READIN: problem reading n-s' )
  !! get poloidal mode number
  ioError = NF_INQ_VARID(ncid, "mpol", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, mPol)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading mpol' )
  !! get toroidal mode number
  ioError = NF_INQ_VARID(ncid, "ntor", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, nTor)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading ntor' )
!  !! get mnmax
!  ioError = NF_INQ_VARID(ncid, "mnmax", id)
!  ioError = ioError + NF_GET_VAR_INT(ncid, id, mnmax)
!  !! get mnmax_nyq
!  ioError = NF_INQ_VARID(ncid, "mnmax_nyq", id)
!  ioError = ioError + NF_GET_VAR_INT(ncid, id, mnmax_nyq)

  !! get iasym
  ioError = NF_INQ_VARID(ncid, "lasym__logical__", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, int_logical) ; lasym=(int_logical.EQ.1)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading lasym' )
  !! get lrfp
  ioError = NF_INQ_VARID(ncid, "lrfp__logical__", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, int_logical); lrfp=(int_logical.EQ.1)
  IF (ioError .NE. 0) THEN
    WRITE(UNIT_stdOut,*) 'INFO VMEC READIN: problem reading lrfp'
  ELSE
    IF (lrfp) THEN
      WRITE(UNIT_stdOut,'(4X,A)') "  VMEC run with lrfp=TRUE !!!"
      STOP
    END IF
  END IF
  !! get B_0
  ioError = NF_INQ_VARID(ncid, "b0", id)
  ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, b0)
  IF (ioError .NE. 0) THEN
    WRITE(UNIT_stdOut,*) 'INFO VMEC READIN: problem reading b0'
  ELSE
    !! check the sign of b0
    IF (b0 < 0) THEN
      WRITE(UNIT_stdOut,'(4X,A)') "  VMEC run with b0 < 0 !!!"
    END IF
  END IF
  !! get signgs
  ioError = NF_INQ_VARID(ncid, "signgs", id)
  ioError = ioError + NF_GET_VAR_INT(ncid, id, signgs)
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading signgs' )
  IF (signgs < 0) THEN
    WRITE(UNIT_stdOut,'(4X,A)') "  VMEC data has sign gs < 0 !!!"
  END IF
  !! get Aminor_p
  ioError = NF_INQ_VARID(ncid, "Aminor_p", id)
  ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, aMinor)
  IF (ioError .NE. 0) WRITE(UNIT_stdOut,*) 'INFO VMEC READIN: problem reading Aminor_p'
  !! get Rmajor_p
  ioError = NF_INQ_VARID(ncid, "Rmajor_p", id)
  ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, rMajor)
  IF (ioError .NE. 0) WRITE(UNIT_stdOut,*) 'INFO VMEC READIN: problem reading Rmajor_p'
  !! get volume_p
  ioError = NF_INQ_VARID(ncid, "volume_p", id)
  ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, volume)
  IF (ioError .NE. 0) WRITE(UNIT_stdOut,*) 'INFO VMEC READIN: problem reading volume_p'

  CALL alloc_all() !needs nFluxVMEC,mn_mode,mn_mode_nyq and lasym

  !! read x_m
  ioError = NF_INQ_VARID(ncid, "xm", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /), (/ mn_mode /),&
       xm(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading xm' )
  !! read x_n
  ioError = NF_INQ_VARID(ncid, "xn", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /), (/ mn_mode /),&
       xn(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading xn' )
  !! read x_m^nyq
  ioError = NF_INQ_VARID(ncid, "xm_nyq", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
       (/ mn_mode_nyq /), xm_nyq(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading xm_nyq' )
  !! read x_n^nyq
  ioError = NF_INQ_VARID(ncid, "xn_nyq", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
       (/ mn_mode_nyq /), xn_nyq(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading xn_nyq' )
  !! read iotaf
  ioError = NF_INQ_VARID(ncid, "iotaf", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
       (/ nFluxVMEC /), iotaf(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading iotaf' )
  !! read presf
  ioError = NF_INQ_VARID(ncid, "presf", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
       (/ nFluxVMEC /), presf(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading presf' )
  !! read phi
  ioError = NF_INQ_VARID(ncid, "phi", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
       (/ nFluxVMEC /), phi(:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading phi' )

  !! scale toroidal flux to get internal VMEC Phi
  phi(:nFluxVMEC) = REAL(signgs, wp) * phi(:nFluxVMEC) / TwoPi

  !! read chi
  ioError = NF_INQ_VARID(ncid, "chi", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
       (/ nFluxVMEC /), chi(:))
  !IF (ioError .NE. 0)  STOP 'VMEC READIN: problem reading chi'
  IF (ioError .NE. 0) THEN
    WRITE(*,'(4X,A)') 'WARNING VMEC READIN: problem reading chi  (not used in GVEC up to now), is set to zero!!!'
    chi=0.
  END IF

  !! scale poloidal flux to get internal VMEC chi (so that iota=chi'/phi' >0 )
  chi(:nFluxVMEC) = chi(:nFluxVMEC) / TwoPi

  !! read phipf
!  ioError = NF_INQ_VARID(ncid, "phipf", id)
!  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!       (/ nFluxVMEC /), phipf(:))
!  !! scale toroidal flux to get internal VMEC Phi
!  phipf(:nFluxVMEC) = REAL(signgs, wp) * phipf(:nFluxVMEC) / TwoPi
  !! read R_mn
  ioError = NF_INQ_VARID(ncid, "rmnc", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
       nFluxVMEC /), rmnc(:, 1:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading Rmnc' )
  !! read Z_mn
  ioError = NF_INQ_VARID(ncid, "zmns", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
       nFluxVMEC /), zmns(:, 1:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading Zmns' )
  !! read lambda_mn on HALF MESH
  lambda_grid="half"
  ioError = NF_INQ_VARID(ncid, "lmns", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
       nFluxVMEC /), lmns(:, 1:))
  IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading lmns' )
  IF(lasym)THEN
    ioError = NF_INQ_VARID(ncid, "rmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), rmns(:, 1:))
    IF (ioError .NE. 0) CALL abort(__STAMP__,&
                              'VMEC READIN: problem reading Rmns' )
    !! read Z_mn
    ioError = NF_INQ_VARID(ncid, "zmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), zmnc(:, 1:))
    IF (ioError .NE. 0) CALL abort(__STAMP__,&
                              'VMEC READIN: problem reading Zmnc' )
    !! read lambda_mn on HALF MESH
    ioError = NF_INQ_VARID(ncid, "lmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), lmnc(:, 1:))
    IF (ioError .NE. 0) CALL abort(__STAMP__,&
                            'VMEC READIN: problem reading lmnc' )
  END IF
  !! read jacobian_mn on HALF MESH!!
  ioError = NF_INQ_VARID(ncid, "gmnc", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
       mn_mode_nyq, nFluxVMEC /), gmnc(:, 1:))
  !! read |B| on HALF MESH!!
  ioError = NF_INQ_VARID(ncid, "bmnc", id)
  ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
       mn_mode_nyq, nFluxVMEC /), bmnc(:, 1:))
  IF(lasym)THEN
    ioError = NF_INQ_VARID(ncid, "gmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), gmns(:, 1:))
    !! read |B| on HALF MESH!!
    ioError = NF_INQ_VARID(ncid, "bmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), bmns(:, 1:))
  END IF

  IF (ioError .NE. 0) CALL abort(__STAMP__,&
               " Cannot read variables from "//TRIM(fileName)//" ! " )

  ioError = NF_CLOSE(ncid)
  WRITE(*,'(4X,A)')'...DONE.'

END SUBROUTINE ReadVMEC_NETCDF
#endif /*NETCDF*/

!===================================================================================================================================
!> allocate all arrays
!!
!===================================================================================================================================
SUBROUTINE alloc_all()
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: aStat, aError
!===================================================================================================================================
  !! allocate memory for fourier arrays
  aError = 0
  ALLOCATE(xm(mn_mode), xn(mn_mode), stat = aStat)
  xm=0.0_wp ; xn =0.0_wp
  aError = aError + aStat
  ALLOCATE(xm_nyq(mn_mode_nyq), xn_nyq(mn_mode_nyq), stat = aStat)
  xm_nyq=0.0_wp ; xn_nyq =0.0_wp
  aError = aError + aStat
  ALLOCATE(iotaf(nFluxVMEC), stat = aStat)
  iotaf=0.0_wp
  aError = aError + aStat
  ALLOCATE(presf(nFluxVMEC), stat = aStat)
  presf=0.0_wp
  aError = aError + aStat
  ALLOCATE(phi(nFluxVMEC), stat = aStat)
  phi=0.0_wp
  aError = aError + aStat
  ALLOCATE(chi(nFluxVMEC), stat = aStat)
  chi=0.0_wp
  aError = aError + aStat
!  ALLOCATE(phipf(nFluxVMEC), stat = aStat)
!  aError = aError + aStat
  ALLOCATE(rmnc(mn_mode, nFluxVMEC), stat = aStat)
  rmnc=0.0_wp
  aError = aError + aStat
  ALLOCATE(zmns(mn_mode, nFluxVMEC), stat = aStat)
  zmns=0.0_wp
  aError = aError + aStat
  ALLOCATE(lmns(mn_mode, nFluxVMEC), stat = aStat)
  lmns=0.0_wp
  aError = aError + aStat
  IF(lasym)THEN
    ALLOCATE(rmns(mn_mode, nFluxVMEC), stat = aStat)
    rmns=0.0_wp
    aError = aError + aStat
    ALLOCATE(zmnc(mn_mode, nFluxVMEC), stat = aStat)
    zmnc=0.0_wp
    aError = aError + aStat
    ALLOCATE(lmnc(mn_mode, nFluxVMEC), stat = aStat)
    lmnc=0.0_wp
    aError = aError + aStat
    ALLOCATE(gmns(mn_mode_nyq, nFluxVMEC), stat = aStat)
    gmns=0.0_wp
    aError = aError + aStat
    ALLOCATE(bmns(mn_mode_nyq, nFluxVMEC), stat = aStat)
    bmns=0.0_wp
    aError = aError + aStat
  END IF
  ALLOCATE(gmnc(mn_mode_nyq, nFluxVMEC), stat = aStat)
  gmnc=0.0_wp
  aError = aError + aStat
  ALLOCATE(bmnc(mn_mode_nyq, nFluxVMEC), stat = aStat)
  aError = aError + aStat
  bmnc=0.0_wp

  IF (aError .NE. 0) THEN
    CALL abort(__STAMP__, &
              "Allocation failure in subroutine ReadVmecOutput!" )
  END IF
END SUBROUTINE alloc_all

!===================================================================================================================================
!> Finalize: Deallocate module variables
!!
!===================================================================================================================================
SUBROUTINE FinalizeReadVMEC()
  USE MODgvec_Globals,ONLY:MPIroot
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT. MPIroot)RETURN
  SDEALLOCATE( xm       )
  SDEALLOCATE( xn       )
  SDEALLOCATE( xm_nyq   )
  SDEALLOCATE( xn_nyq   )
  SDEALLOCATE( iotaf   )
  SDEALLOCATE( presf   )
  SDEALLOCATE( phi     )
  SDEALLOCATE( chi     )
  SDEALLOCATE( rmnc    )
  SDEALLOCATE( rmns    )
  SDEALLOCATE( zmnc    )
  SDEALLOCATE( zmns    )
  SDEALLOCATE( lmnc    )
  SDEALLOCATE( lmns    )

  SDEALLOCATE(gmnc)
  SDEALLOCATE(bmnc)
  SDEALLOCATE(gmns)
  SDEALLOCATE(bmns)
END SUBROUTINE FinalizeReadVMEC

END MODULE MODgvec_VMEC_Readin
