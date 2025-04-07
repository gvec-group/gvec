!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module **VMEC**
!!
!! Initializes variables to evaluate a VMEC dataset.
!! In radial direction, a cubic spline is used to interpolate the data.
!! Calls readin of VMEC "wout" datafile (netcdf format).
!!
!===================================================================================================================================
MODULE MODgvec_VMEC
! MODULES
USE MODgvec_Globals,ONLY:wp,MPIroot
USE MODgvec_cubic_spline, ONLY: t_cubspl
IMPLICIT NONE
PRIVATE

INTERFACE InitVMEC
  MODULE PROCEDURE InitVMEC
END INTERFACE
INTERFACE VMEC_EvalSplMode
  MODULE PROCEDURE VMEC_EvalSplMode
END INTERFACE
INTERFACE FinalizeVMEC
  MODULE PROCEDURE FinalizeVMEC
END INTERFACE

PUBLIC::InitVMEC
PUBLIC::VMEC_EvalSplMode
PUBLIC::FinalizeVMEC
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Initialize VMEC module
!!
!===================================================================================================================================
SUBROUTINE InitVMEC
! MODULES
USE MODgvec_Globals,ONLY:UNIT_stdOut,abort,fmt_sep,GETFREEUNIT
USE MODgvec_rProfile_bspl, ONLY: t_rProfile_bspl
USE MODgvec_cubic_spline, ONLY: interpolate_cubic_spline
USE MODgvec_ReadInTools
USE MODgvec_VMEC_Vars
USE MODgvec_VMEC_Readin
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iMode,bc_unit
LOGICAL              :: useFilter
REAL(wp),ALLOCATABLE :: c_iota(:),knots_iota(:), c_pres(:), knots_pres(:)
REAL(wp),ALLOCATABLE :: c_phi(:),knots_phi(:), c_chi(:), knots_chi(:)
!===================================================================================================================================
IF(.NOT.MPIroot) RETURN
WRITE(UNIT_stdOut,'(A)')'  INIT VMEC INPUT ...'

!VMEC "wout*.nc"  file
VMECdataFile   = GETSTR("VMECwoutfile")
VMECFile_Format= GETINT("VMECwoutfile_format",Proposal=0)

CALL ReadVmec(VMECdataFile,VMECfile_format)

switchzeta=.TRUE. !always switch zeta=-phi_vmec
switchtheta=(signgs>0) !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! transform VMEC data (R,phi=zeta,Z) to GVEC right hand side system (R,Z,phi), swap sign of zeta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(UNIT_stdOut,'(A)')'  ... switching from VMECs (R,phi,Z) to gvecs (R,Z,phi) coordinate system ...'
IF(.NOT.switchtheta)THEN
  WRITE(UNIT_stdOut,'(A)')'  ... VMEC signgs<0, so zeta direction is switched.'
  DO iMode=1,mn_mode
    IF(NINT(xm(iMode)).EQ.0)THEN
      !xn for m=0 are only positive, so we need to swap the sign of sinus coefficients
      ! -coef_{0n} sin(-n*(-zeta))= coef_{0n} sin(-n*zeta)
      !( cosine cos(-n*(-zeta))=cos(-n*zeta), symmetric in zeta for m=0)
      zmns(iMode,:)=-zmns(iMode,:)
      lmns(iMode,:)=-lmns(iMode,:)
      IF(lasym) THEN
        rmns(iMode,:)=-rmns(iMode,:)
      END IF
    ELSE
      !for m>0 , xn are always pairs of negative and positive modes,
      !so here we can simply swap the sign of the mode number
      xn(iMode)=-xn(iMode)
    END IF
  END DO !iMode=1,mn_mode
  !additional nyq data (bmnc,gmnc
  DO iMode=1,mn_mode_nyq
    ! nyq data is only cosine (bmnc,gmnc)
    IF(NINT(xm_nyq(iMode)).NE.0)THEN
      xn_nyq(iMode)=-xn_nyq(iMode)
    END If
  END DO !iMode=1,mn_mode_nyq
  !since sign of zeta has changed, swap sign of Jacobian, too.
  gmnc=-gmnc
  IF(lasym) gmns=-gmns

  ! also iota must change sign, since its sign depend on the coordinate system
  iotaf=-iotaf
ELSE
  WRITE(UNIT_stdOut,'(A)')'  ... VMEC signgs>0, so zeta and theta direction are switched.'
  !cos(-m*theta - (-n*zeta))=cos(m*theta-n*zeta)  !same sign for coefficients
  !only sin(-m*theta - (-n*zeta))=-sin(m*theta-n*zeta)  ! switch sign of coefficients
  zmns=-zmns
  lmns=-lmns
  IF(lasym) THEN
      rmns=-rmns
  END IF
END IF !not switchtheta

!WRITE VMEC_TO_GVEC BOUNDARY AND AXIS file
bc_Unit=GETFREEUNIT()
OPEN(UNIT     = bc_Unit       ,&
     FILE     = "vmec_to_gvec_boundary_and_axis.txt" ,&
     STATUS   = 'REPLACE'   ,&
     ACCESS   = 'SEQUENTIAL' )
  WRITE(bc_unit,'(A)')'! vmec boundary and axis written in the gvec coordinates (right-handed in rho,theta,zeta)'
  WRITE(bc_unit,'(A)')'!'
  WRITE(bc_unit,'(A)')'! non-zero boundary values of rmnc vmec -> X1_b_cos'
  DO iMode=1,mn_mode
    IF(ABS(rmnc(iMode,nFluxVMEC)).GT.1e-12)THEN
      WRITE(bc_unit,'("X1_b_cos(",I5,";",I5,")",X,"=",X,E22.15)')NINT(xm(iMode)),NINT(xn(iMode)/nfp),rmnc(iMode,nFluxVMEC)
    END IF
  END DO !imode
  IF(lasym)THEN
    WRITE(bc_unit,'(A)')'! non-zero boundary values of rmns vmec -> X1_b_sin'
    DO iMode=1,mn_mode
      IF(ABS(rmns(iMode,nFluxVMEC)).GT.1e-12)THEN
        WRITE(bc_unit,'("X1_b_sin(",I5,";",I5,")",X,"=",X,E22.15)')NINT(xm(iMode)),NINT(xn(iMode)/nfp),rmns(iMode,nFluxVMEC)
      END IF
    END DO !imode
  END IF
  IF(lasym)THEN
    WRITE(bc_unit,'(A)')'! non-zero boundary values of zmnc vmec -> X2_b_cos'
    DO iMode=1,mn_mode
      IF(ABS(zmnc(iMode,nFluxVMEC)).GT.1e-12)THEN
        WRITE(bc_unit,'("X2_b_cos(",I5,";",I5,")",X,"=",X,E22.15)')NINT(xm(iMode)),NINT(xn(iMode)/nfp),zmnc(iMode,nFluxVMEC)
      END IF
    END DO !imode
  END IF
  WRITE(bc_unit,'(A)')'! non-zero boundary values of zmns vmec -> X2_b_sin'
  DO iMode=1,mn_mode
    IF(ABS(zmns(iMode,nFluxVMEC)).GT.1e-12)THEN
      WRITE(bc_unit,'("X2_b_sin(",I5,";",I5,")",X,"=",X,E22.15)') NINT(xm(iMode)),NINT(xn(iMode)/nfp),zmns(iMode,nFluxVMEC)
    END IF
  END DO !imode
  !axis
  WRITE(bc_unit,'(A)')'! non-zero axis values of rmnc vmec -> X1_a_cos'
  DO iMode=1,mn_mode
    IF((NINT(xm(imode)).EQ.0).AND.(ABS(rmnc(iMode,nFluxVMEC)).GT.1e-12))THEN
      WRITE(bc_unit,'("X1_a_cos(",I1,";",I5,")",X,"=",X,E22.15)')0,NINT(xn(iMode)/nfp),rmnc(iMode,1)
    END IF
  END DO !imode
  IF(lasym)THEN
    WRITE(bc_unit,'(A)')'! non-zero axis values of rmns vmec -> X1_a_sin'
    DO iMode=1,mn_mode
      IF((NINT(xm(imode)).EQ.0).AND.(ABS(rmns(iMode,nFluxVMEC)).GT.1e-12))THEN
        WRITE(bc_unit,'("X1_a_sin(",I1,";",I5,")",X,"=",X,E22.15)')0,NINT(xn(iMode)/nfp),rmns(iMode,1)
      END IF
    END DO !imode
  END IF
  IF(lasym)THEN
    WRITE(bc_unit,'(A)')'! non-zero axis values of zmnc vmec -> X2_a_cos'
    DO iMode=1,mn_mode
      IF((NINT(xm(imode)).EQ.0).AND.(ABS(zmnc(iMode,nFluxVMEC)).GT.1e-12))THEN
        WRITE(bc_unit,'("X2_a_cos(",I1,";",I5,")",X,"=",X,E22.15)')0,NINT(xn(iMode)/nfp),zmnc(iMode,1)
      END IF
    END DO !imode
  END IF
  WRITE(bc_unit,'(A)')'! non-zero axis values of zmns vmec -> X2_a_sin'
  DO iMode=1,mn_mode
    IF((NINT(xm(imode)).EQ.0).AND.(ABS(zmns(iMode,nFluxVMEC)).GT.1e-12))THEN
      WRITE(bc_unit,'("X2_a_sin(",I1,";",I5,")",X,"=",X,E22.15)')0,NINT(xn(iMode)/nfp),zmns(iMode,1)
    END IF
  END DO !imode
CLOSE(bc_unit)


!toroidal flux from VMEC, now called PHI!
ALLOCATE(Phi_prof(nFluxVMEC))
Phi_prof = phi

!normalized toroidal flux (=flux variable s [0;1] in VMEC)
ALLOCATE(NormFlux_prof(nFluxVMEC))
NormFlux_prof=(Phi_prof-Phi_prof(1))/(Phi_prof(nFluxVMEC)-Phi_prof(1))
WRITE(UNIT_stdOut,'(4X,A,3F10.4)')'normalized toroidal flux of first three flux surfaces',NormFlux_prof(2:4)
!poloidal flux from VMEC
ALLOCATE(chi_prof(nFluxVMEC))
chi_prof=chi
WRITE(UNIT_stdOut,'(4X,A,3F10.4)')'min/max toroidal flux',MINVAL(phi)*TwoPi,MAXVAL(phi)*TwoPi
WRITE(UNIT_stdOut,'(4X,A,3F10.4)')'min/max poloidal flux',MINVAL(chi)*TwoPi,MAXVAL(chi)*TwoPi

WRITE(UNIT_stdOut,'(4X,A, I6)')'Total Number of flux surfaces: ',nFluxVMEC
WRITE(UNIT_stdOut,'(4X,A, I6)')'Total Number of mn-modes     : ',mn_mode
WRITE(UNIT_stdOut,'(4X,A,3I6)')'Max Mode m,n,nfp             : ',NINT(MAXVAL(xm)),NINT(MAXVAL(xn)),nfp
IF(lasym)THEN
  WRITE(UNIT_stdOut,'(6X,A,I6)')   '  lasym=T... ASYMMETRIC'
ELSE
  WRITE(UNIT_stdOut,'(6X,A,I6)')   '  lasym=F... SYMMETRIC (half fourier modes)'
END IF

useFilter=.TRUE. !GETLOGICAL('VMECuseFilter',Proposal=.TRUE.) !SHOULD BE ALWAYS TRUE...

ALLOCATE(xmabs(mn_mode))
DO iMode=1,mn_mode
  xmabs(iMode)=ABS(NINT(xm(iMode)))
  IF(useFilter)THEN
    IF(xmabs(iMode) > 3) THEN !Filtering for |m| > 3
      IF(MOD(xmabs(iMode),2) == 0) THEN
        xmabs(iMode)=2 !Even mode, remove rho**2
      ELSE
        xmabs(iMode)=3 !Odd mode, remove rho**3
      END IF
    END IF
  END IF !usefilter
END DO !iMode=1,mn_mode

!prepare Spline interpolation
ALLOCATE(rho(1:nFluxVMEC))
rho(:)=SQRT(NormFlux_prof(:))


ALLOCATE(Rmnc_Spl(mn_mode))
CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Rmnc,Rmnc_Spl)

ALLOCATE(Zmns_Spl(mn_mode))
CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Zmns,Zmns_Spl)

IF(lasym)THEN
  WRITE(Unit_stdOut,'(4X,A)')'LASYM=TRUE : R,Z,lambda in cos and sin!'
  ALLOCATE(Rmns_Spl(mn_mode))
  CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Rmns,Rmns_Spl)

  ALLOCATE(Zmnc_Spl(mn_mode))
  CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Zmnc,Zmnc_Spl)

END IF


ALLOCATE(lmns_Spl(mn_mode))
IF(lasym) ALLOCATE(lmnc_Spl(mn_mode))
!WRITE(*,*)'DEBUG:lambda_grid:',lambda_grid
IF(lambda_grid.EQ."half")THEN
  !lambda given on half grid
  CALL           FitSplineHalf(mn_mode,nFluxVMEC,xmAbs,lmns,lmns_Spl)
  IF(lasym) CALL FitSplineHalf(mn_mode,nFluxVMEC,xmAbs,lmnc,lmnc_Spl)
ELSEIF(lambda_grid.EQ."full")THEN
  CALL           FitSpline(mn_mode,nFluxVMEC,xmAbs,lmns,lmns_Spl)
  IF(lasym) CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,lmnc,lmnc_Spl)
ELSE
  CALL abort(__STAMP__, &
             'no lambda_grid found!!!! lambda_grid='//TRIM(lambda_grid) )
END IF


CALL interpolate_cubic_spline(rho**2,presf,c_pres,knots_pres, (/0,0/))
vmec_pres_profile = t_rProfile_bspl(knots_pres, c_pres)
SDEALLOCATE(knots_pres)
SDEALLOCATE(c_pres)

CALL interpolate_cubic_spline(rho**2,Phi_prof,c_phi,knots_phi,(/0,0/))
vmec_phi_profile = t_rProfile_bspl(knots_phi, c_phi)
SDEALLOCATE(knots_phi)
SDEALLOCATE(c_phi)

CALL interpolate_cubic_spline(rho**2,chi_prof,c_chi,knots_chi,(/0,0/))
vmec_chi_profile = t_rProfile_bspl(knots_chi, c_chi)
SDEALLOCATE(knots_chi)
SDEALLOCATE(c_chi)

CALL interpolate_cubic_spline(rho**2,iotaf,c_iota,knots_iota,(/0,0/))
vmec_iota_profile = t_rProfile_bspl(knots_iota, c_iota)
SDEALLOCATE(knots_iota)
SDEALLOCATE(c_iota)

WRITE(Unit_stdOut,'(4X,A,3F12.4)')'tor. flux Phi  axis/middle/edge',Phi(1)*TwoPi,Phi(nFluxVMEC/2)*TwoPi,Phi(nFluxVMEC)*TwoPi
WRITE(Unit_stdOut,'(4X,A,3F12.4)')'pol. flux chi  axis/middle/edge',chi(1)*TwoPi,chi(nFluxVMEC/2)*TwoPi,chi(nFluxVMEC)*TwoPi
WRITE(Unit_stdOut,'(4X,A,3F12.4)')'   iota        axis/middle/edge',vmec_iota_profile%eval_at_rho(0.0_wp),&
vmec_iota_profile%eval_at_rho(rho(nFluxVMEC/2)),vmec_iota_profile%eval_at_rho(1.0_wp)
WRITE(Unit_stdOut,'(4X,A,3F12.4)')'  pressure     axis/middle/edge',vmec_pres_profile%eval_at_rho(0.0_wp),&
vmec_pres_profile%eval_at_rho(rho(nFluxVMEC/2)),vmec_pres_profile%eval_at_rho(1.0_wp)


WRITE(UNIT_stdOut,'(A)')'  ... INIT VMEC DONE.'
WRITE(UNIT_stdOut,fmt_sep)

END SUBROUTINE InitVMEC


!===================================================================================================================================
!> Fit disrete data along flux surfaces as spline for each fourier mode
!!
!===================================================================================================================================
SUBROUTINE FitSpline(modes,nFlux,mAbs,Xmn,Xmn_Spl)
! MODULES
USE MODgvec_VMEC_Vars, ONLY: rho
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)  :: modes             !! number of modes
INTEGER, INTENT(IN)  :: nFlux             !! number of flux surfaces
INTEGER, INTENT(IN)  :: mabs(modes)       !! filtered m-mode value
REAL(wp), INTENT(IN) :: Xmn(modes,nFlux)  !! fourier coefficients at all flux surfaces
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(t_cubspl),INTENT(OUT):: Xmn_Spl(modes)  !!  spline fitted fourier coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
REAL(wp)          :: Xmn_val(nFlux)  !
!===================================================================================================================================
  DO iMode=1,modes
    !scaling with rho^|m|
    DO iFlux=2,nFlux
      IF(mabs(iMode).EQ.0)THEN
        Xmn_val(iFlux)=Xmn(iMode,iFlux)
      ELSE
        Xmn_val(iFlux)=Xmn(iMode,iFlux) /(rho(iFlux)**mabs(iMode))
      END IF
    END DO !i
    !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
    Xmn_val(1)=(Xmn_val(2)*rho(3)**2-Xmn_val(3)*rho(2)**2) /(rho(3)**2-rho(2)**2)

    Xmn_spl(iMode)=t_cubspl(rho,Xmn_val, BC=(/1,0/))
  END DO !iMode

END SUBROUTINE FitSpline


!===================================================================================================================================
!> Fit disrete data along flux surfaces as spline for each fourier mode
!! input is given on the half mesh 2:nFluxVMEC
!!
!===================================================================================================================================
SUBROUTINE FitSplineHalf(modes,nFlux,mabs,Xmn_half,Xmn_Spl)
! MODULES
USE MODgvec_VMEC_Vars, ONLY: rho,NormFlux_prof
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: modes                  !! number of modes
INTEGER, INTENT(IN) :: nFlux                  !! number of flux surfaces
INTEGER, INTENT(IN) :: mabs(modes)            !! filtered m-mode value
REAL(wp),INTENT(IN) :: Xmn_half(modes,nFlux)  !! fourier coefficients at all flux surfaces
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(t_cubspl),INTENT(OUT):: Xmn_Spl(1:modes) !!  spline fitted fourier coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
REAL(wp)          :: Xmn_val(nFlux+1)  ! spline fitted fourier coefficients
REAL(wp)          :: rho_half(nFlux+1)
TYPE(t_cubspl),ALLOCATABLE :: spl_half
!===================================================================================================================================
DO iFlux=1,nFlux-1
  rho_half(iFlux+1)=SQRT(0.5_wp*(NormFlux_prof(iFlux+1)+NormFlux_prof(iFlux))) !0.5*(rho(iFlux)+rho(iFlux+1))
END DO
!add end points
rho_half(1)=0.0_wp
rho_half(nFlux+1)=1.0_wp

 DO iMode=1,modes
   !scaling with rho^|m|
   DO iFlux=2,nFlux
     IF(mabs(iMode).EQ.0)THEN
       Xmn_val(iFlux)=Xmn_half(iMode,iFlux)
     ELSE
       Xmn_val(iFlux)=Xmn_half(iMode,iFlux) /(rho_half(iFlux)**mabs(iMode))
     END IF
   END DO !i
   !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
   Xmn_val(1)=(Xmn_val(2)*rho_half(3)**2-Xmn_val(3)*rho_half(2)**2) /(rho_half(3)**2-rho_half(2)**2)
   !Extrapolate to Edge
   Xmn_val(nFlux+1)= ( Xmn_val(nFlux  )*(rho_half(nFlux+1)-rho_half(nFlux-1))     &
                      -Xmn_val(nFlux-1)*(rho_half(nFlux+1)-rho_half(nFlux  )) )   &
                    /(  rho_half(nFlux)   -rho_half(nFlux-1) )

   spl_half=t_cubspl(rho_half,Xmn_val, BC=(/1,0/))
   Xmn_val(1:nFlux) = spl_half%eval(rho,0)
   !respline
   Xmn_val(1)  = ( Xmn_val(2)*rho(3)**2-Xmn_val(3)*rho(2)**2) /(rho(3)**2-rho(2)**2)
   Xmn_Spl(iMode)=t_cubspl(rho,Xmn_val(1:nFlux), BC=(/1,0/))
 END DO !iMode
!
END SUBROUTINE FitSplineHalf


!===================================================================================================================================
!> evaluate spline for specific mode at position s
!!
!===================================================================================================================================
FUNCTION VMEC_EvalSplMode(mn_in,rderiv,rho_in,xx_Spl)
  ! MODULES
  USE MODgvec_VMEC_Readin
  USE MODgvec_VMEC_Vars
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    INTEGER,INTENT(IN)         :: mn_in(:) !of size 1: =jmode, of size 2: find jmode to mn
    INTEGER,INTENT(IN)         :: rderiv !0: eval spl, 1: eval spl deriv
    REAL(wp),INTENT(IN)        :: rho_in(:) !! position to evaluate rho=[0,1], rho=sqrt(phi_norm)
    TYPE(t_cubspl),INTENT(IN)  :: xx_Spl(:)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)                   :: VMEC_EvalSplMode(SIZE(rho_in))
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER                    :: jMode,modefound
    REAL(wp),DIMENSION(SIZE(rho_in))  :: rhom,drhom,xx_eval !for weighted spline interpolation
  !===================================================================================================================================
    IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                          'EvalSpl called from non-MPIroot process, but VMEC data only on root!')
    IF(size(mn_in,1).EQ.2)THEN
      modefound=0
      DO jMode=1,mn_mode
        IF((NINT(xm(jMode)).EQ.mn_in(1)).AND.(NINT(xn(jMode)).EQ.mn_in(2)))THEN
          modefound=jMode
          EXIT
        END IF
      END DO
      IF(modefound.NE.0) THEN
        jMode=modefound
      ELSE
        WRITE(*,*)'Remark: mode m= ',mn_in(1),' n= ',mn_in(2),'not found in VMEC solution, setting to zero!'
        VMEC_EvalSplMode=0.0_wp
        RETURN
      END IF
    ELSEIF(size(mn_in,1).EQ.1)THEN
      jMode=mn_in(1)
    ELSE
      CALL abort(__STAMP__, &
       'mn_in should have size 1 or 2')
    END IF

    SELECT CASE(xmabs(jMode))
    CASE(0)
      rhom=1.0_wp
      drhom=0.0_wp
    CASE(1)
      rhom=rho_in
      drhom=1.0_wp
    CASE(2)
      rhom=rho_in*rho_in
      drhom=2.0_wp*rho_in
    CASE DEFAULT
      rhom=rho_in**xmabs(jMode)
      drhom=REAL(xmabs(jMode),wp)*rho_in**(xmabs(jMode)-1)
    END SELECT
    xx_eval=xx_Spl(jMode)%eval(rho_in,0)  ! includes weight 1/rhom
    IF(rderiv.EQ.0) THEN
      VMEC_EvalSplMode=rhom*xx_eval
    ELSEIF(rderiv.EQ.1) THEN
      VMEC_EvalSplMode=rhom*xx_Spl(jMode)%eval(rho_in,1) + drhom*xx_eval
    ELSE
      CALL abort(__STAMP__, &
       'rderiv should be 0 or 1')
    END IF
  END FUNCTION VMEC_EvalSplMode

!===================================================================================================================================
!> Finalize VMEC module
!!
!===================================================================================================================================
SUBROUTINE FinalizeVMEC
! MODULES
USE MODgvec_VMEC_Vars
USE MODgvec_VMEC_Readin,ONLY:FinalizeReadVMEC
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.MPIroot) RETURN

  CALL FinalizeReadVmec()


  SDEALLOCATE(Phi_prof)
  SDEALLOCATE(NormFlux_prof)
  SDEALLOCATE(chi_prof)
  SDEALLOCATE(xmabs)
  SDEALLOCATE(rho)
  SDEALLOCATE(Rmnc_Spl)
  SDEALLOCATE(Zmns_Spl)
  SDEALLOCATE(Rmns_Spl)
  SDEALLOCATE(Zmnc_Spl)
  SDEALLOCATE(lmns_Spl)
  SDEALLOCATE(lmnc_Spl)
  SDEALLOCATE(vmec_Phi_profile)
  SDEALLOCATE(vmec_chi_profile)
  SDEALLOCATE(vmec_iota_profile)
  SDEALLOCATE(vmec_pres_profile)

END SUBROUTINE FinalizeVMEC

END MODULE MODgvec_VMEC
