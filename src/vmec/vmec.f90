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
IMPLICIT NONE
PRIVATE

INTERFACE InitVMEC 
  MODULE PROCEDURE InitVMEC 
END INTERFACE

INTERFACE VMEC_EvalSpl
  MODULE PROCEDURE VMEC_EvalSpl  
END INTERFACE

INTERFACE VMEC_EvalSplMode
  MODULE PROCEDURE VMEC_EvalSplMode
END INTERFACE

INTERFACE FinalizeVMEC 
  MODULE PROCEDURE FinalizeVMEC 
END INTERFACE

PUBLIC::InitVMEC
PUBLIC::VMEC_EvalSpl
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
USE MODgvec_Globals,ONLY:UNIT_stdOut,abort,fmt_sep
USE MODgvec_ReadInTools
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
USE MODgvec_VMEC_Vars
USE MODgvec_VMEC_lambda, ONLY:recomputeLambda
USE MODgvec_VMEC_Readin
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iMode,nyq,np_m,np_n
LOGICAL              :: useFilter
!===================================================================================================================================
IF(.NOT.MPIroot) RETURN
WRITE(UNIT_stdOut,'(A)')'  INIT VMEC INPUT ...'

!VMEC "wout*.nc"  file
VMECdataFile   = GETSTR("VMECwoutfile")
VMECFile_Format= GETINT("VMECwoutfile_format",Proposal=0)
switchZeta     = GETLOGICAL("VMEC_switchZeta",Proposal=.TRUE.)
relambda       = GETLOGICAL("VMEC_relambda",Proposal=.FALSE.)
IF(relambda) nyq=GETINT("VMEC_Lam_nyq",Proposal=4)

CALL ReadVmec(VMECdataFile,VMECfile_format)

IF(switchZeta)THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! transform VMEC data (R,phi=zeta,Z) to GVEC right hand side system (R,Z,phi), swap sign of zeta  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(UNIT_stdOut,'(A)')'  ... switching from VMECs (R,phi,Z) to gvecs (R,Z,phi) coordinate system ...'
  DO iMode=1,mn_mode
    IF(xm(iMode).EQ.0)THEN
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
    IF(xm_nyq(iMode).NE.0)THEN
      xn_nyq(iMode)=-xn_nyq(iMode)
    END If  
  END DO !iMode=1,mn_mode_nyq
  !since sign of zeta has changed, swap sign of Jacobian, too.
  gmnc=-gmnc
  IF(lasym) gmns=-gmns

  ! also iota must change sign, since its sign depend on the coordinate system
  iotaf=-iotaf
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END IF !switchZeta

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


ALLOCATE(Rmnc_Spl(4,1:nFluxVMEC,mn_mode)) !first dim is for spline interpolation
CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Rmnc,Rmnc_Spl)

ALLOCATE(Zmns_Spl(4,1:nFluxVMEC,mn_mode))
CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Zmns,Zmns_Spl)

IF(lasym)THEN
  WRITE(Unit_stdOut,'(4X,A)')'LASYM=TRUE : R,Z,lambda in cos and sin!'
  ALLOCATE(Rmns_Spl(4,1:nFluxVMEC,mn_mode)) 
  CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Rmns,Rmns_Spl)
  
  ALLOCATE(Zmnc_Spl(4,1:nFluxVMEC,mn_mode))
  CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,Zmnc,Zmnc_Spl)
  
END IF


ALLOCATE(lmns_Spl(4,1:nFluxVMEC,mn_mode))
IF(lasym) ALLOCATE(lmnc_Spl(4,1:nFluxVMEC,mn_mode))
IF(reLambda)THEN
  np_m=1+nyq*2*(NINT(MAXVAL(ABS(xm)))/2) !m_points [0,2pi]
  np_n=1+nyq*2*(NINT(MAXVAL(ABS(xn)))/(2*nfp)) !n_points [0,2pi/nfs] ->  
  !recompute lambda on FULL GRID
  CALL RecomputeLambda(np_m,np_n) 
  lambda_grid="full"
  CALL           FitSpline(mn_mode,nFluxVMEC,xmAbs,lmns,lmns_Spl)
  IF(lasym) CALL FitSpline(mn_mode,nFluxVMEC,xmAbs,lmnc,lmnc_Spl)
ELSE
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
END IF


ALLOCATE(pres_spl(4,1:nFluxVMEC))
pres_spl(1,:)=presf(:)
CALL SPLINE1_FIT(nFluxVMEC,rho,pres_Spl(:,:), K_BC1=3, K_BCN=0)

ALLOCATE(Phi_spl(4,1:nFluxVMEC))
Phi_spl(1,:)=Phi_Prof(:)
CALL SPLINE1_FIT(nFluxVMEC,rho,Phi_Spl(:,:), K_BC1=0, K_BCN=0)

ALLOCATE(chi_spl(4,1:nFluxVMEC))
chi_spl(1,:)=chi_Prof(:)
CALL SPLINE1_FIT(nFluxVMEC,rho,chi_Spl(:,:), K_BC1=0, K_BCN=0)

ALLOCATE(iota_spl(4,1:nFluxVMEC))
iota_spl(1,:)=iotaf(:)
CALL SPLINE1_FIT(nFluxVMEC,rho,iota_Spl(:,:), K_BC1=3, K_BCN=0)

WRITE(Unit_stdOut,'(4X,A,3F12.4)')'tor. flux Phi  axis/middle/edge',Phi(1)*TwoPi,Phi(nFluxVMEC/2)*TwoPi,Phi(nFluxVMEC)*TwoPi
WRITE(Unit_stdOut,'(4X,A,3F12.4)')'pol. flux chi  axis/middle/edge',chi(1)*TwoPi,chi(nFluxVMEC/2)*TwoPi,chi(nFluxVMEC)*TwoPi
WRITE(Unit_stdOut,'(4X,A,3F12.4)')'   iota        axis/middle/edge',iota_Spl(1,1),iota_Spl(1,nFluxVMEC/2),iota_Spl(1,nFluxVMEC)
WRITE(Unit_stdOut,'(4X,A,3F12.4)')'  pressure     axis/middle/edge',pres_Spl(1,1),pres_Spl(1,nFluxVMEC/2),pres_Spl(1,nFluxVMEC)


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
USE SPLINE1_MOD,       ONLY: SPLINE1_FIT 
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
REAL(wp), INTENT(OUT):: Xmn_Spl(4,nFlux,modes)  !!  spline fitted fourier coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
!===================================================================================================================================
Xmn_Spl=0.0_wp
DO iMode=1,modes
  !scaling with rho^|m|
  DO iFlux=2,nFlux
    IF(mabs(iMode).EQ.0)THEN
      Xmn_Spl(1,iFlux,iMode)=Xmn(iMode,iFlux)
    ELSE
      Xmn_Spl(1,iFlux,iMode)=Xmn(iMode,iFlux) /(rho(iFlux)**mabs(iMode))
    END IF
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
  Xmn_Spl(1,1,iMode)=(Xmn_Spl(1,2,iMode)*rho(3)**2-Xmn_Spl(1,3,iMode)*rho(2)**2) /(rho(3)**2-rho(2)**2)
!  !Quadratic extrapolation to axis with dx'(rho=0)=0.0_wp
!  r1=rho(2)**2*rho(3)**4-rho(2)**4*rho(3)**2
!  r2=rho(2)**2*rho(4)**4-rho(2)**4*rho(4)**2
!  Xmn_Spl(1,1,iMode)= ( r1*(Xmn_Spl(1,2,iMode)*rho(4)**4-Xmn_Spl(1,4,iMode)*rho(2)**4) &
!                       -r2*(Xmn_Spl(1,2,iMode)*rho(3)**4-Xmn_Spl(1,3,iMode)*rho(2)**4)) &
!                     /( r1*(rho(4)**4-rho(2)**4)-r2*(rho(3)**4-rho(2)**4))
  CALL SPLINE1_FIT(nFlux,rho,Xmn_Spl(:,:,iMode), K_BC1=3, K_BCN=0)
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
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
USE SPLINE1_MOD,       ONLY:SPLINE1_INTERP 
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
REAL(wp),INTENT(OUT):: Xmn_Spl(4,nFlux,modes) !!  spline fitted fourier coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
REAL(wp)          :: Xmn_half_Spl(4,nFlux+1)  ! spline fitted fourier coefficients 
REAL(wp)          :: rho_half(1:nFlux+1)
INTEGER           :: iFlag
CHARACTER(len=100):: message
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
      Xmn_half_Spl(1,iFlux)=Xmn_half(iMode,iFlux)
    ELSE
      Xmn_half_Spl(1,iFlux)=Xmn_half(iMode,iFlux) /(rho_half(iFlux)**mabs(iMode))
    END IF
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.0_wp
  Xmn_Half_Spl(1,1)=(Xmn_Half_Spl(1,2)*rho_half(3)**2-Xmn_Half_Spl(1,3)*rho_half(2)**2) /(rho_half(3)**2-rho_half(2)**2)
  !Extrapolate to Edge 
  Xmn_Half_Spl(1,nFlux+1)= ( Xmn_half_Spl(1,nFlux  )*(rho_half(nFlux+1)-rho_half(nFlux-1))     &
                            -Xmn_half_Spl(1,nFlux-1)*(rho_half(nFlux+1)-rho_half(nFlux  )) )   &
                               /(  rho_half(nFlux)   -rho_half(nFlux-1) )
  CALL SPLINE1_FIT(nFlux+1,rho_half,Xmn_half_Spl(:,:), K_BC1=3, K_BCN=0)
  iflag=0
  message=''
  CALL SPLINE1_INTERP((/1,0,0/),nFlux+1,rho_half,Xmn_half_Spl, &
                                nFlux  ,rho     ,Xmn_Spl(:,:,iMode),       &
                          iflag,message, K_BC1=3,K_BCN=0)
  !respline
  Xmn_Spl(2:4,:,iMode)=0.0_wp
  Xmn_Spl(1,1,iMode)  =(Xmn_Spl(1,2,iMode)*rho(3)**2-Xmn_Spl(1,3,iMode)*rho(2)**2) /(rho(3)**2-rho(2)**2)
  CALL SPLINE1_FIT(nFlux,rho,Xmn_Spl(:,:,iMode), K_BC1=3, K_BCN=0)
END DO !iMode 

END SUBROUTINE FitSplineHalf

!===================================================================================================================================
!> evaluate 1d spline at position s
!!
!===================================================================================================================================
FUNCTION VMEC_EvalSpl(rderiv,rho_in,xx_spl)
! MODULES
USE MODgvec_VMEC_Readin
USE MODgvec_VMEC_Vars
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER      , INTENT(IN   ) :: rderiv 
  REAL(wp)     , INTENT(IN   ) :: rho_in !! position to evaluate rho=[0,1], rho=sqrt(phi_norm)
  REAL(wp)     , INTENT(IN   ) :: xx_Spl(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)           :: VMEC_EvalSpl
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iGuess
  REAL(wp)           :: splout(3)
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        'EvalSpl called from non-MPIroot process, but VMEC data only on root!')
  CALL SPLINE1_EVAL((/1,rderiv,0/), nFluxVMEC,rho_in,rho,xx_Spl(:,:),iGuess,splout) 
  VMEC_EvalSpl=splout(1+rderiv)
END FUNCTION VMEC_EvalSpl

!===================================================================================================================================
!> evaluate spline for specific mode at position s
!!
!===================================================================================================================================
FUNCTION VMEC_EvalSplMode(mn_in,rderiv,rho_in,xx_Spl)
! MODULES
USE MODgvec_VMEC_Readin
USE MODgvec_VMEC_Vars
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,INTENT(IN)         :: mn_in(:) !of size 1: =jmode, of size 2: find jmode to mn
  INTEGER,INTENT(IN)         :: rderiv !0: eval spl, 1: eval spl deriv
  REAL(wp),INTENT(IN)        :: rho_in !! position to evaluate rho=[0,1], rho=sqrt(phi_norm)
  REAL(wp),INTENT(IN)        :: xx_Spl(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                   :: VMEC_EvalSplMode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                    :: iGuess,jMode,modefound
  REAL(wp)                   :: rhom,drhom,splOut(3) !for weighted spline interpolation
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
    STOP 'mn_in should have size 1 or 2'
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
  CALL SPLINE1_EVAL((/1,rderiv,0/), nFluxVMEC,rho_in,rho,xx_Spl(:,:,jMode),iGuess,splout) 
  VMEC_EvalSplMode=rhom*splout(1+rderiv)+REAL(rderiv,wp)*(drhom*splout(1))
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
  SDEALLOCATE(pres_spl)
  SDEALLOCATE(Phi_spl)
  SDEALLOCATE(chi_spl)
  SDEALLOCATE(iota_spl)

END SUBROUTINE FinalizeVMEC

END MODULE MODgvec_VMEC
