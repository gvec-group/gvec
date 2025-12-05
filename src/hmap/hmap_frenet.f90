!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** hmap_frenet **
!!
!! This map uses the Frenet frame of a given periodic input curve X0(zeta) along the curve parameter zeta in [0,2pi].
!! It uses the signed orthonormal Frenet-Serret frame (TNB frame) that can be computed from derivatives of X0  in zeta.
!! h:  X_0(\zeta) + q_1 \sigma N(\zeta) + q_2 \sigma B(\zeta)
!! with a sign switching function \sigma(\zeta), that accounts for points of zero curvature.
!! the tangent is T=X_0' / |X_0'|, the bi-normal is B= (X_0' x X_0'') / |X_0' x X_0''|, and the normal N= B X T
!! Derivatives use the Frenet-Serret formulas:
!!
!! dT/dl = k N
!! dN/dl = -kappa T + tau B
!! dB/dl = -tau N
!!
!! With  l(\zeta) being the arc-length, and l' = |X_0'|.
!! the curvature is kappa=  |X_0' x  X_0''| / (l')^3,
!! and the torsion tau= (X_0' x X_0'').X_0''' /  |X_0' x X_0''|^2
!!
!! As a first representation of the curve X0(\zeta), we choose zeta to be the geometric toroidal angle zeta=phi, such that
!!             R0(zeta)*cos(zeta)
!!  X0(zeta)=( R0(zeta)*sin(zeta)  )
!!             Z0(zeta)
!! and both R0,Z0 are represented as a real Fourier series with modes 0... n_max and number of Field periods Nfp
!! R0(zeta) = sum_{n=0}^{n_{max}} rc(n)*cos(n*Nfp*zeta) + rs(n)*sin(n*Nfp*zeta)
!===================================================================================================================================
MODULE MODgvec_hmap_frenet
! MODULES
USE MODgvec_Globals, ONLY:PI,TWOPI,CROSS,wp,Unit_stdOut,abort,MPIroot
USE MODgvec_c_hmap,    ONLY:c_hmap, c_hmap_auxvar
IMPLICIT NONE

PUBLIC

!> Store data that can be precomputed on a set ot zeta points
!> depends on hmap_frenet, but could be used for different point sets in zeta
!
TYPE,EXTENDS(c_hmap_auxvar) :: t_hmap_frenet_auxvar
  REAL(wp)  :: lp,kappa,tau,sigma,lp_p,kappa_p,tau_p
  REAL(wp),DIMENSION(3)::X0,T,N,B
END TYPE t_hmap_frenet_auxvar

TYPE,EXTENDS(c_hmap) :: t_hmap_frenet
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL  :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  ! parameters for hmap_frenet:
  !curve description
  !INTEGER             :: nfp  !already part of c_hmap. Is overwritten in init
  REAL(wp),ALLOCATABLE :: rc(:)  !! input cosine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp
  REAL(wp),ALLOCATABLE :: rs(:)  !! input   sine coefficients of R0 as array (0:n_max) of modes (0,1,...,n_max)*nfp
  REAL(wp),ALLOCATABLE :: zc(:)  !! input cosine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp
  REAL(wp),ALLOCATABLE :: zs(:)  !! input   sine coefficients of Z0 as array (0:n_max) of modes (0,1,...,n_max)*nfp
  INTEGER,ALLOCATABLE  :: Xn(:)   !! array of mode numbers,  local variable =(0,1,...,n_max)*nfp
  LOGICAL              :: omnig=.FALSE.   !! omnigenity. True: sign change of frame at pi/nfp , False: no sign change
  !---------------------------------------------------------------------------------------------------------------------------------

  CONTAINS

  FINAL :: hmap_frenet_free
  PROCEDURE :: eval_all        => hmap_frenet_eval_all
  PROCEDURE :: eval            => hmap_frenet_eval
  PROCEDURE :: eval_aux        => hmap_frenet_eval_aux
  PROCEDURE :: eval_dxdq       => hmap_frenet_eval_dxdq
  PROCEDURE :: eval_dxdq_aux   => hmap_frenet_eval_dxdq_aux
  PROCEDURE :: eval_Jh         => hmap_frenet_eval_Jh
  PROCEDURE :: eval_Jh_aux     => hmap_frenet_eval_Jh_aux
  PROCEDURE :: eval_Jh_dq      => hmap_frenet_eval_Jh_dq
  PROCEDURE :: eval_Jh_dq_aux  => hmap_frenet_eval_Jh_dq_aux
  PROCEDURE :: eval_gij        => hmap_frenet_eval_gij
  PROCEDURE :: eval_gij_aux    => hmap_frenet_eval_gij_aux
  PROCEDURE :: eval_gij_dq     => hmap_frenet_eval_gij_dq
  PROCEDURE :: eval_gij_dq_aux => hmap_frenet_eval_gij_dq_aux
  PROCEDURE :: get_dx_dqi       => hmap_frenet_get_dx_dqi
  PROCEDURE :: get_dx_dqi_aux   => hmap_frenet_get_dx_dqi_aux
  PROCEDURE :: get_ddx_dqij     => hmap_frenet_get_ddx_dqij
  PROCEDURE :: get_ddx_dqij_aux => hmap_frenet_get_ddx_dqij_aux

  !---------------------------------------------------------------------------------------------------------------------------------
  ! procedures for hmap_frenet:
  PROCEDURE :: eval_X0       => hmap_frenet_eval_X0_fromRZ
  PROCEDURE :: sigma         => hmap_frenet_sigma
END TYPE t_hmap_frenet

!INITIALIZATION FUNCTION:
INTERFACE t_hmap_frenet
  MODULE PROCEDURE hmap_frenet_init,hmap_frenet_init_params
END INTERFACE t_hmap_frenet

INTERFACE t_hmap_frenet_auxvar
  MODULE PROCEDURE hmap_frenet_init_aux
END INTERFACE t_hmap_frenet_auxvar

LOGICAL :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize the type hmap_frenet, reading from parameterfile and then call init_params
!===================================================================================================================================
FUNCTION hmap_frenet_init() RESULT(sf)
  ! MODULES
  USE MODgvec_ReadInTools, ONLY: GETLOGICAL,GETINT, GETREALARRAY
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    TYPE(t_hmap_frenet) :: sf !! self
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: nfp,n_max,nvisu
    REAL(wp),ALLOCATABLE :: Rc(:),Rs(:),Zc(:),Zs(:)
    LOGICAL :: omnig
  !===================================================================================================================================
    SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: FRENET FRAME OF A CLOSED CURVE . GET PARAMETERS:'

    nfp   = GETINT("hmap_nfp")
    n_max = GETINT("hmap_n_max")
    nvisu = GETINT("hmap_nvisu",-1)

    ALLOCATE(Rc(0:n_max)) ; Rc=0.0_wp ; Rc=GETREALARRAY("hmap_rc",n_max+1,Rc)
    ALLOCATE(Rs(0:n_max)) ; Rs=0.0_wp ; Rs=GETREALARRAY("hmap_rs",n_max+1,Rs)
    ALLOCATE(Zc(0:n_max)) ; Zc=0.0_wp ; Zc=GETREALARRAY("hmap_zc",n_max+1,Zc)
    ALLOCATE(Zs(0:n_max)) ; Zs=0.0_wp ; Zs=GETREALARRAY("hmap_zs",n_max+1,Zs)

    omnig=GETLOGICAL("hmap_omnig",.FALSE.) !omnigenity


    sf=hmap_frenet_init_params(nfp,n_max,nvisu,Rc,Rs,Zc,Zs,omnig)
    DEALLOCATE(rc,rs,zc,zs)
  END FUNCTION hmap_frenet_init

!===================================================================================================================================
!> initialize the type hmap_frenet with number of elements
!!
!===================================================================================================================================
FUNCTION hmap_frenet_init_params(nfp,n_max,nvisu,Rc,Rs,Zc,Zs,omnig) RESULT(sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: nfp         !! number of field periods
  INTEGER, INTENT(IN) :: n_max       !! maximum mode number of guiding curve
  INTEGER, INTENT(IN) :: nvisu       !! number of visualization points per field period (-1: no visualization)
  REAL(wp),INTENT(IN) :: Rc(0:n_max) !! R cos(-n*zeta) modes of guiding curve, 0..n_max
  REAL(wp),INTENT(IN) :: Rs(0:n_max) !! R sin(-n*zeta) modes of guiding curve, 0..n_max
  REAL(wp),INTENT(IN) :: Zc(0:n_max) !! Z cos(-n*zeta) modes of guiding curve, 0..n_max
  REAL(wp),INTENT(IN) :: Zs(0:n_max) !! Z sin(-n*zeta) modes of guiding curve, 0..n_max
  LOGICAL ,INTENT(IN) :: omnig       !! omnigeneity, gives sign function of Frenet frame. False: sigma=1,
                                     !! True: sigma=+1 for 0<=zeta<=pi/nfp, and -1 for pi/nfp<zeta<2pi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_frenet) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: n
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A)')'INIT HMAP :: FRENET FRAME OF A CLOSED CURVE ...'

  sf%nfp=nfp
  IF(sf%nfp.LE.0) &
     CALL abort(__STAMP__, &
          "hmap_frenet init: nfp > 0 not fulfilled!",TypeInfo="MissingParameterError")

  sf%n_max=n_max
  ALLOCATE(sf%Xn(0:sf%n_max))
  DO n=0,sf%n_max
    sf%Xn(n)=n*sf%nfp
  END DO
  ALLOCATE(sf%rc(0:sf%n_max));sf%rc=Rc
  ALLOCATE(sf%rs(0:sf%n_max));sf%rs=Rs
  ALLOCATE(sf%zc(0:sf%n_max));sf%zc=Zc
  ALLOCATE(sf%zs(0:sf%n_max));sf%zs=Zs
  sf%omnig=omnig

  IF (.NOT.(sf%rc(0) > 0.0_wp)) THEN
     CALL abort(__STAMP__, &
          "hmap_frenet init: condition rc(n=0) > 0 not fulfilled!",TypeInfo="InitializationError")
  END IF

  IF(MPIroot)THEN
    IF(nvisu.GT.0) CALL VisuFrenet(sf,nvisu*sf%nfp)

    CALL CheckZeroCurvature(sf)
  END IF

  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'...DONE.'
  IF(.NOT.test_called) CALL hmap_frenet_test(sf)

END FUNCTION hmap_frenet_init_params

!===================================================================================================================================
!> finalize the type hmap_frenet
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_frenet), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN
  DEALLOCATE(sf%rc)
  DEALLOCATE(sf%rs)
  DEALLOCATE(sf%zc)
  DEALLOCATE(sf%zs)

  sf%initialized=.FALSE.

END SUBROUTINE hmap_frenet_free


!===================================================================================================================================
!> Sample axis and check for zero (<1.e-12) curvature
!!
!===================================================================================================================================
SUBROUTINE checkZeroCurvature( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iz,nz
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)              :: lp,absB
  REAL(wp),DIMENSION((sf%n_max+1)*8) :: zeta,kappa
  LOGICAL ,DIMENSION((sf%n_max+1)*8) :: checkzero
!===================================================================================================================================
  nz=(sf%n_max+1)*8
  DO iz=1,nz
    zeta(iz)=REAL(iz-1,wp)/REAL(nz,wp)*TWOPI/sf%nfp  !0...2pi/nfp without endpoint
    CALL sf%eval_X0(zeta(iz),X0,X0p,X0pp,X0ppp)
    lp=SQRT(SUM(X0p*X0p))
    B=CROSS(X0p,X0pp)
    absB=SQRT(SUM(B*B))
    kappa(iz)=absB/(lp**3)
  END DO !iz
  checkzero=(kappa.LT.1.0e-8)
  IF(ANY(checkzero))THEN
    IF(sf%omnig)THEN
      !omnig=True: kappa can only be zero once, at 0,pi/nfp,[2pi/nfp...]
      IF(.NOT.(checkzero(1).AND.checkzero(nz/2+1).AND.(COUNT(checkzero).EQ.2)))THEN
        DO iz=1,nz
          IF(checkzero(iz)) WRITE(UNIT_StdOut,'(A,E15.5)')'         ...curvature <1e-8 at zeta/(2pi/nfp)=',zeta(iz)*sf%nfp/TWOPI
        END DO
        CALL abort(__STAMP__, &
             "hmap_frenet checkZeroCurvature with omnig=True: found additional points with zero curvature",&
             TypeInfo="InitializationError")
      END IF
    ELSE
      DO iz=1,nz
        IF(checkzero(iz)) WRITE(UNIT_StdOut,'(A,E15.5)')'         ...curvature <1e-8 at zeta/(2pi/nfp)=',zeta(iz)*sf%nfp/TWOPI
      END DO
      CALL abort(__STAMP__, &
           "hmap_frenet checkZeroCurvature with omnig=False: found points with zero curvature",&
           TypeInfo="InitializationError")
    END IF
  END IF
END SUBROUTINE CheckZeroCurvature

!===================================================================================================================================
!> Write evaluation of the axis and signed frenet frame to file
!!
!===================================================================================================================================
SUBROUTINE VisuFrenet( sf ,nvisu)
! MODULES
USE MODgvec_Output_CSV,     ONLY: WriteDataToCSV
USE MODgvec_Output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_Output_netcdf,     ONLY: WriteDataToNETCDF
USE MODgvec_Analyze_vars,     ONLY: outfileType
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  INTEGER             , INTENT(IN) :: nvisu     !!
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)              :: zeta,sigma,lp,absB,kappa,tau,eps
  INTEGER               :: iVar,ivisu,itest
  INTEGER,PARAMETER     :: nVars=26
  CHARACTER(LEN=20)     :: VarNames(1:nVars)
  REAL(wp)              :: values(1:nVars,1:nvisu*sf%nfp+1)
!===================================================================================================================================
  IF(nvisu.LE.0) RETURN
  iVar=0
  VarNames(ivar+1:iVar+3)=(/ "x", "y", "z"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"TX","TY","TZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"NX","NY","NZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/"BX","BY","BZ"/);iVar=iVar+3
  VarNames(iVar+1       )="zeta_norm"       ;iVar=iVar+1
  VarNames(iVar+1       )="sigma_sign"      ;iVar=iVar+1
  VarNames(iVar+1       )="lprime"          ;iVar=iVar+1
  VarNames(iVar+1       )="kappa"           ;iVar=iVar+1
  VarNames(iVar+1       )="tau"             ;iVar=iVar+1
  VarNames(ivar+1:iVar+3)=(/ "X0pX", "X0pY", "X0pZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/ "X0ppX", "X0ppY", "X0ppZ"/);iVar=iVar+3
  VarNames(ivar+1:iVar+3)=(/ "X0pppX", "X0pppY", "X0pppZ"/);iVar=iVar+3

!  values=0.
  DO ivisu=1,nvisu*sf%nfp+1
    eps=0.
    kappa=0.
    itest=0
    DO WHILE(kappa.LT.1.0e-6)! for tau being meaningful
      zeta=(REAL(ivisu-1,wp)+eps)/REAL(nvisu*sf%nfp,wp)*TWOPI
      CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
      lp=SQRT(SUM(X0p*X0p))
      T=X0p/lp
      B=CROSS(X0p,X0pp)
      absB=SQRT(SUM(B*B))
      kappa=absB/(lp**3)
      itest=itest+1
      eps=10**REAL(-16+itest,wp)
      IF(itest.EQ.15)THEN !meaningful kappa not found
        B=0.
        kappa=1.0e-6 !-12
        absB=1.
      END IF
    END DO
    tau=SUM(X0ppp*B)/(absB**2)
    B=B/absB
    N=CROSS(B,T)
    sigma=sf%sigma(zeta)
    iVar=0
    values(ivar+1:iVar+3,ivisu)=X0                ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=T                 ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=N*sigma           ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=B*sigma           ;iVar=iVar+3
    values(iVar+1       ,ivisu)=zeta*sf%nfp/TWOPI ;iVar=iVar+1
    values(iVar+1       ,ivisu)=sigma             ;iVar=iVar+1
    values(iVar+1       ,ivisu)=lp                ;iVar=iVar+1
    values(iVar+1       ,ivisu)=kappa             ;iVar=iVar+1
    values(iVar+1       ,ivisu)=tau               ;iVar=iVar+1
    values(ivar+1:iVar+3,ivisu)=X0p               ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=X0pp              ;iVar=iVar+3
    values(ivar+1:iVar+3,ivisu)=X0ppp             ;iVar=iVar+3
  END DO !ivisu
  IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
    CALL WriteDataToVTK(1,3,nVars-3,(/nvisu*sf%nfp/),1,VarNames(4:nVars),values(1:3,:),values(4:nVars,:),"visu_hmap_frenet.vtu")
  END IF
  IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
#if NETCDF
    CALL WriteDataToNETCDF(1,3,nVars-3,(/nvisu*sf%nfp/),(/"dim_zeta"/),VarNames(4:nVars),values(1:3,:),values(4:nVars,:), &
         "visu_hmap_frenet")
#else
    CALL WriteDataToCSV(VarNames(:) ,values, TRIM("out_visu_hmap_frenet.csv") ,append_in=.FALSE.)
#endif
  END IF
END SUBROUTINE VisuFrenet


!===================================================================================================================================
!> initialize the aux variable
!!
!===================================================================================================================================
FUNCTION hmap_frenet_init_aux( sf ,zeta,do_2nd_der) RESULT(xv)
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet),INTENT(IN) :: sf !! self
  REAL(wp)            ,INTENT(IN) :: zeta
  LOGICAL             ,INTENT(IN) :: do_2nd_der !! compute second derivative and store second derivative terms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  TYPE(t_hmap_frenet_auxvar)      :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0p,X0pp,X0ppp,X0p4,Bloc
  REAL(wp)  :: absB,absB_p
!===================================================================================================================================
  xv%do_2nd_der=do_2nd_der
  xv%zeta = zeta
  IF(xv%do_2nd_der)THEN
    CALL sf%eval_X0(zeta, xv%X0, X0p, X0pp, X0ppp,X0p4=X0p4)
  ELSE
    CALL sf%eval_X0(zeta, xv%X0, X0p, X0pp, X0ppp)
  END IF

  xv%lp = SQRT(SUM(X0p*X0p))
  xv%T = X0p / xv%lp
  Bloc = CROSS(X0p, X0pp)
  absB = SQRT(SUM(Bloc*Bloc))
  xv%kappa = absB / (xv%lp**3)
  IF(xv%kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
          "hmap_frenet cannot evaluate frame at curvature < 1e-8 !", RealInfo=zeta)
  xv%sigma = sf%sigma(zeta)
  xv%tau = SUM(X0ppp*Bloc) / (absB**2)
  xv%B = Bloc / absB
  xv%N = CROSS(xv%B, xv%T)

  IF(xv%do_2nd_der)THEN
    xv%lp_p = SUM(X0pp*X0p) / xv%lp
    absB_p = SUM(Bloc*CROSS(X0p, X0ppp)) / absB
    xv%kappa_p = (absB_p*xv%lp -3*absB * xv%lp_p) / (xv%lp**4)
    xv%tau_p   = (SUM(X0p4*Bloc)*absB -2*SUM(X0ppp*Bloc)*absB_p) / (absB**3)
  END IF

  END FUNCTION hmap_frenet_init_aux


!===================================================================================================================================
!> evaluate all metrics necessary for optimizer
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_eval_all(sf,ndims,dim_zeta,xv,&
                                q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                Jh,    g_tt,    g_tz,    g_zz,&
                                Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN)   :: sf
  INTEGER             , INTENT(IN)   :: ndims(3)    !! 3D dimensions of input arrays
  INTEGER             , INTENT(IN)   :: dim_zeta    !! which dimension is zeta dependent
  CLASS(c_hmap_auxvar), INTENT(IN)   :: xv(ndims(dim_zeta))  !! zeta point positions
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(IN) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),DIMENSION(ndims(1),ndims(2),ndims(3)),INTENT(OUT):: Jh,g_tt    ,g_tz    ,g_zz    , &
                                                               Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                                               Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                                               g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,j,k
  !===================================================================================================================================
  SELECT TYPE(xv)
  TYPE IS(t_hmap_frenet_auxvar)
    SELECT CASE(dim_zeta)
    CASE(1)
      !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        CALL hmap_frenet_eval_all_e(xv(i), &
                 q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
                 Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
                 Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
                 Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
                 g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
      END DO; END DO; END DO
      !$OMP END PARALLEL DO
    CASE(2)
       !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        CALL hmap_frenet_eval_all_e(xv(j), &
                 q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
                 Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
                 Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
                 Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
                 g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
      END DO; END DO; END DO
       !$OMP END PARALLEL DO
    CASE(3)
       !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k)
      DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
        CALL hmap_frenet_eval_all_e(xv(k), &
                 q1(i,j,k),q2(i,j,k),dX1_dt(i,j,k),dX2_dt(i,j,k),dX1_dz(i,j,k),dX2_dz(i,j,k), &
                 Jh(i,j,k)    ,g_tt(i,j,k)    ,g_tz(i,j,k)    ,g_zz(i,j,k), &
                 Jh_dq1(i,j,k),g_tt_dq1(i,j,k),g_tz_dq1(i,j,k),g_zz_dq1(i,j,k), &
                 Jh_dq2(i,j,k),g_tt_dq2(i,j,k),g_tz_dq2(i,j,k),g_zz_dq2(i,j,k), &
                 g_t1(i,j,k),g_t2(i,j,k),g_z1(i,j,k),g_z2(i,j,k),Gh11(i,j,k),Gh22(i,j,k) )
      END DO; END DO; END DO
       !$OMP END PARALLEL DO
    END SELECT
  END SELECT !type(xv)

END SUBROUTINE hmap_frenet_eval_all

!===================================================================================================================================
!> evaluate all quantities at one given point (elemental)
!!
!===================================================================================================================================
PURE SUBROUTINE hmap_frenet_eval_all_e(xv,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                       Jh,    g_tt,    g_tz,    g_zz,     &
                                       Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                       Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                       g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  )
! MODULES
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  TYPE(t_hmap_frenet_auxvar),INTENT(IN) :: xv    !! precomputed auxiliary variables
  REAL(wp),INTENT(IN)  :: q1,q2       !! solution variables q1,q2
  REAL(wp),INTENT(IN)  :: dX1_dt,dX2_dt  !! theta derivative of solution variables q1,q2
  REAL(wp),INTENT(IN)  :: dX1_dz,dX2_dz  !!  zeta derivative of solution variables q1,q2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),INTENT(OUT) :: Jh,g_tt,g_tz,g_zz              !! Jac,1/Jac,g_{ab} with a=theta/zeta b=theta/zeta
  REAL(wp),INTENT(OUT) :: Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1  !! and their variation vs q1
  REAL(wp),INTENT(OUT) :: Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2  !! and their variation vs q2
  REAL(wp),INTENT(OUT) :: g_t1,g_t2,g_z1,g_z2,Gh11,Gh22  !! dq^{i}/dtheta*G^{i1}, dq^{i}/dtheta*G^{i2}, and dq^{i}/dzeta*G^{i1}, dq^{i}/dzeta*G^{i2} and G^{11},G^{22}
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: Gh31,Gh32,Gh33
!===================================================================================================================================
  ASSOCIATE(  lp=>xv%lp, tau=>xv%tau, sigma=>xv%sigma, kappa=>xv%kappa)
  Gh11=1.0_wp
  !Gh21=0.0_wp
  Gh22=1.0_wp
  Gh31 =-lp*tau*q2
  Gh32 = lp*tau*q1

  !Jh=lp*(1.0_wp-sigma*kappa*q1)
  Jh_dq1=-lp*sigma*kappa
  Jh_dq2=0.0_wp

  Jh=lp+Jh_dq1*q1
  ! Gh33 = (lp**2)*((1.0_wp-sigma*kappa*q1)**2+tau**2*(q1**2+q2**2))
  Gh33 = Jh*Jh + Gh31*Gh31 + Gh32*Gh32

  g_t1 = dX1_dt
  g_t2 = dX2_dt
  g_z1 = dX1_dz + Gh31
  g_z2 = dX2_dz + Gh32

  g_tt =   dX1_dt *  g_t1         +  dX2_dt *  g_t2
  g_tz =   dX1_dt *  g_z1         +  dX2_dt *  g_z2
  g_zz =   dX1_dz * (g_z1 + Gh31) +  dX2_dz * (g_z2 + Gh32)  + Gh33

  !Gh11/dq1 =0 Gh12/dq1 =0 Gh13/dq1 = 0
  !            Gh22/dq1 =0 Gh23/dq1 = lp*tau
  !                        Gh33/dq1 = 2*(lp**2)*((1.0_wp-sigma*kappa*q1)*(-sigma*kappa)+tau**2 *(q1))
  !Gh11/dq2 =0 Gh12/dq2 =0 Gh13/dq2 = -lp*tau
  !            Gh22/dq2 =0 Gh23/dq2 = 0
  !                        Gh33/dq2 = 2*(lp*tau)**2*(q2)
  ! => g_t1 /dq1 =0, g_t1/dq2 =0, g_t2/dq1 =0, g_t2/dq2 =0
  ! => g_z1 /dq1 = Gh31/dq1, g_z1/dq2 =Gh31/dq2, g_z2/dq1 =Gh32/dq1, g_z2/dq2 =Gh32/dq2
  g_tt_dq1 = 0.0_wp
  g_tt_dq2 = 0.0_wp

  g_tz_dq1 =  lp*tau*dX2_dt
  g_tz_dq2 = -lp*tau*dX1_dt

  g_zz_dq1 =  2.0_wp*(lp*tau*(dX2_dz + Gh32)+Jh*Jh_dq1)
  g_zz_dq2 = -2.0_wp*lp*tau*(dX1_dz + Gh31)
  END ASSOCIATE
END SUBROUTINE hmap_frenet_eval_all_e


!===================================================================================================================================
!> sign function depending on zeta,
!! if omnig=False, sigma=1
!! if omnig=True, sigma=+1 for 0<=zeta<=pi/nfp, and -1 for pi/nfp<zeta<2pi
!!
!===================================================================================================================================
FUNCTION hmap_frenet_sigma(sf,zeta) RESULT(sigma)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: sigma
!===================================================================================================================================
  sigma=MERGE(SIGN(1.0_wp,SIN(sf%nfp*zeta)),1.0_wp,sf%omnig)
END FUNCTION hmap_frenet_sigma

!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z)
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval( sf ,q_in) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: x_out(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)          :: lp,absB,kappa,sigma,Jh
!===================================================================================================================================
  ! q(:) = (q1,q2,zeta) are the variables in the domain of the map
  ! X(:) = (x,y,z) are the variables in the range of the map
  !
  !  |x |
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |

  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)
  sigma=sf%sigma(zeta)
  Jh=lp*(1.0_wp-sigma*kappa*q1)
  IF(Jh.LT.1.0e-12) &
      CALL abort(__STAMP__, &
           "hmap_frenet, evaluation outside curvature radius (sigma*q1 >= 1./(kappa))",RealInfo=zeta*sf%nfp/TWOPI)
  B=B/absB
  N=CROSS(B,T)
  x_out=X0 +sigma*(q1*N + q2*B)
  END ASSOCIATE
END FUNCTION hmap_frenet_eval


!===================================================================================================================================
!> evaluate the mapping h (q1,q2,zeta) -> (x,y,z)
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_aux( sf ,q1,q2,xv) RESULT(x_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: x_out(3)
!===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  x_out=xv%X0 +xv%sigma*(q1*xv%N + q2*xv%B)
  END SELECT !type(xv)
END FUNCTION hmap_frenet_eval_aux


!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_dxdq( sf ,q_in,q_vec) RESULT(dxdq_qvec)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q_in(3)
  REAL(wp)            , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)          :: lp,absB,kappa,tau,sigma,Jh
!===================================================================================================================================
  !  |x |
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |
  !  dh/dq1 =sigma*N , dh/dq2=sigma*B
  !  dh/dq3 = l' [(1-sigma*kappa*q1)T + sigma*tau*(B*q1-N*q2) ]
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)

  sigma=sf%sigma(zeta)
  Jh=lp*(1.0_wp-sigma*kappa*q1)
  IF(Jh.LT.1.0e-12) &
      CALL abort(__STAMP__, &
           "hmap_frenet, evaluation outside curvature radius (sigma*q1 >= 1/(kappa))",RealInfo=zeta*sf%nfp/TWOPI)

  tau=SUM(X0ppp*B)/(absB**2)
  B=B/absB
  N=CROSS(B,T)
  dxdq_qvec(1:3)= sigma*(N*q_vec(1)+B*q_vec(2))+(Jh*T +sigma*lp*tau*(B*q1-N*q2))*q_vec(3)

  END ASSOCIATE !zeta
END FUNCTION hmap_frenet_eval_dxdq


!===================================================================================================================================
!> evaluate total derivative of the mapping  sum k=1,3 (dx(1:3)/dq^k) q_vec^k,
!! where dx(1:3)/dq^k, k=1,2,3 is evaluated at q_in=(X^1,X^2,zeta) ,
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_dxdq_aux( sf ,q1,q2,q1_vec,q2_vec,q3_vec,xv) RESULT(dxdq_qvec)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q1,q2
  REAL(wp)            , INTENT(IN) :: q1_vec,q2_vec,q3_vec
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: dxdq_qvec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)          :: Jh
!===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  Jh=xv%lp*(1.0_wp-xv%sigma*xv%kappa*q1)
  dxdq_qvec(1:3)= xv%sigma*(xv%N*q1_vec+xv%B*q2_vec) &
                  +(Jh*xv%T +xv%sigma*xv%lp*xv%tau*(xv%B*q1-xv%N*q2))*q3_vec
  END SELECT !type(xv)
END FUNCTION hmap_frenet_eval_dxdq_aux

!===============================================================================================================================
!> evaluate all first derivatives dx(1:3)/dq^i, i=1,2,3 , at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_frenet_get_dx_dqi( sf ,q_in,dx_dq1,dx_dq2,dx_dq3)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q_in(3)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT) :: dx_dq1,dx_dq2,dx_dq3
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,T,N,B
  REAL(wp)          :: lp,absB,kappa,tau,sigma
  !===================================================================================================================================
  !  |x |
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |
  !  dh/dq1 =sigma*N , dh/dq2=sigma*B
  !  dh/dq3 = l' [(1-sigma*kappa*q1)T + sigma*tau*(B*q1-N*q2) ]
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  tau=SUM(X0ppp*B)/(absB**2)
  B=B/absB
  N=CROSS(B,T)
  dx_dq1(1:3)= sigma*N
  dx_dq2(1:3)= sigma*B
  dx_dq3(1:3)=lp*((1.0_wp-sigma*kappa*q1)*T +sigma*tau*(B*q1-N*q2))
  END ASSOCIATE !zeta
END SUBROUTINE hmap_frenet_get_dx_dqi

!===============================================================================================================================
!> evaluate all first derivatives dx(1:3)/dq^i, i=1,2,3 , at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_frenet_get_dx_dqi_aux( sf ,q1,q2,xv,dx_dq1,dx_dq2,dx_dq3)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT) :: dx_dq1,dx_dq2,dx_dq3
  !=================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  dx_dq1(1:3)= xv%sigma*xv%N
  dx_dq2(1:3)= xv%sigma*xv%B
  dx_dq3(1:3)=xv%lp*((1.0_wp-xv%sigma*xv%kappa*q1)*xv%T +xv%sigma*xv%tau*(xv%B*q1-xv%N*q2))
  END SELECT !type(xv)
END SUBROUTINE hmap_frenet_get_dx_dqi_aux

!=================================================================================================================================
!> evaluate all second derivatives d^2x(1:3)/(dq^i dq^j), i,j=1,2,3 is evaluated at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_frenet_get_ddx_dqij( sf ,q_in,ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q_in(3)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT) :: ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33
  !-----------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,X0p4,Bloc,T,N,B
  REAL(wp)          :: lp,absB,kappa,tau,sigma
  REAL(wp)          :: lp_p,absB_p,kappa_p,tau_p
!===================================================================================================================================
  !  |x |
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp,X0p4=X0p4)
  lp=SQRT(SUM(X0p*X0p))
  T=X0p/lp
  Bloc=CROSS(X0p,X0pp)
  absB=SQRT(SUM(Bloc*Bloc))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  tau  = SUM(X0ppp*Bloc)/(absB**2)
  lp_p = SUM(X0pp*X0p) / lp
  absB_p = SUM(Bloc*CROSS(X0p, X0ppp)) / absB
  kappa_p = (absB_p*lp -3*absB * lp_p) / (lp**4)
  tau_p   = (SUM(X0p4*Bloc)*absB -2*SUM(X0ppp*Bloc)*absB_p) / (absB**3)
  B=Bloc/absB
  N=CROSS(B,T)
  ddx_dq11=0.0_wp
  ddx_dq12=0.0_wp
  ddx_dq13=sigma*lp*(-kappa*T +tau*B)

  ddx_dq22=0.0_wp
  ddx_dq23=-sigma*lp*tau*N
  ddx_dq33(1:3)= lp_p*((1.0_wp-sigma*kappa*q1)*T +sigma*tau*(B*q1-N*q2)) &
                 +lp*sigma*( -kappa_p*q1*T +tau_p*(B*q1-N*q2)                     &
                            +lp*( (1.0_wp-sigma*kappa*q1)*(kappa*N)             &
                                 +sigma*tau*((-tau*N)*q1-((-kappa*T +tau*B))*q2)) )
  END ASSOCIATE
END SUBROUTINE hmap_frenet_get_ddx_dqij

!=================================================================================================================================
!> evaluate all second derivatives d^2x(1:3)/(dq^i dq^j), i,j=1,2,3 is evaluated at q_in=(X^1,X^2,zeta),
!!
!===============================================================================================================================
SUBROUTINE hmap_frenet_get_ddx_dqij_aux( sf ,q1,q2,xv,ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33)
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL(wp),DIMENSION(3),INTENT(OUT) :: ddx_dq11,ddx_dq12,ddx_dq13,ddx_dq22,ddx_dq23,ddx_dq33
  !===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  ddx_dq11=0.0_wp
  ddx_dq12=0.0_wp
  ddx_dq13=xv%sigma*xv%lp*(-xv%kappa*xv%T +xv%tau*xv%B)
  ddx_dq22=0.0_wp
  ddx_dq23=-xv%sigma*xv%lp*xv%tau*xv%N
  ddx_dq33(1:3)= xv%lp_p*((1.0_wp-xv%sigma*xv%kappa*q1)*xv%T +xv%sigma*xv%tau*(xv%B*q1-xv%N*q2)) &
                 +xv%lp*xv%sigma*( -xv%kappa_p*q1*xv%T +xv%tau_p*(xv%B*q1-xv%N*q2)                             &
                                  +xv%lp*( (1.0_wp-xv%sigma*xv%kappa*q1)*(xv%kappa*xv%N)                      &
                                          +xv%sigma*xv%tau*((-xv%tau*xv%N)*q1-((-xv%kappa*xv%T +xv%tau*xv%B))*q2)) )
  END SELECT !type(xv)
END SUBROUTINE hmap_frenet_get_ddx_dqij_aux

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(q^1,q^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh( sf ,q_in) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: Jh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,B
  REAL(wp)          :: lp,absB,kappa,sigma
!===================================================================================================================================
  ASSOCIATE(q1=>q_in(1),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
  lp=SQRT(SUM(X0p*X0p))
  B=CROSS(X0p,X0pp)
  absB=SQRT(SUM(B*B))
  kappa=absB/(lp**3)
  IF(kappa.LT.1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet cannot evaluate frame at curvature < 1e-8 !",RealInfo=zeta*sf%nfp/TWOPI)
  sigma=sf%sigma(zeta)

  Jh=lp*(1.0_wp-sigma*kappa*q1)
  IF(Jh .LT. 1.0e-8) &
      CALL abort(__STAMP__, &
           "hmap_frenet, evaluation outside curvature radius, Jh<0",RealInfo=zeta*sf%nfp/TWOPI)

  END ASSOCIATE !zeta
END FUNCTION hmap_frenet_eval_Jh

!===================================================================================================================================
!> evaluate Jacobian of mapping h: J_h=sqrt(det(G)) at q=(q^1,q^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh_aux( sf ,q1,q2,xv) RESULT(Jh)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q1,q2
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: Jh
!===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  Jh=xv%lp*(1.0_wp-xv%sigma*xv%kappa*q1)
  END SELECT !type(xv)
END FUNCTION hmap_frenet_eval_Jh_aux

!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: sum_k q_vec^k * dJ_h/dq^k, k=1,2,3 at q=(q^1,q^2,zeta)
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh_dq( sf ,q_in,q_vec) RESULT(Jh_dq)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q_in(3)
  REAL(wp)            , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: Jh_dq
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,Bloc
  REAL(wp)          :: lp,absB,kappa,sigma
  REAL(wp)          :: lp_p,absB_p,kappa_p
!===================================================================================================================================
  !  |x |
  !  |y |=  X0(zeta) + sigma*(N(zeta)*q1 + B(zeta)*q2)
  !  |z |
  ASSOCIATE(q1=>q_in(1),q2=>q_in(2),zeta=>q_in(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
  lp=SQRT(SUM(X0p*X0p))
  Bloc=CROSS(X0p,X0pp)
  absB=SQRT(SUM(Bloc*Bloc))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  lp_p = SUM(X0pp*X0p) / lp
  absB_p = SUM(Bloc*CROSS(X0p, X0ppp)) / absB
  kappa_p = (absB_p*lp -3*absB * lp_p) / (lp**4)
  !tau_p   = (SUM(X0p4*Bloc)*absB -2*SUM(X0ppp*Bloc)*absB_p) / (absB**3)

  Jh_dq=-lp*sigma*kappa*q_vec(1) + (lp_p-sigma*q1*(lp_p*kappa+lp*kappa_p))*q_vec(3) !dsigma/dzeta is a dirac at kappa=0, so it is not evaluated
  END ASSOCIATE !zeta
END FUNCTION hmap_frenet_eval_Jh_dq

!===================================================================================================================================
!> evaluate derivative of Jacobian of mapping h: sum_k q_vec^k * dJ_h/dq^k, k=1,2,3 at q=(q^1,q^2,zeta)
!!
!! NOTE: needs auxvar with do_2nd_der=.TRUE.!! not checked for performance reasons.
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_Jh_dq_aux( sf ,q1,q2,q1_vec,q2_vec,q3_vec,xv) RESULT(Jh_dq)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: q1,q2
  REAL(wp)            , INTENT(IN) :: q1_vec,q2_vec,q3_vec
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: Jh_dq
!===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)

  Jh_dq=-xv%lp*xv%sigma*xv%kappa*q1_vec + (xv%lp_p-xv%sigma*q1*(xv%lp_p*xv%kappa+xv%lp*xv%kappa_p))*q3_vec
  END SELECT !type(xv)
END FUNCTION hmap_frenet_eval_Jh_dq_aux


!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij( sf ,qL_in,q_G,qR_in) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: qL_in(3)
  REAL(wp)            , INTENT(IN) :: q_G(3)
  REAL(wp)            , INTENT(IN) :: qR_in(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,Bloc
  REAL(wp)              :: lp,absB,kappa,tau,sigma
  REAL(wp)              :: Ga, Gb, Gc
!===================================================================================================================================
  ! A = -q2*l' * tau
  ! B =  q1*l' * tau
  ! C = Jh^2 + (l'*tau)^2(q1^2+q2^2)
  !                       |q1  |   |1   0   Ga|        |q1  |
  !q_i G_ij q_j = (dalpha |q2  | ) |0   1   Gb| (dbeta |q2  | )
  !                       |q3  |   |Ga  Gb  Gc|        |q3  |
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp)
  lp=SQRT(SUM(X0p*X0p))
  Bloc=CROSS(X0p,X0pp)
  absB=SQRT(SUM(Bloc*Bloc))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  tau=SUM(X0ppp*Bloc)/(absB**2)

  Ga = -lp*tau*q2
  Gb =  lp*tau*q1
  Gc = (lp**2)*((1.0_wp-sigma*kappa*q1)**2+tau**2*(q1**2+q2**2))
  g_ab=      qL_in(1)*qR_in(1) &
            +qL_in(2)*qR_in(2) &
       + Gc* qL_in(3)*qR_in(3) &
       + Ga*(qL_in(1)*qR_in(3)+qL_in(3)*qR_in(1)) &
       + Gb*(qL_in(2)*qR_in(3)+qL_in(3)*qR_in(2))
  END ASSOCIATE
END FUNCTION hmap_frenet_eval_gij

!===================================================================================================================================
!>  evaluate sum_ij (qL_i (G_ij(q_G)) qR_j) ,,
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij_aux( sf ,qL1,qL2,qL3,q1,q2,qR1,qR2,qR3,xv) RESULT(g_ab)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: qL1,qL2,qL3
  REAL(wp)            , INTENT(IN) :: q1,q2
  REAL(wp)            , INTENT(IN) :: qR1,qR2,qR3
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                        :: g_ab
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)              :: Ga, Gb, Gc
!===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  Ga = -xv%lp*xv%tau*q2
  Gb =  xv%lp*xv%tau*q1
  !Gc = (xv%lp**2)*((1.0_wp-xv%sigma*xv%kappa*q1)**2+xv%tau**2*(q1**2+q2**2))
  Gc = (xv%lp*(1.0_wp-xv%sigma*xv%kappa*q1))**2+ Ga*Ga+ Gb*Gb
  g_ab=      qL1*qR1 &
            +qL2*qR2 &
       + Gc* qL3*qR3 &
       + Ga*(qL1*qR3+qL3*qR1) &
       + Gb*(qL2*qR3+qL3*qR2)
  END SELECT !type(xv)
END FUNCTION hmap_frenet_eval_gij_aux

!===================================================================================================================================
!>  evaluate sum_k=1,3 q_vec^k * sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij_dq( sf ,qL_in,q_G,qR_in,q_vec) RESULT(g_ab_dq)
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: qL_in(3)
  REAL(wp)            , INTENT(IN) :: q_G(3)
  REAL(wp)            , INTENT(IN) :: qR_in(3)
  REAL(wp)            , INTENT(IN) :: q_vec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: g_ab_dq
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(3) :: X0,X0p,X0pp,X0ppp,X0p4,Bloc
  REAL(wp)              :: lp,absB,kappa,tau,sigma
  REAL(wp)              :: lp_p,absB_p,kappa_p,tau_p
!===================================================================================================================================
  !                       |q1  |   |0  0        0           |        |q1  |
  !q_i G_ij q_j = (dalpha |q2  | ) |0  0      l'*tau        | (dbeta |q2  | )
  !                       |q3  |   |0  l'*tau  dG33/dq1     |        |q3  |
  ASSOCIATE(q1=>q_G(1),q2=>q_G(2),zeta=>q_G(3))
  CALL sf%eval_X0(zeta,X0,X0p,X0pp,X0ppp,X0p4=X0p4)
  lp=SQRT(SUM(X0p*X0p))
  Bloc=CROSS(X0p,X0pp)
  absB=SQRT(SUM(Bloc*Bloc))
  kappa=absB/(lp**3)
  sigma=sf%sigma(zeta)
  tau=SUM(X0ppp*Bloc)/(absB**2)
  lp_p = SUM(X0pp*X0p) / lp
  absB_p = SUM(Bloc*CROSS(X0p, X0ppp)) / absB
  kappa_p = (absB_p*lp -3*absB * lp_p) / (lp**4)
  tau_p   = (SUM(X0p4*Bloc)*absB -2*SUM(X0ppp*Bloc)*absB_p) / (absB**3)

  !G13 = -lp*tau*q2
  !G23 =  lp*tau*q1
  !G33 = (lp**2)*((1.0_wp-sigma*kappa*q1)**2+tau**2*(q1**2+q2**2))

  g_ab_dq =-(lp*tau*q_vec(2)+(lp_p*tau+lp*tau_p)*q2*q_vec(3))*(qL_in(1)*qR_in(3)+ qL_in(3)*qR_in(1)) &
           +(lp*tau*q_vec(1)+(lp_p*tau+lp*tau_p)*q1*q_vec(3))*(qL_in(2)*qR_in(3)+ qL_in(3)*qR_in(2)) &
           +2.0_wp*(qL_in(3)*qR_in(3))*( (lp**2)*( q_vec(1)*((tau**2+kappa**2)*q1-sigma*kappa)       &
                                                  +q_vec(2)*  tau**2*q2                        )     &
                                        +q_vec(3)*( lp_p*lp* ( (1.0_wp-sigma*kappa*q1)**2+tau**2*(q1**2+q2**2)) &
                                                   +(lp**2)* ( (1.0_wp-sigma*kappa*q1)*(-sigma*kappa_p*q1)      &
                                                              +tau*tau_p*(q1**2+q2**2) ) ) )
  END ASSOCIATE
END FUNCTION hmap_frenet_eval_gij_dq

!===================================================================================================================================
!>  evaluate sum_ij (qL_i d/dq^k(G_ij(q_G)) qR_j) , k=1,2
!! where qL=(dX^1/dalpha,dX^2/dalpha [,dzeta/dalpha]) and qR=(dX^1/dbeta,dX^2/dbeta [,dzeta/dbeta]) and
!! where qL=(dX^1/dalpha,dX^2/dalpha ,dzeta/dalpha) and qR=(dX^1/dbeta,dX^2/dbeta ,dzeta/dbeta) and
!! dzeta_dalpha then known to be either 0 of ds and dtheta and 1 for dzeta
!!
!! NOTE: needs auxvar with do_2nd_der=.TRUE.!! not checked for performance reasons.
!!
!===================================================================================================================================
FUNCTION hmap_frenet_eval_gij_dq_aux( sf ,qL1,qL2,qL3,q1,q2,qR1,qR2,qR3,q1_vec,q2_vec,q3_vec,xv) RESULT(g_ab_dq)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN) :: sf
  REAL(wp)            , INTENT(IN) :: qL1,qL2,qL3
  REAL(wp)            , INTENT(IN) :: q1,q2
  REAL(wp)            , INTENT(IN) :: qR1,qR2,qR3
  REAL(wp)            , INTENT(IN) :: q1_vec,q2_vec,q3_vec
  CLASS(c_hmap_auxvar), INTENT(IN) :: xv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                         :: g_ab_dq

!===================================================================================================================================
  SELECT TYPE(xv); TYPE IS(t_hmap_frenet_auxvar)
  g_ab_dq =-(xv%lp*xv%tau*q2_vec+(xv%lp_p*xv%tau+xv%lp*xv%tau_p)*q2*q3_vec)*(qL1*qR3+ qL3*qR1) &
           +(xv%lp*xv%tau*q1_vec+(xv%lp_p*xv%tau+xv%lp*xv%tau_p)*q1*q3_vec)*(qL2*qR3+ qL3*qR2) &
          +2.0_wp*(qL3*qR3)*( (xv%lp**2)*( q1_vec*((xv%tau**2+xv%kappa**2)*q1-xv%sigma*xv%kappa)           &
                                          +q2_vec*  xv%tau**2*q2                        )                  &
                             +q3_vec*( xv%lp_p*xv%lp* ( (1.0_wp-xv%sigma*xv%kappa*q1)**2+xv%tau**2*(q1**2+q2**2))   &
                                      +(xv%lp**2)*    ( (1.0_wp-xv%sigma*xv%kappa*q1)*(-xv%sigma*xv%kappa_p*q1)  &
                                                       +xv%tau*xv%tau_p*(q1**2+q2**2) ) ) )
  END SELECT ! TYPE(xv)
END FUNCTION hmap_frenet_eval_gij_dq_aux


!===================================================================================================================================
!> evaluate curve X0(zeta), position and first three derivatives, from given R0,Z0 Fourier
!!
!===================================================================================================================================
PURE SUBROUTINE hmap_frenet_eval_X0_fromRZ( sf,zeta,X0,X0p,X0pp,X0ppp,X0p4)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(IN ) :: sf
  REAL(wp)            , INTENT(IN ) :: zeta       !! position along closed curve parametrized in [0,2pi]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)            , INTENT(OUT) :: X0(1:3)      !! curve position in cartesian coordinates
  REAL(wp)            , INTENT(OUT) :: X0p(1:3)     !! 1st derivative in zeta
  REAL(wp)            , INTENT(OUT) :: X0pp(1:3)    !! 2nd derivative in zeta
  REAL(wp)            , INTENT(OUT) :: X0ppp(1:3)   !! 3rd derivative in zeta
  REAL(wp)            , INTENT(OUT),OPTIONAL :: X0p4(1:3)  !! 4th derivative in zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: R0,R0p,R0pp,R0ppp,R0p4,Z0p4
  REAL(wp) :: coszeta,sinzeta
!===================================================================================================================================
  CALL eval_fourier1d(sf%n_max,sf%Xn,sf%rc,sf%rs,zeta,R0,R0p,R0pp,R0ppp,R0p4)
  CALL eval_fourier1d(sf%n_max,sf%Xn,sf%zc,sf%zs,zeta,X0(3),X0p(3),X0pp(3),X0ppp(3),Z0p4) !=Z0,Z0p,Z0pp,Z0ppp
  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  ASSOCIATE(x   =>X0(1)   ,y   =>X0(2)   , &
            xp  =>X0p(1)  ,yp  =>X0p(2)  , &
            xpp =>X0pp(1) ,ypp =>X0pp(2) , &
            xppp=>X0ppp(1),yppp=>X0ppp(2))
    !! angle zeta=geometric toroidal angle phi=atan(y/x)
    x=R0*coszeta
    y=R0*sinzeta

    xp = R0p*coszeta  - R0*sinzeta
    yp = R0p*sinzeta  + R0*coszeta
    !xp  = R0p*coszeta  -y
    !yp  = R0p*sinzeta  +x

    xpp = R0pp*coszeta - 2*R0p*sinzeta - R0*coszeta
    ypp = R0pp*sinzeta + 2*R0p*coszeta - R0*sinzeta
    !xpp  = R0pp*coszeta -2.0_wp*yp + x
    !ypp  = R0pp*sinzeta +2.0_wp*xp + y

    xppp = R0ppp*coszeta - 3*R0pp*sinzeta - 3*R0p*coszeta + R0*sinzeta
    yppp = R0ppp*sinzeta + 3*R0pp*coszeta - 3*R0p*sinzeta - R0*coszeta
    !xppp  = R0ppp*coszeta +3.0_wp*(xp-ypp) + y
    !yppp  = R0ppp*sinzeta +3.0_wp*(yp+xpp) + x
  IF(PRESENT(X0p4))THEN
    X0p4(1)  = R0p4*coszeta - 4*R0ppp*sinzeta - 6*R0pp*coszeta  + 4*R0p*sinzeta + R0*coszeta
    X0p4(2)  = R0p4*sinzeta + 4*R0ppp*coszeta - 6*R0pp*sinzeta  - 4*R0p*coszeta + R0*sinzeta
    X0p4(3)  = Z0p4
  END IF
  END ASSOCIATE !x,y,xp,yp,...

END SUBROUTINE hmap_frenet_eval_X0_fromRZ


!===================================================================================================================================
!> evaluate 1d fourier series from given cos/sin coefficients and mode numbers xn
!! SUM(xc(0:n_max)*COS(xn(0:n_max)*zeta)+xs(0:n_max)*SIN(xn(0:n_max)*zeta)
!! evaluate all derivatives 1,2,3 alongside
!!
!===================================================================================================================================
PURE SUBROUTINE eval_fourier1d(n_max,xn,xc,xs,zeta,x,xp,xpp,xppp,xp4)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER  , INTENT(IN ) :: n_max        !! number of modes is n_max+1  (0...n_max)
  INTEGER  , INTENT(IN ) :: xn(0:n_max)  !! array of mode numbers
  REAL(wp) , INTENT(IN ) :: xc(0:n_max)  !! cosine coefficients
  REAL(wp) , INTENT(IN ) :: xs(0:n_max)  !!   sine coefficients
  REAL(wp) , INTENT(IN ) :: zeta         !! angular position [0,2pi]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp) , INTENT(OUT) :: x      !! value at zeta
  REAL(wp) , INTENT(OUT) :: xp     !! 1st derivative in zeta
  REAL(wp) , INTENT(OUT) :: xpp    !! 2nd derivative in zeta
  REAL(wp) , INTENT(OUT) :: xppp   !! 3rd derivative in zeta
  REAL(wp) , INTENT(OUT) :: xp4    !! 4th derivative in zeta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp),DIMENSION(0:n_max) :: cos_nzeta,sin_nzeta,xtmp,xptmp
!===================================================================================================================================
  cos_nzeta=COS(REAL(xn,wp)*zeta)
  sin_nzeta=SIN(REAL(xn,wp)*zeta)
  xtmp = xc*cos_nzeta+xs*sin_nzeta
  xptmp= REAL(xn,wp)*(-xc*sin_nzeta+xs*cos_nzeta)
  x    = SUM(xtmp)
  xp   = SUM(xptmp)
  xpp  = SUM(REAL(-xn*xn,wp)*xtmp)
  xppp = SUM(REAL(-xn*xn,wp)*xptmp)
  xp4  = SUM(REAL(xn**4,wp)*xtmp)

END SUBROUTINE eval_fourier1d


!===================================================================================================================================
!> test hmap_frenet - evaluation of the map
!!
!===================================================================================================================================
SUBROUTINE hmap_frenet_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_hmap_frenet), INTENT(INOUT) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,idir,jdir,kdir,qdir,izeta,i,j,k,ndims(1:3),ijk(3)
  INTEGER,PARAMETER  :: nzeta=5
  INTEGER,PARAMETER  :: ns=2
  INTEGER,PARAMETER  :: nthet=3
  REAL(wp)           :: refreal,checkreal,x(3),q_in(3),q_test(3,3),x_eps(3),dxdq(3),gij,gij_eps,Jh_0,Jh_eps
  REAL(wp),ALLOCATABLE :: zeta(:)
  REAL(wp)           :: qloc(3),q_thet(3),q_zeta(3)
  REAL(wp)           :: dxdq_eps(3),dx_dqi(3,3),ddx_dqij(3,3,3)
  REAL(wp),ALLOCATABLE,DIMENSION(:,:,:) :: q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
                                     Jh,g_tt,    g_tz,    g_zz,     &
                                     Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1, &
                                     Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2, &
                                     g_t1,g_t2,g_z1,g_z2,Gh11,Gh22
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  REAL(wp),PARAMETER :: epsFD=1.0e-8
  REAL(wp),PARAMETER :: realtolFD=1.0e-4
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: R0, Z0
  TYPE(t_hmap_frenet_auxvar),ALLOCATABLE :: xv(:)
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN hmap_frenet TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN

    !evaluate on the axis q1=q2=0
    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    q_in=(/0.0_wp, 0.0_wp, 0.335_wp*PI/)
    R0 = SUM(sf%rc(:)*COS(sf%Xn(:)*q_in(3)) + sf%rs(:)*SIN(sf%Xn(:)*q_in(3)))
    Z0 = SUM(sf%zc(:)*COS(sf%Xn(:)*q_in(3)) + sf%zs(:)*SIN(sf%Xn(:)*q_in(3)))
    x = sf%eval(q_in )
    checkreal=SUM((x-(/R0*COS(q_in(3)),R0*SIN(q_in(3)),Z0/))**2)
    refreal = 0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
            '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
       nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3))') &
     '\n =>  should be ', refreal,' : |y-eval_map(x)|^2= ', checkreal
    END IF !TEST

    q_test(1,:)=(/1.0_wp, 0.0_wp, 0.0_wp/)
    q_test(2,:)=(/0.0_wp, 1.0_wp, 0.0_wp/)
    q_test(3,:)=(/0.0_wp, 0.0_wp, 1.0_wp/)
    DO qdir=1,3
      !check dx/dq^i with FD
      iTest=101+qdir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      q_in=(/0.0_wp, 0.0_wp, 0.335_wp*PI/)
      x = sf%eval(q_in )
      x_eps = sf%eval(q_in+epsFD*q_test(qdir,:))
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      checkreal=SQRT(SUM((dxdq - (x_eps-x)/epsFD)**2)/SUM(x*x))
      refreal = 0.0_wp

      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be <',realtolFD,' : |dxdqFD-eval_dxdq|= ', checkreal,", dq=",qdir
      END IF !TEST
    END DO

    !! TEST G_ij
    DO idir=1,3; DO jdir=idir,3
      iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      checkreal= SUM(sf%eval_dxdq(q_in,q_test(idir,:))*sf%eval_dxdq(q_in,q_test(jdir,:)))
      refreal  =sf%eval_gij(q_test(idir,:),q_in,q_test(jdir,:))
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),2(A,I3))') &
       '\n =>  should be ', refreal,' : sum|Gij-eval_gij|= ', checkreal,', i=',idir,', j=',jdir
      END IF !TEST
    END DO; END DO
    !! TEST dJh_dqk with FD
    DO qdir=1,3
      iTest=iTest+1; IF(testdbg)WRITE(*,*)'iTest=',iTest
      Jh_0    = sf%eval_Jh(   q_in                     )
      Jh_eps  = sf%eval_Jh(   q_in+epsFD*q_test(qdir,:))
      refreal = sf%eval_Jh_dq(q_in                     ,q_test(qdir,:))
      checkreal=(Jh_eps-Jh_0)/epsFD-refreal
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be <', realtolFD,' : |dJh_dqFD-eval_Jh_dq|= ', checkreal,', dq=',qdir
      END IF !TEST
    END DO !qdir
    !! TEST dG_ij_dqk with FD
    DO qdir=1,3
    DO idir=1,3; DO jdir=1,3
      iTest=iTest+1; IF(testdbg)WRITE(*,*)'iTest=',iTest
      gij     = sf%eval_gij(   q_test(idir,:),q_in                     ,q_test(jdir,:))
      gij_eps = sf%eval_gij(   q_test(idir,:),q_in+epsFD*q_test(qdir,:),q_test(jdir,:))
      refreal = sf%eval_gij_dq(q_test(idir,:),q_in                     ,q_test(jdir,:),q_test(qdir,:))
      checkreal=(gij_eps-gij)/epsFD-refreal
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),3(A,I3))') &
       '\n =>  should be <', realtolFD,' : |dGij_dqFD-eval_gij_dq|= ', checkreal,', i=',idir,', j=',jdir,', dq=',qdir
      END IF !TEST
    END DO; END DO
    END DO !qdir

    CALL sf%get_dx_dqi(q_in,dx_dqi(1,:),dx_dqi(2,:),dx_dqi(3,:))
    DO qdir=1,3
      !check dx/dq^i with FD
      iTest=10+qdir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      checkreal=SUM(ABS(dxdq-dx_dqi(qdir,:)))
      refreal=0.0_wp
      IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
              '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
         nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),(A,I3))') &
       '\n =>  should be ', refreal,' : |dxdq-eval_dxdq|= ', checkreal,", dq=",qdir
      END IF !TEST
    END DO !qdir
    CALL sf%get_ddx_dqij(q_in,ddx_dqij(1,1,:),ddx_dqij(1,2,:),ddx_dqij(1,3,:), &
                              ddx_dqij(2,2,:),ddx_dqij(2,3,:),ddx_dqij(3,3,:))
    ddx_dqij(2,1,:)=ddx_dqij(1,2,:)
    ddx_dqij(3,1,:)=ddx_dqij(1,3,:)
    ddx_dqij(3,2,:)=ddx_dqij(2,3,:)
    DO qdir=1,3
      dxdq = sf%eval_dxdq(q_in,q_test(qdir,:))
      DO idir=1,3
        iTest=20+idir+3*(qdir-1) ; IF(testdbg)WRITE(*,*)'iTest=',iTest
        dxdq_eps = sf%eval_dxdq(q_in+epsFD*q_test(idir,:),q_test(qdir,:))
        checkreal=SUM(ABS( (dxdq_eps-dxdq)/epsFD-ddx_dqij(qdir,idir,:)))
        refreal=0.0_wp
        IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtolFD))) THEN
           nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
                '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
           nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),2(A,I3))') &
         '\n =>  should be <', realtolFD,' : |ddx_dqijFD-eval_ddx_dqij|= ', checkreal,", dqi=",qdir,", dqj=",idir
        END IF !TEST
      END DO !idir
    END DO !qdir

 END IF !testlevel >=1
 IF (testlevel .GE. 2) THEN
  DO idir=1,3
    SELECT CASE(idir)
    CASE(1)
      jdir=2; kdir=3
    CASE(2)
      jdir=1; kdir=3
    CASE(3)
      jdir=1; kdir=2
    END SELECT
    ndims(idir)=nzeta+idir
    ndims(jdir)=ns
    ndims(kdir)=nthet
    ALLOCATE(zeta(ndims(idir)),xv(ndims(idir)))
    DO izeta=1,ndims(idir)
      zeta(izeta)=0.333_wp+REAL(izeta-1,wp)/REAL(ndims(idir)-1,wp)*0.221_wp
      xv(izeta)=hmap_frenet_init_aux(sf,zeta(izeta),.TRUE.)
    END DO
    ALLOCATE(q1(ndims(1),ndims(2),ndims(3)))
    ALLOCATE(q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz,Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2,g_t1,g_t2,g_z1,g_z2,Gh11,Gh22, &
             mold=q1)
    !assign somewhat randomly
    DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
      q1(i,j,k) = 0.11_wp -0.22_wp *REAL((i+j)*k,wp)/REAL((ndims(idir)+ndims(jdir))*ndims(kdir),wp)
      q2(i,j,k) = 0.15_wp -0.231_wp*REAL((i+k)*j,wp)/REAL((ndims(idir)+ndims(kdir))*ndims(jdir),wp)
      dX1_dt(i,j,k)=-0.1_wp  +0.211_wp*REAL((i+2*j)*k,wp)/REAL((ndims(idir)+2*ndims(jdir))*ndims(kdir),wp)
      dX2_dt(i,j,k)= 0.231_wp-0.116_wp*REAL((2*i+k)*j,wp)/REAL((2*ndims(idir)+ndims(kdir))*ndims(jdir),wp)
      dX1_dz(i,j,k)=-0.024_wp+0.013_wp*REAL((3*i+2*j)*k,wp)/REAL((3*ndims(idir)+2*ndims(jdir))*ndims(kdir),wp)
      dX2_dz(i,j,k)=-0.06_wp +0.031_wp*REAL((2*k+3*k)*i,wp)/REAL((2*ndims(kdir)+3*ndims(kdir))*ndims(idir),wp)
    END DO; END DO; END DO
    CALL sf%eval_all(ndims,idir,xv,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
         Jh,g_tt,g_tz,g_zz,&
         Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,&
         Jh_dq2,g_tt_dq2,g_tz_dq2,g_zz_dq2,&
         g_t1,g_t2,g_z1,g_z2,Gh11,Gh22)
    DO k=1,ndims(3); DO j=1,ndims(2); DO i=1,ndims(1)
      ijk=(/i,j,k/)
      izeta=ijk(idir)
      qloc=(/q1(i,j,k),q2(i,j,k),zeta(izeta)/)
      q_thet=(/dX1_dt(i,j,k),dX2_dt(i,j,k),0.0_wp/)
      q_zeta=(/dX1_dz(i,j,k),dX2_dz(i,j,k),1.0_wp/)
      Jh(i,j,k)       =Jh(i,j,k)       - sf%eval_Jh(qloc)
      g_tt(i,j,k)     =g_tt(i,j,k)     - sf%eval_gij(q_thet,qloc,q_thet)
      g_tz(i,j,k)     =g_tz(i,j,k)     - sf%eval_gij(q_thet,qloc,q_zeta)
      g_zz(i,j,k)     =g_zz(i,j,k)     - sf%eval_gij(q_zeta,qloc,q_zeta)
      Jh_dq1(i,j,k)   =Jh_dq1(i,j,k)   - sf%eval_Jh_dq(qloc,(/1.0_wp,0.0_wp,0.0_wp/))
      Jh_dq2(i,j,k)   =Jh_dq2(i,j,k)   - sf%eval_Jh_dq(qloc,(/0.0_wp,1.0_wp,0.0_wp/))
      g_tt_dq1(i,j,k) =g_tt_dq1(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_thet,(/1.0_wp,0.0_wp,0.0_wp/))
      g_tt_dq2(i,j,k) =g_tt_dq2(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_thet,(/0.0_wp,1.0_wp,0.0_wp/))
      g_tz_dq1(i,j,k) =g_tz_dq1(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_zeta,(/1.0_wp,0.0_wp,0.0_wp/))
      g_tz_dq2(i,j,k) =g_tz_dq2(i,j,k) - sf%eval_gij_dq(q_thet,qloc,q_zeta,(/0.0_wp,1.0_wp,0.0_wp/))
      g_zz_dq1(i,j,k) =g_zz_dq1(i,j,k) - sf%eval_gij_dq(q_zeta,qloc,q_zeta,(/1.0_wp,0.0_wp,0.0_wp/))
      g_zz_dq2(i,j,k) =g_zz_dq2(i,j,k) - sf%eval_gij_dq(q_zeta,qloc,q_zeta,(/0.0_wp,1.0_wp,0.0_wp/))
      g_t1(i,j,k)     =g_t1(i,j,k)     - sf%eval_gij(q_thet,qloc,(/1.0_wp,0.0_wp,0.0_wp/))
      g_t2(i,j,k)     =g_t2(i,j,k)     - sf%eval_gij(q_thet,qloc,(/0.0_wp,1.0_wp,0.0_wp/))
      g_z1(i,j,k)     =g_z1(i,j,k)     - sf%eval_gij(q_zeta,qloc,(/1.0_wp,0.0_wp,0.0_wp/))
      g_z2(i,j,k)     =g_z2(i,j,k)     - sf%eval_gij(q_zeta,qloc,(/0.0_wp,1.0_wp,0.0_wp/))
      Gh11(i,j,k)     =Gh11(i,j,k)     - sf%eval_gij((/1.0_wp,0.0_wp,0.0_wp/),qloc,(/1.0_wp,0.0_wp,0.0_wp/))
      Gh22(i,j,k)     =Gh22(i,j,k)     - sf%eval_gij((/0.0_wp,1.0_wp,0.0_wp/),qloc,(/0.0_wp,1.0_wp,0.0_wp/))
    END DO; END DO; END DO

    iTest=201+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(Jh))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|Jh_all-eval_Jh(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=202+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_tt))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_tt_all-eval_g_tt(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=203+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_tz))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_tz_all-eval_g_tz(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=203+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_zz))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_zz_all-eval_g_zz(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=204+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(Jh_dq1))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|Jh_dq1_all-eval_Jh_dq1(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=205+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(Jh_dq2))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|Jh_dq2_all-eval_Jh_dq2(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=206+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_tt_dq1))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_tt_dq1_all-eval_g_tt_dq1(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=207+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_tz_dq1))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_tz_dq1_all-eval_g_tz_dq1(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=208+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_zz_dq1))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_zz_dq1_all-eval_g_zz_dq1(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=209+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_tt_dq2))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_tt_dq2_all-eval_g_tt_dq2(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=210+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_tz_dq2))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_tz_dq2_all-eval_g_tz_dq2(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=211+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_zz_dq2))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_zz_dq2_all-eval_g_zz_dq2(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=212+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_t1))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_t1_all-eval_g_t1(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=213+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_t2))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_t2_all-eval_g_t2(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=214+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_z1))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_z1_all-eval_g_z1(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=215+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(g_z2))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|g_z2_all-eval_g_z2(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=216+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(Gh11))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|Gh11_all-eval_Gh11(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    iTest=217+20*idir ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal=SUM(ABS(Gh22))/REAL(PRODUCT(ndims),wp)
    refreal=0.0_wp
    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
           '\n!! hmap_frenet TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,E11.3),A,I4)') &
    '\n =>  should be ', refreal,' : |sum(|Gh22_all-eval_Gh22(xall)|)|= ', checkreal, " ,idir=",idir
    END IF

    DEALLOCATE(zeta,q1,q2,dX1_dt,dX2_dt,dX1_dz,dX2_dz, &
               Jh,g_tt,g_tz,g_zz,Jh_dq1,g_tt_dq1,g_tz_dq1,g_zz_dq1,Jh_dq2,&
               g_tt_dq2,g_tz_dq2,g_zz_dq2,g_t1,g_t2,g_z1,g_z2,Gh11,Gh22)
    DEALLOCATE(xv)
  END DO !idir
END IF

 test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE hmap_frenet_test

END MODULE MODgvec_hmap_frenet
