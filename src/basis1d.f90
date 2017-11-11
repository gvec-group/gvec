!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/hopr)
!
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

!==================================================================================================================================
!>
!!# Module ** Basis 1D **
!!
!! Routines to provide and evaluate 1D polynomial Lagrange basis functions, interpolation and integration points
!!
!==================================================================================================================================
MODULE MOD_Basis1D
! MODULES
USE MOD_Globals, ONLY: wp
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE BuildLegendreVdm
   MODULE PROCEDURE BuildLegendreVdm
END INTERFACE

INTERFACE InitializeVandermonde
   MODULE PROCEDURE InitializeVandermonde
END INTERFACE

INTERFACE ChebyshevGaussNodesAndWeights
   MODULE PROCEDURE ChebyshevGaussNodesAndWeights
END INTERFACE

INTERFACE ChebyGaussLobNodesAndWeights
   MODULE PROCEDURE ChebyGaussLobNodesAndWeights
END INTERFACE

INTERFACE ClenshawCurtisNodesAndWeights
   MODULE PROCEDURE ClenshawCurtisNodesAndWeights
END INTERFACE

INTERFACE LegendreGaussNodesAndWeights
   MODULE PROCEDURE LegendreGaussNodesAndWeights
END INTERFACE

INTERFACE LegGaussLobNodesAndWeights
   MODULE PROCEDURE LegGaussLobNodesAndWeights
END INTERFACE

INTERFACE LegendrePolynomialAndDerivative
   MODULE PROCEDURE LegendrePolynomialAndDerivative
END INTERFACE

INTERFACE PolynomialDerivativeMatrix
   MODULE PROCEDURE PolynomialDerivativeMatrix
END INTERFACE

INTERFACE BarycentricWeights
   MODULE PROCEDURE BarycentricWeights
END INTERFACE

INTERFACE LagrangeInterpolationPolys
   MODULE PROCEDURE LagrangeInterpolationPolys
END INTERFACE

INTERFACE EQUALTOTOLERANCE
   MODULE PROCEDURE EQUALTOTOLERANCE
END INTERFACE

PUBLIC::INV
PUBLIC::INV33
PUBLIC::BuildLegendreVdm
PUBLIC::InitializeVandermonde
PUBLIC::LegGaussLobNodesAndWeights
PUBLIC::LegendreGaussNodesAndWeights
PUBLIC::ChebyshevGaussNodesAndWeights
PUBLIC::ChebyGaussLobNodesAndWeights
PUBLIC::ClenshawCurtisNodesAndWeights
PUBLIC::LegendrePolynomialAndDerivative
PUBLIC::PolynomialDerivativeMatrix
PUBLIC::BarycentricWeights
PUBLIC::LagrangeInterpolationPolys
PUBLIC::EQUALTOTOLERANCE

!==================================================================================================================================
REAL(wp),PARAMETER        :: PP_RealTolerance = EPSILON(1.0_wp) !! machine precision
REAL(wp),PARAMETER        :: PP_Pi = ACOS(-1.0_wp)              !! Pi up to machine accuracy


CONTAINS


!==================================================================================================================================
!> Computes matrix inverse using LAPACK
!> Input matrix should be a square matrix
!==================================================================================================================================
FUNCTION INV(A) RESULT(AINV)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(wp),INTENT(IN)  :: A(:,:)                      !! input matrix
REAL(wp)             :: AINV(SIZE(A,1),SIZE(A,2))   !! result: inverse of A
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI
! LOCAL VARIABLES
REAL(wp)    :: work(SIZE(A,1))  ! work array for LAPACK
INTEGER :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER :: n,info
!==================================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
CALL DGETRF(n, n, Ainv, n, ipiv, info)

IF(info.NE.0)THEN
   STOP 'Matrix is numerically singular!'
END IF

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

IF(info.NE.0)THEN
   STOP 'Matrix inversion failed!'
END IF
END FUNCTION INV

!===================================================================================================================================
!> Computes the inverse of a 3x3 matrix
!===================================================================================================================================
SUBROUTINE INV33(M,MInv,detM_out)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)     :: M(3,3)  !! input 3x3 matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp),INTENT(OUT)    :: MInv(3,3) !!inverse 3x3
REAL(wp),INTENT(OUT),OPTIONAL ::detM_out !!determinant of matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL(wp) ::detM !!determinant of matrix
!===================================================================================================================================
detM =   M(1,1)*M(2,2)*M(3,3)  &
       - M(1,1)*M(2,3)*M(3,2)  &
       - M(1,2)*M(2,1)*M(3,3)  &
       + M(1,2)*M(2,3)*M(3,1)  &
       + M(1,3)*M(2,1)*M(3,2)  &
       - M(1,3)*M(2,2)*M(3,1)

IF(PRESENT(detM_out)) detM_out=detM
IF(ABS(detM).LE.1.0E-12_wp)THEN
   MInv = 0.0_wp
   RETURN
END IF

MInv(1,1) =  (M(2,2)*M(3,3)-M(2,3)*M(3,2))
MInv(2,1) = -(M(2,1)*M(3,3)-M(2,3)*M(3,1))
MInv(3,1) =  (M(2,1)*M(3,2)-M(2,2)*M(3,1))
MInv(1,2) = -(M(1,2)*M(3,3)-M(1,3)*M(3,2))
MInv(2,2) =  (M(1,1)*M(3,3)-M(1,3)*M(3,1))
MInv(3,2) = -(M(1,1)*M(3,2)-M(1,2)*M(3,1))
MInv(1,3) =  (M(1,2)*M(2,3)-M(1,3)*M(2,2))
MInv(2,3) = -(M(1,1)*M(2,3)-M(1,3)*M(2,1))
MInv(3,3) =  (M(1,1)*M(2,2)-M(1,2)*M(2,1))
MInv=MInv/detM

END SUBROUTINE INV33


!==================================================================================================================================
!> Build a 1D Vandermonde matrix from an orthonormal Legendre basis to a nodal basis and reverse
!==================================================================================================================================
SUBROUTINE buildLegendreVdm(N_In,xi_In,Vdm_Leg,sVdm_Leg)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_In                    !! input polynomial degree
REAL(wp),INTENT(IN)    :: xi_In(0:N_In)           !! nodal positions [-1,1]
REAL(wp),INTENT(OUT)   ::  Vdm_Leg(0:N_In,0:N_In) !! Vandermonde from Legendre to nodal basis
REAL(wp),INTENT(OUT)   :: sVdm_Leg(0:N_In,0:N_In) !! Vandermonde from nodal basis to Legendre
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j
REAL(wp)               :: dummy
!REAL(wp)               :: wBary_Loc(0:N_In)
!REAL(wp)               :: xGauss(0:N_In),wGauss(0:N_In)
!==================================================================================================================================
! Alternative to matrix inversion: Compute inverse Vandermonde directly
! Direct inversion seems to be more accurate

!CALL BarycentricWeights(N_In,xi_in,wBary_loc)
!! Compute first the inverse (by projection)
!CALL LegendreGaussNodesAndWeights(N_In,xGauss,wGauss)
!!Vandermonde on xGauss
!DO i=0,N_In
!  DO j=0,N_In
!    CALL LegendrePolynomialAndDerivative(j,xGauss(i),Vdm_Leg(i,j),dummy)
!  END DO !i
!END DO !j
!Vdm_Leg=TRANSPOSE(Vdm_Leg)
!DO j=0,N_In
!  Vdm_Leg(:,j)=Vdm_Leg(:,j)*wGauss(j)
!END DO
!!evaluate nodal basis (depends on NodeType, for Gauss: unity matrix)
!CALL InitializeVandermonde(N_In,N_In,wBary_Loc,xi_In,xGauss,sVdm_Leg)
!sVdm_Leg=MATMUL(Vdm_Leg,sVdm_Leg)

!compute the Vandermonde on xGP (Depends on NodeType)
DO i=0,N_In; DO j=0,N_In
  CALL LegendrePolynomialAndDerivative(j,xi_In(i),Vdm_Leg(i,j),dummy)
END DO; END DO !j
sVdm_Leg=INV(Vdm_Leg)
!check (Vdm_Leg)^(-1)*Vdm_Leg := I
dummy=ABS(SUM(ABS(MATMUL(sVdm_Leg,Vdm_Leg)))/REAL(N_In+1,wp)-1.0_wp)
IF(dummy.GT.10.0_wp*PP_RealTolerance) STOP &
                                         'problems in MODAL<->NODAL Vandermonde '
END SUBROUTINE buildLegendreVdm



!===================================================================================================================================
!> Build a 1D Vandermonde matrix using the Lagrange basis functions of degree
!> N_In, evaluated at the interpolation points xi_Out
!===================================================================================================================================
SUBROUTINE InitializeVandermonde(N_In,N_Out,wBary_In,xi_In,xi_Out,Vdm)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_In                       !! (IN)  input polynomial degree
INTEGER,INTENT(IN) :: N_Out                      !! (IN)  output polynomial degree
REAL(wp),INTENT(IN)    :: xi_In(0:N_In)              !! (IN)  input nodal positions [-1,1]
REAL(wp),INTENT(IN)    :: xi_Out(0:N_Out)            !! (IN)  outout nodal positions [-1,1]
REAL(wp),INTENT(IN)    :: wBary_In(0:N_In)           !! (IN)  input interpolation weights
REAL(wp),INTENT(OUT)   :: Vdm(0:N_Out,0:N_In)        !! (OUT) nodal Vandermonde from N_In to N_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iXi
!==================================================================================================================================
DO iXi=0,N_Out
  CALL LagrangeInterpolationPolys(xi_Out(iXi),N_In,xi_In,wBary_In,Vdm(iXi,:)) !l(0:N_In)
END DO
END SUBROUTINE InitializeVandermonde



!===================================================================================================================================
!> Evaluate the Legendre polynomial L_N and its derivative at position x[-1,1]
!> recursive algorithm using the N_in-1 N_in-2 Legendre polynomials
!> algorithm 22, Kopriva book
!===================================================================================================================================
SUBROUTINE LegendrePolynomialAndDerivative(N_in,x,L,Lder)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in   !! (IN)  polynomial degree, (N+1) CLpoints
REAL(wp),INTENT(IN)    :: x      !! (IN)  coordinate value in the interval [-1,1]
REAL(wp),INTENT(OUT)   :: L      !! (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
REAL(wp),INTENT(OUT)   :: Lder   !! (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iLegendre
REAL(wp)    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
REAL(wp)    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
IF(N_in .EQ. 0)THEN
  L=1.0_wp
  Lder=0.0_wp
ELSEIF(N_in .EQ. 1) THEN
  L=x
  Lder=1.0_wp
ELSE ! N_in > 1
  L_Nm2=1.0_wp
  L_Nm1=x
  Lder_Nm2=0.0_wp
  Lder_Nm1=1.0_wp
  DO iLegendre=2,N_in
    L=(REAL(2*iLegendre-1,wp)*x*L_Nm1 - REAL(iLegendre-1,wp)*L_Nm2)/REAL(iLegendre,wp)
    Lder=Lder_Nm2 + REAL(2*iLegendre-1,wp)*L_Nm1
    L_Nm2=L_Nm1
    L_Nm1=L
    Lder_Nm2=Lder_Nm1
    Lder_Nm1=Lder
  END DO !iLegendre=2,N_in
END IF ! N_in
!normalize
L=L*SQRT(REAL(N_in,wp)+0.5_wp)
Lder=Lder*SQRT(REAL(N_in,wp)+0.5_wp)
END SUBROUTINE LegendrePolynomialAndDerivative




!==================================================================================================================================
!> Compute Chebychev-Gauss nodes and integration weights (algorithm 27, Kopriva book)
!==================================================================================================================================
SUBROUTINE ChebyshevGaussNodesAndWeights(N_in,xGP,wGP)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !! polynomial degree, (N_in+1) CLpoints
REAL(wp),INTENT(OUT)          :: xGP(0:N_in)  !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !! Gauss point integration weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iGP
!==================================================================================================================================
DO iGP=0,N_in
  xGP(iGP)=-cos(REAL(2*iGP+1,wp)/REAL(2*N_in+2,wp)*PP_Pi)
END DO
IF(PRESENT(wGP))THEN
  DO iGP=0,N_in
    wGP(iGP)=PP_Pi/REAL(N_in+1,wp)
  END DO
END IF
END SUBROUTINE ChebyshevGaussNodesAndWeights



!==================================================================================================================================
!> Compute Chebychev-Gauss-Lobatto nodes and integration weights (algorithm 27, Kopriva book)
!==================================================================================================================================
SUBROUTINE ChebyGaussLobNodesAndWeights(N_in,xGP,wGP)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !! polynomial degree, (N_in+1) CLpoints
REAL(wp),INTENT(OUT)          :: xGP(0:N_in)  !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !! Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP
!==================================================================================================================================
DO iGP=0,N_in
  xGP(iGP)=-COS(REAL(iGP,wp)/REAL(N_in,wp)*PP_Pi)
END DO
IF(PRESENT(wGP))THEN
  DO iGP=0,N_in
    wGP(iGP)=PP_Pi/REAL(N_in,wp)
  END DO
  wGP(0)=wGP(0)*0.5_wp
  wGP(N_in)=wGP(N_in)*0.5_wp
END IF
END SUBROUTINE ChebyGaussLobNodesAndWeights


!==================================================================================================================================
!> Compute Clenshaw-Curtis nodes and integration weights
!==================================================================================================================================
SUBROUTINE ClenshawCurtisNodesAndWeights(N_in,xGP,wGP)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !! polynomial degree, (N_in+1) CLpoints
REAL(wp),INTENT(OUT)          :: xGP(0:N_in)  !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !! Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iGP,j
REAL(wp)                      :: b,theta
!==================================================================================================================================
IF(N_in.EQ.0)THEN
  xGP(0) = 0.0_wp
  wGP(0) = 2.0_wp
ELSE
  DO iGP=0,N_in
    xGP(iGP) = -COS(REAL(iGP,wp)*PP_Pi/REAL(N_in,wp))
  END DO
  xGP(0)   =-1.0_wp
  IF(MOD(N_in+1,2).EQ.1)THEN
    xGP(N_in/2)=0.0_wp
  END IF
  xGP(N_in)= 1.0_wp
  IF(PRESENT(wGP))THEN
    wGP=1.0_wp
    DO iGP=0,N_in
      theta = REAL(iGP,wp)*PP_Pi/REAL(N_in,wp)
      DO j=1,N_in/2
        b=MERGE(1.0_wp,2.0_wp,2*j.EQ.N_in)
        wGP(iGP) = wGP(iGP) - b * COS(REAL(2*j,wp)*theta) / REAL(4*j*j-1,wp)
      END DO
    END DO
    wGP(1:N_in-1)=2.0_wp*wGP(1:N_in-1)
    wGP=wGP/REAL(N_in,wp)
  END IF
END IF
END SUBROUTINE ClenshawCurtisNodesAndWeights

!==================================================================================================================================
!> @brief Compute Legendre-Gauss nodes and integration weights (algorithm 23, Kopriva book)
!>
!> Starting with Chebychev point positions, a Newton method is used to find the roots
!> of the Legendre Polynomial L_(N_in+1), which are the positions of Gausspoints
!> uses LegendrePolynomialAndDerivative subroutine
!==================================================================================================================================
SUBROUTINE LegendreGaussNodesAndWeights(N_in,xGP,wGP)
!MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in              !! polynomial degree, (N_in+1) Gausspoints
REAL(wp),INTENT(OUT)          :: xGP(0:N_in)       !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)       !! Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nIter = 10        ! max. number of newton iterations
REAL(wp)                  :: Tol   = 1.E-15_wp    ! tolerance for Newton iteration: TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
REAL(wp)                  :: L_Np1,Lder_Np1    ! L_{N_in+1},Lder_{N_in+1}
REAL(wp)                  :: dx                ! Newton step
REAL(wp)                  :: cheb_tmp          ! temporary variable for evaluation of chebychev node positions
!==================================================================================================================================
IF(N_in .EQ. 0) THEN
  xGP=0.0_wp
  IF(PRESENT(wGP))wGP=2.0_wp
  RETURN
ELSEIF(N_in.EQ.1)THEN
  xGP(0)=-sqrt(1.0_wp/3.0_wp)
  xGP(N_in)=-xGP(0)
  IF(PRESENT(wGP))wGP=1.0_wp
  RETURN
ELSE ! N_in>1
  cheb_tmp=2.0_wp*atan(1.0_wp)/REAL(N_in+1,wp) ! pi/(2N+2)
  DO iGP=0,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cheb_tmp*REAL(2*iGP+1,wp)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
      dx=-L_Np1/Lder_Np1
      xGP(iGP)=xGP(iGP)+dx
      IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      WRITE(*,*) 'maximum iteration steps >10 in Newton iteration for Legendre Gausspoint'
      xGP(iGP)=-cos(cheb_tmp*REAL(2*iGP+1)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        WRITE(*,*)iter,xGP(iGP)    !DEBUG
        CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
        dx=-L_Np1/Lder_Np1
        xGP(iGP)=xGP(iGP)+dx
        IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
      END DO !iter
      STOP &
                 'ERROR: Legendre Gauss nodes could not be computed up to desired precision. Code stopped!'
    END IF ! (iter.GT.nIter)
    CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      !wGP(iGP)=2./((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1) !if Legendre not normalized
      wGP(iGP)=REAL(2*N_in+3,wp)/((1.0_wp-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO !iGP
END IF ! N_in
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.0_wp
  CALL LegendrePolynomialAndDerivative(N_in+1,xGP(N_in/2),L_Np1,Lder_Np1)
  !IF(PRESENT(wGP))wGP(N_in/2)=2./(Lder_Np1*Lder_Np1) !if Legendre not normalized
  IF(PRESENT(wGP))wGP(N_in/2)=(REAL(2*N_in+3,wp))/(Lder_Np1*Lder_Np1)
END IF ! (mod(N_in,2) .EQ. 0)
END SUBROUTINE LegendreGaussNodesAndWeights



!==================================================================================================================================
!> Evaluate the polynomial q=L_{N_in+1}-L_{N_in-1} and its derivative at position x in [-1,1]
!> Recursive algorithm using the N_in-1 N_in-2 Legendre polynomials. (Algorithm 24, Kopriva book)
!==================================================================================================================================
SUBROUTINE qAndLEvaluation(N_in,x,q,qder,L)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in                            !! polynomial degree
REAL(wp),INTENT(IN)    :: x                               !! coordinate value in the interval [-1,1]
REAL(wp),INTENT(OUT)   :: L                               !! \f$ L_N(\xi) \f$
REAL(wp),INTENT(OUT)   :: q                               !! \f$ q_N(\xi) \f$
REAL(wp),INTENT(OUT)   :: qder                            !! \f$ \partial/\partial\xi \; L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iLegendre
REAL(wp)               :: L_Nm1,L_Nm2                     ! L_{N_in-2},L_{N_in-1}
REAL(wp)               :: Lder,Lder_Nm1,Lder_Nm2          ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
L_Nm2=1.0_wp
L_Nm1=x
Lder_Nm2=0.0_wp
Lder_Nm1=1.0_wp
DO iLegendre=2,N_in
  L=(REAL(2*iLegendre-1,wp)*x*L_Nm1 - REAL(iLegendre-1,wp)*L_Nm2)/REAL(iLegendre,wp)
  Lder=Lder_Nm2 + REAL(2*iLegendre-1,wp)*L_Nm1
  L_Nm2=L_Nm1
  L_Nm1=L
  Lder_Nm2=Lder_Nm1
  Lder_Nm1=Lder
END DO ! iLegendre
q=REAL(2*N_in+1,wp)/REAL(N_in+1,wp)*(x*L -L_Nm2) !L_{N_in+1}-L_{N_in-1} !L_Nm2 is L_Nm1, L_Nm1 was overwritten!
qder= REAL(2*N_in+1,wp)*L             !Lder_{N_in+1}-Lder_{N_in-1}
END SUBROUTINE qAndLEvaluation



!==================================================================================================================================
!> Starting with initial guess by Parter Relation, a Newton method is used to find the roots
!> of the Legendre Polynomial Lder_(N_in), which are the positions of Gauss-Lobatto points.
!> Uses qAndLEvaluation subroutine.
!> algorithm 25, Kopriva
!==================================================================================================================================
SUBROUTINE LegGaussLobNodesAndWeights(N_in,xGP,wGP)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in             !! polynomial degree (N_in+1) Gausspoints
REAL(wp),INTENT(OUT)          :: xGP(0:N_in)      !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)      !! Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nIter = 10       ! max. number of newton iterations
REAL(wp)                      :: Tol   = 1.E-15_wp   ! tolerance for Newton iteration : TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
REAL(wp)                      :: q,qder,L         ! \f$ q=L_{N_in+1}-L_{N_in-1} \f$ ,qder is derivative, \f$ L=L_{N_in} \f$
REAL(wp)                      :: dx               ! Newton step
REAL(wp)                      :: cont1,cont2      !temporary variable for evaluation of parter nodes positions
!==================================================================================================================================
xGP(0)=-1.0_wp
xGP(N_in)= 1.0_wp
IF(PRESENT(wGP))THEN
  wGP(0)= 2.0_wp/REAL(N_in*(N_in+1),wp)
  wGP(N_in)=wGP(0)
END IF
IF(N_in.GT.1)THEN
  cont1=PP_Pi/REAL(N_in,wp) ! pi/N_in
  cont2=3.0_wp/(REAL(8*N_in,wp)*PP_Pi) ! 3/(8*N_in*pi)
  DO iGP=1,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cont1*(REAL(iGP,wp)+0.25_wp)-cont2/(REAL(iGP,wp)+0.25_wp)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL qAndLEvaluation(N_in,xGP(iGP),q,qder,L)
      dx=-q/qder
      xGP(iGP)=xGP(iGP)+dx
      IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      WRITE(*,*) 'maximum iteration steps >10 in Newton iteration for LGL point:'
      xGP(iGP)=-cos(cont1*(REAL(iGP,wp)+0.25)-cont2/(REAL(iGP,wp)+0.25_wp)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        WRITE(*,*)'iter,x^i',iter,xGP(iGP)     !DEBUG
        CALL qAndLEvaluation(N_in,xGP(iGP),q,qder,L)
        dx=-q/qder
        xGP(iGP)=xGP(iGP)+dx
        IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
      END DO ! iter
      STOP &
                 'ERROR: Legendre Gauss Lobatto nodes could not be computed up to desired precision. Code stopped!'
    END IF ! (iter.GT.nIter)
    CALL qAndLEvaluation(N_in,xGP(iGP),q,qder,L)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      wGP(iGP)=wGP(0)/(L*L)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO ! iGP
END IF !(N_in.GT.1)
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.0_wp
  CALL qAndLEvaluation(N_in,xGP(N_in/2),q,qder,L)
  IF(PRESENT(wGP))wGP(N_in/2)=wGP(0)/(L*L)
END IF ! (mod(N_in,2) .EQ. 0)
END SUBROUTINE LegGaussLobNodesAndWeights



!==================================================================================================================================
!> Computes barycentric (interpolation) weights for interpolation polynomial given by set of nodes. (Algorithm 30, Kopriva book)
!==================================================================================================================================
SUBROUTINE BarycentricWeights(N_in,xGP,wBary)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in               !! polynomial degree
REAL(wp),INTENT(IN)    :: xGP(0:N_in)        !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT)   :: wBary(0:N_in)      !! barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,jGP
!==================================================================================================================================
wBary(:)=1.
DO iGP=1,N_in
  DO jGP=0,iGP-1
    wBary(jGP)=wBary(jGP)*(xGP(jGP)-xGP(iGP))
    wBary(iGP)=wBary(iGP)*(xGP(iGP)-xGP(jGP))
  END DO ! jGP
END DO ! iGP
wBary(:)=1./wBary(:)
END SUBROUTINE BarycentricWeights



!==================================================================================================================================
!> Computes polynomial differentiation matrix for interpolation polynomial given by set of nodes. (Algorithm 37, Kopriva book)
!==================================================================================================================================
SUBROUTINE PolynomialDerivativeMatrix(N_in,xGP,D)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in              !! polynomial degree
REAL(wp),INTENT(IN)    :: xGP(0:N_in)       !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(OUT)   :: D(0:N_in,0:N_in)  !! differentiation Matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,iLagrange
REAL(wp)               :: wBary(0:N_in)
!==================================================================================================================================
CALL BarycentricWeights(N_in,xGP,wBary)
D(:,:)=0.0_wp
DO iLagrange=0,N_in
  DO iGP=0,N_in
    IF(iLagrange.NE.iGP)THEN
      D(iGP,iLagrange)=wBary(iLagrange)/(wBary(iGP)*(xGP(iGP)-xGP(iLagrange)))
      D(iGP,iGP)=D(iGP,iGP)-D(iGP,iLagrange)
    END IF ! (iLagrange.NE.iGP)
  END DO ! iGP
END DO ! iLagrange
END SUBROUTINE PolynomialDerivativeMatrix



!==================================================================================================================================
!> Determines if two REAL(wp) numbers are equal up to a specified tolerance (=PP_RealTolerance, normaly set to machine precision)
!> Takes into account that x,y are located in-between [-1;1]
!> Based on Algorithm 139, Kopriva
!==================================================================================================================================
FUNCTION ALMOSTEQUAL(x,y)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(wp),INTENT(IN) :: x                !! (IN)  first scalar to be compared
REAL(wp),INTENT(IN) :: y                !! (IN)  second scalar to be compared
LOGICAL         :: AlmostEqual      !! (OUT) TRUE if |x-y| < 2*PP_RealTolerance
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
AlmostEqual=.FALSE.
IF((x.EQ.0.0_wp).OR.(y.EQ.0.0_wp)) THEN
  IF(ABS(x-y).LE.2.0_wp*PP_RealTolerance) AlmostEqual=.TRUE.
ELSE ! x, y not zero
  IF((ABS(x-y).LE.PP_RealTolerance*ABS(x)).AND.((ABS(x-y).LE.PP_RealTolerance*ABS(y)))) AlmostEqual=.TRUE.
END IF ! x,y zero
END FUNCTION ALMOSTEQUAL


!==================================================================================================================================
!> Determines if two REAL(wp) numbers are equal up to a given tolerance.
!> Routine requires: x,y > tolerance
!==================================================================================================================================
FUNCTION EQUALTOTOLERANCE(x,y,tolerance) 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(wp),INTENT(IN) :: x                !! (IN)  first scalar to be compared
REAL(wp),INTENT(IN) :: y                !! (IN)  second scalar to be compared
REAL(wp),INTENT(IN) :: tolerance        !! (IN)  Tolerance to be checked against
LOGICAL         :: EqualToTolerance !! (OUT) TRUE if x and y are closer than tolerance 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)            :: diff,maxInput
!===================================================================================================================================
EqualToTolerance = .FALSE.

maxInput = MAX(ABS(x),ABS(y))
diff = ABS(x-y)

! Test absolute error
IF (diff.LE.tolerance) THEN
  EqualToTolerance=.TRUE.
  RETURN
END IF

! Test relative error
IF(diff.LT.maxInput*tolerance) EqualToTolerance=.TRUE.

END FUNCTION EQUALTOTOLERANCE



!============================================================================================================================
!> Computes all Lagrange functions evaluated at position x in [-1,1]
!> For details see paper Barycentric Lagrange Interpolation by Berrut and Trefethen (SIAM 2004)
!> Uses function ALMOSTEQUAL
!> Algorithm 34, Kopriva book
!============================================================================================================================
SUBROUTINE LagrangeInterpolationPolys(x,N_in,xGP,wBary,L)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(wp), INTENT(IN)   :: x                !! Coordinate
INTEGER ,INTENT(IN)    :: N_in             !! polynomial degree
REAL(wp),INTENT(IN)    :: xGP(0:N_in)      !! Gauss point positions for the reference interval [-1,1]
REAL(wp),INTENT(IN)    :: wBary(0:N_in)    !! Barycentric weights
REAL(wp),INTENT(OUT)   :: L(0:N_in)        !! Lagrange basis functions evaluated at x
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP
LOGICAL            :: xEqualGP         ! is x equal to a Gauss Point
REAL(wp)               :: DummySum
!============================================================================================================================
xEqualGP=.FALSE.
DO iGP=0,N_in
  L(iGP)=0.0_wp
  IF(ALMOSTEQUAL(x,xGP(iGP))) THEN
    L(iGP)=1.0_wp
    xEqualGP=.TRUE.
  END IF
END DO

! if x is equal to a Gauss point, L=(0,....,1,....0)
IF(xEqualGP) RETURN
DummySum=0.0_wp
DO iGP=0, N_in
  L(iGP)=wBary(iGP)/(x-xGP(iGP))
  DummySum=DummySum+L(iGP)
END DO

DO iGP=0,N_in
  L(iGP)=L(iGP)/DummySum
END DO
END SUBROUTINE LagrangeInterpolationPolys


END MODULE MOD_Basis1D
