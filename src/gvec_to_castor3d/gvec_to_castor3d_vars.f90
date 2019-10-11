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

!===================================================================================================================================
!>
!!# Module ** gvec_to_castor3d Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_gvec_to_castor3d_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
!INPUT VARIABLES
CHARACTER(LEN=255) :: fileName      !< name of GVEC file
INTEGER            :: Ns_out        !< number of equidistant points in radial s-direction (includes axis and edge!)
INTEGER            :: SFLcoord      !< which angular coordinates to choose: =0: GVEC coord. (no SFL), =1: PEST SFL, =2: BOOZER SFL
INTEGER            :: factorFourier !< factor theta,zeta resolution Ntheta=Factor*m_max, Nzeta=MAX(1,Factor*n_max)
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER               :: nfp_out            !< number of field periods
INTEGER               :: asym_out           !< =0: symmetric configuration (R~cos,Z~sin,lambda~sin), =1 asymmetric 
INTEGER               :: mn_max_out(2)      !< maximum number of modes in m,n
INTEGER               :: Nthet_out          !< total number of points in theta direction theta[0,2pi (
INTEGER               :: Nzeta_out          !< total number of points in zeta direction zeta[0,-2pi/NFP( opposite sign compared to GVEC!!
REAL(wp)              :: PhiEdge            !< total toroidal flux at the last flux surface (0 at axis), *2Pi and opposite compared GVEC!!
REAL(wp)              :: ChiEdge            !< total poloidal flux at the last flux surface (0 at axis), *2pi compared to GVEC!!!
REAL(wp),ALLOCATABLE  :: s_pos(:)           !< positions in s for evaluation s=sqrt(phi/phiEdge), size (Ns_out)
REAL(wp),ALLOCATABLE  :: thet_pos(:)        !< positions in theta for evaluation, size (Nthet_out)
REAL(wp),ALLOCATABLE  :: zeta_pos(:)        !< positions in zeta for evaluation , size (Nzeta_out)

!1D data
INTEGER,PARAMETER     :: nVar1D = 12         !< number of variables in 1d profiles
INTEGER,PARAMETER     :: SPOS__    = 1
INTEGER,PARAMETER     :: PHI__     = 2
INTEGER,PARAMETER     :: DPHIDS__  = 3
INTEGER,PARAMETER     :: CHI__     = 4
INTEGER,PARAMETER     :: DCHIDS__  = 5
INTEGER,PARAMETER     :: IOTA__    = 6
INTEGER,PARAMETER     :: PRESSURE__= 7
INTEGER,PARAMETER     :: ITOR__    = 8
INTEGER,PARAMETER     :: IPOL__    = 9
INTEGER,PARAMETER     :: FAVG__    =10
INTEGER,PARAMETER     :: FMIN__    =11
INTEGER,PARAMETER     :: FMAX__    =12
CHARACTER(LEN=50),DIMENSION(nVar1D),PARAMETER :: StrVarNames1D(nVar1D)=(/ CHARACTER(LEN=50) :: &
                           's'            & ! 1 : position s =sqrt(phi/phiEdge) [0,1]
                          ,'Phi'          & ! 2 : toroidal flux 
                          ,'dPhi_ds'      & ! 3 : derivative of toroidal flux to s coordinate
                          ,'Chi'          & ! 4 : poloidal flux
                          ,'dChi_ds'      & ! 5 : derivative of poloidal flux to s coordinate
                          ,'iota'         & ! 6 : iota profile
                          ,'Pressure'     & ! 7 : pressure 
                          ,'Itor'         & ! 8 : Toroidal current 
                          ,'Ipol'         & ! 9 : Poloidal current 
                          ,'Favg'         & !10 : Only tokamaks(n=0!), toroidal magnetic field strength is F/R (averaged over theta)
                          ,'Fmin'         & !11 : F(s) is averaged over theta, Fmin(s) = min(F(s,theta))
                          ,'Fmax'         & !12 : F(s) is averaged over theta, Fmax(s) = max(F(s,theta))
                                    /)
REAL(wp),ALLOCATABLE  :: data_1D(:,:)        !< 1D profiles size (nVar1D,Ns_out)

!3D scalar data 
INTEGER,PARAMETER     :: nVarScalar3D=5           !< number of variabels in 3D data
INTEGER,PARAMETER     :: X1__     = 1
INTEGER,PARAMETER     :: X2__     = 2
INTEGER,PARAMETER     :: GZETA__  = 3
INTEGER,PARAMETER     :: BSUPT__  = 4
INTEGER,PARAMETER     :: BSUPZ__  = 5
CHARACTER(LEN=50),DIMENSION(nVarScalar3D),PARAMETER :: StrVarNamesScalar3D(nVarScalar3D)=(/ CHARACTER(LEN=50) :: &
                           'X1(R)'       & ! 1 : for Torus map (hmap=1), R=X1
                          ,'X2(Z)'       & ! 2 : for Torus map (hmap=1), Z=X2
                          ,'Gzeta'       & ! 3 : map to geometric toroidal angle, phi = zeta+Gzeta
                          ,'B^theta'     & ! 4 : theta component of magnetic field B^theta = B.grad(theta) 
                          ,'B^zeta'      & ! 5 : zeta component of magnetic field B^theta = B.grad(zeta) 
                                    /)
REAL(wp),ALLOCATABLE  :: data_scalar3D(:,:,:,:)    !< Size (Nthet_out,Nzeta_out,Ns_out,nVar3D)

!3D vector data 
INTEGER,PARAMETER     :: nVarVector3D=4           !< number of variabels in 3D data
INTEGER,PARAMETER     :: BFIELD__     = 1
INTEGER,PARAMETER     :: ECOV_S__     = 2
INTEGER,PARAMETER     :: ECOV_THETA__ = 3
INTEGER,PARAMETER     :: ECOV_ZETA__  = 4
CHARACTER(LEN=50),DIMENSION(nVarVector3D),PARAMETER :: StrVarNamesVector3D(nVarVector3D)=(/ CHARACTER(LEN=50) :: &
                           'Bfield'           & ! 1 : magnetic field vector    , (x,y,z) cartesian components
                          ,'ecov_s'           & ! 2 : covariant vector in s    , (x,y,z) cartesian components
                          ,'ecov_theta'       & ! 3 : covariant vector in theta, (x,y,z) cartesian components
                          ,'ecov_zeta'        & ! 4 : covariant vector in zeta , (x,y,z) cartesian components
                                    /)
REAL(wp),ALLOCATABLE  :: data_vector3D(:,:,:,:,:)    !< Size (3,Nthet_out,Nzeta_out,Ns_out,nVarVector3D)
                                         
!===================================================================================================================================



END MODULE MODgvec_gvec_to_castor3d_Vars

