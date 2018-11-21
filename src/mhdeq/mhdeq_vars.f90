!===================================================================================================================================
! Copyright (C) 2017 - 2018  Florian Hindenlang <hindenlang@gmail.com>
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

!===================================================================================================================================
!>
!!# Module **MHDEQ Variables**
!!
!! Variables for the MHD equilibrium mapping / data
!!
!===================================================================================================================================
MODULE MODgvec_MHDEQ_Vars
! MODULES
USE MODgvec_Globals,ONLY:wp
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER              :: WhichInitEquilibrium      !! 0: no input read (default), 1: use VMEC data, 2: use Solov'ev equilibrium
INTEGER              :: InputCoordSys             !! 0: x_in(1:3)=(x,y,z) with r^2=x^2+y^2,r=[0,1], z=[0,1]
                                                  !! 1: x_in(1:3)= (rho,zeta,theta) rho[0,1],zeta[0,1],theta[0,1] 
LOGICAL              :: useMHDEQ                  !! =(whichEquilibrium>0)
INTEGER              :: nVarMHDEQ=10               
REAL(wp),ALLOCATABLE :: MHDEQoutdataGL(:,:,:,:,:) !! MHD equilibrium data to be written to hdf5 file, on Gauss-Lobatto nodes
REAL(wp),ALLOCATABLE :: MHDEQdataEq(:,:,:,:,:)    !! VMEC data on equidistant nodes (forvisualization) 
INTEGER              :: nRhoCoefs                 !! number of density coefficients 
INTEGER              :: RhoFluxVar                !! Dependent variable for evaluation of Density polynomial rho(x)
                                                  !! =0: x=psinorm Normalized TOROIDAL flux
                                                  !! =1: x=chinorm Normalized POLOIDAL flux 
                                                  !! =2: x=sqrt(psinorm), =3: x=sqrt(chinorm)
REAL(wp),ALLOCATABLE :: RhoCoefs(:)               !! density coefficients of the polynomial coefficients:
                                                  !! rho_1+rho_2*x + rho_3*x^2 ...
                                                
CHARACTER(LEN=255),DIMENSION(10),PARAMETER :: MHDEQvarNames(10)=(/ CHARACTER(LEN=255) :: &
                      'MHDEQ-Density'     & ! 1 
                     ,'MHDEQ-Pressure'    & ! 2
                     ,'MHDEQ-BX'          & ! 3
                     ,'MHDEQ-BY'          & ! 4
                     ,'MHDEQ-BZ'          & ! 5
                     ,'MHDEQ-polflux'     & ! 6
                     ,'MHDEQ-torflux'     & ! 7
                     ,'MHDEQ-AX'          & ! 8
                     ,'MHDEQ-AY'          & ! 9
                     ,'MHDEQ-AZ'          & !10    
                                         /)

!===================================================================================================================================
END MODULE MODgvec_MHDEQ_Vars

