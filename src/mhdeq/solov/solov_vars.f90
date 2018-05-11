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
!!# Module **Soloviev variables**
!!
!!
!===================================================================================================================================
MODULE MODgvec_Solov_Vars
! MODULES
USE MODgvec_Globals, ONLY:wp
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER             :: setup         !! use given parameter setup: 10: circular, 20: iterlike
REAL(wp)            :: p_R0          !! major radius
REAL(wp)            :: p_eps         !! minor/major radius ratio
REAL(wp)            :: p_kappa       !! elliticity
REAL(wp)            :: p_delta       !! triangularity
REAL(wp)            :: asin_delta    !! ASIN(p_delta)
REAL(wp)            :: p_A           !! pA <1 (else deltap>0)
REAL(wp)            :: p_B0          !! toroidal magnetic field strength at magn. axis
REAL(wp)            :: p_qaxis       !! q-factor on axis
REAL(wp)            :: p_paxis       !! pressure on axis
REAL(wp)            :: xaxis(2)      !! xy-position o of magnetic axis
REAL(wp)            :: psi_scale     !! psi scaling between Soloviev and physical psiReal=psi_scale*psi
REAL(wp)            :: psi_axis      !! psi value at the magnetic axis
REAL(wp)            :: psi_edge      !! psi value at the edge
REAL(wp)            :: presEdge      !! pressure on edge, soloviev pressure profile linear in psinorm
REAL(wp)            :: F_axis        !! F^2 on axis, f-profile, dF^2/dpsi=const for soloviev 
REAL(wp)            :: deltaF2       !! F=SQRT(F_axis^2+deltaF2*psinorm)
REAL(wp)            :: psiCoefs(0:7) !! coefficients for representing flux variable psi

!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
END MODULE MODgvec_Solov_Vars

