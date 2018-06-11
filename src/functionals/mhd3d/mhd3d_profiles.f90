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
!!# Module **MHD3D profiles**
!!
!! CONTAINS  all functions to evaluate 1D profiles
!!
!===================================================================================================================================
MODULE MODgvec_MHD3D_profiles
! MODULES
USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep
IMPLICIT NONE
PUBLIC

!===================================================================================================================================

CONTAINS



!===================================================================================================================================
!> evaluate rotational transform, iota(phi_norm(s)) =chi'/phi'
!! NOTE that since VMEC has a definition of a positive iota with respect to R,phi,Z coordinate system, but GVEC is in (R,Z,phi)
!! the sign of iota in GVEC is opposite to the sign in VMEC
!!
!===================================================================================================================================
FUNCTION Eval_iota(spos)
! MODULES
USE MODgvec_Globals    ,ONLY: EVAL1DPOLY
USE MODgvec_MHD3D_Vars ,ONLY: which_init,n_iota_coefs,iota_coefs
USE MODgvec_VMEC       ,ONLY: VMEC_EvalSpl
USE MODgvec_VMEC_vars  ,ONLY: iota_spl
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_iota
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm=Eval_PhiNorm(spos)
  SELECT CASE(which_init)
  CASE(0)
    eval_iota=Eval1DPoly(n_iota_coefs,iota_coefs,phi_norm)
  CASE(1)
    eval_iota= VMEC_EvalSpl(0,SQRT(phi_norm),iota_Spl) !variable rho in vmec evaluations is sqrt(phi/phi_edge)
  END SELECT
END FUNCTION Eval_iota

!===================================================================================================================================
!> evaluate pressure profile p(phi_norm(s))
!!
!===================================================================================================================================
FUNCTION Eval_pres(spos)
! MODULES
USE MODgvec_Globals    ,ONLY: EVAL1DPOLY
USE MODgvec_MHD3D_Vars ,ONLY: which_init,n_pres_coefs,pres_coefs
USE MODgvec_VMEC       ,ONLY: VMEC_EvalSpl
USE MODgvec_VMEC_vars  ,ONLY: pres_spl
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_pres
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm=Eval_PhiNorm(spos)
  SELECT CASE(which_init)
  CASE(0)
      eval_pres=Eval1DPoly(n_pres_coefs,pres_coefs,phi_norm) 
  CASE(1)
    eval_pres=VMEC_EvalSpl(0,SQRT(phi_norm),pres_Spl) !variable rho in vmec evaluations is sqrt(phi/phi_edge)
  END SELECT
END FUNCTION Eval_pres

!===================================================================================================================================
!> evaluate mass profile mass(phi_norm(s))
!!
!===================================================================================================================================
FUNCTION Eval_mass(spos)
! MODULES
USE MODgvec_MHD3D_vars, ONLY:gamm
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_mass
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(gamm).LT.1.0E-12)THEN
  eval_mass=eval_pres(spos)
ELSE
  !TODO would need to multiply by ( V'(s) ) ^gamma
  CALL abort(__STAMP__, &
      ' gamma>0: evaluation of mass profile not yet implemented!')
END IF
END FUNCTION Eval_mass

!===================================================================================================================================
!> evaluate poloidal flux chi(phi_norm), currently only from interpolated VMEC data
!> should be computed as an integral over chiPrime instead!! 
!!
!===================================================================================================================================
FUNCTION Eval_chi(spos)
! MODULES
USE MODgvec_MHD3D_Vars ,ONLY: which_init
USE MODgvec_VMEC       ,ONLY: VMEC_EvalSpl
USE MODgvec_VMEC_vars  ,ONLY: chi_spl
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_chi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: phi_norm
!===================================================================================================================================
  phi_norm=Eval_PhiNorm(spos)
  SELECT CASE(which_init)
  CASE(0)
    WRITE(*,*) 'WARNING: eval_chi not implemented yet for which_init=0'  
  CASE(1)
    eval_chi=VMEC_EvalSpl(0,SQRT(phi_norm),chi_Spl) !variable rho in vmec evaluations is sqrt(phi/phi_edge)
  END SELECT
END FUNCTION Eval_chi

!===================================================================================================================================
!> evaluate s derivative of poloidal flux via iota and phiPrime: chiPrime=-iota*PhiPrime
!!
!===================================================================================================================================
FUNCTION Eval_chiPrime(spos)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_chiPrime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Eval_ChiPrime = Eval_PhiPrime(spos)*Eval_iota(spos)
END FUNCTION Eval_chiPrime

!===================================================================================================================================
!> evaluate toroidal flux Phi 
!!
!===================================================================================================================================
FUNCTION Eval_Phi(spos)
! MODULES
USE MODgvec_MHD3D_Vars,ONLY:Phi_edge
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_Phi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Eval_Phi=Phi_edge*Eval_PhiNorm(spos)

END FUNCTION Eval_Phi

!===================================================================================================================================
!> evaluate s-derivative of toroidal flux Phi
!!
!===================================================================================================================================
FUNCTION Eval_PhiPrime(spos)
! MODULES
USE MODgvec_MHD3D_Vars,ONLY:Phi_edge
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_PhiPrime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Eval_PhiPrime=Phi_edge*Eval_PhiNormPrime(spos)

END FUNCTION Eval_PhiPrime

!===================================================================================================================================
!> evaluate normalized toroidal flux Phi/Phi_edge=s^2, this variable is used for the input profiles of iota and pressure!!!
!!
!===================================================================================================================================
FUNCTION Eval_PhiNorm(spos)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_PhiNorm
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Eval_PhiNorm=(spos**2)

END FUNCTION Eval_PhiNorm

!===================================================================================================================================
!> evaluate s derivative of normalized toroidal flux Phi/Phi_edge =s^2
!!
!===================================================================================================================================
FUNCTION Eval_PhiNormPrime(spos)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp), INTENT(IN   ) :: spos !! s position to evaluate s=[0,1] 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)               :: Eval_PhiNormPrime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Eval_PhiNormPrime=2.0_wp*spos

END FUNCTION Eval_PhiNormPrime

END MODULE MODgvec_MHD3D_profiles
