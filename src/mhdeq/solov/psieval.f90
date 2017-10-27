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
!!# Module **Psi Evaluation**
!!
!! Evaluate psi from a grad-shavranov solution
!!
!===================================================================================================================================
MODULE MOD_PsiEval
! MODULES
USE MOD_Globals, ONLY:wp
IMPLICIT NONE
PUBLIC
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Evaluate Psi functions
!!
!===================================================================================================================================
FUNCTION EvalPsi(x,y)
! MODULES
USE MOD_Solov_Vars,ONLY:psiCoefs
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: EvalPsi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
EvalPsi=SUM(EvalPsiVec(x,y)*PsiCoefs)

END FUNCTION EvalPsi


!===================================================================================================================================
!> Evaluate derivative of Psi functions
!!
!===================================================================================================================================
FUNCTION EvaldPsi(dir,nn,x,y)
! MODULES
USE MOD_Solov_Vars,ONLY:psiCoefs
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,  INTENT(IN) :: dir !! 1:x 2:y
INTEGER,  INTENT(IN) :: nn  !! nth derivative
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: EvaldPsi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(dir)
CASE(1)
  SELECT CASE(nn)
  CASE(1)
    EvaldPsi=SUM(EvaldPsidxVec(x,y)*PsiCoefs)
  CASE(2)
    EvaldPsi=SUM(Evald2PsidxVec(x,y)*PsiCoefs)
  END SELECT
CASE(2)
  SELECT CASE(nn)
  CASE(1)
    EvaldPsi=SUM(EvaldPsidyVec(x,y)*PsiCoefs)
  CASE(2)
    EvaldPsi=SUM(Evald2PsidyVec(x,y)*PsiCoefs)
  END SELECT
CASE(3) !mixed d^2/(dxdy)
  EvaldPsi=SUM(Evald2PsidxdyVec(x,y)*PsiCoefs)
END SELECT

END FUNCTION EvaldPsi


!===================================================================================================================================
!> Evaluate Psi functions
!!
!===================================================================================================================================
FUNCTION EvalPsiVec(x,y)
! MODULES
USE MOD_Solov_Vars,ONLY:p_A
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: EvalPsiVec(0:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)               :: x2,y2,lgx
!===================================================================================================================================
x2=x*x
y2=y*y
lgx=LOG(x)
EvalPsiVec(0) = 0.125_wp*x2*(x2 + p_A * (4.0_wp * lgx - x2 ))
EvalPsiVec(1) = 1.0_wp
EvalPsiVec(2) = x2
EvalPsiVec(3) = y2 - x2 * lgx
EvalPsiVec(4) = x2*(x2 - 4.0_wp * y2)
EvalPsiVec(5) = y2*( 2.0_wp * y2 - 9.0_wp * x2)   + 3.0_wp*lgx*x2*( x2 - 4.0_wp * y2 )
EvalPsiVec(6) = x2*(x2*(x2 - 12.0_wp * y2)) + 8.0_wp*x2*y2*y2
EvalPsiVec(7) = y2*(y2*(8.0_wp*y2 - (140.0_wp+120.0_wp*lgx)* x2 )) + x2*( (75 +180.0_wp*lgx)*x2*y2  - 15.0_wp*x2*x2*lgx)
         
END FUNCTION EvalPsiVec


!===================================================================================================================================
!> Evaluate d/dx (Psi) functions
!!
!===================================================================================================================================
FUNCTION EvaldPsidxVec(x,y)
! MODULES
USE MOD_Solov_Vars,ONLY:p_A
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: EvaldPsidxVec(0:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)               :: x2,y2,lgx
!===================================================================================================================================
x2=x*x
y2=y*y
lgx=LOG(x)
EvaldPsidxVec(0) = x*(p_A*(1.0_wp*lgx + 0.5_wp) + 0.5_wp*x2*(1.0_wp-p_A))
EvaldPsidxVec(1) = 0.0_wp
EvaldPsidxVec(2) = 2.0_wp*x
EvaldPsidxVec(3) = -x*(2.0_wp*lgx + 1.0_wp)
EvaldPsidxVec(4) = x*4.*(x2 - 2.0_wp*y2)
EvaldPsidxVec(5) = x*(x2*(12.0_wp*lgx + 3.0_wp) - y2*(24.*lgx + 30.0_wp))
EvaldPsidxVec(6) = x*(x2*6.*(x2 - 8.0_wp*y2) + 16.0_wp*y2*y2)
EvaldPsidxVec(7) = x*(x2*(y2*(720.0_wp*lgx + 480.0_wp)- x2*(90.0_wp*lgx + 15.0_wp) ) - y2*y2*(240.0_wp*lgx + 400.0_wp)) 
         
END FUNCTION EvaldPsidxVec


!===================================================================================================================================
!> Evaluate d^2/dx^2 (Psi) functions
!!
!===================================================================================================================================
FUNCTION Evald2PsidxVec(x,y)
! MODULES
USE MOD_Solov_Vars,ONLY:p_A
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: Evald2PsidxVec(0:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)               :: x2,y2,lgx
!===================================================================================================================================
x2=x*x
y2=y*y
lgx=LOG(x)
Evald2PsidxVec(0) = p_A*(lgx + 1.5_wp) + 1.5_wp*x2*(1.0_wp-p_A)
Evald2PsidxVec(1) = 0.0_wp
Evald2PsidxVec(2) = 2.0_wp
Evald2PsidxVec(3) = -(2.0_wp*lgx + 3.0_wp)
Evald2PsidxVec(4) = 12.0_wp*x2 - 8.0_wp*y2
Evald2PsidxVec(5) = x2*(36.0_wp*lgx + 21.0_wp) - y2*(24.0_wp*lgx + 54.0_wp)
Evald2PsidxVec(6) = x2*(30.0_wp*x2 - 144.0_wp*y2) + 16.0_wp*y2*y2
Evald2PsidxVec(7) = x2*(y2*(2160.0_wp*lgx + 2160.0_wp) - x2*(450.0_wp*lgx + 165.0_wp) ) - y2*y2*(240.0_wp*lgx + 640.0_wp) 
         
END FUNCTION Evald2PsidxVec


!===================================================================================================================================
!> Evaluate d^2/(dxdy) (Psi) functions
!!
!===================================================================================================================================
FUNCTION Evald2PsidxdyVec(x,y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: Evald2PsidxdyVec(0:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)               :: x2,y2,lgx
!===================================================================================================================================
x2=x*x
y2=y*y
lgx=LOG(x)
Evald2PsidxdyVec(0) = 0.0_wp
Evald2PsidxdyVec(1) = 0.0_wp
Evald2PsidxdyVec(2) = 0.0_wp
Evald2PsidxdyVec(3) = 0.0_wp
Evald2PsidxdyVec(4) = -16.0_wp*x*y
Evald2PsidxdyVec(5) = x*y*(-48.0_wp*lgx - 60.0_wp)
Evald2PsidxdyVec(6) = x*y*(-96.0_wp*x2 + 64.0_wp*y2)
Evald2PsidxdyVec(7) = x*y*(x2*(1440.0_wp*lgx + 960.0_wp) + y2*(-960.0_wp*lgx - 1600.0_wp))
END FUNCTION Evald2PsidxdyVec


!===================================================================================================================================
!> Evaluate d/dy (Psi) functions
!!
!===================================================================================================================================
FUNCTION EvaldPsidyVec(x,y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: EvaldPsidyVec(0:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)               :: x2,y2,lgx
!===================================================================================================================================
x2=x*x
y2=y*y
lgx=LOG(x)
EvaldPsidyVec(0) = 0.0_wp 
EvaldPsidyVec(1) = 0.0_wp
EvaldPsidyVec(2) = 0.0_wp
EvaldPsidyVec(3) = 2.0_wp*y
EvaldPsidyVec(4) = -8.0_wp*x2*y
EvaldPsidyVec(5) =  8.0_wp*y*y2 -x2*y*(24.0_wp*lgx + 18.0_wp)
EvaldPsidyVec(6) = x2*(-24.0_wp*x2*y + 32.0_wp*y*y2)
EvaldPsidyVec(7) = x2*(x2*y*(360.0_wp*lgx + 150.0_wp)) + y2*y*(x2*(-480.0_wp*lgx - 560.0_wp) + 48.0_wp*y2)
         
END FUNCTION EvaldPsidyVec


!===================================================================================================================================
!> Evaluate d/dy (Psi) functions
!!
!===================================================================================================================================
FUNCTION Evald2PsidyVec(x,y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp), INTENT(IN)   :: x   !! x coordinate   
REAL(wp), INTENT(IN)   :: y   !! y coordinate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)               :: Evald2PsidyVec(0:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)               :: x2,y2,lgx
!===================================================================================================================================
x2=x*x
y2=y*y
lgx=LOG(x)
Evald2PsidyVec(0) = 0.0_wp 
Evald2PsidyVec(1) = 0.0_wp
Evald2PsidyVec(2) = 0.0_wp
Evald2PsidyVec(3) = 2.0_wp
Evald2PsidyVec(4) = -8.0_wp*x2
Evald2PsidyVec(5) =  24.0_wp*y2 -x2*(24.0_wp*lgx + 18.0_wp)
Evald2PsidyVec(6) = x2*(-24.0_wp*x2 + 96.0_wp*y2)
Evald2PsidyVec(7) = x2*(x2*(360.0_wp*lgx + 150.0_wp)) + y2*(x2*(-1440.0_wp*lgx - 1680.0_wp) + 240.0_wp*y2)

END FUNCTION Evald2PsidyVec

END MODULE MOD_PsiEval
