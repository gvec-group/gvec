!===================================================================================================================================
! Copyright (c) 2017-2018 Florian Hindenlang
!
! This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================


!===================================================================================================================================
!> Here globally used variables /functions are defined 
!===================================================================================================================================
MODULE MOD_Globals

IMPLICIT NONE

PUBLIC 

!-----------------------------------------------------------------------------------------------------------------------------------
! Select here the working precision wp
!INTEGER, PARAMETER :: wp = selected_real_kind(6,35)   !< single precision
INTEGER, PARAMETER :: wp = selected_real_kind(15,307)  !< double precision
!INTEGER, PARAMETER :: wp = selected_real_kind(33,307) !< quadruple precision
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=20) :: fmt_sep ='(100("="))'              !< formatting of separator line: WRITE(*,fmt_sep)
REAL(wp),PARAMETER  :: Pi   =ACOS(-1.0_wp) !< pi parameter
REAL(wp),PARAMETER  :: TwoPi=2.0_wp*Pi
!-----------------------------------------------------------------------------------------------------------------------------------


INTERFACE CROSS
   MODULE PROCEDURE CROSS
END INTERFACE

INTERFACE NORMALIZE
   MODULE PROCEDURE NORMALIZE
END INTERFACE

INTERFACE GetDet3
   MODULE PROCEDURE GetDet3
END INTERFACE

INTERFACE GetInv3
   MODULE PROCEDURE GetInv3
END INTERFACE

CONTAINS

!===================================================================================================================================
!> normalizes a nDim vector with respect to the eucledian norm
!===================================================================================================================================
PURE FUNCTION NORMALIZE(v1,nVal)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nVal     !< vector size
REAL(wp),INTENT(IN)     :: v1(nVal) !< vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)            :: normalize(nVal) ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
normalize=v1/SQRT(SUM(v1*v1))
END FUNCTION NORMALIZE


!===================================================================================================================================
!> computes the cross product of to 3 dimensional vectors: cross=v1 x v2
!===================================================================================================================================
PURE FUNCTION CROSS(v1,v2)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN) :: v1(3)    ! ? 
REAL(wp),INTENT(IN) :: v2(3)    ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)            :: cross(3) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
cross=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS


!===================================================================================================================================
!> compute determinant of 3x3 matrix
!===================================================================================================================================
FUNCTION getDet3(Mat)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)  :: Mat(3,3) !< input matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)             :: getDet3 !< determinant of the input matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
getDet3=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
         + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
         + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION getDet3


!===================================================================================================================================
!> compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
!===================================================================================================================================
FUNCTION getInv3(Mat,sDet_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(wp),INTENT(IN)             :: Mat(3,3) ! input matrix
REAL(wp),INTENT(IN),OPTIONAL    :: sDet_in  ! determinant of input matrix (otherwise computed here)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(wp)             :: getInv3(3,3) !< inverse matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(wp)             :: sDet 
!===================================================================================================================================
IF(PRESENT(sDet_in))THEN
  sDet=sDet_in
ELSE
  sDet=1.0_wp/getDet3(Mat)
END IF
getInv3(1,1) = ( Mat(2,2) * Mat(3,3) - Mat(2,3) * Mat(3,2) ) * sDet
getInv3(1,2) = ( Mat(1,3) * Mat(3,2) - Mat(1,2) * Mat(3,3) ) * sDet
getInv3(1,3) = ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * sDet
getInv3(2,1) = ( Mat(2,3) * Mat(3,1) - Mat(2,1) * Mat(3,3) ) * sDet
getInv3(2,2) = ( Mat(1,1) * Mat(3,3) - Mat(1,3) * Mat(3,1) ) * sDet
getInv3(2,3) = ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * sDet
getInv3(3,1) = ( Mat(2,1) * Mat(3,2) - Mat(2,2) * Mat(3,1) ) * sDet
getInv3(3,2) = ( Mat(1,2) * Mat(3,1) - Mat(1,1) * Mat(3,2) ) * sDet
getInv3(3,3) = ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * sDet
END FUNCTION getInv3


END MODULE MOD_Globals
