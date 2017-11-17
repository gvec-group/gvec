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
!!# Module ** base **
!!
!! 2D Fourier base in the two angular directions: (poloidal,toroidal) ~ (m,n) ~ (theta,zeta) [0,2pi]x[0,2pi/nfp]
!!
!! explicit real fourier basis: sin(x_mn) or cos(x_mn) with x_mn=(m*theta - n*nfp*zeta)  ,
!! with mode numbers m and n 
!!
!===================================================================================================================================
MODULE MOD_base
! MODULES
USE MOD_Globals ,ONLY: wp,Unit_stdOut,abort
USE MOD_sBase   ,ONLY: t_sbase,sbase_new
USE MOD_fBase   ,ONLY: t_fbase,fbase_new
USE MOD_sGrid   ,ONLY: t_sgrid
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

!TYPE, ABSTRACT :: c_base
!  !---------------------------------------------------------------------------------------------------------------------------------
!  CLASS(c_sbase),ALLOCATABLE  :: s  !! container for radial basis
!  CLASS(c_fbase),ALLOCATABLE  :: f  !! container for angular basis
!  CONTAINS
!    PROCEDURE(i_sub_base_init       ),DEFERRED :: init
!    PROCEDURE(i_sub_base_free       ),DEFERRED :: free
!    PROCEDURE(i_sub_base_copy       ),DEFERRED :: copy
!    PROCEDURE(i_fun_base_evalDOF    ),DEFERRED :: evalDOF
!
!END TYPE c_base
!
!ABSTRACT INTERFACE
!  SUBROUTINE i_sub_base_init(sf)
!    IMPORT c_base
!    CLASS(c_base) , INTENT(INOUT) :: sf
!  END SUBROUTINE i_sub_base_init
!
!  SUBROUTINE i_sub_base_free( sf ) 
!    IMPORT c_base
!    CLASS(c_base), INTENT(INOUT) :: sf
!  END SUBROUTINE i_sub_base_free
!
!  SUBROUTINE i_sub_base_copy( sf, tocopy ) 
!    IMPORT c_base
!    CLASS(c_base), INTENT(IN   ) :: tocopy
!    CLASS(c_base), INTENT(INOUT) :: sf
!  END SUBROUTINE i_sub_base_copy
!
!  FUNCTION i_fun_base_evalDOF( sf,DOFs ) RESULT(y_IP_GP)
!    IMPORT wp,c_base
!  CLASS(c_base), INTENT(IN   ) :: sf
!  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%f%modes,1:sf%s%nBase)
!  REAL(wp)                      :: y_IP_GP(sf%f%mn_IP,sf%s%nGP)
!  END FUNCTION i_fun_base_evalDOF
!
!END INTERFACE
 


!TYPE,EXTENDS(c_base) :: t_base
TYPE                 :: t_base
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters sbase
  INTEGER              :: s_deg                      !! input parameter: degree of Spline/polynomial 
  INTEGER              :: s_degGP                    !! number of Gauss-points (degGP+1) per element >= deg
  INTEGER              :: s_continuity               !! input parameter: full spline (=deg-1) or discontinuous (=-1)
  CLASS(t_sgrid),POINTER :: s_grid                   !! pointer to grid 
  !input parameters fbase
  INTEGER              :: f_mn_max(2)   !! input parameter: maximum number of fourier modes: m_max=mn_max(1),n_max=mn_max(2)
  INTEGER              :: f_mn_nyq(2)   !! number of equidistant integration points (trapezoidal rule) in m and n
  INTEGER              :: f_nfp         !! number of field periods (toroidal repetition after 2pi/nfp)
  INTEGER              :: f_sin_cos     !! can be either only sine: _SIN_  or only cosine _COS_ or full: _SINCOS_
  !input parameters
  LOGICAL              :: f_exclude_mn_zero  !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)

  CLASS(t_sbase),ALLOCATABLE  :: s  !! container for radial basis
  CLASS(t_fbase),ALLOCATABLE  :: f  !! container for angular basis
  
  CONTAINS

  PROCEDURE :: free       => base_free
  PROCEDURE :: copy       => base_copy
  PROCEDURE :: evalDOF    => base_evalDOF

END TYPE t_base

LOGICAL  :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> allocate and initialize the type base 
!!
!===================================================================================================================================
SUBROUTINE Base_new( sf,deg_in,continuity_in,grid_in,degGP_in, &
                     mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   )        :: deg_in        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity_in !! continuity: 
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in       !! grid information
  INTEGER       , INTENT(IN   )        :: degGP_in      !! gauss quadrature points: nGP=degGP+1 
  INTEGER        , INTENT(IN   ) :: mn_max_in(2)  !! maximum mode in m and n
  INTEGER        , INTENT(IN   ) :: mn_nyq_in(2)  !! number of integration points 
  INTEGER        , INTENT(IN   ) :: nfp_in        !! number of field periods 
  CHARACTER(LEN=8),INTENT(IN   ) :: sin_cos_in    !! can be either only sine: " _sin_" only cosine: " _cos_" or full: "_sin_cos_"
  LOGICAL         ,INTENT(IN   ) :: exclude_mn_zero_in !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_base), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(t_base :: sf)
  CALL sbase_new(sf%s,deg_in,continuity_in)
  CALL fbase_new(sf%f)
  CALL sf%s%init(grid_in,degGP_in)
  CALL sf%f%init(mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
  sf%initialized=.TRUE.

  IF(.NOT.test_called) CALL Base_test(sf)

END SUBROUTINE base_new


!===================================================================================================================================
!> finalize the type base
!!
!===================================================================================================================================
SUBROUTINE base_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_base), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.sf%initialized)RETURN
CALL sf%s%free()
CALL sf%f%free()

sf%initialized=.FALSE.
END SUBROUTINE base_free


!===================================================================================================================================
!> copy the type base
!!
!===================================================================================================================================
SUBROUTINE base_copy( sf , tocopy)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: tocopy
  CLASS(t_base), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL sf%s%copy(tocopy%s)
CALL sf%f%copy(tocopy%f)
END SUBROUTINE base_copy

!===================================================================================================================================
!> evaluate all degrees of freedom at all Gauss Points (deriv=0 solution, deriv=1 first derivative d/ds)
!!
!===================================================================================================================================
FUNCTION base_evalDOF(sf,deriv,DOFs) RESULT(y_IP_GP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: sf     !! self
  INTEGER      , INTENT(IN   ) :: deriv   !! =0: base, =1: ds, =2: dtheta, =3: dzeta
  REAL(wp)     , INTENT(IN   ) :: DOFs(1:sf%f%modes,1:sf%s%nBase)  !! array of all modes and all radial dofs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                     :: y_IP_GP(sf%f%mn_IP,sf%s%nGP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                      :: iGP,iMode
  INTEGER,PARAMETER            :: deriv_map_GP(0:3)=(/0,DERIV_S,         0,         0/) !map input deriv to ds deriv
  INTEGER,PARAMETER            :: deriv_map_IP(0:3)=(/0,      0,DERIV_THET,DERIV_ZETA/) !map input deriv to dtheta/dzeta deriv
  REAL(wp)                     :: y_tmp(1:sf%f%modes,sf%s%nGP)
!===================================================================================================================================
  DO iMode=1,sf%f%modes
    y_tmp(iMode,:)=sf%s%evalDOF_GP(deriv_map_GP(deriv),DOFs(iMode,:))
  END DO!iMode
  DO iGP=1,sf%s%nGP
    y_IP_GP(:,iGP)=sf%f%evalDOF_IP(deriv_map_IP(deriv),y_tmp(:,iGP))
  END DO !iGP

END FUNCTION base_evalDOF


!===================================================================================================================================
!> test base variable
!!
!===================================================================================================================================
SUBROUTINE Base_test( sf )
! MODULES
USE MOD_GLobals, ONLY: UNIT_StdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,iMode,iGP
  CHARACTER(LEN=10)  :: fail
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  REAL(wp)           :: checkreal 
  REAL(wp)           :: dofs(1:sf%f%modes,1:sf%s%nBase)
  REAL(wp)           :: g_sIP(1:sf%s%nBase)
  REAL(wp)           :: g_IP_GP(1:sf%f%mn_IP,1:sf%s%nGP)
!===================================================================================================================================
  test_called=.TRUE.
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN BASE TEST ID',nTestCalled,'  >>>>>>>>>'
  ASSOCIATE(modes=>sf%f%modes,sin_range=>sf%f%sin_range,cos_range=>sf%f%cos_range, &
            deg=>sf%s%deg,nBase=>sf%s%nBase,sin_cos=>sf%f%sin_cos,nGP=>sf%s%nGP)
  IF(testlevel.LE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    g_IP_GP(:,:)=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      g_sIP(:)=(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(iMode,:)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +SIN(sf%f%Xmn(1,iMode)*sf%f%x_IP(1,:)-sf%f%Xmn(2,iMode)*sf%f%x_IP(2,:))* &
                        (0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      g_sIP(:)=(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(iMode,:)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +COS(sf%f%Xmn(1,iMode)*sf%f%x_IP(1,:)-sf%f%Xmn(2,iMode)*sf%f%x_IP(2,:))* &
                        (0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode

    g_IP_GP=g_IP_GP - sf%evalDOF(0,dofs)

    checkreal=MAXVAL(ABS(g_IP_GP))
    IF(testdbg.OR.(.NOT.(checkreal .LT. realtol) )) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! BASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,E11.3))') &
      '\n =>  should be 0.0 : MAX(|g_IP_exact-g_IP(dofs)|) = ', checkreal
    END IF !TEST

  END IF !testlevel<1
  END ASSOCIATE !sf

  test_called=.FALSE.

END SUBROUTINE Base_test



END MODULE MOD_base

