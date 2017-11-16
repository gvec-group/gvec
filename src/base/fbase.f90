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
!!# Module ** fBase **
!!
!! 2D Fourier base in the two angular directions: (poloidal,toroidal) ~ (m,n) ~ (theta,zeta) [0,2pi]x[0,2pi/nfp]
!!
!! explicit real fourier basis: sin(x_mn) or cos(x_mn) with x_mn=(m*theta - n*nfp*zeta)  ,
!! with mode numbers m and n 
!!
!===================================================================================================================================
MODULE MOD_fBase
! MODULES
USE MOD_Globals                  ,ONLY: TWOPI,wp,Unit_stdOut,abort
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE, ABSTRACT :: c_fBase
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  INTEGER              :: mn_max(2)   !! input parameter: maximum number of fourier modes: m_max=mn_max(1),n_max=mn_max(2)
  INTEGER              :: mn_nyq(2)   !! number of equidistant integration points (trapezoidal rule) in m and n
  INTEGER              :: nfp         !! number of field periods (toroidal repetition after 2pi/nfp)
  INTEGER              :: sin_cos     !! can be either only sine: _SIN_  or only cosine _COS_ or full: _SINCOS_
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: modes       !! total number of modes in basis (depends if only sin/cos or sin & cos are used)
  CONTAINS
    PROCEDURE(i_sub_fBase_init          ),DEFERRED :: init
    PROCEDURE(i_sub_fBase_free          ),DEFERRED :: free
    PROCEDURE(i_sub_fBase_copy          ),DEFERRED :: copy
    PROCEDURE(i_fun_fBase_evalDOF_x     ),DEFERRED :: evalDOF_x
    PROCEDURE(i_fun_fBase_evalDOF_IP    ),DEFERRED :: evalDOF_IP
    PROCEDURE(i_fun_fBase_initDOF       ),DEFERRED :: initDOF

END TYPE c_fBase

ABSTRACT INTERFACE
  SUBROUTINE i_sub_fBase_init(sf, mn_max_in,mn_nyq_in,nfp_in, &
                                  sin_cos_in,exclude_mn_zero_in)
    IMPORT c_fBase
    CLASS(c_fBase) , INTENT(INOUT) :: sf
    INTEGER        , INTENT(IN   ) :: mn_max_in(2)
    INTEGER        , INTENT(IN   ) :: mn_nyq_in(2)
    INTEGER        , INTENT(IN   ) :: nfp_in      
    CHARACTER(LEN=8),INTENT(IN   ) :: sin_cos_in  
    LOGICAL         ,INTENT(IN   ) :: exclude_mn_zero_in
  END SUBROUTINE i_sub_fBase_init

  SUBROUTINE i_sub_fBase_free( sf ) 
    IMPORT c_fBase
    CLASS(c_fBase), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_fBase_free

  SUBROUTINE i_sub_fBase_copy( sf, tocopy ) 
    IMPORT c_fBase
    CLASS(c_fBase), INTENT(IN   ) :: tocopy
    CLASS(c_fBase), INTENT(INOUT) :: sf
  END SUBROUTINE i_sub_fBase_copy

  FUNCTION i_fun_fBase_initDOF( sf, g_IP ) RESULT(DOFs) 
    IMPORT wp,c_fBase
    CLASS(c_fBase), INTENT(IN   ) :: sf
    REAL(wp)      , INTENT(IN   ) :: g_IP(sf%mn_nyq(1),sf%mn_nyq(2))
    REAL(wp)                      :: DOFs(1:sf%modes)
  END FUNCTION i_fun_fBase_initDOF

  FUNCTION i_fun_fBase_evalDOF_x( sf, x,deriv,DOFs ) RESULT(y) 
    IMPORT wp,c_fBase
  CLASS(c_fBase), INTENT(IN   ) :: sf
  REAL(wp)      , INTENT(IN   ) :: x(2)
  INTEGER       , INTENT(IN   ) :: deriv(2) 
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%modes)
  REAL(wp)                      :: y
  END FUNCTION i_fun_fBase_evalDOF_x

  FUNCTION i_fun_fBase_evalDOF_IP( sf,deriv,DOFs ) RESULT(y_IP)
    IMPORT wp,c_fBase
  CLASS(c_fBase), INTENT(IN   ) :: sf
  INTEGER       , INTENT(IN   ) :: deriv(2) 
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%modes)
  REAL(wp)                      :: y_IP(sf%mn_nyq(1),sf%mn_nyq(2))
  END FUNCTION i_fun_fBase_evalDOF_IP

END INTERFACE
 


TYPE,EXTENDS(c_fBase) :: t_fBase
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters
  LOGICAL              :: exclude_mn_zero=.FALSE.  !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp)             :: d_thet              !! integration weight in theta direction: =2pi/mn_nyq(1)
  REAL(wp)             :: d_zeta              !! integration weight in zeta direction : =nfp*(2pi/nfp)/mn_nyq(2)=2*pi/mn_nyq(2)
  INTEGER,ALLOCATABLE  :: Xmn(:,:)            !! mode number m,n for each iMode=1,modes, size(2,modes)
  INTEGER,ALLOCATABLE  :: zero_odd_even(:)    !! =0 for m=n=0 mode, =1 for m= odd mode, =2 for m=even mode size(modes)
  REAL(wp),ALLOCATABLE :: x_IP(:,:)           !! (theta,zeta)position of interpolation points theta [0,2pi    ]x[0,2pi/nfp] size(2,mn_nyq(1)*mn_nyq(2))
  REAL(wp),ALLOCATABLE :: base_IP(:,:)        !! basis functions,                 size(1:mn_nyq(1)*mn_nyq(2),1:mn_modes)
  REAL(wp),ALLOCATABLE :: base_dthet_IP(:,:)  !! dthet derivative of basis functions, (1:mn_nyq(1)*mn_nyq(2),1:mn_modes)
  REAL(wp),ALLOCATABLE :: base_dzeta_IP(:,:)  !! dzeta derivative of basis functions, (1:mn_nyq(1)*mn_nyq(2),1:mn_modes)

  
  CONTAINS

  PROCEDURE :: init          => fBase_init
  PROCEDURE :: free          => fBase_free
  PROCEDURE :: copy          => fBase_copy
  PROCEDURE :: evalDOF_x     => fBase_evalDOF_x
  PROCEDURE :: evalDOF_IP    => fBase_evalDOF_IP
  PROCEDURE :: initDOF       => fBase_initDOF

END TYPE t_fBase

LOGICAL  :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> initialize the type fBase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE fBase_new( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(t_fBase :: sf)
END SUBROUTINE fBase_new

!===================================================================================================================================
!> initialize the type fBase with polynomial degree, continuity ( -1: disc, 1: full)
!! and number of gauss points per element
!!
!===================================================================================================================================
SUBROUTINE fBase_init( sf, mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER        , INTENT(IN   ) :: mn_max_in(2)  !! maximum mode in m and n
  INTEGER        , INTENT(IN   ) :: mn_nyq_in(2)  !! number of integration points 
  INTEGER        , INTENT(IN   ) :: nfp_in        !! number of field periods 
  CHARACTER(LEN=8),INTENT(IN   ) :: sin_cos_in    !! can be either only sine: " _sin_" only cosine: " _cos_" or full: "_sin_cos_"
  LOGICAL         ,INTENT(IN   ) :: exclude_mn_zero_in !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i,iMode,m,n,mIP,nIP,mn_excl
  INTEGER :: modes_sin,modes_cos
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(4X,A,2(A,I6," , ",I6),A,I4,A,L2,A)')'INIT fBase type:', &
       ' mn_max= (',mn_max_in, &
       ' ), mn_nyq = ',mn_nyq_in, &
       ' )\n      nfp    = ',nfp_in, &
       ' exclude_mn_zero = ',exclude_mn_zero_in, &
       ' ,  sin/cos : '//TRIM(sin_cos_in)//' ...'
  IF(sf%initialized) THEN
    CALL abort(__STAMP__, &
        "Trying to reinit fBase!") 
  END IF
  IF((mn_nyq_in(1)/(mn_max_in(1)+1)).LT.1) &
    CALL abort(__STAMP__, &
        "error in fBase: mn_nyq in theta should be > mn_max(1)!") 
  IF((mn_nyq_in(2)/(mn_max_in(2)+1)).LT.1) &
    CALL abort(__STAMP__, &
         "error in fBase: mn_nyq in zeta should be > mn_max(2)!") 

  sf%mn_max(1:2)  = mn_max_in(1:2)
  sf%mn_nyq(1:2)  = mn_nyq_in(1:2)
  sf%nfp          = nfp_in
  sf%exclude_mn_zero  = exclude_mn_zero_in
  IF(INDEX(TRIM(sin_cos_in),"_sin_").NE.0) THEN
    WRITE(*,*)'DEBUG found _sin_' 
    sf%sin_cos  = _SIN_
    
  ELSEIF(INDEX(TRIM(sin_cos_in),"_cos_").NE.0) THEN
    sf%sin_cos  = _COS_
  ELSEIF(INDEX(TRIM(sin_cos_in),"_sincos_").NE.0) THEN
    sf%sin_cos  = _SINCOS_
  ELSE 
    CALL abort(__STAMP__, &
         "error in fBase: sin/cos not correctly specified, should be either _SIN_, _COS_ or _SINCOS_") 
  END IF
  ASSOCIATE(&
              m_max      => sf%mn_max(1)  &
            , n_max      => sf%mn_max(2)  &
            , m_nyq      => sf%mn_nyq(1)  &
            , n_nyq      => sf%mn_nyq(2)  &
            , nfp         => sf%nfp       &
            , sin_cos     => sf%sin_cos   &
            , modes       => sf%modes     &
            )
  mn_excl=MERGE(1,0,sf%exclude_mn_zero) !=1 if exclude=TRUE, =0 if exclude=FALSE
  ! modes_sin :: m=0: n=1...nMax , m=1...m_max: n=-n_max...n_max. REMARK: for sine, m=0,n=0 is automatically excluded
  ! mode_ cos :: m=0: n=0...nMax , m=1...m_max: n=-n_max...n_max. mn_excl=True will exclude m=n=0 
  modes_sin= (n_max  )         + m_max*(2*n_max+1) 
  modes_cos= (n_max+1)-mn_excl + m_max*(2*n_max+1) 
  SELECT CASE(sin_cos)
  CASE(_SIN_)
    modes  = modes_sin 
  CASE(_COS_)
    modes  = modes_cos 
  CASE(_SINCOS_)
    modes  = modes_sin+modes_cos 
  END SELECT

  CALL fbase_alloc(sf)
  
  iMode=0
  IF((sin_cos.EQ._SIN_).OR.(sin_cos.EQ._SINCOS_) )THEN
    !SINE
    m=0
    DO n=1,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/) 
    END DO !n
    DO m=1,m_max; DO n=-n_max,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/) 
    END DO; END DO !m,n
  ELSEIF((sin_cos.EQ._COS_).OR.(sin_cos.EQ._SINCOS_) )THEN
    !COSINE (for _SINCOS_, it comes after sine)
    m=0
    DO n=mn_excl,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/) 
    END DO !n
    DO m=1,m_max; DO n=-n_max,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/) 
    END DO; END DO !m,n
  END IF !_SIN_,_COS_,_SINCOS_

  IF(iMode.NE.modes) STOP' Problem in Xmn '

  DO iMode=1,modes
    IF((sf%Xmn(1,iMode).EQ.0).AND.(sf%Xmn(2,iMode).EQ.0))THEN
      sf%zero_odd_even(iMode)=MN_ZERO
    ELSE
      IF(MOD(sf%Xmn(1,iMode),2).EQ.0)THEN
        sf%zero_odd_even(iMode)=M_EVEN
      ELSE
        sf%zero_odd_even(iMode)=M_ODD
      END IF
    END IF
  END DO !iMode=1,modes

  sf%d_thet = TWOPI/REAL(m_nyq,wp)
  sf%d_zeta = TWOPI/REAL(n_nyq*nfp,wp)

  i=0
  DO nIP=1,n_nyq
    DO mIP=1,m_nyq
      i=i+1
      sf%x_IP(1,i)=(REAL(mIP,wp)-0.5_wp)*sf%d_thet
      sf%x_IP(2,i)=(REAL(nIP,wp)-0.5_wp)*sf%d_zeta
    END DO !m
  END DO !n

  sf%d_zeta = sf%d_zeta*REAL(nfp,wp) ! to get full integral [0,2pi)
 
  iMode=0
  IF((sin_cos.EQ._SIN_).OR.(sin_cos.EQ._SINCOS_) )THEN
    !SINE
    DO i=1,modes_sin
      iMode=iMode+1
      !TODO base_IP(:,:,iMode)=
      !TODO base_dthet_IP(:,:,iMode)=
      !TODO base_dzeta_IP(:,:,iMode)=
    END DO !n
  ELSEIF((sin_cos.EQ._COS_).OR.(sin_cos.EQ._SINCOS_) )THEN
    !COSINE (for _SINCOS_, it comes after sine)
    !SINE
    DO i=1,modes_cos
      iMode=iMode+1
      !TODO base_IP(:,:,iMode)=
      !TODO base_dthet_IP(:,:,iMode)=
      !TODO base_dzeta_IP(:,:,iMode)=
    END DO !n
  END IF !_SIN_,_COS_,_SINCOS_

  END ASSOCIATE !sf


  sf%initialized=.TRUE.
  SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
  IF(.NOT.test_called) CALL fBase_test(sf)

END SUBROUTINE fBase_init


!===================================================================================================================================
!> allocate all variables in  fBase
!!
!===================================================================================================================================
SUBROUTINE fBase_alloc( sf)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ASSOCIATE(&
              m_nyq      => sf%mn_nyq(1)  &
            , n_nyq      => sf%mn_nyq(2)  &
            , modes       => sf%modes     &
            )
  ALLOCATE(sf%Xmn(        2,1:modes))
  ALLOCATE(sf%zero_odd_even(1:modes))
  ALLOCATE(sf%x_IP(       2,1:m_nyq*n_nyq) )
  ALLOCATE(sf%base_IP(      1:m_nyq*n_nyq,1:modes) )
  ALLOCATE(sf%base_dthet_IP(1:m_nyq*n_nyq,1:modes) )
  ALLOCATE(sf%base_dzeta_IP(1:m_nyq*n_nyq,1:modes) )
  END ASSOCIATE !m_nyq,n_nyq,modes
END SUBROUTINE fBase_alloc


!===================================================================================================================================
!> finalize the type fBase
!!
!===================================================================================================================================
SUBROUTINE fBase_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) RETURN
  !allocatables  
  SDEALLOCATE(sf%Xmn)
  SDEALLOCATE(sf%zero_odd_even)
  SDEALLOCATE(sf%x_IP)
  SDEALLOCATE(sf%base_IP)
  SDEALLOCATE(sf%base_dthet_IP)
  SDEALLOCATE(sf%base_dzeta_IP)
  
  sf%mn_max     =-1
  sf%mn_nyq     =-1 
  sf%nfp        =-1        
  sf%modes      =-1
  sf%sin_cos    =-1
  sf%d_thet     =0.0_wp
  sf%d_zeta     =0.0_wp
  sf%exclude_mn_zero=.FALSE.
  sf%initialized=.FALSE.

END SUBROUTINE fBase_free


!===================================================================================================================================
!> copy the type fBase
!!
!===================================================================================================================================
SUBROUTINE fBase_copy( sf , tocopy)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(c_fBase), INTENT(IN   ) :: tocopy
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=8) :: sin_cos
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPEIS(t_fBase)
  IF(.NOT.tocopy%initialized) THEN
    CALL abort(__STAMP__, &
        "fBase_copy: not initialized fBase from which to copy!")
  END IF
  IF(sf%initialized) THEN
    SWRITE(UNIT_stdOut,'(A)')'WARNING!! reinit of fBase in copy!'
    CALL sf%free()
  END IF
  SELECT CASE(tocopy%sin_cos)
  CASE(_SIN_)
    sin_cos  = "_sin_"
  CASE(_COS_)
    sin_cos  = "_cos_"
  CASE(_SINCOS_)
    sin_cos  = "_sincos_"
  END SELECT
 CALL sf%init(tocopy%mn_max         &
             ,tocopy%mn_nyq         &
             ,tocopy%nfp            &
             ,sin_cos               &
             ,tocopy%exclude_mn_zero)

  END SELECT !TYPE
END SUBROUTINE fBase_copy


!===================================================================================================================================
!> simply evaluate function or derivative at point x
!!
!===================================================================================================================================
FUNCTION fBase_evalDOF_x(sf,x,deriv,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  REAL(wp)      , INTENT(IN   ) :: x(2)   !! point position in [0,2pi]x[0,2pi/nfp]
  INTEGER       , INTENT(IN   ) :: deriv(2)  !! derivative in theta/zeta
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%modes)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
y=0.0_wp !TODO

END FUNCTION fBase_evalDOF_x

!===================================================================================================================================
!> evaluate all degrees of freedom at all Gauss Points (deriv=0 solution, deriv=1 first derivative d/ds)
!!
!===================================================================================================================================
FUNCTION fBase_evalDOF_IP(sf,deriv,DOFs) RESULT(y_IP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: deriv(2)  !! only =0 or =1, derivative in theta and zeta 
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%modes)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y_IP(sf%mn_nyq(1),sf%mn_nyq(2)) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
y_IP=0.0_wp !TODO
END FUNCTION fBase_evalDOF_IP

!===================================================================================================================================
!>  take values interpolated at sf%s_IP positions and give back the degrees of freedom 
!!
!===================================================================================================================================
FUNCTION fBase_initDOF( sf , g_IP) RESULT(DOFs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf    !! self
  REAL(wp)      , INTENT(IN   ) :: g_IP(sf%mn_nyq(1),sf%mn_nyq(2))  !!  interpolation values at theta_IP zeta_IP positions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: DOFs(1:sf%modes)  !! projection to fourier base
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END FUNCTION fBase_initDOF

!===================================================================================================================================
!> test fBase variable
!!
!===================================================================================================================================
SUBROUTINE fBase_test( sf)
! MODULES
USE MOD_GLobals, ONLY: testdbg,testlevel,nfailedMsg,nTestCalled,testfailedMsg
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest
  REAL(wp)           :: checkreal,refreal
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
!===================================================================================================================================
  test_called=.TRUE.
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN FBASE TEST ID',nTestCalled,'    >>>>>>>>>'
  ASSOCIATE(&
              m_max      => sf%mn_max(1)  &
            , n_max      => sf%mn_max(2)  &
            , m_nyq      => sf%mn_nyq(1)  &
            , n_nyq      => sf%mn_nyq(2)  &
            , nfp         => sf%nfp       &
            , sin_cos     => sf%sin_cos   &
            , modes       => sf%modes     &
            )
  IF(testlevel.LE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal =SUM(sf%x_IP(1,:)*sf%x_IP(2,:))*sf%d_thet*sf%d_zeta
    refreal   =(0.5_wp*(TWOPI)**2)*REAL(nfp,wp)*(0.5_wp*(TWOPI/REAL(nfp,wp))**2)

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! fBase TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,2(A,I4),2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : ', sin_cos, &
      '\n =>  should be ', refreal,' : nfp*int(int(theta*zeta, 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

  END IF !testlevel <=1
  END ASSOCIATE !sf
   
END SUBROUTINE fBase_test


END MODULE MOD_fBase

