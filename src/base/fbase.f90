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
  INTEGER              :: mn_IP       !! =mn_nyq(1)*mn_nyq(2) 
  INTEGER              :: nfp         !! number of field periods (toroidal repetition after 2pi/nfp)
  INTEGER              :: sin_cos     !! can be either only sine: _SIN_  or only cosine _COS_ or full: _SINCOS_
  !input parameters
  LOGICAL              :: exclude_mn_zero  !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: modes       !! total number of modes in basis (depends if only sin/cos or sin & cos are used)
  CONTAINS
    PROCEDURE(i_sub_fBase_init          ),DEFERRED :: init
    PROCEDURE(i_sub_fBase_free          ),DEFERRED :: free
    PROCEDURE(i_sub_fBase_copy          ),DEFERRED :: copy
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
    REAL(wp)      , INTENT(IN   ) :: g_IP(sf%mn_IP)
    REAL(wp)                      :: DOFs(1:sf%modes)
  END FUNCTION i_fun_fBase_initDOF

  FUNCTION i_fun_fBase_evalDOF_IP( sf,deriv,DOFs ) RESULT(y_IP)
    IMPORT wp,c_fBase
  CLASS(c_fBase), INTENT(IN   ) :: sf
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%modes)
  REAL(wp)                      :: y_IP(sf%mn_IP)
  END FUNCTION i_fun_fBase_evalDOF_IP

END INTERFACE
 


TYPE,EXTENDS(c_fBase) :: t_fBase
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: sin_range(2)       !! sin_range(1)+1:sin_range(2) is range with sine bases 
  INTEGER              :: cos_range(2)       !! sin_range(1)+1:sin_range(2) is range with sine bases 
  REAL(wp)             :: d_thet              !! integration weight in theta direction: =2pi/mn_nyq(1)
  REAL(wp)             :: d_zeta              !! integration weight in zeta direction : =nfp*(2pi/nfp)/mn_nyq(2)=2*pi/mn_nyq(2)
  INTEGER,ALLOCATABLE  :: Xmn(:,:)            !! mode number m,n for each iMode=1,modes, size(2,modes)
  INTEGER,ALLOCATABLE  :: zero_odd_even(:)    !! =0 for m=n=0 mode, =1 for m= odd mode, =2 for m=even mode size(modes)
  REAL(wp),ALLOCATABLE :: x_IP(:,:)           !! (theta,zeta)position of interpolation points theta [0,2pi]x[0,2pi/nfp]size(2,mn_IP)
  REAL(wp),ALLOCATABLE :: base_IP(:,:)        !! basis functions,                 size(1:mn_IP,1:mn_modes)
  REAL(wp),ALLOCATABLE :: base_dthet_IP(:,:)  !! dthet derivative of basis functions, (1:mn_IP,1:mn_modes)
  REAL(wp),ALLOCATABLE :: base_dzeta_IP(:,:)  !! dzeta derivative of basis functions, (1:mn_IP,1:mn_modes)

  
  CONTAINS

  PROCEDURE :: init             => fBase_init
  PROCEDURE :: free             => fBase_free
  PROCEDURE :: copy             => fBase_copy
  PROCEDURE :: evalDOF_IP       => fBase_evalDOF_IP
  PROCEDURE :: initDOF          => fBase_initDOF

END TYPE t_fBase

LOGICAL  :: test_called=.FALSE.

CHARACTER(LEN=8)   :: sin_cos_map(3)=(/"_sin_   ", &
                                       "_cos_   ", &
                                       "_sincos_" /)

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> allocate the type fBase 
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
!> initialize the type fBase maximum mode numbers, number of integration points, type of basis (sin/cos or sin and cos) 
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
  LOGICAL         ,INTENT(IN   ) :: exclude_mn_zero_in !! =true: exclude m=n=0 mode in the basis (only important if cos is in basis)
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
  sf%mn_IP        = sf%mn_nyq(1)*sf%mn_nyq(2)
  sf%nfp          = nfp_in
  sf%exclude_mn_zero  = exclude_mn_zero_in
  IF(INDEX(TRIM(sin_cos_in),TRIM(sin_cos_map(_SIN_))).NE.0) THEN
    WRITE(*,*)'DEBUG found _sin_' 
    sf%sin_cos  = _SIN_
    
  ELSEIF(INDEX(TRIM(sin_cos_in),TRIM(sin_cos_map(_COS_))).NE.0) THEN
    sf%sin_cos  = _COS_
  ELSEIF(INDEX(TRIM(sin_cos_in),TRIM(sin_cos_map(_SINCOS_))).NE.0) THEN
    sf%sin_cos  = _SINCOS_
  ELSE 
    CALL abort(__STAMP__, &
         "error in fBase: sin/cos not correctly specified, should be either _SIN_, _COS_ or _SINCOS_") 
  END IF
  ASSOCIATE(&
              m_max      => sf%mn_max(1)   &
            , n_max      => sf%mn_max(2)   &
            , m_nyq      => sf%mn_nyq(1)   &
            , n_nyq      => sf%mn_nyq(2)   &
            , nfp         => sf%nfp        &
            , sin_cos     => sf%sin_cos    &
            , modes       => sf%modes      &
            , sin_range  => sf%sin_range &
            , cos_range  => sf%cos_range &
            )
  mn_excl=MERGE(1,0,sf%exclude_mn_zero) !=1 if exclude=TRUE, =0 if exclude=FALSE
  ! modes_sin :: m=0: n=1...nMax , m=1...m_max: n=-n_max...n_max. REMARK: for sine, m=0,n=0 is automatically excluded
  ! mode_ cos :: m=0: n=0...nMax , m=1...m_max: n=-n_max...n_max. mn_excl=True will exclude m=n=0 
  modes_sin= (n_max  )         + m_max*(2*n_max+1) 
  modes_cos= (n_max+1)-mn_excl + m_max*(2*n_max+1) 
  SELECT CASE(sin_cos)
  CASE(_SIN_)
    modes        = modes_sin 
    sin_range(:) = (/0,modes_sin/)
    cos_range(:) = (/0,0/) !off
  CASE(_COS_)
    modes         = modes_cos 
    sin_range(:) = (/0,0/) !off
    cos_range(:) = (/0,modes_cos/)
  CASE(_SINCOS_)
    modes        = modes_sin+modes_cos 
    sin_range(:) = (/0,modes_sin/)
    cos_range(:) = (/modes_sin,modes/) 
  END SELECT

  CALL fbase_alloc(sf)
  
  iMode=0 !first sine then cosine
  !SINE
  IF((sin_range(2)-sin_range(1)).EQ.modes_sin)THEN
    m=0
    DO n=1,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/)  !include nfp here 
    END DO !n
    DO m=1,m_max; DO n=-n_max,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/)  !include nfp here 
    END DO; END DO !m,n
  END IF !sin_range>0

  !COSINE (for _SINCOS_, it comes after sine)
  IF((cos_range(2)-cos_range(1)).EQ.modes_cos)THEN
    m=0
    DO n=mn_excl,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/)  !include nfp here 
    END DO !n
    DO m=1,m_max; DO n=-n_max,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/)  !include nfp here  
    END DO; END DO !m,n
  END IF !cos_range>0

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
 
  !SIN(m*theta-(n*nfp)*zeta)
  DO iMode=sin_range(1)+1,sin_range(2)
    sf%base_IP(:,iMode)      =                          SIN(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    sf%base_dthet_IP(:,iMode)= REAL(sf%Xmn(1,iMode),wp)*COS(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    sf%base_dzeta_IP(:,iMode)=-REAL(sf%Xmn(2,iMode),wp)*COS(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
  END DO !iMode

  !COS(m*theta-(n*nfp)*zeta)
  DO iMode=cos_range(1)+1,cos_range(2)
    sf%base_IP(:,iMode)      =                          COS(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    sf%base_dthet_IP(:,iMode)=-REAL(sf%Xmn(1,iMode),wp)*SIN(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    sf%base_dzeta_IP(:,iMode)= REAL(sf%Xmn(2,iMode),wp)*SIN(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
  END DO !iMode

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
              mn_IP     => sf%mn_IP     &
            , modes     => sf%modes     &
            )
  ALLOCATE(sf%Xmn(        2,1:modes))
  ALLOCATE(sf%zero_odd_even(1:modes))
  ALLOCATE(sf%x_IP(       2,1:mn_IP) )
  ALLOCATE(sf%base_IP(      1:mn_IP,1:modes) )
  ALLOCATE(sf%base_dthet_IP(1:mn_IP,1:modes) )
  ALLOCATE(sf%base_dzeta_IP(1:mn_IP,1:modes) )
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
  sf%mn_IP      =-1 
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
!> evaluate  all modes at all interpolation points 
!!
!===================================================================================================================================
FUNCTION fBase_evalDOF_IP(sf,deriv,DOFs) RESULT(y_IP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: DOFs(1:sf%modes)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y_IP(sf%mn_IP) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(deriv)
  CASE(0)
    y_IP=MATMUL(sf%base_IP(:,:),DOFs(:))
  CASE(DERIV_THET)
    y_IP=MATMUL(sf%base_dthet_IP(:,:),DOFs(:))
  CASE(DERIV_ZETA)
    y_IP=MATMUL(sf%base_dzeta_IP(:,:),DOFs(:))
  END SELECT
END FUNCTION fBase_evalDOF_IP

!===================================================================================================================================
!>  take values interpolated at sf%s_IP positions and project onto fourier basis by integration 
!!
!===================================================================================================================================
FUNCTION fBase_initDOF( sf , g_IP) RESULT(DOFs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf    !! self
  REAL(wp)      , INTENT(IN   ) :: g_IP(sf%mn_IP)  !!  interpolation values at theta_IP zeta_IP positions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: DOFs(1:sf%modes)  !! projection to fourier base
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iMode
!===================================================================================================================================
  DO iMode=1,sf%modes
    DOFs(iMode)=(sf%d_thet*sf%d_zeta/((TWOPI)**2))*SUM(sf%base_IP(:,iMode)*g_IP(:))
  END DO
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
  INTEGER            :: iTest,iMode,jMode,ncoszero,nsinzero
  REAL(wp)           :: checkreal,refreal
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: dofs(1:sf%modes)
  REAL(wp)           :: g_IP(1:sf%mn_IP)
!===================================================================================================================================
  test_called=.TRUE. !avoid infinite loop if init is called here
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
            , mn_IP      => sf%mn_IP      &
            , nfp        => sf%nfp        &
            , sin_cos    => sf%sin_cos    &
            , sin_range  => sf%sin_range  &
            , cos_range  => sf%cos_range  &
            , modes      => sf%modes      &
            )
  IF(testlevel.LE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal =SUM(sf%x_IP(1,:)*sf%x_IP(2,:))*sf%d_thet*sf%d_zeta
    refreal   =(0.5_wp*(TWOPI)**2)*REAL(nfp,wp)*(0.5_wp*(TWOPI/REAL(nfp,wp))**2)

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,2(A,I4),2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : ', sin_cos, &
      '\n =>  should be ', refreal,' : nfp*int(int(theta*zeta, 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    checkreal=0.0_wp
    DO iMode=1,modes
      DO jMode=1,modes
        IF(iMode.NE.jMode)THEN
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_IP(:,jMode))))
        END IF !iMode /=jMode
      END DO
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : nfp*int(int(base(imode)*base(jmode), 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

    iTest=103 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    checkreal=0.0_wp
    nsinzero=0
    DO iMode=sin_range(1)+1,sin_range(2)
      checkreal=checkreal+   ((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_IP(:,iMode)))
      IF(sf%zero_odd_even(iMode).EQ.MN_ZERO) nsinzero=nsinzero+1
    END DO
    ncoszero=0
    DO iMode=cos_range(1)+1,cos_range(2)
      checkreal=checkreal+   ((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_IP(:,iMode)))
      IF(sf%zero_odd_even(iMode).EQ.MN_ZERO) ncoszero=ncoszero+1
    END DO
    checkreal=checkreal/REAL(modes,wp)
    refreal=(TWOPI)**2 *( 0.5*(REAL(cos_range(2)-cos_range(1)-ncoszero,wp) + REAL(sin_range(2)-sin_range(1),wp))  &
                         +REAL(ncoszero,wp) )/REAL(modes,wp)

    IF(testdbg.OR.(.NOT.( (ABS(checkreal-refreal).LT. realtol).AND. & 
                          (nsinzero              .EQ. 0      )      ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,(A,I4),A,(A,I4),2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
       '\n =>  should be  0 : nsinzero = ', nsinzero,  &
       '\n =>  should be ', refreal,' : nfp*int(int(base(imode)*base(imode), 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

    iTest=104 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    checkreal=0.0_wp
    DO iMode=1,modes
      DO jMode=1,modes
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dthet_IP(:,jMode))))
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dzeta_IP(:,jMode))))
      END DO
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : nfp*int(int(base(imode)*base_dthet/dzeta(jmode), 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST
  END IF !testlevel <=1
  IF (testlevel .LE.2)THEN

    iTest=201 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*SIN(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    END DO !iMode 
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*COS(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    END DO !iMode 
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(0,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF(dofs)|) ', checkreal
    END IF !TEST

    iTest=202 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(sf%Xmn(1,iMode),wp)*COS(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    END DO !iMode 
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-sf%Xmn(1,iMode),wp)*SIN(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    END DO !iMode 
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_THET,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dthet(dofs)|) ', checkreal
    END IF !TEST

    iTest=203 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-sf%Xmn(2,iMode),wp)*COS(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    END DO !iMode 
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL( sf%Xmn(2,iMode),wp)*SIN(sf%Xmn(1,iMode)*sf%x_IP(1,:)-sf%Xmn(2,iMode)*sf%x_IP(2,:))
    END DO !iMode 
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_ZETA,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testfailedMsg(nfailedMsg),'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dthet(dofs)|) ', checkreal
    END IF !TEST
  END IF !testlevel <=2
  END ASSOCIATE !sf
  test_called=.FALSE.   

END SUBROUTINE fBase_test


END MODULE MOD_fBase

