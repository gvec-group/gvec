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
MODULE MODgvec_fBase
! MODULES
USE MODgvec_Globals                  ,ONLY: TWOPI,wp,Unit_stdOut,abort
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
    PROCEDURE(i_sub_fBase_compare       ),DEFERRED :: compare
    PROCEDURE(i_sub_fBase_change_base   ),DEFERRED :: change_base
    PROCEDURE(i_fun_fBase_eval          ),DEFERRED :: eval
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

  SUBROUTINE i_sub_fBase_compare( sf, tocompare, is_same, cond_out ) 
    IMPORT c_fBase
    CLASS(c_fBase)  , INTENT(IN   ) :: sf
    CLASS(c_fBase)  , INTENT(IN   ) :: tocompare
    LOGICAL,OPTIONAL, INTENT(  OUT) :: is_same
    LOGICAL,OPTIONAL, INTENT(  OUT) :: cond_out(:)
  END SUBROUTINE i_sub_fBase_compare

  SUBROUTINE i_sub_fBase_change_base( sf, old_fBase, iterDim,old_data,sf_data ) 
    IMPORT c_fBase,wp
    CLASS(c_fBase) , INTENT(IN   ) :: sf
    CLASS(c_fBase) , INTENT(IN   ) :: old_fBase
    INTEGER         ,INTENT(IN   ) :: iterDim
    REAL(wp)        ,INTENT(IN   ) :: old_data(:,:)
    REAL(wp)        ,INTENT(  OUT) :: sf_data(:,:)
  END SUBROUTINE i_sub_fBase_change_base

  FUNCTION i_fun_fBase_initDOF( sf, g_IP ) RESULT(DOFs) 
    IMPORT wp,c_fBase
    CLASS(c_fBase), INTENT(IN   ) :: sf
    REAL(wp)      , INTENT(IN   ) :: g_IP(:)
    REAL(wp)                      :: DOFs(1:sf%modes)
  END FUNCTION i_fun_fBase_initDOF

  FUNCTION i_fun_fBase_eval( sf,deriv,x) RESULT(base_x)
    IMPORT wp,c_fBase
  CLASS(c_fBase), INTENT(IN   ) :: sf
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: x(2)
  REAL(wp)                      :: base_x(sf%modes)
  END FUNCTION i_fun_fBase_eval

  FUNCTION i_fun_fBase_evalDOF_x( sf,x,deriv,DOFs ) RESULT(y)
    IMPORT wp,c_fBase
  CLASS(c_fBase), INTENT(IN   ) :: sf
  REAL(wp)      , INTENT(IN   ) :: x(2)
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)
  REAL(wp)                      :: y
  END FUNCTION i_fun_fBase_evalDOF_x

  FUNCTION i_fun_fBase_evalDOF_IP( sf,deriv,DOFs ) RESULT(y_IP)
    IMPORT wp,c_fBase
  CLASS(c_fBase), INTENT(IN   ) :: sf
  INTEGER       , INTENT(IN   ) :: deriv
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)
  REAL(wp)                      :: y_IP(sf%mn_IP)
  END FUNCTION i_fun_fBase_evalDOF_IP

END INTERFACE
 


TYPE,EXTENDS(c_fBase) :: t_fBase
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: sin_range(2)        !! sin_range(1)+1:sin_range(2) is range with sine bases 
  INTEGER              :: cos_range(2)        !! sin_range(1)+1:sin_range(2) is range with sine bases 
  INTEGER              :: mn_zero_mode        !! points to m=0,n=0 mode in mode array (1:mn_modes) (only one can exist for cosine, else =-1)
  REAL(wp)             :: d_thet              !! integration weight in theta direction: =2pi/mn_nyq(1)
  REAL(wp)             :: d_zeta              !! integration weight in zeta direction : =nfp*(2pi/nfp)/mn_nyq(2)=2*pi/mn_nyq(2)
  INTEGER,ALLOCATABLE  :: Xmn(:,:)            !! mode number (m,n*nfp) for each iMode=1,modes, size(2,modes)
  INTEGER,ALLOCATABLE  :: zero_odd_even(:)    !! =0 for m=n=0 mode, =1 for m= odd mode, =2 for m=even mode size(modes)
  REAL(wp),ALLOCATABLE :: x_IP(:,:)           !! (theta,zeta)position of interpolation points theta [0,2pi]x[0,2pi/nfp]size(2,mn_IP)
  REAL(wp),ALLOCATABLE :: base_IP(:,:)        !! basis functions,                 size(1:mn_IP,1:modes)
  REAL(wp),ALLOCATABLE :: base_dthet_IP(:,:)  !! dthet derivative of basis functions, (1:mn_IP,1:modes)
  REAL(wp),ALLOCATABLE :: base_dzeta_IP(:,:)  !! dzeta derivative of basis functions, (1:mn_IP,1:modes)

  REAL(wp),ALLOCATABLE :: snorm_base(:)       !! 1/norm of each basis function, size(1:mn_modes), norm=int_0^2pi int_0^pi (base_mn(thet,zeta))^2 dthet dzeta 
  INTEGER              :: mTotal1D            !! mTotal1D =mn_max(1)+1  for sin or cos base, and mTotal=2*(mn_max(1)+1) for sin&cos base
  REAL(wp),ALLOCATABLE :: base1D_IPthet(:,:,:) !! 1D basis,  size(1:mn_nyq(1),1:2,1:mTotal1D), 
                                               !! if sin(m t-n z):   sin(m t), -cos(m t) and if cos(m t-n z): cos(m t),sin(m t)
  REAL(wp),ALLOCATABLE :: base1D_dthet_IPthet(:,:,:) !! derivative of 1D basis, size(1:mn_nyq(1),1:2,1:mTotal1D)
                                               !! if sin(m t-n z): m cos(m t),m sin(m t) and if cos(m t-n z): -m sin(m t),m cos(m t)
  REAL(wp),ALLOCATABLE :: base1D_IPzeta(:,:,:) !! 1D basis functions, size(1:2,-mn_max(2):mn_max(2),1:mn_nyq(2))
                                               !! for sin/cos(m t-n z): cos(n z),sin(n z)
  REAL(wp),ALLOCATABLE :: base1D_dzeta_IPzeta(:,:,:) !! derivative of 1D basis functions, size(1:2,-mn_max(2):mn_max(2),1:mn_nyq(2))
                                               !! for sin/cos(m t-n z): -n sin(n z),n cos(n z)
  
  CONTAINS

  PROCEDURE :: init             => fBase_init
  PROCEDURE :: free             => fBase_free
  PROCEDURE :: copy             => fBase_copy
  PROCEDURE :: compare          => fBase_compare
  PROCEDURE :: change_base      => fBase_change_base
  PROCEDURE :: eval             => fBase_eval
  PROCEDURE :: evalDOF_x        => fBase_evalDOF_x
! PROCEDURE :: evalDOF_IP       => fBase_evalDOF_IP !use _tens instead!
  PROCEDURE :: evalDOF_IP       => fBase_evalDOF_IP_tens
  PROCEDURE :: initDOF          => fBase_initDOF

END TYPE t_fBase

CHARACTER(LEN=8)   :: sin_cos_map(3)=(/"_sin_   ", &
                                       "_cos_   ", &
                                       "_sincos_" /)

LOGICAL, PRIVATE  :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> allocate the type fBase 
!!
!===================================================================================================================================
SUBROUTINE fBase_new( sf, mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
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
  CLASS(t_fBase), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(t_fBase :: sf)
  call perfon("fbase_new")
  CALL sf%init(mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)

  call perfoff("fbase_new")
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
  REAL(wp):: mm,nn
!===================================================================================================================================
  IF(.NOT.test_called)THEN
    SWRITE(UNIT_stdOut,'(4X,A,2(A,I6," , ",I6),A,I4,A,L2,A)')'INIT fBase type:', &
         ' mn_max= (',mn_max_in, &
         ' ), mn_nyq = ',mn_nyq_in, &
         ' )\n      nfp    = ',nfp_in, &
         ' exclude_mn_zero = ',exclude_mn_zero_in, &
         ' ,  sin/cos : '//TRIM(sin_cos_in)//' ...'
  END IF
  IF(sf%initialized) THEN
    CALL abort(__STAMP__, &
        "Trying to reinit fBase!") 
  END IF
  IF((mn_nyq_in(1)/(mn_max_in(1)+1)).LT.1) &
    CALL abort(__STAMP__, &
        "error in fBase: mn_nyq in theta should be > mn_max(1)!") 
  IF((mn_nyq_in(2)/(mn_max_in(2)+1)).LT.1) &
    CALL abort(__STAMP__, &
         "error in fBase: mn_nyq in zeta should be > mn_max(2)!",mn_nyq_in(2),REAL(mn_max_in(2))) 

  sf%mn_max(1:2)  = mn_max_in(1:2)
  sf%mn_nyq(1:2)  = mn_nyq_in(1:2)
  sf%mn_IP        = sf%mn_nyq(1)*sf%mn_nyq(2)
  sf%nfp          = nfp_in
  sf%exclude_mn_zero  = exclude_mn_zero_in
  IF(INDEX(TRIM(sin_cos_in),TRIM(sin_cos_map(_SIN_))).NE.0) THEN
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
    sf%mTotal1D=m_max+1
  CASE(_COS_)
    modes         = modes_cos 
    sin_range(:) = (/0,0/) !off
    cos_range(:) = (/0,modes_cos/)
    sf%mTotal1D=m_max+1
  CASE(_SINCOS_)
    modes        = modes_sin+modes_cos 
    sin_range(:) = (/0,modes_sin/)
    cos_range(:) = (/modes_sin,modes/) 
    sf%mTotal1D=2*(m_max+1)
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

  sf%mn_zero_mode=-1 !default 

  !COSINE (for _SINCOS_, it comes after sine)
  IF((cos_range(2)-cos_range(1)).EQ.modes_cos)THEN
    m=0
    IF(mn_excl.EQ.0) sf%mn_zero_mode=iMode+1
    DO n=mn_excl,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/)  !include nfp here 
    END DO !n
    DO m=1,m_max; DO n=-n_max,n_max
      iMode=iMode+1
      sf%Xmn(:,iMode)=(/m,n*nfp/)  !include nfp here  
    END DO; END DO !m,n
  END IF !cos_range>0

  IF(iMode.NE.modes) STOP ' Problem in Xmn '

  DO iMode=1,modes
    m=sf%Xmn(1,iMode)
    n=sf%Xmn(2,iMode)
    ! set odd/even/zero for m-modes
    IF((m.EQ.0))THEN  !m=0
      IF((n.EQ.0))THEN !n=0
        sf%zero_odd_even(iMode)=MN_ZERO
      ELSE !n /=0
        sf%zero_odd_even(iMode)=M_ZERO
      END IF
    ELSE  !m /=0
      IF(MOD(m,2).EQ.0)THEN
        sf%zero_odd_even(iMode)=M_EVEN
      ELSE 
        IF(m.EQ.1)THEN
          sf%zero_odd_even(iMode)=M_ODD_FIRST
        ELSE
          sf%zero_odd_even(iMode)=M_ODD
        END IF
      END IF
    END IF
    ! compute 1/norm, with norm=int_0^2pi int_0^2pi (base_mn(thet,zeta))^2 dthet dzeta
    IF((m.EQ.0).AND.(n.EQ.0))THEN !m=n=0 (only needed for cos)
      sf%snorm_base(iMode)=1.0_wp/(TWOPI*TWOPI) !norm=4pi^2 
    ELSE  
      sf%snorm_base(iMode)=2.0_wp/(TWOPI*TWOPI) !norm=2pi^2
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
 
!! !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)  
!!   DO i=1,sf%mn_IP
!!     sf%base_IP(      i,:)=sf%eval(         0,sf%x_IP(:,i))
!!     sf%base_dthet_IP(i,:)=sf%eval(DERIV_THET,sf%x_IP(:,i))
!!     sf%base_dzeta_IP(i,:)=sf%eval(DERIV_ZETA,sf%x_IP(:,i))
!!   END DO
!! !$OMP END PARALLEL DO 

  !SIN(m*theta-(n*nfp)*zeta)
  DO iMode=sin_range(1)+1,sin_range(2)
    mm=REAL(sf%Xmn(1,iMode),wp)
    nn=REAL(sf%Xmn(2,iMode),wp)
    sf%base_IP(:,iMode)      =    SIN(mm*sf%x_IP(1,:)-nn*sf%x_IP(2,:))
    sf%base_dthet_IP(:,iMode)= mm*COS(mm*sf%x_IP(1,:)-nn*sf%x_IP(2,:))
    sf%base_dzeta_IP(:,iMode)=-nn*COS(mm*sf%x_IP(1,:)-nn*sf%x_IP(2,:))
  END DO !iMode

  !COS(m*theta-(n*nfp)*zeta)
  DO iMode=cos_range(1)+1,cos_range(2)
    mm=REAL(sf%Xmn(1,iMode),wp)
    nn=REAL(sf%Xmn(2,iMode),wp)
    sf%base_IP(:,iMode)      =    COS(mm*sf%x_IP(1,:)-nn*sf%x_IP(2,:))
    sf%base_dthet_IP(:,iMode)=-mm*SIN(mm*sf%x_IP(1,:)-nn*sf%x_IP(2,:))
    sf%base_dzeta_IP(:,iMode)= nn*SIN(mm*sf%x_IP(1,:)-nn*sf%x_IP(2,:))
  END DO !iMode

  !!!! 1D BASE
  !sin(mt-nz) = sin(mt)*cos(nz)-cos(mt)*sin(nz) =as(1,mt)*b(1,nz) + as(2,mt)*b(2,nz)
  !cos(mt-nz) = cos(mt)*cos(nz)+sin(mt)*sin(nz) =ac(1,mt)*b(1,nz) + ac(2,mt)*b(2,nz)

  IF((sf%sin_cos.EQ._SIN_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !a1s 
    DO m=0,m_max
      mm=REAL(m,wp)
      DO mIP=1,m_nyq
        ASSOCIATE(xm=>sf%X_IP(1,mIP))
        sf%base1D_IPthet(      mIP,1:2,1+m)  =(/    SIN(mm*xm),   -COS(mm*xm)/)
        sf%base1D_dthet_IPthet(mIP,1:2,1+m)  =(/ mm*COS(mm*xm), mm*SIN(mm*xm)/)
        END ASSOCIATE
      END DO
    END DO
  END IF
  IF((sf%sin_cos.EQ._COS_).OR.(sf%sin_cos.EQ._SINCOS_))THEN 
    i=sf%mTotal1D-(sf%mn_max(1)+1) !=offset, =0 if cos, =m_max+1 if sincos
    DO m=0,m_max
      mm=REAL(m,wp)
      DO mIP=1,m_nyq
        ASSOCIATE(xm=>sf%X_IP(1,mIP))
        sf%base1D_IPthet(      mIP,1:2,i+1+m)  =(/    COS(mm*xm),    SIN(mm*xm)/)
        sf%base1D_dthet_IPthet(mIP,1:2,i+1+m)  =(/-mm*SIN(mm*xm), mm*COS(mm*xm)/)
        END ASSOCIATE
      END DO
    END DO
  END IF

  DO n=-n_max,n_max  
    nn=REAL(n*nfp,wp)
    DO nIP=1,n_nyq
      ASSOCIATE(xn=>sf%X_IP(2,1+m_nyq*(nIP-1)))
      sf%base1D_IPzeta(      1:2,n,nIP)  =(/    COS(nn*xn),    SIN(nn*xn)/)
      sf%base1D_dzeta_IPzeta(1:2,n,nIP)  =(/-nn*SIN(nn*xn), nn*COS(nn*xn)/)
      END ASSOCIATE
    END DO
  END DO

  END ASSOCIATE !sf


  sf%initialized=.TRUE.
  IF(.NOT.test_called) THEN
    SWRITE(UNIT_stdOut,'(4X,A)')'... DONE'
    CALL fBase_test(sf)
  END IF

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
  ALLOCATE(sf%snorm_base(1:modes) )
  ALLOCATE(sf%base1D_IPthet(      1:sf%mn_nyq(1),1:2,1:sf%mTotal1D) )
  ALLOCATE(sf%base1D_dthet_IPthet(1:sf%mn_nyq(1),1:2,1:sf%mTotal1D) )
  ALLOCATE(sf%base1D_IPzeta(      1:2,-sf%mn_max(2):sf%mn_max(2),1:sf%mn_nyq(2)) )
  ALLOCATE(sf%base1D_dzeta_IPzeta(1:2,-sf%mn_max(2):sf%mn_max(2),1:sf%mn_nyq(2)) )
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
  SDEALLOCATE(sf%snorm_base)
  SDEALLOCATE(sf%base1D_IPthet)
  SDEALLOCATE(sf%base1D_dthet_IPthet)
  SDEALLOCATE(sf%base1D_IPzeta)
  SDEALLOCATE(sf%base1D_dzeta_IPzeta)
  
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
  CLASS(c_fBase), INTENT(IN   ) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=8) :: sin_cos
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPE IS(t_fBase)
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
!> compare sf with the input type fBase
!!
!===================================================================================================================================
SUBROUTINE fBase_compare( sf , tocompare,is_same, cond_out)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase),  INTENT(IN   ) :: sf !! self
  CLASS(c_fBase),  INTENT(IN   ) :: tocompare
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  LOGICAL,OPTIONAL,INTENT(  OUT) :: is_same
  LOGICAL,OPTIONAL,INTENT(  OUT) :: cond_out(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL  :: cond(5)
!===================================================================================================================================
  SELECT TYPE(tocompare); TYPE IS(t_fBase)
  IF(.NOT.tocompare%initialized) THEN
    CALL abort(__STAMP__, &
        "fBase_compare: tried to compare with non-initialized fBase!")
  END IF

  cond(1)= ALL( sf%mn_max(:)      .EQ.  tocompare%mn_max(:)       )
  cond(2)=    ( sf%nfp            .EQ.  tocompare%nfp             )
  cond(3)=    ( sf%modes          .EQ.  tocompare%modes           )
  cond(4)=    ( sf%sin_cos        .EQ.  tocompare%sin_cos         )
  cond(5)=    ( sf%exclude_mn_zero.EQV. tocompare%exclude_mn_zero ) 

  IF(PRESENT(is_same)) is_same=ALL(cond)
  IF(PRESENT(cond_out)) cond_out(1:5)=cond

  END SELECT !TYPE
END SUBROUTINE fBase_compare


!===================================================================================================================================
!> change data from oldBase to self. 
!! Forier modes are directly copied so, if new mode space is smaller, its like a Fourier cut-off.
!! if new modes do not match old ones, they are set to zero. 
!! Note that a change of nfp is not possible· as well as a change from sine to cosine
!!
!===================================================================================================================================
SUBROUTINE fBase_change_base( sf,old_fBase,iterDim,old_data,sf_data)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase),  INTENT(IN   ) :: sf !! self
  CLASS(c_fBase),  INTENT(IN   ) :: old_fBase       !! base of old_data
  INTEGER         ,INTENT(IN   ) :: iterDim        !! iterate on first or second dimension or old_data/sf_data
  REAL(wp)        ,INTENT(IN   ) :: old_data(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)        ,INTENT(  OUT) :: sf_data(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL             :: cond(5)
  INTEGER             :: iMode 
  INTEGER,ALLOCATABLE :: modeMapSin(:,:),modeMapCos(:,:)
!===================================================================================================================================
  SELECT TYPE(old_fBase); TYPE IS(t_fBase)
  IF(.NOT.old_fBase%initialized) THEN
    CALL abort(__STAMP__, &
        "fBase_change_base: tried to change base with non-initialized fBase!")
  END IF
  IF((iterDim.LT.1).OR.(iterDim.GT.2))THEN
    CALL abort(__STAMP__, &
        "fBase_change_base: iterDim can only be 1 or 2!")
  END IF
  IF(SIZE(old_data,iterDim).NE.SIZE(sf_data,iterDim)) THEN
    CALL abort(__STAMP__, &
        "fBase_change_base: iteration dimenion of old_data and sf_data have to be the same!")
  END IF
  IF(SIZE(old_data,3-iterDim).NE.old_fBase%modes) THEN
    CALL abort(__STAMP__, &
        "fBase_change_base: old_data size does not match old_fBase!")
  END IF
  IF(SIZE( sf_data,3-iterDim).NE.      sf%modes) THEN
    CALL abort(__STAMP__, &
        "fBase_change_base: sf_data size does not match sf fBase!")
  END IF

  CALL sf%compare(old_fBase,cond_out=cond(1:5))

  IF(ALL(cond))THEN
   !same base
   sf_data=old_data
  ELSE 
    !actually change base
    IF(.NOT.cond(2)) THEN !nfp
      CALL abort(__STAMP__, &
          "fBase_change_base: different nfp found, cannot change base!")
    END IF
    IF(.NOT.cond(4)) THEN !sin_cos /= sin_cos_old
      ! sin <-> cos : not ok
      ! cos <-> sin : not ok
      ! sin <-> sin_cos : ok
      ! cos <-> sin_cos : ok
      IF(.NOT.(ANY((/sf%sin_cos,old_fBase%sin_cos/).EQ._SINCOS_)))THEN
      CALL abort(__STAMP__, &
          "fBase_change_base: cannot change base between sine and cosine!")
      END IF
    END IF
    ASSOCIATE(mn_max    => old_fBase%mn_max   ,&
              nfp       => old_fBase%nfp      ,&
              Xmn       => old_fBase%Xmn      ,&
              sin_range => old_fBase%sin_range,&
              cos_range => old_fBase%cos_range ) 
    ALLOCATE(modeMapSin( 0:mn_max(1),-mn_max(2):mn_max(2)))
    ALLOCATE(modeMapCos( 0:mn_max(1),-mn_max(2):mn_max(2)))
    modeMapSin=-1
    DO iMode=sin_range(1)+1,sin_range(2)
      modeMapSin(Xmn(1,iMode),Xmn(2,iMode)/nfp)=iMode
    END DO
    modeMapCos=-1
    DO iMode=cos_range(1)+1,cos_range(2)
      modeMapCos(Xmn(1,iMode),Xmn(2,iMode)/nfp)=iMode 
    END DO
    END ASSOCIATE !old_fBase%...

    sf_data=0.0_wp
    IF((old_fBase%sin_range(2)-old_fBase%sin_range(1)).GT.0)THEN ! =_SIN_ / _SIN_COS_
      DO iMode=sf%sin_range(1)+1,sf%sin_range(2)
        IF(    sf%Xmn(1,iMode) .GT.old_fBase%mn_max(1))CYCLE ! remains zero
        IF(ABS(sf%Xmn(2,iMode)/sf%nfp).GT.old_fBase%mn_max(2))CYCLE ! remains zero
        SELECT CASE(iterDim)
        CASE(1)
          sf_data(:,iMode)=old_data(:,modeMapSin(sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp))
        CASE(2)
          sf_data(iMode,:)=old_data(modeMapSin(sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp),:)
        END SELECT
      END DO 
    END IF !old_fBase  no sine
    IF((old_fBase%cos_range(2)-old_fBase%cos_range(1)).GT.0)THEN ! =_COS_ / _SIN_COS_
      DO iMode=sf%cos_range(1)+1,sf%cos_range(2)
        IF(    sf%Xmn(1,iMode) .GT.old_fBase%mn_max(1))CYCLE !  m  > m_max_old, remains zero
        IF(ABS(sf%Xmn(2,iMode)/sf%nfp).GT.old_fBase%mn_max(2))CYCLE ! |n| > n_max_old, remains zero
        SELECT CASE(iterDim)
        CASE(1)
          sf_data(:,iMode)=old_data(:,modeMapCos(sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp))
        CASE(2)
          sf_data(iMode,:)=old_data(modeMapCos(sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp),:)
        END SELECT
      END DO 
    END IF !old_fBase  no sine
             
    DEALLOCATE(modeMapSin)
    DEALLOCATE(modeMapCos)
  END IF !same base
     
  END SELECT !TYPE
END SUBROUTINE fBase_change_base

!===================================================================================================================================
!> evaluate  all modes at specific given point
!!
!===================================================================================================================================
FUNCTION fBase_eval(sf,deriv,x) RESULT(base_x)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: x(2) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: base_x(sf%modes) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iMode
!===================================================================================================================================
  ASSOCIATE(sin_range=>sf%sin_range,cos_range=>sf%cos_range,Xmn=>sf%Xmn) 
  SELECT CASE(deriv)
  CASE(0)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_x(iMode)=                       SIN(REAL(Xmn(1,iMode),wp)*x(1)-REAL(Xmn(2,iMode),wp)*x(2))
    END DO !iMode                                                                               
    DO iMode=cos_range(1)+1,cos_range(2)                                                          
      base_x(iMode)=                       COS(REAL(Xmn(1,iMode),wp)*x(1)-REAL(Xmn(2,iMode),wp)*x(2))
    END DO !iMode                                                                                 
  CASE(DERIV_THET)                                                                                
    DO iMode=sin_range(1)+1,sin_range(2)                                                          
      base_x(iMode)= REAL(Xmn(1,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*x(1)-REAL(Xmn(2,iMode),wp)*x(2))
    END DO !iMode                                                                                 
    DO iMode=cos_range(1)+1,cos_range(2)                                                          
      base_x(iMode)=-REAL(Xmn(1,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*x(1)-REAL(Xmn(2,iMode),wp)*x(2))
    END DO !iMode                                                                                 
  CASE(DERIV_ZETA)                                                                                
    DO iMode=sin_range(1)+1,sin_range(2)                                                          
      base_x(iMode)=-REAL(Xmn(2,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*x(1)-REAL(Xmn(2,iMode),wp)*x(2))
    END DO !iMode                                                                                 
    DO iMode=cos_range(1)+1,cos_range(2)                                                          
      base_x(iMode)= REAL(Xmn(2,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*x(1)-REAL(Xmn(2,iMode),wp)*x(2))
    END DO !iMode
  CASE DEFAULT 
    CALL abort(__STAMP__, &
         "fbase_evalDOF_IP: derivative must be 0,DERIV_THET,DERIV_ZETA!")
  END SELECT
  END ASSOCIATE
END FUNCTION fBase_eval

!===================================================================================================================================
!> evaluate  all modes at all interpolation points 
!!
!===================================================================================================================================
FUNCTION fBase_evalDOF_x(sf,x,deriv,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  REAL(wp)      , INTENT(IN   ) :: x(2)   !! input coordinate theta,zeta in [0,2pi]^2
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                      :: base_x(1:sf%modes)
!===================================================================================================================================
IF(SIZE(DOFs,1).NE.sf%modes) CALL abort(__STAMP__, &
       'nDOF not correct when calling fBase_evalDOF_IP' )
  base_x=sf%eval(deriv,x)
  y=DOT_PRODUCT(base_x,DOFs(:))

END FUNCTION fBase_evalDOF_x

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
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y_IP(sf%mn_IP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
call perfon('evaldof_ip')
  IF(SIZE(DOFs,1).NE.sf%modes) CALL abort(__STAMP__, &
       'nDOF not correct when calling fBase_evalDOF_IP' )
  SELECT CASE(deriv)
  CASE(0)
    !y_IP=MATMUL(sf%base_IP(:,:),DOFs(:))
    __MATVEC_N(y_IP,sf%base_IP,DOFs)
  CASE(DERIV_THET)
    !y_IP=MATMUL(sf%base_dthet_IP(:,:),DOFs(:))
    __MATVEC_N(y_IP,sf%base_dthet_IP,DOFs)
  CASE(DERIV_ZETA)
    !y_IP=MATMUL(sf%base_dzeta_IP(:,:),DOFs(:))
    __MATVEC_N(y_IP,sf%base_dzeta_IP,DOFs)
  CASE DEFAULT 
    CALL abort(__STAMP__, &
         "fbase_evalDOF_IP: derivative must be 0,DERIV_THET,DERIV_ZETA!")
  END SELECT
call perfoff('evaldof_ip')
END FUNCTION fBase_evalDOF_IP

!===================================================================================================================================
!> evaluate  all modes at all interpolation points, making use of the tensor product:
!> y_ij = DOFs_mn * SIN(m*t_i - n*z_j ) => SIN(m*t_i) DOFs_mn COS(n*z_j) -COS(m*t_i) DOFs_mn SIN(n*z_j)
!> y_ij = DOFs_mn * COS(m*t_i - n*z_j ) => COS(m*t_i) DOFs_mn COS(n*z_j) +SIN(m*t_i) DOFs_mn SIN(n*z_j)
!>                                     => a1_im DOFs_mn b1_nj + a2_im DOFs_mn b2_nj
!> can be written as 2 SPECIAL MATMAT operations:
!> c(i,1,n)=a1(i,m) DOFs(m,n) , c(i,2,n) = a2(i,m) DOFs(m,n)  => c(i,d,n) = DOT_PROD(a(i,d,1:mmax),DOFs(1:mmax,n))
!> y(i,j) = c(i,1,n) b1(n,j) + c(i,2,n) b2(n,j) 
!>        = DOT_PROD(c(i,1:2,1:nmax),b(1:2,1:nmax,j)
!> the 1D ordering in y does not neead a reshape, y(i,j) => y(1:mn_IP), 1D array data can be kept, 
!> as it is passed (with its start adress) to DGEMM.
!!
!===================================================================================================================================
FUNCTION fBase_evalDOF_IP_tens(sf,deriv,DOFs) RESULT(y_IP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y_IP(sf%mn_IP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iMode,offset,mTotal,nTotal
  REAL(wp)                      :: Amn(1:sf%mTotal1D,-sf%mn_max(2):sf%mn_max(2))
  REAL(wp)                      :: Ctmp(1:sf%mn_nyq(1),1:2,-sf%mn_max(2):sf%mn_max(2))
!===================================================================================================================================
call perfon('evaldof_ip_tens')
  IF(SIZE(DOFs,1).NE.sf%modes) CALL abort(__STAMP__, &
         'nDOF not correct when calling fBase_evalDOF_IP' )

  offset=sf%mTotal1D-(sf%mn_max(1)+1) !=0 if sin or cos, =sf%mn_max(1)+1 if sin+cos
  !initialize non existing modes to zero
  Amn(1,-sf%mn_max(2):0)=0.0_wp 
  IF(offset.GT.0) Amn(offset+1,-sf%mn_max(2):0)=0.0_wp 

  !copy DOFs to  (0:m_max , -n_max:n_max) matrix, careful: Xmn(2,:) has nfp factor!
  DO iMode=sf%sin_range(1)+1,sf%sin_range(2)
    Amn(1+sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp)=DOFs(iMode)
  END DO
  DO iMode=sf%cos_range(1)+1,sf%cos_range(2)
    Amn(offset+1+sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp)=DOFs(iMode)
  END DO

  mTotal=  sf%mTotal1D
  nTotal=2*sf%mn_max(2)+1 !-n_max:n_nax

  SELECT CASE(deriv)
  CASE(0)
!    DO n=-sf%mn_max(2),sf%mn_max(2)
!      DO i=1,sf%mn_nyq(1)
!        Ctmp(i,1,n)=SUM(sf%base1D_IPthet(i,1,:)*Amn(:,n))
!        Ctmp(i,2,n)=SUM(sf%base1D_IPthet(i,2,:)*Amn(:,n))
!      END DO !i
!    END DO !n
!    k=0
!    DO j=1,sf%mn_nyq(2)
!      DO i=1,sf%mn_nyq(1)
!        k=k+1
!        y_IP(k)=SUM(Ctmp(i,1:2,:)*sf%base1D_IPzeta(1:2,:,j))
!      END DO !i
!    END DO !j
    __DGEMM_NN(Ctmp,2*sf%mn_nyq(1),  mTotal,sf%base1D_IPthet,  mTotal,      nTotal,Amn)
    __DGEMM_NN(y_IP,  sf%mn_nyq(1),2*nTotal,            Ctmp,2*nTotal,sf%mn_nyq(2),sf%base1D_IPzeta) 
  CASE(DERIV_THET)
    __DGEMM_NN(Ctmp,2*sf%mn_nyq(1),  mTotal,sf%base1D_dthet_IPthet,  mTotal,      nTotal,Amn)
    __DGEMM_NN(y_IP,  sf%mn_nyq(1),2*nTotal,                  Ctmp,2*nTotal,sf%mn_nyq(2),sf%base1D_IPzeta) 
  CASE(DERIV_ZETA)
    __DGEMM_NN(Ctmp,2*sf%mn_nyq(1),  mTotal,sf%base1D_IPthet,  mTotal,      nTotal,Amn)
    __DGEMM_NN(y_IP,  sf%mn_nyq(1),2*nTotal,            Ctmp,2*nTotal,sf%mn_nyq(2),sf%base1D_dzeta_IPzeta) 
  CASE DEFAULT 
    CALL abort(__STAMP__, &
         "fbase_evalDOF_IP_tens: derivative must be 0,DERIV_THET,DERIV_ZETA!")
  END SELECT
call perfoff('evaldof_ip_tens')
END FUNCTION fBase_evalDOF_IP_tens


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
  REAL(wp)      , INTENT(IN   ) :: g_IP(:)  !!  interpolation values at theta_IP zeta_IP positions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: DOFs(1:sf%modes)  !! projection to fourier base
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(SIZE(g_IP,1).NE.sf%mn_IP) CALL abort(__STAMP__, &
       'nDOF not correct when calling fBase_initDOF' )
  __AMATVEC_T(DOFs,(sf%d_thet*sf%d_zeta),sf%base_IP,g_IP)
  DOFs(:)=sf%snorm_base(:)*DOFs(:)
END FUNCTION fBase_initDOF

!===================================================================================================================================
!> test fBase variable
!!
!===================================================================================================================================
SUBROUTINE fBase_test( sf)
! MODULES
USE MODgvec_GLobals, ONLY: testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,iMode,jMode,ncoszero,nsinzero,i_mn
  REAL(wp)           :: checkreal,refreal
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
  REAL(wp)           :: dofs(1:sf%modes)
  REAL(wp)           :: g_IP(1:sf%mn_IP)
  TYPE(t_fbase)      :: testfBase
  LOGICAL            :: check(5)
  REAL(wp),ALLOCATABLE :: oldDOF(:,:),newDOF(:,:)
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
            , Xmn        => sf%Xmn        &
            )
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    checkreal =SUM(sf%x_IP(1,:)*sf%x_IP(2,:))*sf%d_thet*sf%d_zeta
    refreal   =(0.5_wp*(TWOPI)**2)*REAL(nfp,wp)*(0.5_wp*(TWOPI/REAL(nfp,wp))**2)

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,2(A,I4),2(A,E11.3))') &
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
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
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
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,(A,I4),2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
       '\n =>  should be  0 : nsinzero = ', nsinzero,  &
       '\n =>  should be ', refreal,' : nfp*int(int(base(imode)*base(imode), 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

    iTest=104 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    checkreal=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      DO jMode=sin_range(1)+1,sin_range(2)
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dthet_IP(:,jMode))))
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dzeta_IP(:,jMode))))
      END DO
    END DO
    DO iMode=cos_range(1)+1,cos_range(2)
      DO jMode=cos_range(1)+1,cos_range(2)
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dthet_IP(:,jMode))))
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dzeta_IP(:,jMode))))
      END DO
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : nfp*int(int(base(imode)*base_dthet/dzeta(jmode), 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

    !get new fbase and check compare
    iTest=111 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testfBase%init(sf%mn_max,sf%mn_nyq,sf%nfp,sin_cos_map(sf%sin_cos),sf%exclude_mn_zero)
    CALL testfBase%compare(sf,is_same=check(1))
    CALL testfBase%free()
    IF(.NOT.check(1))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,A)') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be true' 
    END IF !TEST

    !get new fbase and check compare
    iTest=112 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testfBase%init(sf%mn_max,sf%mn_nyq,sf%nfp+1,sin_cos_map(sf%sin_cos),(.NOT.sf%exclude_mn_zero))
    CALL testfBase%compare(sf,cond_out=check(1:5))
    CALL testfBase%free()
    IF(ALL(check))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,A)') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be false' 
    END IF !TEST

    !get new fbase and check compare
    iTest=113 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testfBase%init(2*sf%mn_max,2*sf%mn_nyq,sf%nfp,sin_cos_map(sf%sin_cos),sf%exclude_mn_zero)
    CALL testfBase%compare(sf,cond_out=check)
    CALL testfBase%free()
    IF(ALL(check))THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,A)') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be false' 
    END IF !TEST

    !get new fbase and check change_base execution  (can fail by abort)
    iTest=121 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testfBase%init(2*sf%mn_max,2*sf%mn_nyq,sf%nfp,sin_cos_map(sf%sin_cos),sf%exclude_mn_zero)
    ALLOCATE(oldDOF(1:sf%modes,2),newDOF(1:testfBase%modes,2))
    oldDOF(:,1)=1.1_wp
    oldDOF(:,2)=2.2_wp
    
    CALL testfBase%change_base(sf,2,oldDOF,newDOF)
    checkreal=SUM(newDOF)
    refreal  =SUM(oldDOF)
    CALL testfBase%free()
    DEALLOCATE(oldDOF,newDOF)
    IF(testdbg.OR.(.NOT.( (ABS(checkreal-refreal).LT. realtol) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
       '\n =>  should be ', refreal,' : ', checkreal
    END IF !TEST

    !get new fbase and check change_base execution only (can only fail by abort)
    iTest=122 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testfBase%init((/sf%mn_max(1)/2,sf%mn_max(2)/),(/sf%mn_nyq(1)/2,sf%mn_nyq(2)/),sf%nfp,sin_cos_map(sf%sin_cos),.TRUE.)
    ALLOCATE(oldDOF(3,1:sf%modes),newDOF(3,1:testfBase%modes))
    oldDOF(1,:)=-1.1_wp
    oldDOF(2,:)=-2.2_wp
    oldDOF(3,:)=-3.3_wp
    
    CALL testfBase%change_base(sf,1,oldDOF,newDOF)
    checkreal=SUM(newDOF)/REAL(testfBase%modes,wp)
    refreal  =-6.6_wp
    CALL testfBase%free()
    DEALLOCATE(oldDOF,newDOF)
    IF(testdbg.OR.(.NOT.( (ABS(checkreal-refreal).LT. realtol) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
       '\n =>  should be ', refreal,' : ', checkreal
    END IF !TEST

  END IF !testlevel <=1
  IF (testlevel .GE.2)THEN

    iTest=201 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*SIN(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode 
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*COS(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode 
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(0,dofs))) 
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF(dofs)|) ', checkreal
    END IF !TEST

    iTest=202 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    !use g_IP /dofs from test 201

    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal, ABS(g_IP(i_mn)-sf%evalDOF_x(sf%X_IP(:,i_mn),0,dofs))) 
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(x,(:),dofs)|) ', checkreal
STOP
    END IF !TEST

    iTest=203 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    !use g_IP /dofs from test 201

    checkreal=MAXVAL(ABS(sf%initDOF(g_IP)-dofs))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|initDOF(g_IP)-dofs|) ', checkreal
    END IF !TEST

    iTest=204 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL( Xmn(1,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode 
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-Xmn(1,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode 
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_THET,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dthet(dofs)|) ', checkreal
    END IF !TEST

    iTest=205 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    ! use g_IP and dofs from test 204
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal,ABS(g_IP(i_mn)-sf%evalDOF_x(sf%x_IP(:,i_mn),DERIV_THET,dofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dthet,x(:),dofs)|) ', checkreal
STOP
    END IF !TEST

    iTest=206 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-Xmn(2,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode 
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL( Xmn(2,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode 
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_ZETA,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dthet(dofs)|) ', checkreal
    END IF !TEST

    iTest=207 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    
    !use g_IP / dofs from test 206
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal,ABS(g_IP(i_mn)-sf%evalDOF_x(sf%x_IP(:,i_mn),DERIV_ZETA,dofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dthet,x(:),dofs)|) ', checkreal
STOP
    END IF !TEST
  END IF !testlevel <=2
  END ASSOCIATE !sf
  test_called=.FALSE.   

END SUBROUTINE fBase_test


END MODULE MODgvec_fBase

