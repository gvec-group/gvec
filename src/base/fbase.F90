!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
USE MODgvec_Globals                  ,ONLY: TWOPI,wp,Unit_stdOut,abort,MPIRoot
IMPLICIT NONE

PRIVATE
PUBLIC t_fbase,fbase_new,sin_cos_map

TYPE :: t_fBase
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
  INTEGER              :: modes       !! total(global) number of modes in basis (depends if only sin/cos or sin & cos are used)
  INTEGER              :: modes_str, modes_end   !! local range of modes, when distributed over MPI subdomains
  INTEGER,ALLOCATABLE  :: offset_modes(:)        !! allocated (0:nRanks), gives range on each rank:
                                                 !!   modes_str:modes_end=offset_modes(myRank)+1:offset_modes(myRank+1)
  INTEGER,ALLOCATABLE  :: whichRank(:)           !! know the MPI rank for each mode
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER              :: sin_range(2)        !! sin_range(1)+1:sin_range(2) is range with sine bases
  INTEGER              :: cos_range(2)        !! sin_range(1)+1:sin_range(2) is range with cosine bases
  INTEGER              :: mn_zero_mode        !! points to m=0,n=0 mode in mode array (1:mn_modes) (only one can exist for cosine, else =-1)
  REAL(wp)             :: d_thet              !! integration weight in theta direction: =2pi/mn_nyq(1)
  REAL(wp)             :: d_zeta              !! integration weight in zeta direction : =nfp*(2pi/nfp)/mn_nyq(2)=2*pi/mn_nyq(2)
  INTEGER,ALLOCATABLE  :: Xmn(:,:)            !! mode number (m,n*nfp) for each iMode=1,modes, size(2,modes)
  INTEGER,ALLOCATABLE  :: zero_odd_even(:)    !! =0 for m=n=0 mode, =1 for m= odd mode, =2 for m=even mode size(modes)
  REAL(wp),ALLOCATABLE :: x_IP(:,:)           !! (theta,zeta)position of interpolation points theta [0,2pi]x[0,2pi/nfp]size(2,mn_IP)
  REAL(wp),ALLOCATABLE :: thet_IP(:)           !! 1d theta position of interpolation points theta [0,2pi] size(mn_nyq(1))
  REAL(wp),ALLOCATABLE :: zeta_IP(:)           !! 1d zeta position of interpolation points theta [0,2pi/nfp] size(mn_nyq(2))
  REAL(wp),ALLOCATABLE :: base_IP(:,:)        !! basis functions,                 size(1:mn_IP,1:modes)
  REAL(wp),ALLOCATABLE :: base_dthet_IP(:,:)  !! dthet derivative of basis functions, (1:mn_IP,1:modes)
  REAL(wp),ALLOCATABLE :: base_dzeta_IP(:,:)  !! dzeta derivative of basis functions, (1:mn_IP,1:modes)

  REAL(wp),ALLOCATABLE :: snorm_base(:)       !! 1/norm of each basis function, size(1:modes), norm=int_0^2pi int_0^pi (base_mn(thet,zeta))^2 dthet dzeta
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
  PROCEDURE :: eval_xn          => fBase_eval_xn
  PROCEDURE :: evalDOF_x        => fBase_evalDOF_x
  PROCEDURE :: evalDOF_xn       => fBase_evalDOF_xn
  PROCEDURE :: evalDOF_xn_tens  => fBase_evalDOF_xn_tens
  ! PROCEDURE :: evalDOF_IP       => fBase_evalDOF_IP !use _tens instead!
  PROCEDURE :: evalDOF_IP       => fBase_evalDOF_IP_tens
  !  PROCEDURE :: projectIPtoDOF   => fBase_projectIPtoDOF
  PROCEDURE :: projectIPtoDOF   => fBase_projectIPtoDOF_tens
  PROCEDURE :: projectxntoDOF   => fBase_projectxntoDOF
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
  TYPE(t_fBase), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(sf)
  __PERFON("fbase_new")
  CALL sf%init(mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)

  __PERFOFF("fbase_new")
END SUBROUTINE fBase_new

!===================================================================================================================================
!> initialize the type fBase maximum mode numbers, number of integration points, type of basis (sin/cos or sin and cos)
!!
!===================================================================================================================================
SUBROUTINE fBase_init( sf, mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
! MODULES
USE MODgvec_Globals, ONLY: myRank, nRanks
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
  INTEGER :: i,iMode,m,n,mIP,nIP,mn_excl,iRank
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
  IF(mn_nyq_in(1).LT.(2*mn_max_in(1)+1)) &
    CALL abort(__STAMP__, &
        "error in fBase: mn_nyq in theta should be >= 2*mn_max(1)+1!",mn_nyq_in(1),REAL(mn_max_in(1)),&
        TypeInfo="InvalidParameterError")
  IF(mn_nyq_in(2).LT.(2*mn_max_in(2)+1)) &
    CALL abort(__STAMP__, &
         "error in fBase: mn_nyq in zeta should be >= 2*mn_max(2)+1!",mn_nyq_in(2),REAL(mn_max_in(2)),&
        TypeInfo="InvalidParameterError")

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
         "error in fBase: sin/cos not correctly specified, should be either _SIN_, _COS_ or _SINCOS_",&
        TypeInfo="InvalidParameterError")
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
  ! modes_sin :: m=0: n=1...n_max , m=1...m_max: n=-n_max...n_max. REMARK: for sine, m=0,n=0 is automatically excluded
  ! modes_cos :: m=0: n=0...n_max , m=1...m_max: n=-n_max...n_max. mn_excl=True will exclude m=n=0
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

  !MPI DECOMPOSITION OF MODES INTO EQUAL MODE GROUPS, FOR PARALLELIZING COMPUTATIONS DONE SEPARATELY ON EACH MODE
  ALLOCATE(sf%offset_modes(0:nRanks))
  sf%offset_modes(0)=0
  DO iRank=0,nRanks-1
    sf%offset_modes(iRank+1)=(modes*(iRank+1))/nRanks
  END DO
  sf%modes_str = sf%offset_modes(myRank  )+1
  sf%modes_end = sf%offset_modes(myRank+1)
  IF(MPIroot)THEN
    iRank=COUNT((sf%offset_modes(1:nRanks)-sf%offset_modes(0:nRanks-1)) .EQ. 0)
    IF (iRank.GT.0) THEN
      WRITE(UNIT_stdout,'(5X,A,I4,A,I4,A)') &
      'WARNING: more MPI ranks than number of modes! ', iRank , ' ranks of ' ,nRanks, &
              ' have no modes associated. This only affects MPI scaling.'
    END IF
  END IF

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

  IF(iMode.NE.modes) CALL abort(__STAMP__,&
                                ' Problem in Xmn ')

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

  DO mIP=1,m_nyq
    sf%thet_IP(mIP)=(REAL(mIP,wp)-0.5_wp)*sf%d_thet
  END DO
  DO nIP=1,n_nyq
    sf%zeta_IP(nIP)=(REAL(nIP,wp)-0.5_wp)*sf%d_zeta
  END DO

  i=0
  DO nIP=1,n_nyq
    DO mIP=1,m_nyq
      i=i+1
      sf%x_IP(1,i)=sf%thet_IP(mIP)
      sf%x_IP(2,i)=sf%zeta_IP(nIP)
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
        ASSOCIATE(xm=>sf%thet_IP(mIP))
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
        ASSOCIATE(xm=>sf%thet_IP(mIP))
        sf%base1D_IPthet(      mIP,1:2,i+1+m)  =(/    COS(mm*xm),    SIN(mm*xm)/)
        sf%base1D_dthet_IPthet(mIP,1:2,i+1+m)  =(/-mm*SIN(mm*xm), mm*COS(mm*xm)/)
        END ASSOCIATE
      END DO
    END DO
  END IF

  DO n=-n_max,n_max
    nn=REAL(n*nfp,wp)
    DO nIP=1,n_nyq
      ASSOCIATE(xn=>sf%zeta_IP(nIP))
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
  ALLOCATE(sf%thet_IP(1:sf%mn_nyq(1)))
  ALLOCATE(sf%zeta_IP(1:sf%mn_nyq(2)))
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
  SDEALLOCATE(sf%thet_IP)
  SDEALLOCATE(sf%zeta_IP)
  SDEALLOCATE(sf%base_IP)
  SDEALLOCATE(sf%base_dthet_IP)
  SDEALLOCATE(sf%base_dzeta_IP)
  SDEALLOCATE(sf%snorm_base)
  SDEALLOCATE(sf%base1D_IPthet)
  SDEALLOCATE(sf%base1D_dthet_IPthet)
  SDEALLOCATE(sf%base1D_IPzeta)
  SDEALLOCATE(sf%base1D_dzeta_IPzeta)
  SDEALLOCATE(sf%offset_modes)

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
  CLASS(t_fBase), INTENT(IN   ) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_fBase), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=8) :: sin_cos
!===================================================================================================================================
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
  TYPE(t_fBase),  INTENT(IN   ) :: tocompare
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  LOGICAL,OPTIONAL,INTENT(  OUT) :: is_same
  LOGICAL,OPTIONAL,INTENT(  OUT) :: cond_out(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL  :: cond(5)
!===================================================================================================================================
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
  CLASS(t_fBase),  INTENT(IN   ) :: old_fBase       !! base of old_data
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
END SUBROUTINE fBase_change_base

!===================================================================================================================================
!> evaluate  all modes at specific given point in theta and zeta
!!
!===================================================================================================================================
FUNCTION fBase_eval(sf,deriv,x) RESULT(base_x)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: x(2)   !! theta,zeta point position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: base_x(sf%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
base_x =  RESHAPE(sf%eval_xn(deriv,1,x),(/sf%modes/))
END FUNCTION fbase_eval

!===================================================================================================================================
!> evaluate  all modes at a list of given points in theta and zeta
!!
!===================================================================================================================================
FUNCTION fBase_eval_xn(sf,deriv,np,xn) RESULT(base_xn)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf         !! self
  INTEGER       , INTENT(IN   ) :: deriv      !! =0: base, =2: dthet , =3: dzeta
  INTEGER       , INTENT(IN   ) :: np         !! number of points in xn
  REAL(wp)      , INTENT(IN   ) :: xn(2,1:np) !! theta,zeta point positions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: base_xn(1:np,sf%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iMode
!===================================================================================================================================
  ASSOCIATE(sin_range=>sf%sin_range,cos_range=>sf%cos_range,Xmn=>sf%Xmn)
  SELECT CASE(deriv)
  CASE(0)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_xn(:,iMode)=                       SIN(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      base_xn(:,iMode)=                       COS(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
  CASE(DERIV_THET)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_xn(:,iMode)= REAL(Xmn(1,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      base_xn(:,iMode)=-REAL(Xmn(1,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
  CASE(DERIV_ZETA)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_xn(:,iMode)=-REAL(Xmn(2,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      base_xn(:,iMode)= REAL(Xmn(2,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
  CASE(DERIV_THET_THET)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_xn(:,iMode)=-REAL(Xmn(1,iMode)**2,wp)*SIN(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      base_xn(:,iMode)=-REAL(Xmn(1,iMode)**2,wp)*COS(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
  CASE(DERIV_THET_ZETA)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_xn(:,iMode)= REAL(Xmn(1,iMode)*Xmn(2,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      base_xn(:,iMode)= REAL(Xmn(1,iMode)*Xmn(2,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
  CASE(DERIV_ZETA_ZETA)
    DO iMode=sin_range(1)+1,sin_range(2)
      base_xn(:,iMode)=-REAL(Xmn(2,iMode)**2,wp)*SIN(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      base_xn(:,iMode)=-REAL(Xmn(2,iMode)**2,wp)*COS(REAL(Xmn(1,iMode),wp)*xn(1,:)-REAL(Xmn(2,iMode),wp)*xn(2,:))
    END DO !iMode
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "fbase_evalDOF_IP: derivative must be 0,DERIV_THET,_ZETA,_THET_THET,_THET_ZETA,_ZETA_ZETA!")
  END SELECT
  END ASSOCIATE
END FUNCTION fbase_eval_xn

!===================================================================================================================================
!> evaluate special 1D base in theta direction (cos(m*t_i),sin(m*t_i)) or its derivative(s) on a given set of points
!! for tensor-product evaluation of 2D sin and cos base:
!!   sin(m*thet-n*zeta) = sin(m*thet)*cos(n*zeta)-cos(m*thet)*sin(n*zeta)
!!     == dot_product( (sin(m*thet),-cos(m*thet)) , (cos(n*zeta),sin(n*zeta)))
!!   cos(m*thet-n*zeta) = cos(m*thet)*cos(n*zeta)+sin(m*thet)*sin(n*zeta)
!!     == dot_product( (cos(m*thet), sin(m*thet)) , (cos(n*zeta),sin(n*zeta)))
!! so for the 1D base, mTotal1d depends on using sin/cos/sin+cos base.
!===================================================================================================================================
FUNCTION fBase_eval1d_thet(sf,deriv,nthet,thet) RESULT(base1d_thet)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf         !! self
  INTEGER       , INTENT(IN   ) :: deriv !! =0: base, =1: dthet , =2: dthet^2
  INTEGER       , INTENT(IN   ) :: nthet       !! number of points in theta
  REAL(wp)      , INTENT(IN   ) :: thet(1:nthet)   !! theta 1D point positions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: base1d_thet(1:nthet,1:2,1:sf%mTotal1D)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: m,m_max,i
  REAL(wp):: mm
!===================================================================================================================================
m_max=sf%mn_max(1)

SELECT CASE(deriv)
CASE(0)
  IF((sf%sin_cos.EQ._SIN_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !2D SINE
    DO m=0,m_max
      mm=REAL(m,wp)
      base1d_thet(:,1,1+m)  = SIN(mm*thet(:))
      base1d_thet(:,2,1+m)  =-COS(mm*thet(:))
    END DO
  END IF
  IF((sf%sin_cos.EQ._COS_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !2D cosine
    i=sf%mTotal1D-(sf%mn_max(1)+1) !=offset, =0 if cos, =m_max+1 if sincos
    DO m=0,m_max
      mm=REAL(m,wp)
      base1d_thet(:,1,i+1+m)  =COS(mm*thet(:))
      base1d_thet(:,2,i+1+m)  =SIN(mm*thet(:))
    END DO
  END IF
CASE(1)
  IF((sf%sin_cos.EQ._SIN_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !2D SINE
    DO m=0,m_max
      mm=REAL(m,wp)
      base1d_thet(:,1,1+m)  = mm*COS(mm*thet(:))
      base1d_thet(:,2,1+m)  = mm*SIN(mm*thet(:))
    END DO
  END IF
  IF((sf%sin_cos.EQ._COS_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !2D cosine
    i=sf%mTotal1D-(sf%mn_max(1)+1) !=offset, =0 if cos, =m_max+1 if sincos
    DO m=0,m_max
      mm=REAL(m,wp)
      base1d_thet(:,1,i+1+m)  =-mm*SIN(mm*thet(:))
      base1d_thet(:,2,i+1+m)  = mm*COS(mm*thet(:))
    END DO
  END IF
CASE(2)
  IF((sf%sin_cos.EQ._SIN_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !2D SINE
    DO m=0,m_max
      mm=REAL(m,wp)
      base1d_thet(:,1,1+m)  =-mm*mm*SIN(mm*thet(:))
      base1d_thet(:,2,1+m)  = mm*mm*COS(mm*thet(:))
    END DO
  END IF
  IF((sf%sin_cos.EQ._COS_).OR.(sf%sin_cos.EQ._SINCOS_))THEN !2D cosine
    i=sf%mTotal1D-(sf%mn_max(1)+1) !=offset, =0 if cos, =m_max+1 if sincos
    DO m=0,m_max
      mm=REAL(m,wp)
      base1d_thet(:,1,i+1+m)  =-mm*mm*COS(mm*thet(:))
      base1d_thet(:,2,i+1+m)  =-mm*mm*SIN(mm*thet(:))
    END DO
  END IF
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "fBase_eval1d_thet: derivative must be 0,1,2 !")
  END SELECT
END FUNCTION fBase_eval1d_thet



!===================================================================================================================================
!> evaluate special 1D base in zeta direction (cos(m*t_i),sin(m*t_i)) or its derivative(s) on a given set of points
!! for tensor-product evaluation of 2D sin and cos base:
!!   sin(m*thet-n*zeta) = sin(m*thet)*cos(n*zeta)-cos(m*thet)*sin(n*zeta)
!!     == dot_product( (sin(m*thet),-cos(m*thet)) , (cos(n*zeta),sin(n*zeta)))
!!   cos(m*thet-n*zeta) = cos(m*thet)*cos(n*zeta)+sin(m*thet)*sin(n*zeta)
!!     == dot_product( (cos(m*thet), sin(m*thet)) , (cos(n*zeta),sin(n*zeta)))
!! so for the 1D base, nTotal1d is always 2*n_max+1
!===================================================================================================================================
FUNCTION fBase_eval1d_zeta(sf,deriv,nzeta,zeta) RESULT(base1d_zeta)
  ! MODULES
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_fBase), INTENT(IN   ) :: sf         !! self
    INTEGER       , INTENT(IN   ) :: deriv !! =0: base, =1: dzeta , =2: dzeta^2
    INTEGER       , INTENT(IN   ) :: nzeta       !! number of points in zeta
    REAL(wp)      , INTENT(IN   ) :: zeta(1:nzeta)   !! zeta 1D point positions
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)                      :: base1d_zeta(1:2,-sf%mn_max(2):sf%mn_max(2),1:nzeta)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    INTEGER :: n,n_max,nfp
    REAL(wp):: nn
  !===================================================================================================================================
  n_max=sf%mn_max(2)
  nfp=sf%nfp

  SELECT CASE(deriv)
  CASE(0)
    DO n=-n_max,n_max
      nn=REAL(n*nfp,wp)
      base1D_zeta(      1,n,:)  = COS(nn*zeta(:))
      base1D_zeta(      2,n,:)  = SIN(nn*zeta(:))
    END DO
  CASE(1)  !
    DO n=-n_max,n_max
      nn=REAL(n*nfp,wp)
      base1D_zeta(      1,n,:)  = -nn*SIN(nn*zeta(:))
      base1D_zeta(      2,n,:)  =  nn*COS(nn*zeta(:))
    END DO
  CASE(2)
    DO n=-n_max,n_max
      nn=REAL(n*nfp,wp)
      base1D_zeta(      1,n,:)  = -nn*nn*COS(nn*zeta(:))
      base1D_zeta(      2,n,:)  = -nn*nn*SIN(nn*zeta(:))
    END DO
  CASE DEFAULT
    CALL abort(__STAMP__, &
           "fBase_eval1d_zeta: derivative must be 0,1,2 !")
    END SELECT
  END FUNCTION fBase_eval1d_zeta

!===================================================================================================================================
!> evaluate  all modes at a given interpolation point
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
       'nDOF not correct when calling fBase_evalDOF_x' )
  base_x=sf%eval(deriv,x)
  y=DOT_PRODUCT(base_x,DOFs(:))

END FUNCTION fBase_evalDOF_x

!===================================================================================================================================
!> evaluate  all modes at a list of interpolation points
!!
!===================================================================================================================================
FUNCTION fBase_evalDOF_xn(sf,np,xn,deriv,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: np     !! number of points to be evaluated
  REAL(wp)      , INTENT(IN   ) :: xn(2,1:np)   !! input coordinate theta,zeta in [0,2pi]^2
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y(1:np)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                      :: base_xn(1:np,1:sf%modes)
!===================================================================================================================================
IF(SIZE(DOFs,1).NE.sf%modes) CALL abort(__STAMP__, &
       'nDOF not correct when calling fBase_evalDOF_x' )
  base_xn=sf%eval_xn(deriv,np,xn)
  __MATVEC_N(y,base_xn,DOFs)

END FUNCTION fBase_evalDOF_xn

!===================================================================================================================================
!> evaluate  all modes on a tensor-produc grid (t_i,z_j), making use of the tensor product in the fourier series:
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
FUNCTION fBase_evalDOF_xn_tens(sf,nthet,nzeta,thet,zeta,deriv,DOFs) RESULT(y)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  INTEGER       , INTENT(IN   ) :: nthet  !! number of points in theta
  INTEGER       , INTENT(IN   ) :: nzeta  !! number of points in zeta
  REAL(wp)      , INTENT(IN   ) :: thet(1:nthet) !! theta positions
  REAL(wp)      , INTENT(IN   ) :: zeta(1:nzeta) !! zeta positions
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y(1:nthet*nzeta)   !! DOFS evaluated on tensor-product grid,
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iMode,offset,mTotal,nTotal
  REAL(wp)                      :: Amn(1:sf%mTotal1D,-sf%mn_max(2):sf%mn_max(2))
  REAL(wp)                      :: Ctmp(1:nthet,1:2,-sf%mn_max(2):sf%mn_max(2))
  REAL(wp)                      :: base1D_thet(1:nthet,1:2,1:sf%mTotal1D)
  REAL(wp)                      :: base1D_zeta(1:2,-sf%mn_max(2):sf%mn_max(2),1:nzeta)
!===================================================================================================================================
  IF(SIZE(DOFs,1).NE.sf%modes) CALL abort(__STAMP__, &
         'nDOF not correct when calling fBase_evalDOF_IP_tens' )

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
    base1d_thet=fBase_eval1d_thet(sf,0,nthet,thet)
    base1d_zeta=fBase_eval1d_zeta(sf,0,nzeta,zeta)
  CASE(DERIV_THET)
    base1d_thet=fBase_eval1d_thet(sf,1,nthet,thet)
    base1d_zeta=fBase_eval1d_zeta(sf,0,nzeta,zeta)
  CASE(DERIV_ZETA)
    base1d_thet=fBase_eval1d_thet(sf,0,nthet,thet)
    base1d_zeta=fBase_eval1d_zeta(sf,1,nzeta,zeta)
  CASE(DERIV_THET_THET)
    base1d_thet=fBase_eval1d_thet(sf,2,nthet,thet)
    base1d_zeta=fBase_eval1d_zeta(sf,0,nzeta,zeta)
  CASE(DERIV_THET_ZETA)
    base1d_thet=fBase_eval1d_thet(sf,1,nthet,thet)
    base1d_zeta=fBase_eval1d_zeta(sf,1,nzeta,zeta)
  CASE(DERIV_ZETA_ZETA)
    base1d_thet=fBase_eval1d_thet(sf,0,nthet,thet)
    base1d_zeta=fBase_eval1d_zeta(sf,2,nzeta,zeta)
  CASE DEFAULT  !for other derivatives, resort to not precomputed/ explicit computation:
    CALL abort(__STAMP__, &
         "fbase_evalDOF_xn_tens: derivative must be 0,DERIV_THET,_ZETA,_THET_THET,_THET_ZETA,_ZETA_ZETA!")
  END SELECT
  __DGEMM_NN(Ctmp,2*nthet,  mTotal,base1D_thet,  mTotal, nTotal,Amn)
  __DGEMM_NN(y   ,  nthet,2*nTotal,       Ctmp,2*nTotal, nzeta ,base1D_zeta)
END FUNCTION fBase_evalDOF_xn_tens

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
  CASE DEFAULT  !for other derivatives, resort to not precomputed/ explicit computation:
     y_IP = sf%evalDOF_xn(sf%mn_IP,sf%x_IP,deriv,DOFs)
  END SELECT
END FUNCTION fBase_evalDOF_IP

!===================================================================================================================================
!> project from interpolation points to all modes
!!  DOFs = add*DOFs+ fac *MATMUL(base_IP_DOF,y_IP)
!===================================================================================================================================
SUBROUTINE fBase_projectIPtoDOF(sf,add,factor,deriv,y_IP,DOFs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  LOGICAL       , INTENT(IN   ) :: add    !! =F initialize DOFs , =T add to DOFs
  REAL(wp)      , INTENT(IN   ) :: factor !! scale result by factor, before adding to DOFs (should be =1.0_wp if not needed)
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: y_IP(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)      , INTENT(INOUT) :: DOFs(1:sf%modes)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                      :: radd
!===================================================================================================================================
  IF(SIZE(y_IP,1).NE.sf%mn_IP) CALL abort(__STAMP__, &
       'y_IP not correct when calling fBase_projectIPtoDOF' )
  radd=MERGE(1.0_wp,0.0_wp,add)
  SELECT CASE(deriv)
  CASE(0)
    __PAMATVEC_T(radd,DOFs,factor,sf%base_IP,y_IP)
  CASE(DERIV_THET)
    __PAMATVEC_T(radd,DOFs,factor,sf%base_dthet_IP,y_IP)
  CASE(DERIV_ZETA)
    __PAMATVEC_T(radd,DOFs,factor,sf%base_dzeta_IP,y_IP)
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "fbase_projectIPtoDOF: derivative must be 0,DERIV_THET,DERIV_ZETA!")
  END SELECT
END SUBROUTINE fBase_projectIPtoDOF

!===================================================================================================================================
!> project from any 2D set of interpolation points, at tensor-product of (theta,zeta) positions given by "xn", to all modes
!!  DOFs = add*DOFs+ fac *MATMUL(base_xn,yn)
!===================================================================================================================================
SUBROUTINE fBase_projectxntoDOF(sf,add,factor,deriv,np,xn,yn,DOFs)
  ! MODULES
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
    CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
    LOGICAL       , INTENT(IN   ) :: add    !! =F initialize DOFs , =T add to DOFs
    REAL(wp)      , INTENT(IN   ) :: factor !! scale result by factor, before adding to DOFs (should be =1.0_wp if not needed)
    INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
    INTEGER       , INTENT(IN   ) :: np     !! total number of 2D interpolation points
    REAL(wp)      , INTENT(IN   ) :: xn(2,1:np)  !!  (theta=1,zeta=2) position of tensor-product interpolation points, [0,2pi]x[0,2pi/nfp],size(2,mn_IP)
    REAL(wp)      , INTENT(IN   ) :: yn(1:np)  !! value at interpolation points
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
    REAL(wp)      , INTENT(INOUT) :: DOFs(1:sf%modes)  !! array of all modes
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
    REAL(wp)                      :: radd
    REAL(wp)                      :: base_xn(1:np,1:sf%modes)
  !===================================================================================================================================
    base_xn=sf%eval_xn(deriv,np,xn)
    radd=MERGE(1.0_wp,0.0_wp,add)
    __PAMATVEC_T(radd,DOFs,factor,base_xn,yn)
  END SUBROUTINE fBase_projectxntoDOF

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
  REAL(wp)      , INTENT(IN   ) :: DOFs(:)!! array of all modes (sf%modes)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: y_IP(sf%mn_IP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iMode,offset,mTotal,nTotal
  REAL(wp)                      :: Amn(1:sf%mTotal1D,-sf%mn_max(2):sf%mn_max(2))
  REAL(wp)                      :: Ctmp(1:sf%mn_nyq(1),1:2,-sf%mn_max(2):sf%mn_max(2))
!===================================================================================================================================
  IF(SIZE(DOFs,1).NE.sf%modes) CALL abort(__STAMP__, &
         'nDOF not correct when calling fBase_evalDOF_IP_tens' )

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
  CASE DEFAULT  !for other derivatives, resort to not precomputed/ explicit computation:
     y_IP = sf%evalDOF_xn(sf%mn_IP,sf%x_IP,deriv,DOFs)
  END SELECT
END FUNCTION fBase_evalDOF_IP_tens


!===================================================================================================================================
!> inverse of fBase_evalDOF_IP_tens
!!
!===================================================================================================================================
SUBROUTINE fBase_projectIPtoDOF_tens(sf,add,factor,deriv,y_IP,DOFs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf     !! self
  LOGICAL       , INTENT(IN   ) :: add    !! =F initialize DOFs , =T add to DOFs
  REAL(wp)      , INTENT(IN   ) :: factor !! scale result by factor, before adding to DOFs (should be =1.0_wp if not needed)
  INTEGER       , INTENT(IN   ) :: deriv  !! =0: base, =2: dthet , =3: dzeta
  REAL(wp)      , INTENT(IN   ) :: y_IP(:) !! point values (at sf%x_IP if x_IP_in not given)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)      , INTENT(INOUT) :: DOFs(1:sf%modes)  !! array of all modes
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iMode,offset,mTotal,nTotal
  REAL(wp)                      :: Amn(1:sf%mTotal1D,-sf%mn_max(2):sf%mn_max(2))
  REAL(wp)                      :: Ctmp(1:sf%mn_nyq(1),1:2,-sf%mn_max(2):sf%mn_max(2))
!===================================================================================================================================
  IF(SIZE(y_IP,1).NE.sf%mn_IP) CALL abort(__STAMP__, &
         'y_IP not correct when calling fBase_projectIPtoDOF_tens' )
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

    __DGEMM_NT(Ctmp,  sf%mn_nyq(1),sf%mn_nyq(2),y_IP,  2*nTotal,sf%mn_nyq(2),sf%base1D_IPzeta)
    __ADGEMM_TN(Amn,factor, 2*sf%mn_nyq(1),mTotal,sf%base1D_IPthet,  2*sf%mn_nyq(1),nTotal,Ctmp)

  CASE(DERIV_THET)
    __DGEMM_NT(Ctmp,  sf%mn_nyq(1),sf%mn_nyq(2),y_IP,  2*nTotal,sf%mn_nyq(2),sf%base1D_IPzeta)
    __ADGEMM_TN(Amn,factor, 2*sf%mn_nyq(1),mTotal,sf%base1D_dthet_IPthet,  2*sf%mn_nyq(1),nTotal,Ctmp)
  CASE(DERIV_ZETA)
    __DGEMM_NT(Ctmp,  sf%mn_nyq(1),sf%mn_nyq(2),y_IP,  2*nTotal,sf%mn_nyq(2),sf%base1D_dzeta_IPzeta)
    __ADGEMM_TN(Amn,factor, 2*sf%mn_nyq(1),mTotal,sf%base1D_IPthet,  2*sf%mn_nyq(1),nTotal,Ctmp)
  CASE DEFAULT
    CALL abort(__STAMP__, &
         "fbase_evalDOF_IP_tens: derivative must be 0,DERIV_THET,DERIV_ZETA!")
  END SELECT

  offset=sf%mTotal1D-(sf%mn_max(1)+1) !=0 if sin or cos, =sf%mn_max(1)+1 if sin+cos
  !copy modes back
  IF(add)THEN
    DO iMode=sf%sin_range(1)+1,sf%sin_range(2)
      DOFs(iMode)=DOFs(iMode)+Amn(1+sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp)
    END DO
    DO iMode=sf%cos_range(1)+1,sf%cos_range(2)
      DOFs(iMode)=DOFs(iMode)+Amn(offset+1+sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp)
    END DO
  ELSE
    DO iMode=sf%sin_range(1)+1,sf%sin_range(2)
      DOFs(iMode)=Amn(1+sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp)
    END DO
    DO iMode=sf%cos_range(1)+1,sf%cos_range(2)
      DOFs(iMode)=Amn(offset+1+sf%Xmn(1,iMode),sf%Xmn(2,iMode)/sf%nfp)
    END DO
  END IF !add
END SUBROUTINE fBase_projectIPtoDOF_tens

!===================================================================================================================================
!>  take values interpolated at sf%s_IP positions and project onto fourier basis by integration
!!
!===================================================================================================================================
FUNCTION fBase_initDOF( sf , g_IP,thet_zeta_start) RESULT(DOFs)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_fBase), INTENT(IN   ) :: sf    !! self
  REAL(wp)      , INTENT(IN   ) :: g_IP(:)  !!  interpolation values at theta_IP zeta_IP positions
  REAL(wp),INTENT(IN),OPTIONAL :: thet_zeta_start(2) !theta,zeta value of first point (points must remain equidistant and of size mn_nyq(1),mn_nyq(2))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                      :: DOFs(1:sf%modes)  !! projection to fourier base
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)                      :: x_IP_shift(2,sf%mn_IP)
!===================================================================================================================================
  IF(SIZE(g_IP,1).NE.sf%mn_IP) CALL abort(__STAMP__, &
       'nDOF not correct when calling fBase_initDOF' )
  IF(.NOT.(PRESENT(thet_zeta_start)))THEN
    CALL sf%projectIPtoDOF(.FALSE.,(sf%d_thet*sf%d_zeta),0,g_IP,DOFs)
  ELSE
    x_IP_shift(1,:)=sf%x_IP(1,:)-sf%x_IP(1,1)+thet_zeta_start(1)
    x_IP_shift(2,:)=sf%x_IP(2,:)-sf%x_IP(2,1)+thet_zeta_start(2)
    CALL sf%projectxntoDOF(.FALSE.,(sf%d_thet*sf%d_zeta),0,sf%mn_IP,x_IP_shift,g_IP,DOFs)
  END IF
  DOFs(:)=sf%snorm_base(:)*DOFs(:)  !normalize with inverse mass matrix diagonal
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
  REAL(wp)           :: dofs(1:sf%modes),tmpdofs(1:sf%modes),dangle(2)
  REAL(wp)           :: g_IP(1:sf%mn_IP)
  TYPE(t_fbase)      :: testfBase
  LOGICAL            :: check(5)
  REAL(wp),ALLOCATABLE :: oldDOF(:,:),newDOF(:,:)
!===================================================================================================================================
  test_called=.TRUE. !avoid infinite loop if init is called here
  IF(testlevel.LE.0) RETURN
  IF(.NOT.MPIroot) RETURN
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


    ! check off-diagonals of mass matrix =0
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
      '\n =>  should be ', refreal,' : OFF-DIAGONALS of mass matrix 0=:int(int(base(imode)*base(jmode), 0, 2pi),0,2pi/nfp)= ', checkreal
    END IF !TEST

    ! check off-diagonals of mass matrix =0
    iTest=1021 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    checkreal=0.0_wp
    DO iMode=1,modes
      !DIAGONAL
      checkreal=MAX(checkreal,ABS(1.0_wp-sf%snorm_base(iMode)*(sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_IP(:,iMode))))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : DIAGONAL OF MASS MATRIX 0=:1-snorm(iMode)*int(int(base(imode)*base(imode), 0, 2pi),0,2pi/nfp)= ', checkreal
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

    !test mass matrix of base
    iTest=104 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    checkreal=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      DO jMode=sin_range(1)+1,sin_range(2)
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dthet_IP(:,jMode)))/REAL(1+ABS(sf%Xmn(1,jmode)),wp))
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dzeta_IP(:,jMode)))/REAL(1+ABS(sf%Xmn(2,jmode)),wp))
      END DO
    END DO
    DO iMode=cos_range(1)+1,cos_range(2)
      DO jMode=cos_range(1)+1,cos_range(2)
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dthet_IP(:,jMode)))/REAL(1+ABS(sf%Xmn(1,jmode)),wp))
          checkreal=MAX(checkreal,ABS((sf%d_thet*sf%d_zeta)*SUM(sf%base_IP(:,iMode)*sf%base_dzeta_IP(:,jMode)))/REAL(1+ABS(sf%Xmn(2,jmode)),wp))
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

    IF(sf%mn_max(1).GT.1)THEN
    !get new fbase and check change_base execution only (can only fail by abort)
    iTest=122 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    CALL testfBase%init((/sf%mn_max(1)/2,sf%mn_max(2)/),(/sf%mn_nyq(1)/2+1,sf%mn_nyq(2)/),sf%nfp,sin_cos_map(sf%sin_cos),.TRUE.)
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
    END IF !sf%mn_max>1


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


    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_xn_tens(sf%mn_nyq(1),sf%mn_nyq(2),sf%X_IP(1,1:sf%mn_nyq(1)),sf%X_IP(2,1:PRODUCT(sf%mn_nyq(1:2)):sf%mn_nyq(1)),0,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_xn_tens(dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

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
    END IF !TEST

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP /dofs from test 201
    tmpdofs=sf%initDOF(g_IP)
    checkreal=MAXVAL(ABS(tmpdofs-dofs))
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

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP /dofs from test 201
    tmpdofs=sf%initDOF(g_IP,thet_zeta_start=(/sf%thet_IP(1),sf%zeta_IP(1)/))
    checkreal=MAXVAL(ABS(tmpdofs-dofs))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|initDOF(g_IP)-initDOF(g_IP,x_IP)|) ', checkreal
    END IF !TEST

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP  from test 201
    IF(sin_cos.EQ.3)THEN
      dangle=(/0.333_wp,-0.222_wp/)
    ELSE
      dangle=(/TWOPI,-2*TWOPI/)
    END IF
    tmpdofs=sf%initDOF(g_IP,thet_zeta_start=(/sf%x_IP(1,1),sf%x_IP(2,1)/)+dangle)
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal, ABS(g_IP(i_mn)-sf%evalDOF_x((sf%X_IP(:,i_mn)+dangle),0,tmpdofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(x+delta,initdof(g_IP,xIP+delta)|)', checkreal
    END IF !TEST

  END IF !testlevel <=1
  IF (testlevel .GE.2)THEN

    iTest=2031 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP /dofs from test 201
    tmpdofs=sf%initDOF(g_IP,thet_zeta_start=(/sf%thet_IP(1),sf%zeta_IP(1)/))
    checkreal=MAXVAL(ABS(tmpdofs-dofs))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|initDOF(g_IP)-initDOF(g_IP,x_IP)|) ', checkreal
    END IF !TEST

    iTest=2032 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP  from test 201
    IF(sin_cos.EQ.3)THEN
      dangle=(/0.333_wp,-0.222_wp/)
    ELSE
      dangle=(/TWOPI,-2*TWOPI/)
    END IF
    tmpdofs=sf%initDOF(g_IP,thet_zeta_start=(/sf%x_IP(1,1),sf%x_IP(2,1)/)+dangle)
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal, ABS(g_IP(i_mn)-sf%evalDOF_x((sf%X_IP(:,i_mn)+dangle),0,tmpdofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(x+delta,initdof(g_IP,xIP+delta)|)', checkreal
    END IF !TEST

  END IF !testlevel <=1
  IF (testlevel .GE.2)THEN

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

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

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_xn_tens(sf%mn_nyq(1),sf%mn_nyq(2),sf%X_IP(1,1:sf%mn_nyq(1)),sf%X_IP(2,1:PRODUCT(sf%mn_nyq(1:2)):sf%mn_nyq(1)),DERIV_THET,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_xn_tens_dthet(dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

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
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dzeta(dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1  ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_xn_tens(sf%mn_nyq(1),sf%mn_nyq(2),sf%X_IP(1,1:sf%mn_nyq(1)),sf%X_IP(2,1:PRODUCT(sf%mn_nyq(1:2)):sf%mn_nyq(1)),DERIV_ZETA,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_xn_tens_dzeta(dofs)|) ', checkreal
    END IF !TEST


    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

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
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dzeta,x(:),dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-Xmn(1,iMode)**2,wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-Xmn(1,iMode)**2,wp)*COS(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_THET_THET,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dthet_dthet(dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP / dofs from test 208
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal,ABS(g_IP(i_mn)-sf%evalDOF_x(sf%x_IP(:,i_mn),DERIV_THET_THET,dofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dthet_dthet,x(:),dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(Xmn(1,iMode)*Xmn(2,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(Xmn(1,iMode)*Xmn(2,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_THET_ZETA,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dthet_dzeta(dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP / dofs from test 210
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal,ABS(g_IP(i_mn)-sf%evalDOF_x(sf%x_IP(:,i_mn),DERIV_THET_ZETA,dofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dthet_dzeta,x(:),dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    g_IP=0.
    DO iMode=sin_range(1)+1,sin_range(2)
      dofs(iMode)=0.1_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))/(1.0_wp+SQRT(REAL(Xmn(1,iMode)**2+Xmn(2,iMode)**2,wp)))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-Xmn(2,iMode)**2,wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      dofs(iMode)=0.2_wp*(REAL(iMode-modes/2,wp)/REAL(modes,wp))/(1.0_wp+SQRT(REAL(Xmn(1,iMode)**2+Xmn(2,iMode)**2,wp)))
      g_IP(:) =g_IP(:)+dofs(iMode)*REAL(-Xmn(2,iMode)**2,wp)*COS(REAL(Xmn(1,iMode),wp)*sf%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%x_IP(2,:))
    END DO !iMode
    checkreal=MAXVAL(ABS(g_IP-sf%evalDOF_IP(DERIV_ZETA_ZETA,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP-evalDOF_dzeta_dzeta(dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP / dofs from test 212
    checkreal=0.0_wp
    DO i_mn=1,sf%mn_IP
      checkreal=MAX(checkreal,ABS(g_IP(i_mn)-sf%evalDOF_x(sf%x_IP(:,i_mn),DERIV_ZETA_ZETA,dofs)))
    END DO
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dzeta_dzeta,x(:),dofs)|) ', checkreal
    END IF !TEST

    iTest=iTest+1 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    !use g_IP / dofs from test 212, test evalDOF_xn
    checkreal=MAXVAL(ABS(g_IP(1:sf%mn_IP/2)-sf%evalDOF_xn(sf%mn_IP/2,sf%x_IP(1:2,1:sf%mn_IP/2),DERIV_ZETA_ZETA,dofs)))
    refreal=0.0_wp

    IF(testdbg.OR.(.NOT.( ABS(checkreal-refreal).LT. realtol))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! FBASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,I6," , ",I6,(A,I4),A,2(A,E11.3))') &
       ' mn_max= (',m_max,n_max, &
       ' )  nfp    = ',nfp, &
       ' ,  sin/cos : '//TRIM( sin_cos_map(sin_cos)), &
      '\n =>  should be ', refreal,' : MAX(|g_IP(:)-evalDOF_x(dzeta_dzeta,x(:),dofs)|) ', checkreal
    END IF !TEST

  END IF !testlevel <=2
  END ASSOCIATE !sf
  test_called=.FALSE.

END SUBROUTINE fBase_test


END MODULE MODgvec_fBase
