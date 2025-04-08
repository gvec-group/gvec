!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

!===================================================================================================================================
!>
!!# Module ** sol_var_MHD3D **
!!
!! Solution variable for MHD3D functional
!!
!===================================================================================================================================
MODULE MODgvec_sol_var_MHD3D
! MODULES
USE MODgvec_Globals,ONLY:wp,Unit_stdOut,abort,MPIRoot
USE MODgvec_c_sol_var,ONLY:c_sol_var
IMPLICIT NONE

PRIVATE
PUBLIC t_sol_var_MHD3D

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE,EXTENDS(c_sol_var) :: t_sol_var_MHD3D
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL               :: initialized=.FALSE.
  !---------------------------------------------------------------------------------------------------------------------------------
  REAL(wp)              :: W_MHD3D
  REAL(wp),CONTIGUOUS,POINTER :: X1(:,:)    !! X1 variable, shape=(base_s%nBase,base_f%mn_mode)
  REAL(wp),CONTIGUOUS,POINTER :: X2(:,:)    !! X2 variable
  REAL(wp),CONTIGUOUS,POINTER :: LA(:,:)    !! lambda variable
  REAL(wp),CONTIGUOUS,POINTER :: q(:)       !! 1d array container for all variables (is the one allocated)
  INTEGER, ALLOCATABLE  :: varsize(:,:)     !! varsize(a,b): size of rank a of variable b =1..nvars
  INTEGER, ALLOCATABLE  :: offset(:)        !! offset(0:nvars): size of variable b =offset(b)-offset(b-1)
  !---------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

  PROCEDURE  :: init          => sol_var_MHD3D_init
  PROCEDURE  :: free          => sol_var_MHD3D_free
  PROCEDURE  :: copy          => sol_var_MHD3D_copy
  PROCEDURE  :: set_to_solvar => sol_var_MHD3D_set_to_solvar
  PROCEDURE  :: set_to_scalar => sol_var_MHD3D_set_to_scalar
  GENERIC    :: set_to        => set_to_solvar,set_to_scalar  !chooses right routine depending on input type!
  PROCEDURE  :: norm_2        => sol_var_MHD3D_norm_2
  PROCEDURE  :: AXBY          => sol_var_MHD3D_AXBY
  !---------------------------------------------------------------------------------------------------------------------------------
END TYPE t_sol_var_MHD3D

LOGICAL :: test_called=.FALSE.
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> initialize (=allocate) sf of type t_sol_var
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_init( sf,varsize)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: sf  !!sf
  INTEGER         , INTENT(IN   ) :: varsize(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(sf%initialized)THEN
    WRITE(*,*)'WARNING: REALLOCATION OF SOL_VAR_MHD3D!!'
    CALL sf%free()
  END IF
  IF((size(varsize(:)).NE.6)) THEN
    CALL abort(__STAMP__, &
        "sol_var_MHD3D init: varsize is a vector of 3 entries!")
  END IF
  sf%nvars=3
  ALLOCATE(sf%varSize(2,sf%nvars),sf%offset(0:sf%nvars))
  sf%varSize(1,:)=varSize(1:3)
  sf%varSize(2,:)=varSize(4:6)
  sf%offset(0)=0
  sf%offset(1)=sf%offset(0)+PRODUCT(sf%varsize(:,1))
  sf%offset(2)=sf%offset(1)+PRODUCT(sf%varsize(:,2))
  sf%offset(3)=sf%offset(2)+PRODUCT(sf%varsize(:,3))
  ALLOCATE(sf%q(1:sf%offset(sf%nvars)))

  sf%X1(1:sf%varSize(1,1),1:sf%varSize(2,1))=>sf%q(1+sf%offset(0) : sf%offset(1))
  sf%X2(1:sf%varSize(1,2),1:sf%varSize(2,2))=>sf%q(1+sf%offset(1) : sf%offset(2))
  sf%LA(1:sf%varSize(1,3),1:sf%varSize(2,3))=>sf%q(1+sf%offset(2) : sf%offset(3))
  sf%W_MHD3D=-1.0_wp
  sf%initialized=.TRUE.
  CALL sf%set_to(0.0_wp)

  IF(.NOT.test_called) CALL sol_var_MHD3D_test(sf)
END SUBROUTINE sol_var_MHD3D_init


!===================================================================================================================================
!> free (=deallocate) sf of type t_sol_var_MHD3D
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_free( sf )
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized)RETURN
  sf%nvars=-1
  DEALLOCATE(sf%offset)
  DEALLOCATE(sf%varsize)
  NULLIFY(sf%X1)
  NULLIFY(sf%X2)
  NULLIFY(sf%LA)
  DEALLOCATE(sf%q)
  sf%W_MHD3D=-1.0_wp
  sf%initialized=.FALSE.
END SUBROUTINE sol_var_MHD3D_free

!===================================================================================================================================
!> copy tocopy  => sf
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_copy( sf ,tocopy)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_sol_var), INTENT(IN   ) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  SELECT TYPE(tocopy); TYPE IS(t_sol_var_MHD3D)
  IF(.NOT.tocopy%initialized)THEN
    CALL abort(__STAMP__, &
        "sol_var_MHD3D_copy: not initialized sol_var from which to copy!")
  END IF
  IF(sf%initialized) THEN
    WRITE(*,*)'WARNING!! reinit of sol_var_MHD3D in copy'
    CALL sf%free()
  END IF
  CALL sf%init((/size(tocopy%X1,1),size(tocopy%X2,1),size(tocopy%LA,1),  &
                 size(tocopy%X1,2),size(tocopy%X2,2),size(tocopy%LA,2)/) )
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,sf%offset(3)
    sf%q(i)=tocopy%q(i)
  END DO
!$OMP END PARALLEL DO
  sf%W_MHD3D=tocopy%W_MHD3D

  END SELECT !TYPE
END SUBROUTINE sol_var_MHD3D_copy

!===================================================================================================================================
!> set all variables to scalar
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_set_to_scalar( sf, scalar)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: sf  !!sf
  REAL(wp)              , INTENT(IN   ) :: scalar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  IF(.NOT.sf%initialized)THEN
    CALL abort(__STAMP__, &
        "sol_var_MHD3D not initialized in set_to!")
  END IF
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,sf%offset(3)
    sf%q(i)=scalar
  END DO
!$OMP END PARALLEL DO
  sf%W_MHD3D=0.0_wp
END SUBROUTINE sol_var_MHD3D_set_to_Scalar

!===================================================================================================================================
!> set variabes X1,X2,LA of toset  => sf, optional argument to scale toset with a scalar (for example -1.0_wp)
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_set_to_solvar( sf ,toset,scal_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(c_sol_var), INTENT(IN   ) :: toset
  REAL(wp),INTENT(IN),OPTIONAL    :: scal_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  SELECT TYPE(toset); TYPE IS(t_sol_var_MHD3D)
  IF(.NOT.toset%initialized)THEN
    CALL abort(__STAMP__, &
        "sol_var_MHD3D_set_to: not initialized sol_var from which to set!")
  END IF
  IF(PRESENT(scal_in))THEN
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
    DO i=1,sf%offset(3)
      sf%q(i)=scal_in*toset%q(i)
    END DO
!$OMP END PARALLEL DO
    sf%W_MHD3D=0.0_wp
  ELSE
!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
    DO i=1,sf%offset(3)
      sf%q(i)=toset%q(i)
    END DO
!$OMP END PARALLEL DO
    sf%W_MHD3D=toset%W_MHD3D
  END IF

  END SELECT !TYPE
END SUBROUTINE sol_var_MHD3D_set_to_solvar

!===================================================================================================================================
!> |X|^2, where X is of type t_var_sol, so three values are returned: |X1|^2,|X2|^2,|LA|^2
!!
!===================================================================================================================================
FUNCTION sol_var_MHD3D_norm_2( sf ) RESULT(norm_2)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: sf  !!self
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                       :: norm_2(sf%nvars)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF(.NOT.sf%initialized) CALL abort(__STAMP__,&
                              'taking the norm of a sol_var that is not initialized' )

  norm_2(1)=SUM(sf%X1*sf%X1)
  norm_2(2)=SUM(sf%X2*sf%X2)
  norm_2(3)=SUM(sf%LA*sf%LA)
END FUNCTION sol_var_MHD3D_norm_2

!===================================================================================================================================
!> res=a*X+b*Y , where X,Y,res are of type t_var_sol
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_AXBY(sf,aa,X,bb,Y)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)        , INTENT(IN   ) :: aa
  CLASS(c_sol_var), INTENT(IN   ) :: X
  REAL(wp)        , INTENT(IN   ) :: bb
  CLASS(c_sol_var), INTENT(IN   ) :: Y
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(INOUT) :: sf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: i
!===================================================================================================================================
  SELECT TYPE(X);  TYPE IS(t_sol_var_MHD3D)
  IF(.NOT.X%initialized) CALL abort(__STAMP__,&
                                      'AXBY: X not initialized')

  SELECT TYPE(Y);  TYPE IS(t_sol_var_MHD3D)
  IF(.NOT.Y%initialized) CALL abort(__STAMP__,&
                                       'AXBY: Y not initialized')

!$OMP PARALLEL DO        &
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i)
  DO i=1,sf%offset(3)
    sf%q(i) = aa*X%q(i) + bb*Y%q(i)
  END DO
!$OMP END PARALLEL DO
  END SELECT !Type
  END SELECT !Type
END SUBROUTINE sol_var_MHD3D_AXBY


!===================================================================================================================================
!> test sol_var_MHD3D
!!
!===================================================================================================================================
SUBROUTINE sol_var_MHD3D_test( sf )
USE MODgvec_GLobals, ONLY: UNIT_stdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: sf  !!sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp)           :: refreal,checkreal
  INTEGER            :: iTest
  TYPE(t_sol_var_MHD3D)    :: Utest(3)
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  CHARACTER(LEN=10)  :: fail
!===================================================================================================================================
  test_called=.TRUE. ! to prevent infinite loop in this routine
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN SOL_VAR_MHD3D TEST ID',nTestCalled,'    >>>>>>>>>'
  IF(testlevel.GE.1)THEN
    CALL Utest(1)%copy(sf)

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    CALL Utest(1)%set_to(0.435_wp)
    refreal = (0.435_wp)**2
    checkreal =SUM(Utest(1)%norm_2())/SUM(sf%varsize(1,:)*sf%varsize(2,:))
    IF(testdbg.OR.(.NOT.( (checkreal-refreal).LT. realtol) )) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SOL_VAR_MHD3D TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,3I4),2(A,E11.3))') &
      '   varsize_s= ',sf%varsize(1,:), &
      '   varsize_f= ',sf%varsize(2,:), &
      '\n =>  should be ', refreal, ' : norm_2(U=0)= ',checkreal
    END IF !TEST
    CALL Utest(1)%free()
  END IF !testlevel >=1
  IF(testlevel.GE.2)THEN
    CALL Utest(1)%copy(sf)
    CALL Utest(1)%set_to(0.5_wp)
    CALL Utest(2)%copy(sf)
    CALL Utest(2)%set_to(Utest(1))
    CALL Utest(3)%copy(sf)

    iTest=201 ; IF(testdbg)WRITE(*,*)'iTest=',iTest

    CALL Utest(1)%set_to(0.8_wp)
    CALL Utest(2)%set_to(-0.53_wp)
    CALL Utest(3)%AXBY(1.1_wp,Utest(1),-5.5_wp,Utest(2))
    refreal = (1.1_wp*0.8_wp-5.5_wp*(-0.53_wp))**2
    checkreal=SUM(Utest(3)%norm_2())/SUM(sf%varsize(1,:)*sf%varsize(2,:))
    IF(testdbg.OR.(.NOT.( (ABS(checkreal-refreal).LT. realtol) ))) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! SOL_VAR_MHD3D TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(2(A,3I4),2(A,E11.3))') &
      '   varsize_s= ',sf%varsize(1,:), &
      '   varsize_f= ',sf%varsize(2,:), &
      '\n =>  should be ', refreal, ' : norm_2(U%AXBY)= ',checkreal
    END IF !TEST

    CALL Utest(1)%free()
    CALL Utest(2)%free()
    CALL Utest(3)%free()
  END IF !testlevel >=2

  test_called=.FALSE. ! to prevent infinite loop in this routine


END SUBROUTINE sol_var_MHD3D_test

END MODULE MODgvec_Sol_var_MHD3D
