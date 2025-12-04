!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!>
!!# **TEST GVEC TO GENE** Driver program
!!
!! to test, just execute in ini/toksy
!! ../../build/bin/test_gvec_to_gene TOKSY_State_0000_00000000.dat
!! or ini/w7x
!! ../../build/bin/test_gvec_to_gene W7X_State_0000_00000000.dat
!===================================================================================================================================
PROGRAM TEST_GVEC_TO_GENE
USE MODgvec_Globals
USE MODgvec_gvec_to_gene
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
INTEGER                 :: nArgs
CHARACTER(LEN=255)      :: filename
REAL(wp)                :: StartTime,EndTime
REAL(wp)                :: Fa,minor_r,spos,q,q_prime,p,p_prime,phiPrime_edge,q_edge
INTEGER                 :: n0_global,is,i,j
INTEGER,PARAMETER       :: nthet=11
INTEGER,PARAMETER       :: nzeta=22
INTEGER                 :: SFLcoord_test
REAL(wp),DIMENSION(nthet,nzeta)   :: theta_star,theta,zeta
REAL(wp),DIMENSION(3,nthet,nzeta) :: cart_coords,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB
!===================================================================================================================================
  CALL CPU_TIME(StartTime)
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF(nArgs.GE.1)THEN
    CALL GET_COMMAND_ARGUMENT(1,filename)
  ELSE
    STOP 'GVEC_TO_GENE: gvec filename not given, usage: "./executable gvec_file.dat"'
  END IF


  !header
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(5(("*",A128,2X,"*",:,"\n")))')&
 '  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - -  GVEC ==> GENE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - -  GVEC ==> GENE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - -  GVEC ==> GENE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '&
,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  '
  WRITE(Unit_stdOut,'(132("="))')

  DO SFLcoord_test=0,2
  WRITE(Unit_stdOut,'(A,I4)')'TESTING SFLCOORD= ',SFLcoord_test
  WRITE(Unit_stdOut,'(132("="))')
  !initialization phase
  CALL Init_gvec_to_gene(filename,SFLcoord_in=SFLcoord_test,factorSFL_in=2)  !factorSFL=2 for testing purposes only, safe is 4

  CALL gvec_to_gene_scalars(Fa,minor_r,phiPrime_edge,q_edge,n0_global)
  WRITE(UNIT_stdOut,'(A,g21.13)')'Fa',Fa
  WRITE(UNIT_stdOut,'(A,g21.13)')'minor_r',minor_r
  WRITE(UNIT_stdOut,'(A,g21.13)')"Phi'(s=1)",PhiPrime_edge
  WRITE(UNIT_stdOut,'(A,g21.13)')"q(s=1)=Phi'(s=1)/chi'(s=1)",q_edge
  WRITE(UNIT_stdOut,'(A,g21.13)')'n0_global',n0_global
  WRITE(UNIT_stdOut,'(80("-"))')
  DO is=0,8
    spos=0.01+REAL(is)/REAL(8)*0.98
    CALL gvec_to_gene_profile(spos,q,q_prime,p,p_prime)
    WRITE(UNIT_stdOut,'(4(A,g21.13))') &
               's= ',spos &
              ,', q(s)= ', q &
              ,', q_prime(s)= ',q_prime &
              ,', p(s)= ',p &
              ,', p_prime(s)= ',p_prime
    DO i=1,nthet; DO j=1,nzeta
      zeta(i,j)=-PI + REAL(j-1)/REAL(nthet)*2.*Pi
      theta_star(i,j)=REAL(i-1)/REAL(nthet)*2.*Pi - 1.5*zeta(i,j)
    END DO ; END DO
    WRITE(UNIT_stdOut,*)'TESTING gvec_to_gene_coords...'
    CALL gvec_to_gene_coords( nthet,nzeta,spos,theta_star,zeta,theta,cart_coords)
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN x,y,z   : ',MINVAL(cart_coords(1,:,:)) &
                                          ,MINVAL(cart_coords(2,:,:)) &
                                          ,MINVAL(cart_coords(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX x,y,z   : ',MAXVAL(cart_coords(1,:,:)) &
                                          ,MAXVAL(cart_coords(2,:,:)) &
                                          ,MAXVAL(cart_coords(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN th*,th  : ',MINVAL(theta_star) &
                                          ,MINVAL(theta)
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX th*,th  : ',MAXVAL(theta_star) &
                                          ,MAXVAL(theta)

    WRITE(UNIT_stdOut,*)'TESTING gvec_to_gene_metrics...'
    CALL gvec_to_gene_metrics(nthet,nzeta,spos,theta_star,zeta,grad_s,grad_theta_star,grad_zeta,Bfield,grad_absB)
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN grads   : ',MINVAL(grad_s(1,:,:)) &
                                          ,MINVAL(grad_s(2,:,:)) &
                                          ,MINVAL(grad_s(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX grads   : ',MAXVAL(grad_s(1,:,:)) &
                                          ,MAXVAL(grad_s(2,:,:)) &
                                          ,MAXVAL(grad_s(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN gradth* : ',MINVAL(grad_theta_star(1,:,:)) &
                                          ,MINVAL(grad_theta_star(2,:,:)) &
                                          ,MINVAL(grad_theta_star(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX gradth* : ',MAXVAL(grad_theta_star(1,:,:)) &
                                          ,MAXVAL(grad_theta_star(2,:,:)) &
                                          ,MAXVAL(grad_theta_star(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN grad_zet: ',MINVAL(grad_zeta(1,:,:)) &
                                          ,MINVAL(grad_zeta(2,:,:)) &
                                          ,MINVAL(grad_zeta(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX grad_zet: ',MAXVAL(grad_zeta(1,:,:)) &
                                          ,MAXVAL(grad_zeta(2,:,:)) &
                                          ,MAXVAL(grad_zeta(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN Bfield  : ',MINVAL(Bfield(1,:,:)) &
                                          ,MINVAL(Bfield(2,:,:)) &
                                          ,MINVAL(Bfield(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX Bfield  : ',MAXVAL(Bfield(1,:,:)) &
                                          ,MAXVAL(Bfield(2,:,:)) &
                                          ,MAXVAL(Bfield(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MIN grad|B| : ',MINVAL(grad_absB(1,:,:)) &
                                          ,MINVAL(grad_absB(2,:,:)) &
                                          ,MINVAL(grad_absB(3,:,:))
    WRITE(UNIT_stdOut,'(A,3g21.13)')'MAX grad|B| : ',MAXVAL(grad_absB(1,:,:)) &
                                          ,MAXVAL(grad_absB(2,:,:)) &
                                          ,MAXVAL(grad_absB(3,:,:))
    WRITE(UNIT_stdOut,'(80("-"))')
  END DO !spos



  CALL Finalize_gvec_to_gene()

  CALL CPU_TIME(EndTime)
  WRITE(Unit_stdOut,fmt_sep)
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' TEST GVEC TO GENE FINISHED! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,fmt_sep)
  END DO !SFLcoord_test

END PROGRAM TEST_GVEC_TO_GENE
