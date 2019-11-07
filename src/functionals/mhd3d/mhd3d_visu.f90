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
!!# Module ** MHD3D Variables **
!!
!!
!!
!===================================================================================================================================
MODULE MODgvec_MHD3D_visu
! MODULES
USE MODgvec_Globals,ONLY: wp,Unit_stdOut,abort
IMPLICIT NONE
PUBLIC



!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE visu_BC_face(mn_IP ,minmax,fileID)
! MODULES
USE MODgvec_Globals,    ONLY: TWOPI
USE MODgvec_MHD3D_vars, ONLY: X1_base,X2_base,LA_base,hmap,U
USE MODgvec_output_vtk, ONLY: WriteDataToVTK
USE MODgvec_Output_vars,ONLY: Projectname,OutputLevel
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: mn_IP(2) !! muber of points in theta,zeta direction
  REAL(wp)      , INTENT(IN   ) :: minmax(2:3,0:1) !! min/max of theta,zeta [0,1] 
  INTEGER       , INTENT(IN   ) :: fileID          !! added to file name before the ending
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i_m,i_n,nplot(2),iMode
  REAL(wp) :: xIP(2),q(3)
  REAL(wp) :: X1_visu,X2_visu
  REAL(wp) :: coord_visu(3,mn_IP(1),mn_IP(2),1)
  INTEGER,PARAMETER  :: nVal=1
  REAL(wp) :: var_visu(nVal,mn_IP(1),mn_IP(2),1)
  REAL(wp) :: thet(mn_IP(1)),zeta(mn_IP(2))
  REAL(wp) :: spos
  REAL(wp) :: X1_s(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes)
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
  CHARACTER(LEN=255) :: FileName
!===================================================================================================================================
  IF((minmax(2,1)-minmax(2,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visuBC, nothing to visualize since theta-range is <=0, theta_min= ',minmax(2,0),', theta_max= ',minmax(2,1)
    RETURN
  ELSEIF((minmax(3,1)-minmax(3,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visuBC, nothing to visualize since zeta-range is <=0, zeta_min= ',minmax(3,0),', zeta_max= ',minmax(3,1)
    RETURN
  END IF
  DO i_m=1,mn_IP(1)
    thet(i_m)= TWOPI*(minmax(2,0)+(minmax(2,1)-minmax(2,0))*REAL(i_m-1,wp)/REAL(mn_IP(1)-1,wp)) !repeat point exactly
  END DO
  DO i_n=1,mn_IP(2)
    zeta(i_n)=TWOPI*(minmax(3,0)+(minmax(3,1)-minmax(3,0))*REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp))
  END DO
  IF(hmap%which_hmap.NE.3)THEN !not for cylinder
    IF(ABS((minMax(2,1)-minmax(2,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
      thet(mn_IP(1))=thet(1)
    END IF
    IF(ABS((minMax(3,1)-minmax(3,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
      zeta(mn_IP(2))=zeta(1)
    END IF
  END IF!hmap not cylinder
  spos=0.99999999_wp
  DO iMode=1,X1_base%f%modes
    X1_s( iMode)= X1_base%s%evalDOF_s(spos,      0,U(0)%X1(:,iMode))
  END DO
  DO iMode=1,X2_base%f%modes
    X2_s(iMode) = X2_base%s%evalDOF_s(spos,      0,U(0)%X2(:,iMode))
  END DO
  DO iMode=1,LA_base%f%modes
    LA_s(iMode) = LA_base%s%evalDOF_s(spos,      0,U(0)%LA(:,iMode))
  END DO
  DO i_n=1,mn_IP(2)
    xIP(2)  = zeta(i_n)
    DO i_m=1,mn_IP(1)
      xIP(1)= thet(i_m)
      X1_visu=X1_base%f%evalDOF_x(xIP,0,X1_s)
      X2_visu=X2_base%f%evalDOF_x(xIP,0,X2_s)
      q=(/X1_visu,X2_visu,xIP(2)/)
      coord_visu(  :,i_m,i_n,1)=hmap%eval(q)
      var_visu(1,i_m,i_n,1)=LA_base%f%evalDOF_x(xIP,0,LA_s)
    END DO !i_m
  END DO !i_n
  VarNames(1)="lambda"
  nplot(:)=mn_IP-1
  WRITE(filename,'(A,"_visu_BC_",I4.4,"_",I8.8,".vtu")')TRIM(Projectname),outputLevel,fileID
  CALL WriteDataToVTK(2,3,nVal,nplot,1,VarNames,coord_visu,var_visu,TRIM(filename))

END SUBROUTINE visu_BC_face


!===================================================================================================================================
!> visualize the mapping and additional variables in 3D, either on zeta=const planes or fully 3D 
!!
!===================================================================================================================================
SUBROUTINE visu_3D(np_in,minmax,only_planes,fileID )
! MODULES
USE MODgvec_Globals,        ONLY: TWOPI,PI,CROSS
USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,hmap,sgrid,U,F
USE MODgvec_MHD3D_Profiles, ONLY: Eval_iota,Eval_pres,Eval_Phi,Eval_PhiPrime,Eval_chiPrime
USE MODgvec_output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_Output_vars,    ONLY: Projectname,OutputLevel
USE MODgvec_Analyze_Vars,   ONLY: SFL_theta
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: np_in(3)     !! (1) #points in s & theta per element,(2:3) #elements in  theta,zeta
  REAL(wp)      , INTENT(IN   ) :: minmax(3,0:1)  !! minimum /maximum range in s,theta,zeta [0,1] 
  LOGICAL       , INTENT(IN   ) :: only_planes  !! true: visualize only planes, false:  full 3D
  INTEGER       , INTENT(IN   ) :: fileID          !! added to file name before the ending
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i_s,j_s,i_m,i_n,iElem,nElems,nplot(3),minElem,maxElem
  REAL(wp) :: spos,xIP(2),q(3)
  REAL(wp) :: X1_s(X1_base%f%modes),F_X1_s(X1_base%f%modes),dX1ds(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes),F_X2_s(X2_base%f%modes),dX2ds(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes),F_LA_s(LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu,dX1_ds_visu,dX2_ds_visu,dX1_dthet_visu,dX1_dzeta_visu,dX2_dthet_visu,dX2_dzeta_visu
  REAL(wp) :: dLA_dthet_visu,dLA_dzeta_visu,iota_s,pres_s,chiPrime_s,phiPrime_s,e_s(3),e_thet(3),e_zeta(3)

  INTEGER,PARAMETER  :: nVal=14
  REAL(wp) :: coord_visu( 3,np_in(1),np_in(1),np_in(3),np_in(2),sgrid%nElems)
  REAL(wp) :: var_visu(nVal,np_in(1),np_in(1),np_in(3),np_in(2),sgrid%nElems)
  REAL(wp) :: thet(np_in(1),np_in(2)),zeta(np_in(3))
  REAL(wp) :: theta_star,sqrtG
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
  CHARACTER(LEN=255) :: filename
!===================================================================================================================================
  __PERFON("output_visu")
  IF(only_planes)THEN
    SWRITE(UNIT_stdOut,'(A)') 'Start visu planes...'
  ELSE
    SWRITE(UNIT_stdOut,'(A)') 'Start visu 3D...'
  END IF
  IF((minmax(1,1)-minmax(1,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
     'WARNING visu3D, nothing to visualize since s-range is <=0, s_min= ',minmax(1,0),', s_max= ',minmax(1,1)
    RETURN
  ELSEIF((minmax(2,1)-minmax(2,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visu3D, nothing to visualize since theta-range is <=0, theta_min= ',minmax(2,0),', theta_max= ',minmax(2,1)
    RETURN
  ELSEIF((minmax(3,1)-minmax(3,0)).LE.1e-08)THEN
    SWRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visu3D, nothing to visualize since zeta-range is <=0, zeta_min= ',minmax(3,0),', zeta_max= ',minmax(3,1)
    RETURN
  END IF
  __PERFON("prepare_visu")
  VarNames( 1)="lambda"
  VarNames( 2)="sqrtG"
  VarNames( 3)="Phi"
  VarNames( 4)="iota"
  VarNames( 5)="pressure"
  VarNames( 6)="BvecX"
  VarNames( 7)="BvecY"
  VarNames( 8)="BvecZ"
  VarNames( 9)="F_X1"
  VarNames(10)="F_X2"
  VarNames(11)="F_LA"
  VarNames(12)="s"
  VarNames(13)="theta"
  VarNames(14)="zeta"

  var_visu=0.

  ASSOCIATE(n_s=>np_in(1), mn_IP=>np_in(2:3) )
  DO i_m=1,mn_IP(1)
    DO j_s=1,n_s
      thet(j_s,i_m)=TWOPI*(minmax(2,0)+(minmax(2,1)-minmax(2,0)) &
                            *REAL((j_s-1)+(i_m-1)*(n_s-1),wp)/REAL((np_in(1)-1)*mn_IP(1),wp))
    END DO !j_s
  END DO
  DO i_n=1,mn_IP(2)
    zeta(i_n)=TWOPI*(minmax(3,0)+(minmax(3,1)-minmax(3,0))*REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp))
  END DO

  nElems=sgrid%nElems
  DO iElem=1,nElems
    DO i_s=1,n_s
!      spos=sgrid%sp(iElem-1)+(1.0e-06_wp+REAL(i_s-1,wp))/(2.0e-06_wp+REAL(n_s-1,wp))*sgrid%ds(iElem)
      spos=MAX(1.0e-06,sgrid%sp(iElem-1)+(REAL(i_s-1,wp))/(REAL(n_s-1,wp))*sgrid%ds(iElem))
 
      X1_s(:)   = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,      0,U(0)%X1(:,:))
      dX1ds(:)  = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,DERIV_S,U(0)%X1(:,:))
      F_X1_s(:) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,      0,F(0)%X1(:,:))
      X2_s(:)   = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,      0,U(0)%X2(:,:))
      dX2ds(:)  = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,DERIV_S,U(0)%X2(:,:))
      F_X2_s(:) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,      0,F(0)%X2(:,:))
      LA_s(:)   = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,      0,U(0)%LA(:,:))
      F_LA_s(:) = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,      0,F(0)%LA(:,:))

      iota_s=Eval_iota(spos)
      pres_s=Eval_pres(spos)
      phiPrime_s=Eval_PhiPrime(spos)
      chiPrime_s=Eval_chiPrime(spos)
      var_visu( 3,i_s,:,:,:,iElem) =Eval_Phi(spos)
      var_visu( 4,i_s,:,:,:,iElem) =iota_s
      var_visu( 5,i_s,:,:,:,iElem) =pres_s
      var_visu(12,i_s,:,:,:,iElem) =spos
      !define theta2, which corresponds to the theta angle of a given theta_star=theta

!$OMP PARALLEL DO  COLLAPSE(3)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_m,i_n,j_s,xIP,q,sqrtG,e_s,e_thet,e_zeta,theta_star, &
!$OMP           X1_visu,dX1_ds_visu,dX1_dthet_visu,dX1_dzeta_visu,    &
!$OMP           X2_visu,dX2_ds_visu,dX2_dthet_visu,dX2_dzeta_visu,dLA_dthet_visu,dLA_dzeta_visu )  &
!$OMP   SHARED(np_in,i_s,iElem,thet,zeta,SFL_theta,X1_base,X2_base,LA_base,X1_s,X2_s,LA_s,dX1ds,dX2ds,&
!$OMP          F_X1_s,F_X2_s,F_LA_s,hmap,coord_visu,var_visu,chiPrime_s,phiPrime_s)
      DO i_m=1,mn_IP(1)
        DO i_n=1,mn_IP(2)
          DO j_s=1,n_s
            xIP(2)  = zeta(i_n)
            IF(SFL_theta)THEN
              theta_star=thet(j_s,i_m)
              CALL Get_SFL_theta(theta_star,xIP(2),LA_base,LA_s,xIP(1))
            ELSE
              xIP(1)= thet(j_s,i_m)
            END IF
            var_visu(13,i_s,j_s,i_n,i_m,iElem)=xIP(1) !theta for evaluation of X1,X2,LA
            var_visu(14,i_s,j_s,i_n,i_m,iElem)=xIP(2) !zeta  for evaluation of X1,X2,LA

            X1_visu         =X1_base%f%evalDOF_x(xIP,         0,X1_s )
            X2_visu         =X2_base%f%evalDOF_x(xIP,         0,X2_s )
            
            dX1_ds_visu     =X1_base%f%evalDOF_x(xIP,         0,dX1ds)
            dX2_ds_visu     =X2_base%f%evalDOF_x(xIP,         0,dX2ds)
            
            dX1_dthet_visu  =X1_base%f%evalDOF_x(xIP,DERIV_THET,X1_s )
            dX2_dthet_visu  =X2_base%f%evalDOF_x(xIP,DERIV_THET,X2_s )
            dLA_dthet_visu  =LA_base%f%evalDOF_x(xIP,DERIV_THET,LA_s )
            
            dX1_dzeta_visu  =X1_base%f%evalDOF_x(xIP,DERIV_ZETA,X1_s )
            dX2_dzeta_visu  =X2_base%f%evalDOF_x(xIP,DERIV_ZETA,X2_s )
            dLA_dzeta_visu  =LA_base%f%evalDOF_x(xIP,DERIV_ZETA,LA_s )
            
            q=(/X1_visu,X2_visu,xIP(2)/)
            !x,y,z
            coord_visu(:,i_s,j_s,i_n,i_m,iElem )=hmap%eval(q)

            e_s   =hmap%eval_dxdq(q,(/dX1_ds_visu,dX2_ds_visu,0.0_wp/))
            e_thet=hmap%eval_dxdq(q,(/dX1_dthet_visu,dX2_dthet_visu,0.0_wp/))
            e_zeta=hmap%eval_dxdq(q,(/dX1_dzeta_visu,dX2_dzeta_visu,1.0_wp/))

            sqrtG = hmap%eval_Jh(q)*(dX1_ds_visu*dX2_dthet_visu -dX2_ds_visu*dX1_dthet_visu) 
            !IF(ABS(sqrtG- SUM(e_s*(CROSS(e_thet,e_zeta)))).GT.1.0e-04) STOP 'test sqrtg failed'
            
            !lambda
            var_visu(  1,i_s,j_s,i_n,i_m,iElem) = LA_base%f%evalDOF_x(xIP,0,LA_s)
            !sqrtG
            var_visu(  2,i_s,j_s,i_n,i_m,iElem) = sqrtG
            !F_X1,F_X2,F_LA  
            var_visu(  9,i_s,j_s,i_n,i_m,iElem) = X1_base%f%evalDOF_x(xIP,         0,F_X1_s )
            var_visu( 10,i_s,j_s,i_n,i_m,iElem) = X2_base%f%evalDOF_x(xIP,         0,F_X2_s )
            var_visu( 11,i_s,j_s,i_n,i_m,iElem) = LA_base%f%evalDOF_x(xIP,         0,F_LA_s )
            !Bvec
            var_visu(6:8,i_s,j_s,i_n,i_m,iElem)= & 
                          (  e_thet(:)*(chiPrime_s-PhiPrime_s*dLA_dzeta_visu)  &
                           + e_zeta(:)*PhiPrime_s*(1.0_wp+dLA_dthet_visu) )*(1.0_wp/(sqrtG+sign(sqrtG,1.)*1.0e-12))
          END DO !j_s
        END DO !i_n
      END DO !i_m
!OMP END PARALLEL DO
    END DO !i_s
  END DO !iElem
  
  !make grid exactly periodic
  IF(hmap%which_hmap.NE.3)THEN !not for cylinder
    !make theta direction exactly periodic
    IF(ABS((minMax(2,1)-minmax(2,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
!$OMP PARALLEL DO  COLLAPSE(3)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iElem,i_n,i_s)
      DO iElem=1,nElems; DO i_n=1,mn_IP(2); DO i_s=1,n_s
        coord_visu( :,i_s,n_s,i_n,mn_IP(1),iElem)=coord_visu( :,i_s,1,i_n,1,iElem)
      END DO; END DO; END DO
!$OMP END PARALLEL DO
    END IF
    !make zeta direction exactly periodic, only for 3Dvisu
    IF(.NOT.only_planes)THEN
      IF(ABS((minMax(3,1)-minmax(3,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
!$OMP PARALLEL DO  COLLAPSE(4)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iElem,i_n,i_s,j_s)
      DO iElem=1,nElems; DO i_m=1,mn_IP(1); DO j_s=1,n_s; DO i_s=1,n_s
        coord_visu( :,i_s,j_s,mn_IP(2),i_m,iElem)=coord_visu( :,i_s,j_s,1,i_m,iElem)
      END DO; END DO; END DO; END DO
!$OMP END PARALLEL DO
      END IF
    END IF
  END IF!hmap not cylinder
  __PERFOFF("prepare_visu")
  __PERFON("write_visu")
  !range s: include all elements belonging to [smin,smax]
  minElem=MAX(     1,sgrid%find_elem(minmax(1,0))-1)
  maxElem=MIN(nElems,sgrid%find_elem(minmax(1,1))+1)
  IF(only_planes)THEN
    nplot(1:2)=(/n_s,n_s/)-1
    WRITE(filename,'(A,"_visu_planes_",I4.4,"_",I8.8,".vtu")')TRIM(Projectname),outputLevel,fileID
    CALL WriteDataToVTK(2,3,nVal,nplot(1:2),(mn_IP(1)*mn_IP(2)*(maxElem-minElem+1)),VarNames, &
                        coord_visu(:,:,:,:,:,minElem:maxElem), &
                          var_visu(:,:,:,:,:,minElem:maxElem),TRIM(filename))
  ELSE
    !3D
    nplot(1:3)=(/n_s,n_s,mn_IP(2)/)-1
    WRITE(filename,'(A,"_visu_3D_",I4.4,"_",I8.8,".vtu")')TRIM(Projectname),outputLevel,fileID
    CALL WriteDataToVTK(3,3,nVal,nplot,mn_IP(1)*(maxElem-minElem+1),VarNames, &
                        coord_visu(:,:,:,:,:,minElem:maxElem), &
                          var_visu(:,:,:,:,:,minElem:maxElem),TRIM(filename))
  END IF
  __PERFOFF("write_visu")
  
  END ASSOCIATE!n_s,mn_IP
  SWRITE(UNIT_stdOut,'(A)') '... DONE.'
  __PERFOFF("output_visu")
END SUBROUTINE visu_3D


!===================================================================================================================================
!> Finds theta_out that satisfies nonlinear equation:
!> theta_star=theta_out + lambda(theta_out,zeta)
!!
!===================================================================================================================================
SUBROUTINE Get_SFL_theta(theta_star,zeta,LA_base_in,LA,theta_out) 
! MODULES
USE MODgvec_Globals,ONLY: PI
USE MODgvec_base   ,ONLY: t_base
USE MODgvec_Newton ,ONLY: NewtonRoot1D_FdF
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL(wp)     ,INTENT(IN ) :: theta_star                 !< theta* (=start value for Newton)
  REAL(wp)     ,INTENT(IN ) :: zeta                       !< fixed zeta position
  CLASS(t_Base),INTENT(IN ) :: LA_base_in                 !< basis for lambda
  REAL(wp)     ,INTENT(IN ) :: LA(1:LA_base_in%f%modes)   !< modes of Lambda (for a fixed spos)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)     ,INTENT(OUT) :: theta_out                  !< theta_out => 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
 !                                     a             b             maxstep  , xinit       , F0         ,func
 theta_out=NewtonRoot1D_FdF(1.0e-12_wp,theta_star-PI,theta_star+PI,0.1_wp*PI,theta_star   , theta_star,FRdFR)

!for iteration on theta^*
CONTAINS 

  FUNCTION FRdFR(theta_iter)
    !uses current zeta where newton is called, and LA_s from subroutine above
    IMPLICIT NONE
    REAL(wp) :: theta_iter
    REAL(wp) :: FRdFR(2) !output
    !--------------------------------------------------- 
    FRdFR(1)=theta_iter+LA_base_in%f%evalDOF_x((/theta_iter,zeta/),         0,LA(:))  !theta_iter+lambda
    FRdFR(2)=1.0_wp    +LA_base_in%f%evalDOF_x((/theta_iter,zeta/),DERIV_THET,LA(:)) !1+dlambda/dtheta
  END FUNCTION FRdFR

END SUBROUTINE Get_SFL_theta


!===================================================================================================================================
!> check distance between two solutions, via sampling X1,X2 at theta*=theta+lambda, and comparing the distance of 
!> the sampled x,y,z coordinates 
!!
!===================================================================================================================================
SUBROUTINE CheckDistance(U,V,maxDist,avgDist) 
! MODULES
  USE MODgvec_Globals,        ONLY: TWOPI
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,hmap,sgrid
  USE MODgvec_sol_var_MHD3D,  ONLY: t_sol_var_mhd3d
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN) :: U !U and V must be have the same basis and grid!
  CLASS(t_sol_var_MHD3D), INTENT(IN) :: V
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),INTENT(OUT)    :: maxDist,avgDist
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: n_s,mn_IP(2),nElems
  INTEGER  :: i_s,i_m,i_n,iElem
  REAL(wp) :: spos,zeta,theta,theta0
  REAL(wp) :: UX1_s(1:X1_base%f%modes),VX1_s(1:X1_base%f%modes)
  REAL(wp) :: UX2_s(1:X2_base%f%modes),VX2_s(1:X2_base%f%modes)
  REAL(wp) :: ULA_s(1:LA_base%f%modes),VLA_s(1:LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu,LA_visu
  REAL(wp) :: q(3),xU(3),xV(3),dist,xIP(2)
  REAL(wp),ALLOCATABLE :: theta1D(:),zeta1D(:)
  LOGICAL  :: SFL_theta=.TRUE.
!===================================================================================================================================
  __PERFON("checkDistance")
  n_s=3 !number of points to check per element (1 at the left boundary, 2 inner, none at the right)
  mn_IP(1)   = MAX(1,X1_base%f%mn_nyq(1)/2)
  mn_IP(2)   = MAX(1,X1_base%f%mn_nyq(2)/2)
  nElems=sgrid%nElems
  
  maxDist=0.  
  avgDist=0.

  ALLOCATE(theta1D(1:mn_IP(1)),zeta1D(1:mn_IP(2)))
  DO i_n=1,mn_IP(2)
    zeta1D(i_n)  = TWOPI*REAL(i_n-1,wp)/REAL(mn_IP(2),wp) !do not include periodic point 
  END DO
  DO i_m=1,mn_IP(1)
    theta1D(i_m)= TWOPI*REAL(i_m-1,wp)/REAL(mn_IP(1),wp)  !do not include periodic point
  END DO

  DO iElem=1,nElems
    DO i_s=1,n_s
      spos=MAX(1.0e-06,sgrid%sp(iElem-1)+(REAL(i_s-1,wp))/(REAL(n_s,wp))*sgrid%ds(iElem)) !includes axis but not edge

      UX1_s(:) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,0,U%X1(:,:))
      VX1_s(:) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,0,V%X1(:,:))
      UX2_s(:) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,0,U%X2(:,:))
      VX2_s(:) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,0,V%X2(:,:))
      ULA_s(:) = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,0,U%LA(:,:))
      VLA_s(:) = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,0,V%LA(:,:))

!$OMP PARALLEL DO  COLLAPSE(2)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   REDUCTION(+:avgDist) REDUCTION(max:maxDist) &
!$OMP   PRIVATE(i_m,i_n,xIP,q,theta,zeta,theta0,X1_visu,X2_visu,LA_visu,xU,xV,dist ) &
!$OMP   SHARED(mn_IP,theta1D,zeta1D,SFL_theta,X1_base,X2_base,LA_base,UX1_s,UX2_s,ULA_s,VX1_s,VX2_s,VLA_s,hmap)
      DO i_n=1,mn_IP(2)
          DO i_m=1,mn_IP(1)
            zeta  = zeta1D(i_n)
            theta0= theta1D(i_m)
            !for xU
            IF(SFL_theta)THEN
              CALL Get_SFL_theta(theta0,zeta,LA_base,ULA_s, theta)
            ELSE
              LA_visu = LA_base%f%evalDOF_x((/theta0,zeta/),0,ULA_s(:) )
              theta = theta0 + LA_visu
            END IF
          
            xIP=(/theta,zeta/)
          
            X1_visu    = X1_base%f%evalDOF_x(xIP,0,UX1_s(:) )
            X2_visu    = X2_base%f%evalDOF_x(xIP,0,UX2_s(:) )
            
            q=(/X1_visu,X2_visu,zeta/)
            !x,y,z
            xU(:)=hmap%eval(q)
            
            !for xV
            IF(SFL_theta)THEN
              CALL Get_SFL_theta(theta0,zeta,LA_base,VLA_s, theta)
            ELSE
              LA_visu = LA_base%f%evalDOF_x((/theta0,zeta/),0,VLA_s(:) )
              theta = theta0 + LA_visu
            END IF
          
            xIP=(/theta,zeta/)
          
            X1_visu    = X1_base%f%evalDOF_x(xIP,0,VX1_s(:) )
            X2_visu    = X2_base%f%evalDOF_x(xIP,0,VX2_s(:) )
            
            q=(/X1_visu,X2_visu,zeta/)
            !x,y,z
            xV(:)=hmap%eval(q)
            
            dist=SQRT(SUM((xU(:)-xV(:))**2))
            maxDist = MAX(maxDist,dist)
            avgDist = avgDist+dist
          END DO !i_m
      END DO !i_n
!OMP END PARALLEL DO
    END DO !i_s
  END DO !iElem
  avgDist=avgDist/REAL(nElems*n_s*mn_IP(1)*mn_IP(2),wp)

  DEALLOCATE(theta1D,zeta1D)

  __PERFOFF("checkDistance")
END SUBROUTINE CheckDistance


!===================================================================================================================================
!> Visualize 
!!
!===================================================================================================================================
SUBROUTINE visu_1d_modes(n_s,fileID)
! MODULES
USE MODgvec_Analyze_Vars,  ONLY: visu1D
USE MODgvec_fbase,         ONLY: sin_cos_map
USE MODgvec_MHD3D_Vars,    ONLY: U,X1_base,X2_base,LA_base
USE MODgvec_Output_vars,   ONLY: outputLevel
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN   ) :: n_s    !! number of visualization points per element
  INTEGER, INTENT(IN   ) :: fileID !! added to file name before the ending
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL            :: vcase(4)
  CHARACTER(LEN=4)   :: vstr
  CHARACTER(LEN=80)  :: vname,fname
!===================================================================================================================================
  !visu1D: all possible combinations: 1,2,3,4,12,13,14,23,24,34,123,124,234,1234
  WRITE(vstr,'(I4)')visu1D
  vcase=.FALSE.
  IF(INDEX(vstr,'1').NE.0) vcase(1)=.TRUE.
  IF(INDEX(vstr,'2').NE.0) vcase(2)=.TRUE.
  IF(INDEX(vstr,'3').NE.0) vcase(3)=.TRUE.
  IF(INDEX(vstr,'4').NE.0) vcase(4)=.TRUE.
  IF(.NOT.(ANY(vcase))) THEN
    WRITE(*,*)'visu1D case not found:',visu1D,' nothing visualized...'
    RETURN
  END IF
  
  IF(vcase(1))THEN
    WRITE(*,*)'1.1) Visualize 1d profiles of derived quantities...'
    WRITE(fname,'(A,I4.4,"_",I8.8,A4)')'1Dprofiles_',outputLevel,FileID,'.csv'
    CALL eval_1d_profiles(n_s,fname) 
 
    WRITE(*,*)'1.2) Visualize gvec modes in 1D: R,Z,lambda interpolated...'
    vname="X1"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//TRIM(sin_cos_map(X1_base%f%sin_cos))//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,X1_base,U(0)%X1)
    vname="X2"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//TRIM(sin_cos_map(X1_base%f%sin_cos))//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,X2_base,U(0)%X2)
    vname="LA"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//TRIM(sin_cos_map(X1_base%f%sin_cos))//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,LA_base,U(0)%LA)
  END IF
  IF(vcase(2))THEN
    WRITE(*,*)'2) Visualize gvec modes in 1D: dRrho,dZrho interpolated...'
    vname="dX1"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//TRIM(sin_cos_map(X1_base%f%sin_cos))//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,X1_base,U(0)%X1)
    vname="dX2"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//TRIM(sin_cos_map(X1_base%f%sin_cos))//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,X2_base,U(0)%X2)
  END IF
  
  !
  
END SUBROUTINE visu_1d_modes

!===================================================================================================================================
!> 
!!
!===================================================================================================================================
SUBROUTINE eval_1d_profiles(n_s,fname_in)
! MODULES
  USE MODgvec_MHD3D_Profiles,ONLY: Eval_iota,Eval_pres,Eval_Phi
  USE MODgvec_MHD3D_Vars,    ONLY: sgrid
  USE MODgvec_Output_CSV, ONLY:WriteDataToCSV
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,         INTENT(IN   ) :: n_s    !! number of visualization points per element
  CHARACTER(LEN=*),INTENT(IN   ) :: fname_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: i,i_s,iElem,iVar,nVars,nvisu
  CHARACTER(LEN=120),ALLOCATABLE :: VarNames(:) 
  REAL(wp)          ,ALLOCATABLE :: values_visu(:,:)
!===================================================================================================================================
  nVars = 4
  nvisu = sgrid%nElems*n_s
  ALLOCATE(VarNames(nVars))
  ALLOCATE(values_visu(nVars,nvisu))
  iVar=1
  VarNames(1)='rho'
  DO iElem=1,sgrid%nElems
    DO i_s=1,n_s
      values_visu(1,i_s+(iElem-1)*n_s)=sgrid%sp(iElem-1)+(1.0e-06_wp+REAL(i_s-1,wp))/(2.0e-06_wp+REAL(n_s-1,wp))*sgrid%ds(iElem)
    END DO
  END DO
  ASSOCIATE(s_visu=>values_visu(1,:))
  iVar=iVar+1
  VarNames(iVar)='Phi'
  DO i=1,nvisu
    values_visu( iVar,i)=Eval_Phi(s_visu(i))
  END DO !i

  iVar=iVar+1
  Varnames(iVar)='iota(Phi_norm)'
  
  DO i=1,nvisu
    values_visu(  iVar,i)=Eval_iota(s_visu(i))
  END DO !i
  
  iVar=iVar+1
  Varnames(iVar)='pres(Phi_norm)'
  DO i=1,nvisu
    values_visu(  iVar,i)=Eval_pres(s_visu(i))
  END DO !i

  END ASSOCIATE !s_visu
  CALL WriteDataToCSV(VarNames(:) ,values_visu(:,:) ,TRIM(fname_in)  &
                                  ,append_in=.FALSE.,vfmt_in='E15.5')

END SUBROUTINE eval_1d_profiles

!===================================================================================================================================
!> Write all modes of one variable
!!
!===================================================================================================================================
SUBROUTINE writeDataMN_visu(n_s,fname_in,vname,rderiv,base_in,xx_in)
! MODULES
  USE MODgvec_base,          ONLY: t_base
  USE MODgvec_MHD3D_Profiles,ONLY: Eval_iota,Eval_pres,Eval_Phi
  USE MODgvec_MHD3D_Vars,    ONLY: sgrid
  USE MODgvec_write_modes,   ONLY: write_modes
  USE MODgvec_output_vars,   ONLY: Projectname
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,         INTENT(IN   ) :: n_s    !! number of visualization points per element
  INTEGER         ,INTENT(IN   ) :: rderiv !! 0: eval spl, 1: eval spl deriv
  CHARACTER(LEN=*),INTENT(IN   ) :: fname_in
  CHARACTER(LEN=*),INTENT(IN   ) :: vname
  TYPE(t_base)    ,INTENT(IN   ) :: base_in
  REAL(wp)        ,INTENT(INOUT) :: xx_in(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: i,i_s,iElem,nvisu
  INTEGER                        :: nVal,addval
  INTEGER                        :: iMode,j
  CHARACTER(LEN=80)              :: fname
  CHARACTER(LEN=120),ALLOCATABLE :: varnames(:) 
  REAL(wp)          ,ALLOCATABLE :: values_visu(:,:)
  REAL(wp)          ,ALLOCATABLE :: s_visu(:)
!===================================================================================================================================
  WRITE(fname,'(A,A,".csv")')TRIM(ProjectName)//'_modes_',TRIM(fname_in)
  nvisu   =sgrid%nElems*n_s
  
  addval = 5
  ALLOCATE(varnames(   addval+2*base_in%f%modes+2))
  ALLOCATE(values_visu(addval+2*base_in%f%modes+2,nvisu))
  ALLOCATE(s_visu(nvisu))
  
  DO iElem=1,sgrid%nElems
    DO i_s=1,n_s
      s_visu(i_s+(iElem-1)*n_s)=sgrid%sp(iElem-1)+(1.0e-06_wp+REAL(i_s-1,wp))/(2.0e-06_wp+REAL(n_s-1,wp))*sgrid%ds(iElem)
    END DO
  END DO
  
  nVal=1
  Varnames(   nVal)='Phi'
  DO i=1,nvisu
    values_visu(  nVal,i)=Eval_Phi(s_visu(i))
  END DO !i
  
  nVal=nVal+1
  Varnames(   nVal)='chi'
  values_visu(nVal,:)=0.0_wp !TODO
  
  nVal=nVal+1
  Varnames(   nVal)='rho'
  values_visu(nVal,:)=s_visu(:)
  
  nVal=nVal+1
  Varnames(nVal)='iota(Phi_norm)'
  
  DO i=1,nvisu
    values_visu(  nVal,i)=Eval_iota(s_visu(i))
  END DO !i
  
  nVal=nVal+1
  Varnames(nVal)='pres(Phi_norm)'
  DO i=1,nvisu
    values_visu(  nVal,i)=Eval_pres(s_visu(i))
  END DO !i

  DO iMode=base_in%f%sin_range(1)+1,base_in%f%sin_range(2)
    nVal=nVal+1
    WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname)//"_sin", &
      base_in%f%Xmn(1,iMode),base_in%f%Xmn(2,iMode)/base_in%f%nfp
    DO j=1,nvisu
      values_visu(nVal,j)=base_in%s%evalDOF_s(s_visu(j),rderiv,xx_in(:,iMode))
    END DO !j
  END DO
  DO iMode=base_in%f%cos_range(1)+1,base_in%f%cos_range(2)
    nVal=nVal+1
    WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname)//"_cos", &
      base_in%f%Xmn(1,iMode),base_in%f%Xmn(2,iMode)/base_in%f%nfp
    DO j=1,nvisu
      values_visu(nVal,j)=base_in%s%evalDOF_s(s_visu(j),rderiv,xx_in(:,iMode))
    END DO !j
  END DO

  CALL write_modes(fname,vname,nVal,base_in%f%modes,base_in%f%Xmn(1,:), &
                   base_in%f%Xmn(2,:),s_visu,sgrid%sp(1),values_visu(:,:),VarNames) 

  DEALLOCATE(varnames)
  DEALLOCATE(values_visu)
  DEALLOCATE(s_visu)
END SUBROUTINE writeDataMN_visu

END MODULE MODgvec_MHD3D_visu

