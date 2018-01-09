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
MODULE MOD_MHD3D_visu
! MODULES
USE MOD_Globals,ONLY: wp,Unit_stdOut,abort
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
USE MOD_Globals,    ONLY: TWOPI
USE MOD_MHD3D_vars, ONLY: X1_base,X2_base,LA_base,hmap,X1_b,X2_b,LA_b
USE MOD_output_vtk, ONLY: WriteDataToVTK
USE MOD_Output_vars,ONLY: Projectname,OutputLevel
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
  INTEGER  :: i_m,i_n,nplot(2)
  REAL(wp) :: xIP(2),q(3)
  REAL(wp) :: X1_visu,X2_visu
  REAL(wp) :: coord_visu(3,mn_IP(1),mn_IP(2),1)
  INTEGER,PARAMETER  :: nVal=1
  REAL(wp) :: var_visu(nVal,mn_IP(1),mn_IP(2),1)
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
  DO i_n=1,mn_IP(2)
    xIP(2)  = TWOPI*(minmax(3,0)+(minmax(3,1)-minmax(3,0))*REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp))
    DO i_m=1,mn_IP(1)
      xIP(1)= TWOPI*(minmax(2,0)+(minmax(2,1)-minmax(2,0))*REAL(i_m-1,wp)/REAL(mn_IP(1)-1,wp))
      X1_visu=X1_base%f%evalDOF_x(xIP,0,X1_b)
      X2_visu=X2_base%f%evalDOF_x(xIP,0,X2_b)
      q=(/X1_visu,X2_visu,xIP(2)/)
      coord_visu(  :,i_m,i_n,1)=hmap%eval(q)
      var_visu(1,i_m,i_n,1)=LA_base%f%evalDOF_x(xIP,0,LA_b)
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
USE MOD_Globals,        ONLY: TWOPI
USE MOD_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,hmap,sgrid,U,F
USE MOD_MHD3D_Profiles, ONLY: Eval_iota,Eval_pres,Eval_Phi,Eval_PhiPrime
USE MOD_output_vtk,     ONLY: WriteDataToVTK
USE MOD_Output_vars,    ONLY: Projectname,OutputLevel
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   ) :: np_in(3)     !! number of points in s,theta,zeta direction
  REAL(wp)      , INTENT(IN   ) :: minmax(3,0:1)  !! minimum /maximum range in s,theta,zeta [0,1] 
  LOGICAL       , INTENT(IN   ) :: only_planes  !! true: visualize only planes, false:  full 3D
  INTEGER       , INTENT(IN   ) :: fileID          !! added to file name before the ending
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: iMode,i_s,i_m,i_n,iElem,nElems,nplot(3),minElem,maxElem
  REAL(wp) :: spos,xIP(2),q(3)
  REAL(wp) :: X1_s(X1_base%f%modes),F_X1_s(X1_base%f%modes),dX1ds(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes),F_X2_s(X2_base%f%modes),dX2ds(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes),F_LA_s(LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu,dX1_ds_visu,dX2_ds_visu,dX1_dthet_visu,dX1_dzeta_visu,dX2_dthet_visu,dX2_dzeta_visu
  REAL(wp) :: dLA_dthet_visu,dLA_dzeta_visu,iota_s,pres_s,phiPrime_s,e_thet(3),e_zeta(3)

  INTEGER,PARAMETER  :: nVal=11
  REAL(wp) :: coord_visu(     3,np_in(1),np_in(2),np_in(3),sgrid%nElems)
  REAL(wp) :: var_visu(nVal,np_in(1),np_in(2),np_in(3),sgrid%nElems)
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
  CHARACTER(LEN=255) :: filename
!===================================================================================================================================
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
  VarNames(1)="lambda"
  VarNames(2)="sqrtG"
  VarNames(3)="Phi"
  VarNames(4)="iota"
  VarNames(5)="pressure"
  VarNames(6)="BvecX"
  VarNames(7)="BvecY"
  VarNames(8)="BvecZ"
  VarNames( 9)="F_X1"
  VarNames(10)="F_X2"
  VarNames(11)="F_LA"

  ASSOCIATE(n_s=>np_in(1), mn_IP=>np_in(2:3) )
  nElems=sgrid%nElems
  DO iElem=1,nElems
    DO i_s=1,n_s
      spos=sgrid%sp(iElem-1)+(1.0e-06_wp+REAL(i_s-1,wp))/(2.0e-06_wp+REAL(n_s-1,wp))*sgrid%ds(iElem)
      DO iMode=1,X1_base%f%modes
        X1_s( iMode)= X1_base%s%evalDOF_s(spos,      0,U(0)%X1(:,iMode))
        F_X1_s( iMode)= X1_base%s%evalDOF_s(spos,      0,F(0)%X1(:,iMode))
        dX1ds(iMode)= X1_base%s%evalDOF_s(spos,DERIV_S,U(0)%X1(:,iMode))
      END DO
      DO iMode=1,X2_base%f%modes
        X2_s(iMode) = X2_base%s%evalDOF_s(spos,      0,U(0)%X2(:,iMode))
        F_X2_s( iMode)= X2_base%s%evalDOF_s(spos,      0,F(0)%X2(:,iMode))
        dX2ds(iMode)= X2_base%s%evalDOF_s(spos,DERIV_S,U(0)%X2(:,iMode))
      END DO
      DO iMode=1,LA_base%f%modes
        LA_s(iMode) = LA_base%s%evalDOF_s(spos,      0,U(0)%LA(:,iMode))
        F_LA_s( iMode)= LA_base%s%evalDOF_s(spos,      0,F(0)%LA(:,iMode))
      END DO
      iota_s=Eval_iota(spos)
      pres_s=Eval_pres(spos)
      phiPrime_s=Eval_PhiPrime(spos)
      var_visu(3,i_s,:,:,iElem) =Eval_Phi(spos)
      var_visu(4,i_s,:,:,iElem) =iota_s
      var_visu(5,i_s,:,:,iElem) =pres_s
      DO i_n=1,mn_IP(2)
        xIP(2)  = TWOPI*(minmax(3,0)+(minmax(3,1)-minmax(3,0))*REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp))
        DO i_m=1,mn_IP(1)
          xIP(1)= TWOPI*(minmax(2,0)+(minmax(2,1)-minmax(2,0))*REAL(i_m-1,wp)/REAL(mn_IP(1)-1,wp))

          ASSOCIATE(lambda_visu  => var_visu(  1,i_s,i_m,i_n,iElem), &
                    sqrtG_visu   => var_visu(  2,i_s,i_m,i_n,iElem), &
                    Bvec_visu    => var_visu(6:8,i_s,i_m,i_n,iElem), &
                    F_X1_visu    => var_visu(  9,i_s,i_m,i_n,iElem), &
                    F_X2_visu    => var_visu( 10,i_s,i_m,i_n,iElem), &
                    F_LA_visu    => var_visu( 11,i_s,i_m,i_n,iElem)  )


          X1_visu         =X1_base%f%evalDOF_x(xIP,         0,X1_s )
          X2_visu         =X2_base%f%evalDOF_x(xIP,         0,X2_s )

          F_X1_visu       =X1_base%f%evalDOF_x(xIP,         0,F_X1_s )
          F_X2_visu       =X2_base%f%evalDOF_x(xIP,         0,F_X2_s )
          F_LA_visu       =LA_base%f%evalDOF_x(xIP,         0,F_LA_s )

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
          coord_visu(:,i_s,i_m,i_n,iElem )=hmap%eval(q)
          !lambda
          lambda_visu = LA_base%f%evalDOF_x(xIP,0,LA_s)
          !sqrtG
          sqrtG_visu  = hmap%eval_Jh(q)*(dX1_ds_visu*dX2_dthet_visu -dX2_ds_visu*dX1_dthet_visu) 

          e_thet=hmap%eval_dxdq(q,(/dX1_dthet_visu,dX2_dthet_visu,0.0_wp/))
          e_zeta=hmap%eval_dxdq(q,(/dX1_dzeta_visu,dX2_dzeta_visu,1.0_wp/))
          var_visu(6:8,i_s,i_m,i_n,iElem)= &
!          Bvec_visu(:)= & 
                       (  e_thet(:)*(iota_s-dLA_dzeta_visu)  &
                        + e_zeta(:)*(1.0_wp+dLA_dthet_visu) )*(PhiPrime_s/MAX(1.0e-12_wp,sqrtG_visu))
          END ASSOCIATE !lambda,sqrtG,Bvec,F_X1/X2/LA
        END DO !i_m
      END DO !i_n
    END DO !i_s
  END DO !iElem
  
  !range s: include all elements belonging to [smin,smax]
  minElem=MAX(     1,sgrid%find_elem(minmax(1,0))-1)
  maxElem=MIN(nElems,sgrid%find_elem(minmax(1,1))+1)
  IF(only_planes)THEN
    nplot(1:2)=(/n_s,mn_IP(1)/)-1
    WRITE(filename,'(A,"_visu_planes_",I4.4,"_",I8.8,".vtu")')TRIM(Projectname),outputLevel,fileID
    CALL WriteDataToVTK(2,3,nVal,nplot(1:2),(mn_IP(2)*(maxElem-minElem+1)),VarNames, &
                        coord_visu(:,:,:,:,minElem:maxElem), &
                          var_visu(:,:,:,:,minElem:maxElem),TRIM(filename))
  ELSE
    !3D
    nplot(1:3)=(/n_s,mn_IP(1),mn_IP(2)/)-1
    WRITE(filename,'(A,"_visu_3D_",I4.4,"_",I8.8,".vtu")')TRIM(Projectname),outputLevel,fileID
    CALL WriteDataToVTK(3,3,nVal,nplot,(maxElem-minElem+1),VarNames, &
                        coord_visu(:,:,:,:,minElem:maxElem), &
                          var_visu(:,:,:,:,minElem:maxElem),TRIM(filename))
  END IF
  
  END ASSOCIATE!n_s,mn_IP


END SUBROUTINE visu_3D

!===================================================================================================================================
!> Visualize 
!!
!===================================================================================================================================
SUBROUTINE visu_1d_modes(n_s,fileID)
! MODULES
USE MOD_Analyze_Vars,  ONLY: visu1D
USE MOD_fbase,         ONLY: sin_cos_map
USE MOD_MHD3D_Vars,    ONLY: U,X1_base,X2_base,LA_base
USE MOD_Output_vars,   ONLY: outputLevel
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
    WRITE(*,*)'1) Visualize gvec modes in 1D: R,Z,lambda interpolated...'
    vname="X1"//TRIM(sin_cos_map(X1_base%f%sin_cos))
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,X1_base,U(0)%X1)
    vname="X2"//TRIM(sin_cos_map(X2_base%f%sin_cos))
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,X2_base,U(0)%X2)
    vname="LA"//TRIM(sin_cos_map(LA_base%f%sin_cos))
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,LA_base,U(0)%LA)
  END IF
  IF(vcase(2))THEN
    WRITE(*,*)'2) Visualize gvec modes in 1D: dRrho,dZrho interpolated...'
    vname="dX1"//TRIM(sin_cos_map(X1_base%f%sin_cos))
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,X1_base,U(0)%X1)
    vname="dX2"//TRIM(sin_cos_map(X2_base%f%sin_cos))
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,X2_base,U(0)%X2)
  END IF
  
  !
  
END SUBROUTINE visu_1d_modes


!===================================================================================================================================
!> Write all modes of one variable
!!
!===================================================================================================================================
SUBROUTINE writeDataMN_visu(n_s,fname_in,vname,rderiv,base_in,xx_in)
! MODULES
  USE MOD_base,          ONLY: t_base
  USE MOD_MHD3D_Profiles,ONLY: Eval_iota,Eval_pres,Eval_Phi
  USE MOD_MHD3D_Vars,    ONLY: sgrid
  USE MOD_write_modes,   ONLY: write_modes
  USE MOD_output_vars,   ONLY: Projectname
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

  DO iMode=1,base_in%f%modes
    nVal=nVal+1
    WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname),base_in%f%Xmn(1,iMode),base_in%f%Xmn(2,iMode)/base_in%f%nfp
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

END MODULE MOD_MHD3D_visu

