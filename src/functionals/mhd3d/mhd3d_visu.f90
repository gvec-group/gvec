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
  INTEGER  :: i_m,i_n,nplot(2),iMode,iVal
  REAL(wp) :: xIP(2),q(3)
  REAL(wp) :: X1_visu,X2_visu
  REAL(wp) :: coord_visu(3,mn_IP(1),mn_IP(2),1)
  INTEGER,PARAMETER  :: nVal=3
  INTEGER  :: VP_LAMBDA,VP_theta,VP_zeta
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
  ival=1
  VP_LAMBDA=iVal;iVal=iVal+1; VarNames(VP_LAMBDA)="lambda"
  VP_theta =iVal;iVal=iVal+1; VarNames(VP_theta )="theta"
  VP_zeta  =iVal;iVal=iVal+1; VarNames(VP_zeta  )="zeta"
  DO i_m=1,mn_IP(1)
    thet(i_m)= TWOPI*(minmax(2,0)+(minmax(2,1)-minmax(2,0))*REAL(i_m-1,wp)/REAL(mn_IP(1)-1,wp)) !repeat point exactly
  END DO
  DO i_n=1,mn_IP(2)
    zeta(i_n)=TWOPI*(minmax(3,0)+(minmax(3,1)-minmax(3,0))*REAL(i_n-1,wp)/REAL(mn_IP(2)-1,wp))
  END DO
  !make theta direction fully periodic
    IF(ABS((minMax(2,1)-minmax(2,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
      thet(mn_IP(1))=thet(1)
    END IF
  IF(hmap%which_hmap.NE.3)THEN !not for cylinder
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
      var_visu(VP_LAMBDA,i_m,i_n,1)=LA_base%f%evalDOF_x(xIP,0,LA_s)
      var_visu(VP_theta ,i_m,i_n,1)=xIP(1)
      var_visu(VP_zeta  ,i_m,i_n,1)=xIP(2)
    END DO !i_m
  END DO !i_n
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
USE MODgvec_MHD3D_Profiles, ONLY: Eval_iota,Eval_pres,Eval_Phi,Eval_PhiPrime,Eval_chiPrime,Eval_p_prime
USE MODgvec_MHD3D_Profiles, ONLY: Eval_iota_Prime,Eval_Phi_TwoPrime
USE MODgvec_output_vtk,     ONLY: WriteDataToVTK
USE MODgvec_output_netcdf,  ONLY: WriteDataToNETCDF
USE MODgvec_Output_CSV,     ONLY: WriteDataToCSV
USE MODgvec_Output_vars,    ONLY: Projectname,OutputLevel
USE MODgvec_Analyze_Vars,   ONLY: SFL_theta,outfileType
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
  INTEGER  :: i_s,j_s,i_m,i_n,iElem,nElems,nplot(3),minElem,maxElem,n_s,mn_IP(2),ival,i,j
  REAL(wp) :: spos,xIP(2),q(3),q_thet(3),q_zeta(3)
  REAL(wp) :: X1_s(X1_base%f%modes),F_X1_s(X1_base%f%modes),dX1ds(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes),F_X2_s(X2_base%f%modes),dX2ds(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes),F_LA_s(LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu,dX1_ds,dX2_ds,dX1_dthet,dX1_dzeta,dX2_dthet,dX2_dzeta
  REAL(wp) :: dLA_dthet,dLA_dzeta,iota_s,pres_s,phiPrime_s,e_s(3),e_thet(3),e_zeta(3)
#if (defined(VISU_J_FD) || defined(VISU_J_EXACT))
  INTEGER,PARAMETER  :: nVal=34
  INTEGER            :: VP_J
#else
  INTEGER,PARAMETER  :: nVal=31
#endif
  INTEGER  ::VP_LAMBDA,VP_SQRTG,VP_PHI,VP_IOTA,VP_PRES,VP_dp_ds,VP_B,VP_F_X1,VP_F_X2,VP_F_LA, &
             VP_s,VP_theta,VP_zeta,VP_g_tt,VP_g_tz,VP_g_zz,VP_gr_s,VP_gr_t,VP_gr_z,VP_Mscale ,VP_MscaleF,&
             VP_Ipol,VP_Itor
  REAL(wp) :: coord_visu( 3,np_in(1),np_in(1),np_in(3),np_in(2),sgrid%nElems)
  REAL(wp) :: var_visu(nVal,np_in(1),np_in(1),np_in(3),np_in(2),sgrid%nElems)
  REAL(wp) :: var_visu_1d(nVal+3,(np_in(1)-1)*sgrid%nElems+1)
  REAL(wp) :: thet(np_in(1),np_in(2)),zeta(np_in(3))
  REAL(wp) :: theta_star,sqrtG
  CHARACTER(LEN=40) :: CoordNames(3)
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
  CHARACTER(LEN=255) :: filename
  REAL(wp) :: Bthet, Bzeta,Bcart(3),Ipol_int,Itor_int
  REAL(wp) :: grad_s(3), grad_thet(3),grad_zeta(3)
#ifdef VISU_J_FD
  REAL(wp) :: xIP_eps(2)
  REAL(wp) :: X1_s_eps(X1_base%f%modes),dX1ds_eps(X1_base%f%modes)
  REAL(wp) :: X2_s_eps(X2_base%f%modes),dX2ds_eps(X2_base%f%modes)
  REAL(wp) :: LA_s_eps(LA_base%f%modes)
  REAL(wp) :: X1_eps,X2_eps,dX1_ds_eps,dX2_ds_eps,dX1_dthet_eps,dX1_dzeta_eps,dX2_dthet_eps,dX2_dzeta_eps
  REAL(wp) :: dLA_dthet_eps,dLA_dzeta_eps,iota_s_eps,pres_s_eps,phiPrime_s_eps
  REAL(wp) :: B_ds(3), B_dthet(3), B_dzeta(3), grad_Bcart(3, 3)          !< cartesion current density and gradient of magnetic field components
  INTEGER  :: sgn
  REAL(wp) :: delta_s,delta_thet,delta_zeta
  REAL(wp),PARAMETER :: eps   = 1.0e-8 !theta,zeta
  REAL(wp),PARAMETER :: eps_s   = 1.0e-4 !
#endif
#ifdef VISU_J_EXACT
  REAL(wp) :: dX1ds_ds(X1_base%f%modes),dX2ds_ds(X2_base%f%modes),dLAds(LA_base%f%modes)
  REAL(wp) :: phiPrime_s_s,iota_s_s
  REAL(wp) :: dX1_ds_ds,dX1_ds_dthet,dX1_ds_dzeta,dX1_dthet_dthet,dX1_dthet_dzeta,dX1_dzeta_dzeta
  REAL(wp) :: dX2_ds_ds,dX2_ds_dthet,dX2_ds_dzeta,dX2_dthet_dthet,dX2_dthet_dzeta,dX2_dzeta_dzeta 
  REAL(wp) ::           dLA_ds_dthet,dLA_ds_dzeta,dLA_dthet_dthet,dLA_dthet_dzeta,dLA_dzeta_dzeta 
  REAL(wp) :: dBthet_ds,dBthet_dthet,dBthet_dzeta,dBzeta_ds,dBzeta_dthet,dBzeta_dzeta
  REAL(wp),DIMENSION(3)::q_s,q_s_s,q_s_thet,q_s_zeta,q_thet_thet,q_thet_zeta,q_zeta_zeta
  REAL(wp) :: Jh,Jh_dq1,Jh_dq2,dJh_ds,dJh_dthet,dJh_dzeta
  REAL(wp) :: Jp              ,dJp_ds,dJp_dthet,dJp_dzeta
  REAL(wp) :: dsqrtg_ds,dsqrtg_dthet,dsqrtg_dzeta
  REAL(wp) :: g_st,g_st_dq1,g_st_dq2         ,dg_st_dthet,dg_st_dzeta
  REAL(wp) :: g_sz,g_sz_dq1,g_sz_dq2         ,dg_sz_dthet,dg_sz_dzeta
  REAL(wp) :: g_tt,g_tt_dq1,g_tt_dq2,dg_tt_ds,dg_tt_dzeta
  REAL(wp) :: g_tz,g_tz_dq1,g_tz_dq2,dg_tz_ds,dg_tz_dthet,dg_tz_dzeta
  REAL(wp) :: g_zz,g_zz_dq1,g_zz_dq2,dg_zz_ds            ,dg_zz_dthet
  REAL(wp) :: dBsubs_dthet,dBsubs_dzeta,dBsubthet_ds,dBsubthet_dzeta,dBsubzeta_ds,dBsubzeta_dthet

  REAL(wp) :: Js,Jthet,Jzeta
#endif
  REAL(wp) :: Jcart(3)
  REAL(wp),ALLOCATABLE :: tmpcoord(:,:,:,:),tmpvar(:,:,:,:) 
!===================================================================================================================================
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
  __PERFON("output_visu")
  __PERFON("prepare_visu")
  iVal=1
  VP_s      =iVal;iVal=iVal+1; VarNames(VP_s     )="s"
  VP_theta  =iVal;iVal=iVal+1; VarNames(VP_theta )="theta"
  VP_zeta   =iVal;iVal=iVal+1; VarNames(VP_zeta  )="zeta"
  VP_PHI    =iVal;iVal=iVal+1; VarNames(VP_PHI   )="Phi"
  VP_IOTA   =iVal;iVal=iVal+1; VarNames(VP_IOTA  )="iota"
  VP_PRES   =iVal;iVal=iVal+1; VarNames(VP_PRES  )="pressure"
  VP_DP_DS  =iVal;iVal=iVal+1; VarNames(VP_DP_DS )="dp_ds"
  VP_Mscale =iVal;iVal=iVal+1; VarNames(VP_Mscale)="Mscale"
  VP_MscaleF=iVal;iVal=iVal+1; VarNames(VP_MscaleF)="MscaleForce"
  VP_LAMBDA =iVal;iVal=iVal+1; VarNames(VP_LAMBDA)="lambda" 
  VP_SQRTG  =iVal;iVal=iVal+1; VarNames(VP_SQRTG )="sqrtG"
  VP_g_tt   =iVal;iVal=iVal+1; VarNames(VP_g_tt  )="g_tt"
  VP_g_tz   =iVal;iVal=iVal+1; VarNames(VP_g_tz  )="g_tz"
  VP_g_zz   =iVal;iVal=iVal+1; VarNames(VP_g_zz  )="g_zz"
  VP_B      =iVal;iVal=iVal+3; VarNames(VP_B     )="BvecX"
                               VarNames(VP_B+1   )="BvecY"
                               VarNames(VP_B+2   )="BvecZ"
  VP_F_X1   =iVal;iVal=iVal+1; VarNames(VP_F_X1  )="F_X1"
  VP_F_X2   =iVal;iVal=iVal+1; VarNames(VP_F_X2  )="F_X2"
  VP_F_LA   =iVal;iVal=iVal+1; VarNames(VP_F_LA  )="F_LA"
  VP_gr_s   =iVal;iVal=iVal+3; VarNames(VP_gr_s  )="grad_sX"
                               VarNames(VP_gr_s+1)="grad_sY"
                               VarNames(VP_gr_s+2)="grad_sZ"
  VP_gr_t   =iVal;iVal=iVal+3; VarNames(VP_gr_t  )="grad_tX"
                               VarNames(VP_gr_t+1)="grad_tY"
                               VarNames(VP_gr_t+2)="grad_tZ"
  VP_gr_z   =iVal;iVal=iVal+3; VarNames(VP_gr_z  )="grad_zX"
                               VarNames(VP_gr_z+1)="grad_zY"
                               VarNames(VP_gr_z+2)="grad_zZ"
  VP_Ipol   =iVal;iVal=iVal+1; VarNames(VP_Ipol  )="Ipol"
  VP_Itor   =iVal;iVal=iVal+1; VarNames(VP_Itor  )="Itor"
#if (defined(VISU_J_FD) || defined(VISU_J_EXACT))
  VP_J      =iVal;iVal=iVal+3; VarNames(VP_J   )="JvecX"
                               VarNames(VP_J+1 )="JvecY"
                               VarNames(VP_J+2 )="JvecZ"
#endif

  var_visu=0.

  n_s=np_in(1)
  mn_IP=np_in(2:3)

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
!      spos=MAX(1.0e-04,sgrid%sp(iElem-1)+1e-08+(REAL(i_s-1,wp))/(REAL(n_s-1,wp))*(sgrid%ds(iElem)-2*1e-8)) !for discont. data
      spos=MAX(1.0e-04,sgrid%sp(iElem-1)+(REAL(i_s-1,wp))/(REAL(n_s-1,wp))*(sgrid%ds(iElem)))
 
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
#ifdef VISU_J_EXACT
      dX1ds_ds(:)  = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,DERIV_S_S,U(0)%X1(:,:))
      dX2ds_ds(:)  = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,DERIV_S_S,U(0)%X2(:,:))
      dLAds(:)     = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,DERIV_S  ,U(0)%LA(:,:))
      iota_s_s=Eval_iota_Prime(spos)
      phiPrime_s_s=Eval_Phi_TwoPrime(spos)
#endif
      var_visu(VP_S    ,i_s,:,:,:,iElem) =spos
      var_visu(VP_PHI  ,i_s,:,:,:,iElem) =Eval_Phi(spos)
      var_visu(VP_IOTA ,i_s,:,:,:,iElem) =iota_s
      var_visu(VP_PRES ,i_s,:,:,:,iElem) =pres_s
      var_visu(VP_DP_DS,i_s,:,:,:,iElem) =Eval_p_prime(spos)
      var_visu(VP_Mscale,i_s,:,:,:,iElem) = (SUM(X1_base%f%Xmn(1,:)**(4+1)*X1_s(:)**2)+SUM(X2_base%f%Xmn(1,:)**(4+1)*X2_s(:)**2))/&  !pexp=4, qexp=1
                                            (SUM(X1_base%f%Xmn(1,:)**(4  )*X1_s(:)**2)+SUM(X2_base%f%Xmn(1,:)**(4  )*X2_s(:)**2))
      var_visu(VP_MscaleF,i_s,:,:,:,iElem)= (SUM(X1_base%f%Xmn(1,:)**(4+1)*F_X1_s(:)**2)+SUM(X2_base%f%Xmn(1,:)**(4+1)*F_X2_s(:)**2))/&  !pexp=4, qexp=1
                                            (SUM(X1_base%f%Xmn(1,:)**(4  )*F_X1_s(:)**2)+SUM(X2_base%f%Xmn(1,:)**(4  )*F_X2_s(:)**2)+1.0e-14)
#ifdef VISU_J_FD
      ! for Finite  Difference in s
      if (i_s .ne. n_s) then !switch sign of finite difference at last point
        sgn = 1
      else
        sgn = -1
      endif
      delta_s=sgn*eps_s*sgrid%ds(iElem)
      X1_s_eps(:)   = X1_base%s%evalDOF2D_s(spos+delta_s,X1_base%f%modes,      0,U(0)%X1(:,:))
      dX1ds_eps(:)  = X1_base%s%evalDOF2D_s(spos+delta_s,X1_base%f%modes,DERIV_S,U(0)%X1(:,:))
      X2_s_eps(:)   = X2_base%s%evalDOF2D_s(spos+delta_s,X2_base%f%modes,      0,U(0)%X2(:,:))
      dX2ds_eps(:)  = X2_base%s%evalDOF2D_s(spos+delta_s,X2_base%f%modes,DERIV_S,U(0)%X2(:,:))
      LA_s_eps(:)   = LA_base%s%evalDOF2D_s(spos+delta_s,LA_base%f%modes,      0,U(0)%LA(:,:))
      iota_s_eps=Eval_iota(spos+delta_s)
      pres_s_eps=Eval_pres(spos+delta_s)
      phiPrime_s_eps=Eval_PhiPrime(spos+delta_s)
#endif
      !define theta2, which corresponds to the theta angle of a given theta_star=theta
      Itor_int = 0.
      Ipol_int = 0.

!$OMP PARALLEL DO COLLAPSE(3)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_m,i_n,j_s,xIP,q,q_thet,q_zeta,sqrtG,e_s,e_thet,e_zeta,theta_star, &
!$OMP           X1_visu,dX1_ds,dX1_dthet,dX1_dzeta,    &
!$OMP           X2_visu,dX2_ds,dX2_dthet,dX2_dzeta,dLA_dthet,dLA_dzeta, &
#ifdef VISU_J_FD
!$OMP           X1_eps,dX1_ds_eps,dX1_dthet_eps,dX1_dzeta_eps,    &
!$OMP           X2_eps,dX2_ds_eps,dX2_dthet_eps,dX2_dzeta_eps,dLA_dthet_eps,dLA_dzeta_eps, &
!$OMP           JCart, B_ds, B_dthet, B_dzeta, grad_Bcart, xIP_eps, delta_thet, delta_zeta,&
#endif
#ifdef VISU_J_EXACT
!$OMP           dX1_ds_ds,dX1_ds_dthet,dX1_ds_dzeta,dX1_dthet_dthet,dX1_dthet_dzeta,dX1_dzeta_dzeta, &
!$OMP           dX2_ds_ds,dX2_ds_dthet,dX2_ds_dzeta,dX2_dthet_dthet,dX2_dthet_dzeta,dX2_dzeta_dzeta, &
!$OMP                     dLA_ds_dthet,dLA_ds_dzeta,dLA_dthet_dthet,dLA_dthet_dzeta,dLA_dzeta_dzeta, &
!$OMP           dBthet_ds,dBthet_dthet,dBthet_dzeta,dBzeta_ds,dBzeta_dthet,dBzeta_dzeta, &
!$OMP           q_s,q_s_s,q_s_thet,q_s_zeta,q_thet_thet,q_thet_zeta,q_zeta_zeta, &
!$OMP           Jh,Jh_dq1,Jh_dq2,dJh_ds,dJh_dthet,dJh_dzeta, &
!$OMP           Jp              ,dJp_ds,dJp_dthet,dJp_dzeta, &
!$OMP           dsqrtg_ds,dsqrtg_dthet,dsqrtg_dzeta, &
!$OMP           g_st,g_st_dq1,g_st_dq2         ,dg_st_dthet,dg_st_dzeta, &
!$OMP           g_sz,g_sz_dq1,g_sz_dq2         ,dg_sz_dthet,dg_sz_dzeta, &
!$OMP           g_tt,g_tt_dq1,g_tt_dq2,dg_tt_ds,dg_tt_dzeta, &
!$OMP           g_tz,g_tz_dq1,g_tz_dq2,dg_tz_ds,dg_tz_dthet,dg_tz_dzeta, &
!$OMP           g_zz,g_zz_dq1,g_zz_dq2,dg_zz_ds            ,dg_zz_dthet, &
!$OMP           dBsubs_dthet,dBsubs_dzeta,dBsubthet_ds,dBsubthet_dzeta,dBsubzeta_ds,dBsubzeta_dthet, &
!$OMP           Js,Jthet,Jzeta,Jcart, &
#endif
!$OMP           Bcart, Bthet, Bzeta, grad_s, grad_thet, grad_zeta) &
!$OMP   REDUCTION(+:Itor_int,Ipol_int) &
!$OMP   SHARED(np_in,i_s,iElem,thet,zeta,SFL_theta,X1_base,X2_base,LA_base,X1_s,X2_s,LA_s,dX1ds,dX2ds,&
!$OMP          VP_LAMBDA,VP_SQRTG,VP_B,VP_F_X1,VP_F_X2,VP_F_LA, VP_Ipol,VP_Itor,&
!$OMP          VP_theta,VP_zeta,VP_g_tt,VP_g_tz,VP_g_zz,VP_gr_s,VP_gr_t,VP_gr_z,iota_s, &
#ifdef VISU_J_FD
!$OMP          X1_s_eps,X2_s_eps,LA_s_eps,dX1ds_eps,dX2ds_eps,VP_J,iota_s_eps,PhiPrime_s_eps,delta_s,&
#endif
#ifdef VISU_J_EXACT
!$OMP          dX1ds_ds,dX2ds_ds,dLAds,iota_s_s,phiPrime_s_s,VP_J, &
#endif
!$OMP          F_X1_s,F_X2_s,F_LA_s,hmap,coord_visu,var_visu,phiPrime_s,mn_IP,n_s)
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
            var_visu(VP_theta,i_s,j_s,i_n,i_m,iElem)=xIP(1) !theta for evaluation of X1,X2,LA
            var_visu(VP_zeta ,i_s,j_s,i_n,i_m,iElem)=xIP(2) !zeta  for evaluation of X1,X2,LA

            X1_visu         =X1_base%f%evalDOF_x(xIP,         0,X1_s )
            X2_visu         =X2_base%f%evalDOF_x(xIP,         0,X2_s )
            
            dX1_ds     =X1_base%f%evalDOF_x(xIP,         0,dX1ds)
            dX2_ds     =X2_base%f%evalDOF_x(xIP,         0,dX2ds)
            
            dX1_dthet  =X1_base%f%evalDOF_x(xIP,DERIV_THET,X1_s )
            dX2_dthet  =X2_base%f%evalDOF_x(xIP,DERIV_THET,X2_s )
            dLA_dthet  =LA_base%f%evalDOF_x(xIP,DERIV_THET,LA_s )
            
            dX1_dzeta  =X1_base%f%evalDOF_x(xIP,DERIV_ZETA,X1_s )
            dX2_dzeta  =X2_base%f%evalDOF_x(xIP,DERIV_ZETA,X2_s )
            dLA_dzeta  =LA_base%f%evalDOF_x(xIP,DERIV_ZETA,LA_s )
            
            q=(/X1_visu,X2_visu,xIP(2)/)
            q_thet=(/dX1_dthet,dX2_dthet,0.0_wp/)
            q_zeta=(/dX1_dzeta,dX2_dzeta,1.0_wp/)
            coord_visu(:,i_s,j_s,i_n,i_m,iElem )=hmap%eval(q)

            e_s   =hmap%eval_dxdq(q,(/dX1_ds,dX2_ds,0.0_wp/))
            e_thet=hmap%eval_dxdq(q,q_thet)
            e_zeta=hmap%eval_dxdq(q,q_zeta)

           !sqrtG = hmap%eval_Jh(q)*(dX1_ds*dX2_dthet -dX2_ds*dX1_dthet) 
            sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))
            !IF(ABS(sqtG- SUM(e_s*(CROSS(e_thet,e_zeta)))).GT.1.0e-04) STOP 'test sqrtg failed'
            
            ! Get contra-variant basis vectors
            grad_s    = CROSS(e_thet,e_zeta) /sqrtG
            grad_thet = CROSS(e_zeta,e_s   ) /sqrtG
            grad_zeta = CROSS(e_s   ,e_thet) /sqrtG
            var_visu(VP_gr_s:VP_gr_s+2,i_s,j_s,i_n,i_m,iElem) = grad_s
            var_visu(VP_gr_t:VP_gr_t+2,i_s,j_s,i_n,i_m,iElem) = grad_thet
            var_visu(VP_gr_z:VP_gr_z+2,i_s,j_s,i_n,i_m,iElem) = grad_zeta

            var_visu(VP_g_tt ,i_s,j_s,i_n,i_m,iElem) = hmap%eval_gij(q_thet,q,q_thet)   !g_theta,theta
            var_visu(VP_g_tz ,i_s,j_s,i_n,i_m,iElem) = hmap%eval_gij(q_thet,q,q_zeta)   !g_theta,zeta =g_zeta,theta
            var_visu(VP_g_zz ,i_s,j_s,i_n,i_m,iElem) = hmap%eval_gij(q_zeta,q,q_zeta)   !g_zeta,zeta
      
            !lambda
            var_visu(VP_LAMBDA,i_s,j_s,i_n,i_m,iElem) = LA_base%f%evalDOF_x(xIP,0,LA_s)
            !sqrtG
            var_visu(VP_SQRTG,i_s,j_s,i_n,i_m,iElem) = sqrtG
            !F_X1,F_X2,F_LA  
            var_visu(VP_F_X1,i_s,j_s,i_n,i_m,iElem) = X1_base%f%evalDOF_x(xIP,         0,F_X1_s )
            var_visu(VP_F_X2,i_s,j_s,i_n,i_m,iElem) = X2_base%f%evalDOF_x(xIP,         0,F_X2_s )
            var_visu(VP_F_LA,i_s,j_s,i_n,i_m,iElem) = LA_base%f%evalDOF_x(xIP,         0,F_LA_s )
            !Bvec
            Bthet   = (iota_s - dLA_dzeta ) * phiPrime_s   !/sqrtG
            Bzeta   = (1.0_wp + dLA_dthet ) * phiPrime_s       !/sqrtG
            Bcart(:) =  ( e_thet(:) * Bthet + e_zeta(:) * Bzeta) /sqrtG

            var_visu(VP_B:VP_B+2,i_s,j_s,i_n,i_m,iElem)= Bcart(:)
            !poloidal and toroidal current profiles, line integral: integration over one angle /average over other...
            !Itor= int_0^2pi B_theta dtheta = (nfp/2pi) int_0^2pi int_0^(2pi/nfp) B_theta dtheta dzeta
            !Ipol= int_0^2pi B_zeta  dzeta  = nfp* int_0^(2pi/nfp) B_zeta dzeta = (nfp/2pi) int_0^2pi int_0^(2pi/nfp) B_zeta  dtheta dzeta
            Itor_int = Itor_int+ SUM(Bcart(:)*e_thet(:))   !B_theta=B.e_thet 
            Ipol_int = Ipol_int+ SUM(Bcart(:)*e_zeta(:))   !B_zeta =B.e_zeta

           ! Get J components:
            
#ifdef VISU_J_EXACT
            dX1_ds_ds       = X1_base%f%evalDOF_x(xIP,         0,dX1ds_ds)
            dX2_ds_ds       = X2_base%f%evalDOF_x(xIP,         0,dX2ds_ds)

            dX1_ds_dthet    = X1_base%f%evalDOF_x(xIP,DERIV_THET,dX1ds )
            dX2_ds_dthet    = X2_base%f%evalDOF_x(xIP,DERIV_THET,dX2ds )
            dLA_ds_dthet    = LA_base%f%evalDOF_x(xIP,DERIV_THET,dLAds )

            dX1_ds_dzeta    = X1_base%f%evalDOF_x(xIP,DERIV_ZETA,dX1ds )
            dX2_ds_dzeta    = X2_base%f%evalDOF_x(xIP,DERIV_ZETA,dX2ds )
            dLA_ds_dzeta    = LA_base%f%evalDOF_x(xIP,DERIV_ZETA,dLAds )

            dX1_dthet_dthet = X1_base%f%evalDOF_x(xIP,DERIV_THET_THET,X1_s )
            dX2_dthet_dthet = X2_base%f%evalDOF_x(xIP,DERIV_THET_THET,X2_s )
            dLA_dthet_dthet = LA_base%f%evalDOF_x(xIP,DERIV_THET_THET,LA_s )

            dX1_dthet_dzeta = X1_base%f%evalDOF_x(xIP,DERIV_THET_ZETA,X1_s )
            dX2_dthet_dzeta = X2_base%f%evalDOF_x(xIP,DERIV_THET_ZETA,X2_s )
            dLA_dthet_dzeta = LA_base%f%evalDOF_x(xIP,DERIV_THET_ZETA,LA_s )

            dX1_dzeta_dzeta = X1_base%f%evalDOF_x(xIP,DERIV_ZETA_ZETA,X1_s )
            dX2_dzeta_dzeta = X2_base%f%evalDOF_x(xIP,DERIV_ZETA_ZETA,X2_s )
            dLA_dzeta_dzeta = LA_base%f%evalDOF_x(xIP,DERIV_ZETA_ZETA,LA_s )


            dBthet_ds    =  (iota_s - dLA_dzeta ) * phiPrime_s_s + (iota_s_s - dLA_ds_dzeta ) * phiPrime_s 
            dBthet_dthet =                                         (      - dLA_dthet_dzeta ) * phiPrime_s 
            dBthet_dzeta =                                         (      - dLA_dzeta_dzeta ) * phiPrime_s 
            dBzeta_ds    =  (1.0_wp + dLA_dthet ) * phiPrime_s_s +             dLA_ds_dthet   * phiPrime_s 
            dBzeta_dthet =                                                  dLA_dthet_dthet   * phiPrime_s 
            dBzeta_dzeta =                                                  dLA_dthet_dzeta   * phiPrime_s 


            
            q_s          = (/dX1_ds         ,dX2_ds         ,0.0_wp/)
            q_s_s        = (/dX1_ds_ds      ,dX2_ds_ds      ,0.0_wp/)
            q_s_thet     = (/dX1_ds_dthet   ,dX2_ds_dthet   ,0.0_wp/)
            q_s_zeta     = (/dX1_ds_dzeta   ,dX2_ds_dzeta   ,0.0_wp/)
            q_thet_thet  = (/dX1_dthet_dthet,dX2_dthet_dthet,0.0_wp/)
            q_thet_zeta  = (/dX1_dthet_dzeta,dX2_dthet_dzeta,0.0_wp/)
            q_zeta_zeta  = (/dX1_dzeta_dzeta,dX2_dzeta_dzeta,0.0_wp/)

            Jh           = hmap%eval_Jh(q)
            Jh_dq1       = hmap%eval_Jh_dq1(q)
            Jh_dq2       = hmap%eval_Jh_dq2(q)
            dJh_ds       = q_s(   1)*Jh_dq1 + q_s(   2)*Jh_dq2
            dJh_dthet    = q_thet(1)*Jh_dq1 + q_thet(2)*Jh_dq2
            dJh_dzeta    = q_zeta(1)*Jh_dq1 + q_zeta(2)*Jh_dq2

            Jp           = q_s(     1)*q_thet(     2)-q_s(     2)*q_thet(     1)
            dJp_ds       = q_s_s(   1)*q_thet(     2)-q_s_s(   2)*q_thet(     1) &
                          +q_s(     1)*q_s_thet(   2)-q_s(     2)*q_s_thet(   1)
            dJp_dthet    = q_s_thet(1)*q_thet(     2)-q_s_thet(2)*q_thet(     1) &
                          +q_s(     1)*q_thet_thet(2)-q_s(     2)*q_thet_thet(1)
            dJp_dzeta    = q_s_zeta(1)*q_thet(     2)-q_s_zeta(2)*q_thet(     1) &
                          +q_s(     1)*q_thet_zeta(2)-q_s(     2)*q_thet_zeta(1) 
             
            sqrtg        = Jh*Jp
            dsqrtg_ds    = Jh*dJp_ds    + dJh_ds   *Jp 
            dsqrtg_dthet = Jh*dJp_dthet + dJh_dthet*Jp
            dsqrtg_dzeta = Jh*dJp_dzeta + dJh_dzeta*Jp

            g_st         = hmap%eval_gij(     q_s,q,q_thet)
            g_st_dq1     = hmap%eval_gij_dq1( q_s,q,q_thet)
            g_st_dq2     = hmap%eval_gij_dq2( q_s,q,q_thet)
            dg_st_dthet  = hmap%eval_gij(q_s_thet,q,q_thet)     &
                          +hmap%eval_gij(     q_s,q,q_thet_thet) + q_thet(1)*g_st_dq1 + q_thet(2)*g_st_dq2 
            dg_st_dzeta  = hmap%eval_gij(q_s_zeta,q,q_thet)     &
                          +hmap%eval_gij(     q_s,q,q_thet_zeta) + q_zeta(1)*g_st_dq1 + q_zeta(2)*g_st_dq2 

            g_sz         = hmap%eval_gij(     q_s,q,q_zeta)
            g_sz_dq1     = hmap%eval_gij_dq1( q_s,q,q_zeta)
            g_sz_dq2     = hmap%eval_gij_dq2( q_s,q,q_zeta)
            dg_sz_dthet  = hmap%eval_gij(q_s_thet,q,q_zeta)     &
                          +hmap%eval_gij(     q_s,q,q_thet_zeta) + q_thet(1)*g_sz_dq1 + q_thet(2)*g_sz_dq2 
            dg_sz_dzeta  = hmap%eval_gij(q_s_zeta,q,q_zeta)     & 
                          +hmap%eval_gij(     q_s,q,q_zeta_zeta) + q_zeta(1)*g_sz_dq1 + q_zeta(2)*g_sz_dq2 

            g_tt         = hmap%eval_gij(            q_thet,q,q_thet)
            g_tt_dq1     = hmap%eval_gij_dq1(        q_thet,q,q_thet)
            g_tt_dq2     = hmap%eval_gij_dq2(        q_thet,q,q_thet)
            dg_tt_ds     = 2.0_wp*hmap%eval_gij(q_s_thet   ,q,q_thet) + q_s(   1)*g_tt_dq1 + q_s(   2)*g_tt_dq2
            dg_tt_dzeta  = 2.0_wp*hmap%eval_gij(q_thet_zeta,q,q_thet) + q_zeta(1)*g_tt_dq1 + q_zeta(2)*g_tt_dq2

            g_tz         = hmap%eval_gij(     q_thet,q,q_zeta)
            g_tz_dq1     = hmap%eval_gij_dq1( q_thet,q,q_zeta)
            g_tz_dq2     = hmap%eval_gij_dq2( q_thet,q,q_zeta)
            dg_tz_ds     = hmap%eval_gij(q_s_thet   ,q,q_zeta) + hmap%eval_gij(q_thet,q,q_s_zeta   ) + q_s(   1) * g_tz_dq1 + q_s(   2) * g_tz_dq2 
            dg_tz_dthet  = hmap%eval_gij(q_thet_thet,q,q_zeta) + hmap%eval_gij(q_thet,q,q_thet_zeta) + q_thet(1) * g_tz_dq1 + q_thet(2) * g_tz_dq2 
            dg_tz_dzeta  = hmap%eval_gij(q_thet_zeta,q,q_zeta) + hmap%eval_gij(q_thet,q,q_zeta_zeta) + q_zeta(1) * g_tz_dq1 + q_zeta(2) * g_tz_dq2 

            g_zz         = hmap%eval_gij(            q_zeta,q,q_zeta)
            g_zz_dq1     = hmap%eval_gij_dq1(        q_zeta,q,q_zeta)
            g_zz_dq2     = hmap%eval_gij_dq2(        q_zeta,q,q_zeta)
            dg_zz_ds     = 2.0_wp*hmap%eval_gij(q_s_zeta   ,q,q_zeta) + q_s(   1)*g_zz_dq1 + q_s(   2)*g_zz_dq2
            dg_zz_dthet  = 2.0_wp*hmap%eval_gij(q_thet_zeta,q,q_zeta) + q_thet(1)*g_zz_dq1 + q_thet(2)*g_zz_dq2

            dBsubs_dthet     = (-(Bthet*g_st+Bzeta*g_sz)/sqrtg*dsqrtg_dthet +Bthet*dg_st_dthet + dBthet_dthet*g_st + Bzeta*dg_sz_dthet + dBzeta_dthet*g_sz)/sqrtg
            dBsubs_dzeta     = (-(Bthet*g_st+Bzeta*g_sz)/sqrtg*dsqrtg_dzeta +Bthet*dg_st_dzeta + dBthet_dzeta*g_st + Bzeta*dg_sz_dzeta + dBzeta_dzeta*g_sz)/sqrtg

            dBsubthet_ds     = (-(Bthet*g_tt+Bzeta*g_tz)/sqrtg*dsqrtg_ds    +Bthet*dg_tt_ds    + dBthet_ds   *g_tt + Bzeta*dg_tz_ds    + dBzeta_ds   *g_tz)/sqrtg
            dBsubthet_dzeta  = (-(Bthet*g_tt+Bzeta*g_tz)/sqrtg*dsqrtg_dzeta +Bthet*dg_tt_dzeta + dBthet_dzeta*g_tt + Bzeta*dg_tz_dzeta + dBzeta_dzeta*g_tz)/sqrtg

            dBsubzeta_ds     = (-(Bthet*g_tz+Bzeta*g_zz)/sqrtg*dsqrtg_ds    +Bthet*dg_tz_ds    + dBthet_ds   *g_tz + Bzeta*dg_zz_ds    + dBzeta_ds   *g_zz)/sqrtg
            dBsubzeta_dthet  = (-(Bthet*g_tz+Bzeta*g_zz)/sqrtg*dsqrtg_dthet +Bthet*dg_tz_dthet + dBthet_dthet*g_tz + Bzeta*dg_zz_dthet + dBzeta_dthet*g_zz)/sqrtg 
           
            
            Js     = (dBsubzeta_dthet - dBsubthet_dzeta)
            Jthet  = (dBsubs_dzeta    - dBsubzeta_ds   )
            Jzeta  = (dBsubthet_ds    - dBsubs_dthet   )

            Jcart(:) = (e_s(:)*Js+ e_thet(:) * Jthet + e_zeta(:) * Jzeta)/sqrtg
            var_visu(VP_J:VP_J+2,i_s,j_s,i_n,i_m,iElem) = Jcart(:)/(2.0e-7_wp*TWOPI)  !*1/mu_0
#endif /*VISU_J_EXACT*/


#ifdef VISU_J_FD
            ! Get J components - finite difference bases
  
            
            ! Calculate ds derivative of B
            X1_eps         = X1_base%f%evalDOF_x(xIP, 0, X1_s_eps )
            X2_eps         = X2_base%f%evalDOF_x(xIP, 0, X2_s_eps )
            dX1_ds_eps     = X1_base%f%evalDOF_x(xIP, 0, dX1ds_eps)
            dX2_ds_eps     = X2_base%f%evalDOF_x(xIP, 0, dX2ds_eps)
            
            dX1_dthet_eps  = X1_base%f%evalDOF_x(xIP, DERIV_THET,X1_s_eps )
            dX2_dthet_eps  = X2_base%f%evalDOF_x(xIP, DERIV_THET,X2_s_eps )
            dLA_dthet_eps  = LA_base%f%evalDOF_x(xIP, DERIV_THET,LA_s_eps )
            
            dX1_dzeta_eps  = X1_base%f%evalDOF_x(xIP,DERIV_ZETA,X1_s_eps )
            dX2_dzeta_eps  = X2_base%f%evalDOF_x(xIP,DERIV_ZETA,X2_s_eps )
            dLA_dzeta_eps  = LA_base%f%evalDOF_x(xIP,DERIV_ZETA,LA_s_eps )

            q        = (/ X1_eps, X2_eps, xIP(2) /) !(X1,X2,zeta)
            e_s      = hmap%eval_dxdq(q,(/dX1_ds_eps   ,dX2_ds_eps   , 0.0    /)) !dxvec/ds
            e_thet   = hmap%eval_dxdq(q,(/dX1_dthet_eps,dX2_dthet_eps, 0.0    /)) !dxvec/dthet
            e_zeta   = hmap%eval_dxdq(q,(/dX1_dzeta_eps,dX2_dzeta_eps, 1.0_wp /)) !dxvec/dzeta
            sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))

            Bthet   = (iota_s_eps - dLA_dzeta_eps ) * phiPrime_s_eps   !/sqrtG
            Bzeta   = (1.0_wp  + dLA_dthet_eps ) * phiPrime_s_eps       !/sqrtG
            B_ds(:) =  (( e_thet(:) * Bthet + e_zeta(:) * Bzeta) /sqrtG - Bcart(:)) / (delta_s)      ! calculating dBx_ds, dBy_ds, dBz_ds
            
            ! Calculate dtheta derivative of B
            delta_thet = eps*SQRT(SUM(grad_thet*grad_thet))
            xIP_eps        = (/xIP(1)+delta_thet, xIP(2)/)
            X1_eps         = X1_base%f%evalDOF_x(xIP_eps, 0, X1_s )
            X2_eps         = X2_base%f%evalDOF_x(xIP_eps, 0, X2_s )
            dX1_ds_eps     = X1_base%f%evalDOF_x(xIP_eps, 0, dX1ds)
            dX2_ds_eps     = X2_base%f%evalDOF_x(xIP_eps, 0, dX2ds)
            
            dX1_dthet_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_THET, X1_s )
            dX2_dthet_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_THET, X2_s )
            dLA_dthet_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_THET, LA_s )
            
            dX1_dzeta_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X1_s )
            dX2_dzeta_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X2_s )
            dLA_dzeta_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, LA_s )

            q        = (/ X1_eps, X2_eps, xIP_eps(2) /) !(X1,X2,zeta)
            e_s      = hmap%eval_dxdq(q,(/dX1_ds_eps   ,dX2_ds_eps   , 0.0    /)) !dxvec/ds
            e_thet   = hmap%eval_dxdq(q,(/dX1_dthet_eps,dX2_dthet_eps, 0.0    /)) !dxvec/dthet
            e_zeta   = hmap%eval_dxdq(q,(/dX1_dzeta_eps,dX2_dzeta_eps, 1.0_wp /)) !dxvec/dzeta
            sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))

            Bthet   = (iota_s - dLA_dzeta_eps ) * phiPrime_s   !/sqrtG
            Bzeta   = (1.0_wp  + dLA_dthet_eps ) * phiPrime_s       !/sqrtG
            B_dthet(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / (delta_thet)      ! calculating dBx_dtheta, dBy_dtheta, dBz_dtheta
!
!           ! Calculate dzeta derivative of B
            delta_zeta = eps*SQRT(SUM(grad_zeta*grad_zeta))
            xIP_eps = (/xIP(1), xIP(2)+delta_zeta/)
            X1_eps         = X1_base%f%evalDOF_x(xIP_eps, 0, X1_s )
            X2_eps         = X2_base%f%evalDOF_x(xIP_eps, 0, X2_s )
            dX1_ds_eps     = X1_base%f%evalDOF_x(xIP_eps, 0, dX1ds)
            dX2_ds_eps     = X2_base%f%evalDOF_x(xIP_eps, 0, dX2ds)
            
            dX1_dthet_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_THET, X1_s )
            dX2_dthet_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_THET, X2_s )
            dLA_dthet_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_THET, LA_s )
            
            dX1_dzeta_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X1_s )
            dX2_dzeta_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X2_s )
            dLA_dzeta_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, LA_s )

            q        = (/ X1_eps, X2_eps, xIP_eps(2) /) !(X1,X2,zeta)
            e_s      = hmap%eval_dxdq(q,(/dX1_ds_eps   ,dX2_ds_eps   , 0.0    /)) !dxvec/ds
            e_thet   = hmap%eval_dxdq(q,(/dX1_dthet_eps,dX2_dthet_eps, 0.0    /)) !dxvec/dthet
            e_zeta   = hmap%eval_dxdq(q,(/dX1_dzeta_eps,dX2_dzeta_eps, 1.0_wp /)) !dxvec/dzeta
            sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))

            Bthet   = (iota_s - dLA_dzeta_eps ) * phiPrime_s   !/sqrtG
            Bzeta   = (1.0_wp  + dLA_dthet_eps ) * phiPrime_s       !/sqrtG
            B_dzeta(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / (delta_zeta)    ! calculating dBx_dzeta, dBy_dzeta, dBz_dzeta

            ! Calculate B derivatives by finite difference
            grad_Bcart(1, :) = B_ds(1) * grad_s(:) + B_dthet(1) * grad_thet(:) + B_dzeta(1) * grad_zeta(:)   ! grad_BX
            grad_Bcart(2, :) = B_ds(2) * grad_s(:) + B_dthet(2) * grad_thet(:) + B_dzeta(2) * grad_zeta(:)   ! grad_BY
            grad_Bcart(3, :) = B_ds(3) * grad_s(:) + B_dthet(3) * grad_thet(:) + B_dzeta(3) * grad_zeta(:)   ! grad_BZ 

            ! Calculate current cartesian components
            Jcart(1) = grad_Bcart(3, 2) - grad_Bcart(2, 3)   ! dBZ_dY - dBY_dZ
            Jcart(2) = grad_Bcart(1, 3) - grad_Bcart(3, 1)   ! dBX_dZ - dBZ_dX
            Jcart(3) = grad_Bcart(2, 1) - grad_Bcart(1, 2)   ! dBY_dX - dBX_dY
            var_visu(VP_J:VP_J+2,i_s,j_s,i_n,i_m,iElem) = Jcart(:)/(2.0e-7_wp*TWOPI)  !*1/mu_0
#endif /*VISU_J_FD*/
          END DO !j_s
        END DO !i_n
      END DO !i_m
!OMP END PARALLEL DO
      Itor_int = Itor_int*TWOPI/(REAL((mn_IP(1)*mn_IP(2)*n_s),wp)) !(2pi)^2/nfp /(Nt*Nz) * nfp/(2pi)
      Ipol_int = Ipol_int*TWOPI/(REAL((mn_IP(1)*mn_IP(2)*n_s),wp))
      var_visu(VP_Itor,i_s,:,:,:,iElem) = Itor_int/(2.0e-7_wp*TWOPI) !*1/mu_0
      var_visu(VP_Ipol,i_s,:,:,:,iElem) = Ipol_int/(2.0e-7_wp*TWOPI) !*1/mu_0
    END DO !i_s
  END DO !iElem

  ! average data in theta at the axis:
  !IF(minMax(1,0).LE.1e-4)THEN
  !  DO i_n=1,mn_IP(2)
  !    DO iVal=1,nVal
  !      var_visu(iVal,1,:,i_n,:,1)=SUM(var_visu(iVal,1,:,i_n,:,1))/REAL(mn_IP(1)*n_s,wp)
  !    END DO !iVal
  !  END DO !i_n
  !END IF

  !make grid exactly periodic
    !make theta direction exactly periodic
    IF(ABS((minMax(2,1)-minmax(2,0))-1.0_wp).LT.1.0e-04)THEN !fully periodic
!$OMP PARALLEL DO  COLLAPSE(3)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iElem,i_n,i_s)
    DO iElem=1,nElems; DO i_n=1,mn_IP(2); DO i_s=1,n_s
      coord_visu( :,i_s,n_s,i_n,mn_IP(1),iElem)=coord_visu( :,i_s,1,i_n,1,iElem)
    END DO; END DO; END DO
!$OMP END PARALLEL DO
    END IF
  IF(hmap%which_hmap.NE.3)THEN !not for cylinder
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
  var_visu(VP_MscaleF,n_s,:,:,:,nElems)= var_visu(VP_MscaleF,n_s-1,:,:,:,nElems) !boundary force=0
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
    WRITE(filename,'(A,"_visu_3D_",I4.4,"_",I8.8)')TRIM(Projectname),outputLevel,fileID
    IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
    CALL WriteDataToVTK(3,3,nVal,nplot,mn_IP(1)*(maxElem-minElem+1),VarNames, &
                        coord_visu(:,:,:,:,:,minElem:maxElem), &
                            var_visu(:,:,:,:,:,minElem:maxElem),TRIM(filename)//".vtu")
    END IF
    IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
      ALLOCATE(tmpcoord(1:3,1:(n_s-1)*(maxElem-minElem+1)+1,1:(n_s-1)*mn_IP(1)+1,mn_IP(2)))
      ALLOCATE(tmpvar(1:nVal,1:(n_s-1)*(maxElem-minElem+1)+1,1:(n_s-1)*mn_IP(1)+1,mn_IP(2)))
      DO i_n=1,mn_IP(2)
        j=1
        DO i_m=1,mn_IP(1); DO j_s=1,MERGE(n_s-1,n_s,i_m.LT.mn_IP(1))
           i=1
           DO iElem=minElem,maxElem;   DO i_s=1,MERGE(n_s-1,n_s,iElem.LT.maxElem)
             tmpcoord(:,i,j,i_n)=coord_visu( :,i_s,j_s,i_n,i_m,iElem)
             tmpvar(  :,i,j,i_n)=var_visu(   :,i_s,j_s,i_n,i_m,iElem)
             i=i+1
           END DO; END DO
           j=j+1
        END DO; END DO
      END DO
      CALL WriteDataToNETCDF(3,3,nVal,(/(maxElem-minElem+1)*(n_s-1)+1,mn_IP(1)*(n_s-1)+1,mn_IP(2)/),&
                          (/"dim_rho  ","dim_theta","dim_zeta "/),VarNames, &
                          tmpcoord,tmpvar, TRIM(filename))
      DEALLOCATE(tmpcoord,tmpvar)
    END IF !outfileType
  END IF
  __PERFOFF("write_visu")
  WRITE(filename,'(A,"_visu_1D_",I4.4,"_",I8.8)') &
    TRIM(Projectname),outputLevel,fileID
  CoordNames(1)="X"
  CoordNames(2)="Y"
  CoordNames(3)="Z"
  i=1 
  DO iElem=1,nElems;   DO i_s=1,MERGE(n_s-1,n_s,iElem.LT.nElems)
    var_visu_1d(1:3,i)     =coord_visu(:,i_s,1,1,1,iElem)
    var_visu_1d(4:3+nval,i)=var_visu(  :,i_s,1,1,1,iElem)
    i=i+1
  END DO; END DO
#if NETCDF
  CALL WriteDataToNETCDF(1,3,nVal-3,(/(n_s-1)*nElems+1/),(/"dim_rho"/), &
       VarNames(4:nVal),var_visu_1d(1:3,:),var_visu_1d(4:3+nVal,:), TRIM(filename))
#else
  CALL WriteDataToCSV((/CoordNames,VarNames(:)/) ,var_visu_1d,TRIM(filename)//".csv"  &
                                  ,append_in=.FALSE.,vfmt_in='E15.5')
#endif
  
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
!> convert solution Uin to straight-field line coordinates, and then write to visualization/netcdf file.
!! evaluation at given SFLout_radialpos. Passed to a grid, then a deg=1 spline is used, which is interpolatory at the grid points.
!!
!===================================================================================================================================
SUBROUTINE WriteSFLoutfile(Uin,fileID)
! MODULES
  USE MODgvec_MHD3D_Vars,     ONLY: hmap,X1_base,X2_base,LA_base
  USE MODgvec_MHD3D_Profiles, ONLY: Eval_iota,Eval_PhiPrime
  USE MODgvec_Transform_SFL,  ONLY: t_transform_sfl,transform_sfl_new
  USE MODgvec_output_netcdf,  ONLY: WriteDataToNETCDF
  USE MODgvec_output_vtk,     ONLY: WriteDataToVTK
  USE MODgvec_Output_vars,    ONLY: ProjectName,outputLevel
  USE MODgvec_Analyze_Vars,   ONLY: outfileType,SFLout,SFLout_nrp,SFLout_mn_pts,SFLout_mn_max,SFLout_radialpos
  USE MODgvec_sol_var_MHD3D,  ONLY: t_sol_var_mhd3d
  USE MODgvec_Globals,        ONLY: TWOPI,CROSS
  IMPLICIT NONE 
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin      !! input solution 
  INTEGER , INTENT(IN   ) :: fileID          !! added to file name before the ending
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  CLASS(t_transform_sfl),ALLOCATABLE :: trafoSFL
  REAL(wp),ALLOCATABLE       :: coord_out(:,:,:,:),var_out(:,:,:,:),thetstar_pos(:),zetastar_pos(:)
  INTEGER                    :: i_rp,izeta,ithet,nthet_out,nzeta_out,i
  INTEGER                    :: mn_max(2),factorSFL,iVal
  REAL(wp)                   :: spos,xp(2),sqrtG
  REAL(wp)                   :: dX1ds,dX1dthetstar,dX1dzetastar
  REAL(wp)                   :: dX2ds,dX2dthetstar,dX2dzetastar
  REAL(wp)                   :: phiPrime_int,iota_int,Itor_int,Ipol_int
  REAL(wp)                   :: X1_int,X2_int,GZ_int,dGZds,dGZdthetstar,dGZdzetastar,LA_int,dLA_dthet,dLA_dzeta
  REAL(wp)                   :: Gt_int,dGZdthet,dGZdzeta,dX1dthet,dX1dzeta,dX2dthet,dX2dzeta
  REAL(wp)                   :: dthetstar_dthet ,dthetstar_dzeta ,dzetastar_dthet ,dzetastar_dzeta,Jstar
  REAL(wp)                   :: dthet_dthetstarJ,dthet_dzetastarJ,dzeta_dthetstarJ,dzeta_dzetastarJ
  REAL(wp)                   :: Bthet,Bzeta,Bthetstar,Bzetastar
  REAL(wp),DIMENSION(3)      :: qvec,e_s,e_thet,e_zeta,e_thetstar,e_zetastar,Bfield
  REAL(wp),ALLOCATABLE       :: X1_s(:),dX1ds_s(:),X2_s(:),dX2ds_s(:),Gz_s(:),dGZds_s(:),LA_s(:),Gt_s(:)
  INTEGER                    :: VP_rho,VP_iota,VP_thetastar,VP_zetastar,VP_zeta,VP_X1sfl,VP_X2sfl,VP_SQRTG,VP_SQRTGstar,&
                                VP_B,VP_modB,VP_Bsubthet,VP_Bsubzeta,VP_Bthet,VP_Bzeta,VP_grads,VP_theta,&
                                VP_Bsubthetstar,VP_Bsubzetastar,VP_Bthetstar,VP_Bzetastar,VP_Itor,VP_Ipol
  INTEGER,PARAMETER          :: nVal=25
  CHARACTER(LEN=40)          :: VarNames(nval)  
  CHARACTER(LEN=10)          :: sfltype 
  CHARACTER(LEN=255)         :: filename
  LOGICAL                    :: useSFLcoords
  INTEGER                    :: dbg
  !=================================================================================================================================
  IF(SFLout.EQ.0) RETURN
  sfltype=MERGE("_boozer","_pest  ",SFLout.EQ.2)
  WRITE(filename,'(A,"_",I4.4,"_",I8.8,"")') & 
  TRIM(Projectname)//TRIM(sfltype),outputLevel,fileID
  SWRITE(UNIT_stdOut,'(A,A,A)') 'WRITING SFL output: ',TRIM(filename),' ...'
  __PERFON("output_sfl")
  iVal=1
  VP_rho        =iVal;iVal=iVal+1; VarNames(VP_rho      )="rho"
  VP_iota       =iVal;iVal=iVal+1; VarNames(VP_iota     )="iota"
  VP_Itor       =iVal;iVal=iVal+1; VarNames(VP_Itor     )="Itor"
  VP_Ipol       =iVal;iVal=iVal+1; VarNames(VP_Ipol     )="Ipol"
  VP_thetastar  =iVal;iVal=iVal+1; VarNames(VP_thetastar)="thetastar"
  VP_zetastar   =iVal;iVal=iVal+1; VarNames(VP_zetastar )="zetastar"
  VP_theta      =iVal;iVal=iVal+1; VarNames(VP_theta    )="theta"
  VP_zeta       =iVal;iVal=iVal+1; VarNames(VP_zeta     )="zeta"
  VP_SQRTG      =iVal;iVal=iVal+1; VarNames(VP_SQRTG    )="sqrtG"
  VP_SQRTGstar  =iVal;iVal=iVal+1; VarNames(VP_SQRTGstar)="sqrtGstar"
  VP_modB       =iVal;iVal=iVal+1; VarNames(VP_modB     )="modB"
  VP_B          =iVal;iVal=iVal+3; VarNames(VP_B        )="BvecX"
                                   VarNames(VP_B+1      )="BvecY"
                                   VarNames(VP_B+2      )="BvecZ"
  VP_grads      =iVal;iVal=iVal+3; VarNames(VP_grads    )="grad_sX"
                                   VarNames(VP_grads+1  )="grad_sY"
                                   VarNames(VP_grads+2  )="grad_sZ"
  VP_Bsubthet   =iVal;iVal=iVal+1; VarNames(VP_Bsubthet )="Bsubtheta"
  VP_Bsubzeta   =iVal;iVal=iVal+1; VarNames(VP_Bsubzeta )="Bsubzeta"
  VP_Bthet      =iVal;iVal=iVal+1; VarNames(VP_Bthet    )="Btheta"
  VP_Bzeta      =iVal;iVal=iVal+1; VarNames(VP_Bzeta    )="Bzeta"
  VP_Bsubthetstar   =iVal;iVal=iVal+1; VarNames(VP_Bsubthetstar )="Bsubthetastar"
  VP_Bsubzetastar   =iVal;iVal=iVal+1; VarNames(VP_Bsubzetastar )="Bsubzetastar"
  VP_Bthetstar      =iVal;iVal=iVal+1; VarNames(VP_Bthetstar    )="Bthetastar"
  VP_Bzetastar      =iVal;iVal=iVal+1; VarNames(VP_Bzetastar    )="Bzetastar"


  IF(iVal.NE.Nval+1) CALL abort(__STAMP__,"nVal parameter not correctly set")


  factorSFL=4
  DO i=1,2
    IF(SFLout_mn_max(i).EQ.-1)THEN !input =-1, automatic
      mn_max(i) = factorSFL*MAXVAL((/X1_base%f%mn_max(i),X2_base%f%mn_max(i),LA_base%f%mn_max(i)/))
    ELSE 
      mn_max(i) = SFLout_mn_max(i) !user defined
    END IF
  END DO

  CALL transform_sfl_new(trafoSFL,mn_max,SFLout,.false.,&  ! relambda=false
                         X1_base%s%deg,X1_base%s%continuity,X1_base%s%degGP,X1_base%s%grid,&
                         hmap,X1_base,X2_base,LA_base,Eval_PhiPrime,Eval_iota)  !same grid and degree as variable X1.
  CALL trafoSFL%buildTransform(X1_base,X2_base,LA_base,Uin%X1,Uin%X2,Uin%LA)

  Nthet_out=MERGE(2*mn_max(1)+1,SFLout_mn_pts(1),SFLout_mn_pts(1).EQ.-1) !if input =-1, automatically 2*m_max+1, else user defined
  Nzeta_out=MERGE(2*mn_max(2)+1,SFLout_mn_pts(2),SFLout_mn_pts(2).EQ.-1) !if input =-1, automatically 2*n_max+1
  
  DO dbg=1,2
  ASSOCIATE(n_rp=>SFLout_nrp, rp=>SFLout_radialpos,SFLcoord=>trafoSFL%whichSFLcoord)

  ALLOCATE(thetstar_pos(Nthet_out))
  ALLOCATE(zetastar_pos(Nzeta_out))
  DO ithet=1,Nthet_out
    thetstar_pos(ithet)=(TWOPI*(REAL(ithet,wp)-0.5_wp))/REAL(Nthet_out) 
  END DO
  DO izeta=1,Nzeta_out
    zetastar_pos(izeta)=(TWOPI*(REAL(izeta,wp)-0.5_wp))/REAL((Nzeta_out*trafoSFL%nfp),wp)
  END DO
  ALLOCATE(coord_out(3,Nthet_out,Nzeta_out,n_rp),var_out(nVal,Nthet_out,Nzeta_out,n_rp))
  var_out=0.
  
  useSFLcoords=(dbg.EQ.1)
  WRITE(*,*)'USESFLCOORDS=',useSFLcoords

  IF(.NOT.useSFLcoords)THEN !use quantities given in GVEC theta and zeta:
    ASSOCIATE(G_base=>trafoSFL%GZ_base,Gt=>trafoSFL%Gthet,Gz=>trafoSFL%Gz)
    ALLOCATE(X1_s(X1_base%f%modes),dX1ds_s(X1_base%f%modes))
    ALLOCATE(X2_s(X2_base%f%modes),dX2ds_s(X2_base%f%modes))
    ALLOCATE(LA_s(LA_base%f%modes))
    IF(SFLcoord.EQ.2) ALLOCATE(Gt_s(G_base%f%modes),Gz_s(G_base%f%modes))
    DO i_rp=1,n_rp
      Itor_int = 0.
      Ipol_int = 0.
      spos=rp(i_rp)
      iota_int=Eval_iota(spos)
      phiPrime_int=Eval_PhiPrime(spos)
      var_out(VP_rho ,:,:,i_rp)=spos
      var_out(VP_iota,:,:,i_rp)=iota_int
      GZ_int   = 0.0_wp !only changed for SFLcoords=2
      dGZdthetstar = 0.0_wp !only changed for SFLcoords=2
      dGZdzetastar = 0.0_wp !only changed for SFLcoords=2
      !interpolate radially
      X1_s(   :) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,       0,Uin%X1(:,:))
      dX1ds_s(:) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes, DERIV_S,Uin%X1(:,:))
    
      X2_s(   :) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,       0,Uin%X2(:,:))
      dX2ds_s(:) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes, DERIV_S,Uin%X2(:,:))
      !IF(SFLcoord.EQ.1) THEN
      LA_s=LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,0,Uin%LA)
      IF(SFLcoord.EQ.2)THEN
        Gt_s = G_base%s%evalDOF2D_s(spos,G_base%f%modes, 0,Gt)
        Gz_s = G_base%s%evalDOF2D_s(spos,G_base%f%modes, 0,Gz)
      END IF
      DO izeta=1,Nzeta_out; DO ithet=1,Nthet_out
        xp=(/thetstar_pos(ithet),zetastar_pos(izeta)/) !=theta,zeta GVEC !!!
        IF(SFLcoord.EQ.1)THEN
          var_out(VP_thetastar,ithet,izeta,i_rp)=thetstar_pos(ithet)+LA_base%f%evalDOF_x(xp,0,LA_s)
          var_out(VP_zetastar ,ithet,izeta,i_rp)=zetastar_pos(izeta)
          dthetstar_dthet=1.+LA_base%f%evalDOF_x(xp,DERIV_THET, LA_s)
          dthetstar_dzeta=1.+LA_base%f%evalDOF_x(xp,DERIV_ZETA, LA_s)
          dzetastar_dthet=0.
          dzetastar_dzeta=1.
        ELSEIF(SFLcoord.EQ.2)THEN
          var_out(VP_thetastar,ithet,izeta,i_rp)=thetstar_pos(ithet)+G_base%f%evalDOF_x(xp, 0, Gt_s)
          var_out(VP_zetastar ,ithet,izeta,i_rp)=zetastar_pos(izeta)+G_base%f%evalDOF_x(xp, 0, Gz_s)
          dGZdthet=G_base%f%evalDOF_x(xp,DERIV_THET, Gz_s)  !d/dthet, but evaluated at theta*,zeta*!
          dGZdzeta=G_base%f%evalDOF_x(xp,DERIV_ZETA, Gz_s)
          dthetstar_dthet=1.+G_base%f%evalDOF_x(xp,DERIV_THET, Gt_s)
          dthetstar_dzeta=   G_base%f%evalDOF_x(xp,DERIV_ZETA, Gt_s)
          dzetastar_dthet=   dGZdthet
          dzetastar_dzeta=1.+dGZdzeta
        END IF
        !inverse:
        Jstar=dthetstar_dthet*dzetastar_dzeta-dthetstar_dzeta*dzetastar_dthet
        dthet_dthetstarJ= dzetastar_dzeta !/Jstar*Jstar
        dzeta_dzetastarJ= dthetstar_dthet !/Jstar*Jstar
        dthet_dzetastarJ=-dthetstar_dzeta !/Jstar*Jstar
        dzeta_dthetstarJ=-dzetastar_dthet !/Jstar*Jstar

        X1_int   = X1_base%f%evalDOF_x(xp,          0, X1_s  )
        dX1ds    = X1_base%f%evalDOF_x(xp,          0,dX1ds_s)
        dX1dthet = X1_base%f%evalDOF_x(xp, DERIV_THET, X1_s  ) 
        dX1dzeta = X1_base%f%evalDOF_x(xp, DERIV_ZETA, X1_s  ) 
        
        X2_int   = X2_base%f%evalDOF_x(xp,          0, X2_s  )
        dX2ds    = X2_base%f%evalDOF_x(xp,          0,dX2ds_s)
        dX2dthet = X2_base%f%evalDOF_x(xp, DERIV_THET, X2_s  )
        dX2dzeta = X2_base%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )
        
        dLA_dthet= LA_base%f%evalDOF_x(xp, DERIV_THET, LA_s  )
        dLA_dzeta= LA_base%f%evalDOF_x(xp, DERIV_ZETA, LA_s  )
        
        ! !transform derivative from dthet,dzeta=>dthet*,dzeta*
        ! dX1dthetstar = (dX1dthet*dthet_dthetstarJ+dX1dzeta*dzeta_dthetstarJ)/Jstar
        ! dX2dthetstar = (dX2dthet*dthet_dthetstarJ+dX2dzeta*dzeta_dthetstarJ)/Jstar

        ! dX1dzetastar = (dX1dthet*dthet_dzetastarJ+dX1dzeta*dzeta_dzetastarJ)/Jstar
        ! dX2dzetastar = (dX2dthet*dthet_dzetastarJ+dX2dzeta*dzeta_dzetastarJ)/Jstar
        ! IF(SFLcoord.EQ.2)THEN
        !   dGZdthetstar=(dGZdthet*dthet_dthetstarJ+dGZdzeta*dzeta_dthetstarJ)/Jstar
        !   dGZdzetastar=(dGZdthet*dthet_dzetastarJ+dGZdzeta*dzeta_dzetastarJ)/Jstar
        ! END IF

        qvec=(/X1_int,X2_int,xp(2)/)
        coord_out(:,ithet,izeta,i_rp)=hmap%eval(qvec)
        e_s   =hmap%eval_dxdq(qvec,(/dX1ds   ,dX2ds   ,0.0_wp/))
        e_thet=hmap%eval_dxdq(qvec,(/dX1dthet,dX2dthet,0.0_wp/))
        e_zeta=hmap%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta,1.0_wp/))
        sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))
        e_thetstar=(e_thet*dthet_dthetstarJ+e_zeta*dzeta_dthetstarJ)/Jstar
        e_zetastar=(e_thet*dthet_dzetastarJ+e_zeta*dzeta_dzetastarJ)/Jstar

        Bthet   = (iota_int - dLA_dzeta ) * phiPrime_int   !/sqrtG
        Bzeta   = (1.0_wp + dLA_dthet ) * phiPrime_int       !/sqrtG
        Bfield(:) =  ( e_thet(:) * Bthet + e_zeta(:) * Bzeta) /sqrtG

        Itor_int = Itor_int+ SUM(Bfield(:)*e_thet(:))   !B_theta=B.e_thet 
        Ipol_int = Ipol_int+ SUM(Bfield(:)*e_zeta(:))   !B_zeta =B.e_zeta
        ! !e_s          = hmap%eval_dxdq(qvec,(/dX1ds       ,dX2ds       ,       -dGZds       /)) !dxvec/ds
        ! e_thetstar   = hmap%eval_dxdq(qvec,(/dX1dthetstar,dX2dthetstar,       -dGZdthetstar/)) !dxvec/dthetstar
        ! e_zetastar   = hmap%eval_dxdq(qvec,(/dX1dzetastar,dX2dzetastar,1.0_wp -dGZdzetastar/)) !dxvec/dzetastar
        ! !sqrtG        = SUM(e_s*(CROSS(e_thetstar,e_zetastar)))
        ! sqrtG        = hmap%eval_Jh(qvec)*(dX1ds*dX2dthetstar-dX1dthetstar*dX2ds)
        ! Bthetstar    = iota_int*PhiPrime_int   !/sqrtG
        ! Bzetastar    =          PhiPrime_int   !/sqrtG
        ! Bfield(:)    =  ( e_thetstar(:)*Bthetstar+e_zetastar(:)*Bzetastar) /sqrtG
        
        var_out(VP_theta    ,ithet,izeta,i_rp)=xp(1)
        var_out(VP_zeta     ,ithet,izeta,i_rp)=xp(2)
        var_out(VP_SQRTG    ,ithet,izeta,i_rp)=sqrtG
        var_out(VP_SQRTGstar,ithet,izeta,i_rp)=sqrtG/Jstar !=sqrtGstar
        var_out(VP_B:VP_B+2 ,ithet,izeta,i_rp)=Bfield
        var_out(VP_grads:VP_grads+2 ,ithet,izeta,i_rp)=CROSS(e_thet,e_zeta)/sqrtG
        var_out(VP_modB     ,ithet,izeta,i_rp)=SQRT(SUM(Bfield*Bfield))
        var_out(VP_Bthet    ,ithet,izeta,i_rp)=SUM(Bfield*CROSS(e_zeta,e_s   ))/sqrtG
        var_out(VP_Bzeta    ,ithet,izeta,i_rp)=SUM(Bfield*CROSS(e_thet,e_zeta))/sqrtG
        var_out(VP_Bthetstar    ,ithet,izeta,i_rp)=SUM(Bfield*CROSS(e_zetastar,e_s       ))*Jstar/sqrtG
        var_out(VP_Bzetastar    ,ithet,izeta,i_rp)=SUM(Bfield*CROSS(e_thetstar,e_zetastar))
        var_out(VP_Bsubthet ,ithet,izeta,i_rp)=SUM(Bfield*e_thet)
        var_out(VP_Bsubzeta ,ithet,izeta,i_rp)=SUM(Bfield*e_zeta)
        var_out(VP_Bsubthetstar ,ithet,izeta,i_rp)=SUM(Bfield*e_thetstar)
        var_out(VP_Bsubzetastar ,ithet,izeta,i_rp)=SUM(Bfield*e_zetastar)
      END DO; END DO !izeta,ithet
      var_out(VP_Itor ,:,:,i_rp)= Itor_int*TWOPI/(REAL((nthet_out*nzeta_out),wp)) !(2pi)^2/nfp /(Nt*Nz) * nfp/(2pi)
      var_out(VP_Itor ,:,:,i_rp)= Ipol_int*TWOPI/(REAL((nthet_out*nzeta_out),wp))
    END DO !i_rp=1,n_rp

    END ASSOCIATE
  ELSE !use quantities given in SFL coords 
    ASSOCIATE(X1sfl_base=>trafoSFL%X1sfl_base,X1sfl=>trafoSFL%X1sfl,&
              X2sfl_base=>trafoSFL%X2sfl_base,X2sfl=>trafoSFL%X2sfl,&
              Gtsfl_base=>trafoSFL%GZsfl_base,Gtsfl=>trafoSFL%Gtsfl,&
              GZsfl_base=>trafoSFL%GZsfl_base,GZsfl=>trafoSFL%GZsfl )
    ALLOCATE(X1_s(X1sfl_base%f%modes),dX1ds_s(X1sfl_base%f%modes))
    ALLOCATE(X2_s(X2sfl_base%f%modes),dX2ds_s(X2sfl_base%f%modes)) 
    IF(SFLcoord.EQ.2) ALLOCATE(Gt_s(GZsfl_base%f%modes),GZ_s(GZsfl_base%f%modes),dGZds_s(GZsfl_base%f%modes))
    DO i_rp=1,n_rp
      spos=rp(i_rp)
      iota_int=Eval_iota(spos)
      Itor_int=0.
      Ipol_int=0.
      phiPrime_int=Eval_PhiPrime(spos)
      var_out(VP_rho ,:,:,i_rp)=spos
      var_out(VP_iota,:,:,i_rp)=iota_int
      GZ_int   = 0.0_wp !only changed for SFLcoords=2
      !dGZds    = 0.0_wp !only changed for SFLcoords=2
      dGZdthetstar = 0.0_wp !only changed for SFLcoords=2
      dGZdzetastar = 0.0_wp !only changed for SFLcoords=2
        !interpolate radially
      X1_s(   :) = X1sfl_base%s%evalDOF2D_s(spos,X1sfl_base%f%modes,       0,X1sfl(:,:))
      dX1ds_s(:) = X1sfl_base%s%evalDOF2D_s(spos,X1sfl_base%f%modes, DERIV_S,X1sfl(:,:))
    
      X2_s(   :) = X2sfl_base%s%evalDOF2D_s(spos,X2sfl_base%f%modes,       0,X2sfl(:,:))
      dX2ds_s(:) = X2sfl_base%s%evalDOF2D_s(spos,X2sfl_base%f%modes, DERIV_S,X2sfl(:,:))
      IF(SFLcoord.EQ.2)THEN !BOOZER
        Gt_s(   :) = GZsfl_base%s%evalDOF2D_s(spos,GZsfl_base%f%modes,      0,Gtsfl(:,:))
        GZ_s(   :) = GZsfl_base%s%evalDOF2D_s(spos,GZsfl_base%f%modes,      0,GZsfl(:,:))
        dGZds_s(:) = GZsfl_base%s%evalDOF2D_s(spos,GZsfl_base%f%modes,DERIV_S,GZsfl(:,:))
      END IF
      DO izeta=1,Nzeta_out; DO ithet=1,Nthet_out
        xp=(/thetstar_pos(ithet),zetastar_pos(izeta)/)
        X1_int       = X1sfl_base%f%evalDOF_x(xp,          0, X1_s  )
        dX1ds        = X1sfl_base%f%evalDOF_x(xp,          0,dX1ds_s)
        dX1dthetstar = X1sfl_base%f%evalDOF_x(xp, DERIV_THET, X1_s  )
        dX1dzetastar = X1sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )
        
        X2_int       = X2sfl_base%f%evalDOF_x(xp,          0, X2_s  )
        dX2ds        = X2sfl_base%f%evalDOF_x(xp,          0,dX2ds_s)
        dX2dthetstar = X2sfl_base%f%evalDOF_x(xp, DERIV_THET, X2_s  )
        dX2dzetastar = X2sfl_base%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )
  
        
        IF(SFLcoord.EQ.2)THEN !BOOZER coordinates (else=0)
          GZ_int       = GZsfl_base%f%evalDOF_x(xp,         0, GZ_s)
          dGZds        = GZsfl_base%f%evalDOF_x(xp,         0, dGZds_s)
          dGZdthetstar = GZsfl_base%f%evalDOF_x(xp,DERIV_THET, GZ_s)
          dGZdzetastar = GZsfl_base%f%evalDOF_x(xp,DERIV_ZETA, GZ_s)
        END IF
        
        qvec=(/X1_int,X2_int,zetastar_pos(izeta)-GZ_int/) !(X1,X2,"zeta=zetastar-GZ(spos,thetstar,zetastar)")
        coord_out(:,ithet,izeta,i_rp)=hmap%eval(qvec)
        e_s          = hmap%eval_dxdq(qvec,(/dX1ds       ,dX2ds       ,       -dGZds       /)) !dxvec/ds
        e_thetstar   = hmap%eval_dxdq(qvec,(/dX1dthetstar,dX2dthetstar,       -dGZdthetstar/)) !dxvec/dthetstar
        e_zetastar   = hmap%eval_dxdq(qvec,(/dX1dzetastar,dX2dzetastar,1.0_wp -dGZdzetastar/)) !dxvec/dzetastar
        sqrtG        = SUM(e_s*(CROSS(e_thetstar,e_zetastar)))
        !sqrtG        = hmap%eval_Jh(qvec)*(dX1ds*dX2dthetstar-dX1dthetstar*dX2ds)
        Bthetstar    = iota_int*PhiPrime_int   !/sqrtG
        Bzetastar    =          PhiPrime_int   !/sqrtG
        Bfield(:)    =  ( e_thetstar(:)*Bthetstar+e_zetastar(:)*Bzetastar) /sqrtG
        Itor_int = Itor_int+ SUM(Bfield(:)*e_thetstar(:))   !B_theta=B.e_thet 
        Ipol_int = Ipol_int+ SUM(Bfield(:)*e_zetastar(:))   !B_zeta =B.e_zeta
        var_out(VP_thetastar,ithet,izeta,i_rp)=thetstar_pos(ithet)
        var_out(VP_zetastar ,ithet,izeta,i_rp)=zetastar_pos(izeta)
        var_out(VP_theta    ,ithet,izeta,i_rp)=thetstar_pos(ithet)-GZsfl_base%f%evalDOF_x(xp,         0, Gt_s)
        var_out(VP_zeta     ,ithet,izeta,i_rp)=zetastar_pos(izeta)-GZ_int
        var_out(VP_SQRTGstar,ithet,izeta,i_rp)=sqrtG
        var_out(VP_B:VP_B+2 ,ithet,izeta,i_rp)=Bfield
        var_out(VP_grads:VP_grads+2,ithet,izeta,i_rp)=CROSS(e_thetstar,e_zetastar)/sqrtG
        var_out(VP_modB     ,ithet,izeta,i_rp)=SQRT(SUM(Bfield*Bfield))
        var_out(VP_Bthet    ,ithet,izeta,i_rp)=Bthetstar
        var_out(VP_Bzeta    ,ithet,izeta,i_rp)=Bzetastar
        var_out(VP_Bsubthetstar,ithet,izeta,i_rp)=SUM(Bfield*e_thetstar)
        var_out(VP_Bsubzetastar,ithet,izeta,i_rp)=SUM(Bfield*e_zetastar)
      END DO; END DO !ithet,izeta
      var_out(VP_Itor ,:,:,i_rp)= Itor_int*TWOPI/(REAL((nthet_out*nzeta_out),wp)) !(2pi)^2/nfp /(Nt*Nz) * nfp/(2pi)
      var_out(VP_Itor ,:,:,i_rp)= Ipol_int*TWOPI/(REAL((nthet_out*nzeta_out),wp))
    END DO !i_rp=1,,n_rp
    END ASSOCIATE
  END IF !useSFLcoords
  DEALLOCATE(X1_s,dX1ds_s,X2_s,dX2ds_s,thetstar_pos,zetastar_pos)
  SDEALLOCATE(LA_s)
  SDEALLOCATE(Gt_s)
  SDEALLOCATE(dGzds_s)
  SDEALLOCATE(Gz_s)


  IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
   CALL WriteDataToVTK(3,3,nVal,(/Nthet_out-1,Nzeta_out-1,n_rp-1/),1,VarNames, &
                      coord_out(1:3 ,1:Nthet_out,1:Nzeta_out,1:n_rp), &
                      var_out(1:nval,1:Nthet_out,1:Nzeta_out,1:n_rp),TRIM(filename)//TRIM(MERGE('_full  ','_direct',dbg.EQ.1))//".vtu")
  END IF
  IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
    CALL WriteDataToNETCDF(3,3,nVal,(/Nthet_out,Nzeta_out,n_rp/),&
                           (/"dim_theta","dim_zeta ","dim_rho  "/),VarNames, &
                           coord_out(1:3 ,1:Nthet_out,1:Nzeta_out,1:n_rp), &
                           var_out(1:nval,1:Nthet_out,1:Nzeta_out,1:n_rp), TRIM(filename)//TRIM(MERGE('_full  ','_direct',dbg.EQ.1)))
  END IF!outfileType
  END ASSOCIATE !n_rp,rp
  DEALLOCATE(coord_out,var_out)
END DO !dbg
  CALL trafoSFL%free()
  SWRITE(UNIT_stdOut,'(A)') '... DONE.'
  __PERFOFF("output_sfl")
END SUBROUTINE WriteSFLoutfile


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
    zeta1D(i_n)  = TWOPI*REAL(i_n-1,wp)/REAL(mn_IP(2)*X1_base%f%nfp,wp) !do not include periodic point 
  END DO
  DO i_m=1,mn_IP(1)
    theta1D(i_m)= TWOPI*REAL(i_m-1,wp)/REAL(mn_IP(1),wp)  !do not include periodic point
  END DO


!$OMP PARALLEL DO &
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   REDUCTION(+:avgDist) REDUCTION(max:maxDist) &
!$OMP   PRIVATE(i_m,i_n,xIP,q,theta,zeta,theta0,X1_visu,X2_visu,LA_visu,xU,xV,dist, &
!$OMP           UX1_s,UX2_s,ULA_s,VX1_s,VX2_s,VLA_s,spos,iElem,i_s) &
!$OMP   SHARED(nElems,n_s,mn_IP,theta1D,zeta1D,SFL_theta,X1_base,X2_base,LA_base,hmap,U,V,sgrid)
  DO iElem=1,nElems
    DO i_s=1,n_s
      spos=MAX(1.0e-06,sgrid%sp(iElem-1)+(REAL(i_s-1,wp))/(REAL(n_s,wp))*sgrid%ds(iElem)) !includes axis but not edge

      UX1_s(:) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,0,U%X1(:,:))
      VX1_s(:) = X1_base%s%evalDOF2D_s(spos,X1_base%f%modes,0,V%X1(:,:))
      UX2_s(:) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,0,U%X2(:,:))
      VX2_s(:) = X2_base%s%evalDOF2D_s(spos,X2_base%f%modes,0,V%X2(:,:))
      ULA_s(:) = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,0,U%LA(:,:))
      VLA_s(:) = LA_base%s%evalDOF2D_s(spos,LA_base%f%modes,0,V%LA(:,:))

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
    END DO !i_s
  END DO !iElem
!OMP$ END PARALLEL DO
  avgDist=avgDist/REAL(nElems*n_s*mn_IP(1)*mn_IP(2),wp)

  DEALLOCATE(theta1D,zeta1D)

  __PERFOFF("checkDistance")
END SUBROUTINE CheckDistance


!===================================================================================================================================
!> check distance between two solutions, via sampling X1,X2 at theta*=theta+lambda, and comparing the distance of 
!> the sampled x,y,z coordinates 
!!
!===================================================================================================================================
SUBROUTINE CheckAxis(U,n_zeta,AxisPos) 
! MODULES
  USE MODgvec_Globals,        ONLY: TWOPI
  USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base
  USE MODgvec_sol_var_MHD3D,  ONLY: t_sol_var_mhd3d
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN) :: U !U 
  INTEGER               , INTENT(IN) :: n_zeta  !! number of points checked along axis
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp),INTENT(OUT)    :: AxisPos(1:2,n_zeta)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i_n
  REAL(wp) :: zeta,UX1_s(1:X1_base%f%modes),UX2_s(1:X2_base%f%modes)
!===================================================================================================================================
  UX1_s(:) = X1_base%s%evalDOF2D_s(0.0_wp,X1_base%f%modes,0,U%X1(:,:))
  UX2_s(:) = X2_base%s%evalDOF2D_s(0.0_wp,X2_base%f%modes,0,U%X2(:,:))

  DO i_n=1,n_zeta
    zeta = TWOPI*REAL(i_n-1,wp)/REAL(n_zeta*X1_base%f%nfp,wp) !do not include periodic point 
    AxisPos(1,i_n) = X1_base%f%evalDOF_x((/0.0_wp,zeta/),0,UX1_s(:) )
    AxisPos(2,i_n) = X2_base%f%evalDOF_x((/0.0_wp,zeta/),0,UX2_s(:) )
  END DO !i_n
END SUBROUTINE CheckAxis

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
  LOGICAL            :: vcase(5)
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
  IF(INDEX(vstr,'5').NE.0) vcase(5)=.TRUE.
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
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,X1_base,U(0)%X1)
    vname="X2"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,X2_base,U(0)%X2)
    vname="LA"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,0,LA_base,U(0)%LA)
  END IF
  IF(vcase(2))THEN
    WRITE(*,*)'2) Visualize gvec modes in 1D: dX1rho,dX2rho,dLAdrho interpolated...'
    vname="dX1ds"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,X1_base,U(0)%X1)
    vname="dX2ds"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,X2_base,U(0)%X2)
    vname="dLAds"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,DERIV_S,LA_base,U(0)%LA)
  END IF
  IF(vcase(3))THEN
    WRITE(*,*)'3) Visualize gvec modes in 1D: (d/drho)^2 X1/X2/LA interpolated...'
    vname="dX1dss"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,2,X1_base,U(0)%X1)
    vname="dX2dss"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,2,X2_base,U(0)%X2)
    vname="dLAdss"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,2,LA_base,U(0)%LA)
  END IF
  IF(vcase(4))THEN
    WRITE(*,*)'4) Visualize gvec modes in 1D:  |X1|/|X2|/|LA| interpolated...'
    vname="absX1"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,-4,X1_base,U(0)%X1)
    vname="absX2"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,-4,X2_base,U(0)%X2)
    vname="absLA"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,-4,LA_base,U(0)%LA)
  END IF
  IF(vcase(5))THEN
    WRITE(*,*)'5) Visualize gvec modes in 1D:  |X1|/|X2|/|LA| / rho^m interpolated...'
    vname="absX1orhom"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,-5,X1_base,U(0)%X1)
    vname="absX2orhom"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,-5,X2_base,U(0)%X2)
    vname="absLAorhom"
    WRITE(fname,'(A,I4.4,"_",I8.8)')'U0_'//TRIM(vname)//'_',outputLevel,FileID
    CALL writeDataMN_visu(n_s,fname,vname,-5,LA_base,U(0)%LA)
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
  !first element blending of  logarithmic*(1-xi^2) / linear *xi^2
  DO i_s=1,n_s
    values_visu(1,i_s) =values_visu(1,i_s)*(REAL(i_s-1,wp)/REAL(n_s-1,wp))**2 + (1.-(REAL(i_s-1,wp)/REAL(n_s-1,wp))**2) * &
                        (sgrid%sp(0)+ (10**(8.0_wp*(-1.0_wp+REAL(i_s-1,wp)/REAL(n_s-1,wp)))) &
                                                 *(1.0e-06_wp+REAL(n_s+1,wp))/(2.0e-06_wp+REAL(n_s+1,wp))*sgrid%ds(1))
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
                                  ,append_in=.FALSE.)

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
  INTEGER         ,INTENT(IN   ) :: rderiv !! 0: eval spl, 1: eval spl deriv, (negative used as flag)
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
  INTEGER                        :: iMode,j,m
  CHARACTER(LEN=255)             :: fname
  CHARACTER(LEN=120),ALLOCATABLE :: varnames(:) 
  REAL(wp)          ,ALLOCATABLE :: values_visu(:,:)
  REAL(wp)          ,ALLOCATABLE :: s_visu(:)
  REAL(wp)                       :: rhom,val
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
  !first element blending of  logarithmic*(1-xi^2) / linear *xi^2
  DO i_s=1,n_s
    s_visu(i_s) =s_visu(i_s)*(REAL(i_s-1,wp)/REAL(n_s-1,wp))**2 + (1.-(REAL(i_s-1,wp)/REAL(n_s-1,wp))**2) * &
                        (sgrid%sp(0)+ (10**(10.0_wp*(-1.0_wp+REAL(i_s-1,wp)/REAL(n_s-1,wp)))) &
                                                   *(1.0e-06_wp+REAL(n_s+1,wp))/(2.0e-06_wp+REAL(n_s+1,wp))*sgrid%ds(1))
  END DO
  
  nVal=1
  Varnames(   nVal)='rho'
  values_visu(nVal,:)=s_visu(:)
  
  nVal=nVal+1
  Varnames(   nVal)='Phi'
  DO i=1,nvisu
    values_visu(  nVal,i)=Eval_Phi(s_visu(i))
  END DO !i
  
  !nVal=nVal+1
  !Varnames(   nVal)='chi'
  !values_visu(nVal,:)=0.0_wp !TODO
  
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
    IF((iMode.GE.base_in%f%sin_range(1)+1).AND.(iMode.LE.base_in%f%sin_range(2)))THEN
    WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname)//"_sin", &
      base_in%f%Xmn(1,iMode),base_in%f%Xmn(2,iMode)/base_in%f%nfp
    ELSE
    WRITE(VarNames(nVal),'(A,", m=",I4.3,", n=",I4.3)')TRIM(vname)//"_cos", &
      base_in%f%Xmn(1,iMode),base_in%f%Xmn(2,iMode)/base_in%f%nfp
    END IF
    DO j=1,nvisu
      val=base_in%s%evalDOF_s(s_visu(j),MAX(0,rderiv),xx_in(:,iMode))
      IF(rderiv.EQ.-5)THEN !visualize with 1/rho^m factor
        rhom=1.0_wp
        DO m=1,base_in%f%Xmn(1,iMode)
          rhom=rhom*s_visu(j)
        END DO
        values_visu(nVal,j)=ABS(val)/rhom
        !IF(ABS(val).GE.1e-18)THEN
        !  values_visu(nVal,j)=ABS(val)/rhom
        !ELSE
        !  values_visu(nVal,j)=0.0_wp
        !END IF
        !values_visu(nVal,j)=ABS(val)/(rhom+1.0e-16) + 1.0e-16
        !values_visu(nVal,j)=ABS(val)/(s_visu(j)**REAL(base_in%f%Xmn(1,iMode),wp))+1.0e-15
        !rhom=val
        !DO m=1,base_in%f%Xmn(1,iMode)
        !  rhom=rhom/s_visu(j)
        !END DO
        !values_visu(nVal,j)=rhom+1.0e-15
      ELSEIF(rderiv.EQ.-4)THEN !visualize with ABS
        values_visu(nVal,j)=ABS(val)
      ELSE
        values_visu(nVal,j)=val
      END IF
    END DO !j
  END DO

  CALL write_modes(fname,vname,nVal,base_in%f%modes,base_in%f%Xmn(1,:), &
                   base_in%f%Xmn(2,:),s_visu,sgrid%sp(1),values_visu(:,:),VarNames) 

  DEALLOCATE(varnames)
  DEALLOCATE(values_visu)
  DEALLOCATE(s_visu)
END SUBROUTINE writeDataMN_visu

END MODULE MODgvec_MHD3D_visu

