!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
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
USE MODgvec_Globals,ONLY: wp,Unit_stdOut,abort,MPIroot
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
USE MODgvec_output_netcdf,  ONLY: WriteDataToNETCDF
USE MODgvec_Output_vars,ONLY: Projectname,OutputLevel
USE MODgvec_Analyze_Vars,ONLY: outfileType
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
  INTEGER,PARAMETER  :: nVal=5
  INTEGER  :: VP_LAMBDA,VP_theta,VP_zeta,VP_X1,VP_X2
  REAL(wp) :: var_visu(nVal,mn_IP(1),mn_IP(2),1)
  REAL(wp) :: thet(mn_IP(1)),zeta(mn_IP(2))
  REAL(wp) :: rhopos
  REAL(wp) :: X1_s(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes)
  CHARACTER(LEN=40) :: VarNames(nVal)          !! Names of all variables that will be written out
  CHARACTER(LEN=255) :: FileName
  CHARACTER(LEN=255) :: var_visu_attr(nVal,2)
  CHARACTER(LEN=255) :: coord_attr(2,2)
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        "visu_BC_face should only be called by MPIroot")
  IF((minmax(2,1)-minmax(2,0)).LE.1e-08)THEN
    WRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visuBC, nothing to visualize since theta-range is <=0, theta_min= ',minmax(2,0),', theta_max= ',minmax(2,1)
    RETURN
  ELSEIF((minmax(3,1)-minmax(3,0)).LE.1e-08)THEN
    WRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visuBC, nothing to visualize since zeta-range is <=0, zeta_min= ',minmax(3,0),', zeta_max= ',minmax(3,1)
    RETURN
  END IF
  ival=1
  VP_LAMBDA=iVal;iVal=iVal+1; VarNames(VP_LAMBDA)="LA"
  VP_theta =iVal;iVal=iVal+1; VarNames(VP_theta )="theta_grid"
  VP_zeta  =iVal;iVal=iVal+1; VarNames(VP_zeta  )="zeta_grid"
  VP_X1    =iVal;iVal=iVal+1; VarNames(VP_X1  )="X1"
  VP_X2    =iVal;iVal=iVal+1; VarNames(VP_X2  )="X2"
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
  rhopos=0.99999999_wp
  DO iMode=1,X1_base%f%modes
    X1_s( iMode)= X1_base%s%evalDOF_s(rhopos,      0,U(0)%X1(:,iMode))
  END DO
  DO iMode=1,X2_base%f%modes
    X2_s(iMode) = X2_base%s%evalDOF_s(rhopos,      0,U(0)%X2(:,iMode))
  END DO
  DO iMode=1,LA_base%f%modes
    LA_s(iMode) = LA_base%s%evalDOF_s(rhopos,      0,U(0)%LA(:,iMode))
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
      var_visu(VP_X1    ,i_m,i_n,1)=X1_visu
      var_visu(VP_X2    ,i_m,i_n,1)=X2_visu
    END DO !i_m
  END DO !i_n
  nplot(:)=mn_IP-1
  WRITE(filename,'(A,"_visu_BC_",I4.4,"_",I8.8)')TRIM(Projectname),outputLevel,fileID
  IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
    CALL WriteDataToVTK(2,3,nVal,nplot,1,VarNames,coord_visu,var_visu,TRIM(filename)//".vtu")
  END IF
  IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
    var_visu_attr(VP_LAMBDA,1) = "straight field line potential";              var_visu_attr(VP_LAMBDA,2) = "\lambda"
    var_visu_attr(VP_theta,1)  = "Logical poloidal angle on the pol-tor grid"; var_visu_attr(VP_theta,2)  = "\\theta"
    var_visu_attr(VP_zeta,1)   = "Logical toroidal angle on the pol-tor grid"; var_visu_attr(VP_zeta,2)   = "\zeta"
    var_visu_attr(VP_X1,1)     = "first reference coordinate";                 var_visu_attr(VP_X1,2)     = "X^1"
    var_visu_attr(VP_X2,1)     = "second reference coordinate";                var_visu_attr(VP_X2,2)     = "X^2"

    coord_attr(1,1) = "Logical poloidal angle"; coord_attr(1,2) = "\\theta"
    coord_attr(2,1) = "Logical toroidal angle"; coord_attr(2,2) = "\zeta"
    CALL WriteDataToNETCDF(2,3,nVal,mn_IP,(/"pol","tor"/),VarNames,&
    coord_visu,var_visu, TRIM(filename), coord1=thet, coord2=zeta, attr_values=var_visu_attr, CoordNames=(/"theta","zeta "/), attr_coords=coord_attr)
  END IF

END SUBROUTINE visu_BC_face


!===================================================================================================================================
!> visualize the mapping and additional variables in 3D, either on zeta=const planes or fully 3D 
!!
!===================================================================================================================================
SUBROUTINE visu_3D(np_in,minmax,only_planes,fileID )
! MODULES
USE MODgvec_Globals,        ONLY: TWOPI,CROSS
USE MODgvec_MHD3D_vars,     ONLY: X1_base,X2_base,LA_base,hmap,sgrid,U,F
USE MODgvec_MHD3D_vars,     ONLY: Phi_profile, chi_profile, iota_profile, pres_profile
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
  REAL(wp) :: rhopos,xIP(2),q(3),q_thet(3),q_zeta(3)
  REAL(wp) :: X1_s(X1_base%f%modes),F_X1_s(X1_base%f%modes),dX1ds(X1_base%f%modes)
  REAL(wp) :: X2_s(X2_base%f%modes),F_X2_s(X2_base%f%modes),dX2ds(X2_base%f%modes)
  REAL(wp) :: LA_s(LA_base%f%modes),F_LA_s(LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu,dX1_dr,dX2_dr,dX1_dthet,dX1_dzeta,dX2_dthet,dX2_dzeta
  REAL(wp) :: dLA_dthet,dLA_dzeta,iota_s,pres_s,phiPrime_s,e_s(3),e_thet(3),e_zeta(3)
#if (defined(VISU_J_FD) || defined(VISU_J_EXACT))
  INTEGER,PARAMETER  :: nVal=36
  INTEGER            :: VP_J
  REAL(wp)           :: Jcart(3)
#else
  INTEGER,PARAMETER  :: nVal=33
#endif
  INTEGER  ::VP_LAMBDA,VP_SQRTG,VP_PHI,VP_IOTA,VP_PRES,VP_dp_dr,VP_B,VP_F_X1,VP_F_X2,VP_F_LA, &
             VP_s,VP_theta,VP_zeta,VP_g_tt,VP_g_tz,VP_g_zz,VP_gr_s,VP_gr_t,VP_gr_z,VP_Mscale ,VP_MscaleF,&
             VP_Ipol,VP_Itor,VP_X1,VP_X2
  REAL(wp) :: coord_visu( 3,np_in(1),np_in(1),np_in(3),np_in(2),sgrid%nElems)
  REAL(wp) :: var_visu(nVal,np_in(1),np_in(1),np_in(3),np_in(2),sgrid%nElems)
  REAL(wp) :: var_visu_1d(nVal+3,(np_in(1)-1)*sgrid%nElems+1)
  CHARACTER(LEN=255) :: var_visu_attr(nVal,2) !! attributes of all variables that will be written out
  CHARACTER(LEN=255) :: coord_attr(3,2)
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
  REAL(wp) :: X1_eps,X2_eps,dX1_dr_eps,dX2_dr_eps,dX1_dthet_eps,dX1_dzeta_eps,dX2_dthet_eps,dX2_dzeta_eps
  REAL(wp) :: dLA_dthet_eps,dLA_dzeta_eps,iota_s_eps,pres_s_eps,phiPrime_s_eps
  REAL(wp) :: B_dr(3), B_dthet(3), B_dzeta(3), grad_Bcart(3, 3)          !< cartesion current density and gradient of magnetic field components
  INTEGER  :: sgn
  REAL(wp) :: delta_s,delta_thet,delta_zeta
  REAL(wp),PARAMETER :: eps   = 1.0e-8 !theta,zeta
  REAL(wp),PARAMETER :: eps_s   = 1.0e-4 !
#endif
#ifdef VISU_J_EXACT
  REAL(wp) :: dX1ds_dr(X1_base%f%modes),dX2ds_dr(X2_base%f%modes),dLAds(LA_base%f%modes)
  REAL(wp) :: phiPrime_dr,iota_dr
  REAL(wp) :: dX1_dr_dr,dX1_dr_dthet,dX1_dr_dzeta,dX1_dthet_dthet,dX1_dthet_dzeta,dX1_dzeta_dzeta
  REAL(wp) :: dX2_dr_dr,dX2_dr_dthet,dX2_dr_dzeta,dX2_dthet_dthet,dX2_dthet_dzeta,dX2_dzeta_dzeta 
  REAL(wp) ::           dLA_dr_dthet,dLA_dr_dzeta,dLA_dthet_dthet,dLA_dthet_dzeta,dLA_dzeta_dzeta 
  REAL(wp) :: dBthet_dr,dBthet_dthet,dBthet_dzeta,dBzeta_dr,dBzeta_dthet,dBzeta_dzeta
  REAL(wp),DIMENSION(3)::q_s,q_s_s,q_s_thet,q_s_zeta,q_thet_thet,q_thet_zeta,q_zeta_zeta
  REAL(wp) :: Jh,Jh_dq1,Jh_dq2,dJh_dr,dJh_dthet,dJh_dzeta
  REAL(wp) :: Jp              ,dJp_dr,dJp_dthet,dJp_dzeta
  REAL(wp) :: dsqrtg_dr,dsqrtg_dthet,dsqrtg_dzeta
  REAL(wp) :: g_st,g_st_dq1,g_st_dq2         ,dg_st_dthet,dg_st_dzeta
  REAL(wp) :: g_sz,g_sz_dq1,g_sz_dq2         ,dg_sz_dthet,dg_sz_dzeta
  REAL(wp) :: g_tt,g_tt_dq1,g_tt_dq2,dg_tt_dr,dg_tt_dzeta
  REAL(wp) :: g_tz,g_tz_dq1,g_tz_dq2,dg_tz_dr,dg_tz_dthet,dg_tz_dzeta
  REAL(wp) :: g_zz,g_zz_dq1,g_zz_dq2,dg_zz_dr            ,dg_zz_dthet
  REAL(wp) :: dBsubs_dthet,dBsubs_dzeta,dBsubthet_dr,dBsubthet_dzeta,dBsubzeta_dr,dBsubzeta_dthet

  REAL(wp) :: Js,Jthet,Jzeta
#endif
  REAL(wp),ALLOCATABLE :: tmpcoord(:,:,:,:),tmpvar(:,:,:,:), coord1(:), coord2(:), coord3(:)
  INTEGER :: tmp_nrho, tmp_ntheta, tmp_nzeta
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        "visu_3D should only be called by MPIroot")
  IF(only_planes)THEN
    WRITE(UNIT_stdOut,'(A)') 'Start visu planes...'
  ELSE
    WRITE(UNIT_stdOut,'(A)') 'Start visu 3D...'
  END IF
  IF((minmax(1,1)-minmax(1,0)).LE.1e-08)THEN
    WRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
     'WARNING visu3D, nothing to visualize since s-range is <=0, s_min= ',minmax(1,0),', s_max= ',minmax(1,1)
    RETURN
  ELSEIF((minmax(2,1)-minmax(2,0)).LE.1e-08)THEN
    WRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visu3D, nothing to visualize since theta-range is <=0, theta_min= ',minmax(2,0),', theta_max= ',minmax(2,1)
    RETURN
  ELSEIF((minmax(3,1)-minmax(3,0)).LE.1e-08)THEN
    WRITE(UNIT_stdOut,'(A,F6.3,A,F6.3)') &
      'WARNING visu3D, nothing to visualize since zeta-range is <=0, zeta_min= ',minmax(3,0),', zeta_max= ',minmax(3,1)
    RETURN
  END IF
  __PERFON("output_visu")
  __PERFON("prepare_visu")
  iVal=1
  VP_s      =iVal;iVal=iVal+1; VarNames(VP_s     )="rho"
  VP_theta  =iVal;iVal=iVal+1; VarNames(VP_theta )="theta"
  VP_zeta   =iVal;iVal=iVal+1; VarNames(VP_zeta  )="zeta"
  VP_PHI    =iVal;iVal=iVal+1; VarNames(VP_PHI   )="Phi"
  VP_IOTA   =iVal;iVal=iVal+1; VarNames(VP_IOTA  )="iota"
  VP_PRES   =iVal;iVal=iVal+1; VarNames(VP_PRES  )="p"
  VP_DP_DR  =iVal;iVal=iVal+1; VarNames(VP_DP_DR )="dp_dr"
  VP_Mscale =iVal;iVal=iVal+1; VarNames(VP_Mscale)="Mscale"
  VP_MscaleF=iVal;iVal=iVal+1; VarNames(VP_MscaleF)="MscaleForce"
  VP_LAMBDA =iVal;iVal=iVal+1; VarNames(VP_LAMBDA)="LA"
  VP_SQRTG  =iVal;iVal=iVal+1; VarNames(VP_SQRTG )="sqrtG"
  VP_g_tt   =iVal;iVal=iVal+1; VarNames(VP_g_tt  )="g_tt"
  VP_g_tz   =iVal;iVal=iVal+1; VarNames(VP_g_tz  )="g_tz"
  VP_g_zz   =iVal;iVal=iVal+1; VarNames(VP_g_zz  )="g_zz"
  VP_B      =iVal;iVal=iVal+3; VarNames(VP_B     )="BX"
                               VarNames(VP_B+1   )="BY"
                               VarNames(VP_B+2   )="BZ"
  VP_X1     =iVal;iVal=iVal+1; VarNames(VP_X1    )="X1"
  VP_X2     =iVal;iVal=iVal+1; VarNames(VP_X2    )="X2"
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

  ! Set NETCDF attributes
  var_visu_attr = ""

  var_visu_attr(VP_Itor,1)   = "toroidal current";                                  var_visu_attr(VP_Itor,2)   = "I_{tor}"
  var_visu_attr(VP_Ipol,1)   = "poloidal current";                                  var_visu_attr(VP_Ipol,2)   = "I_{pol}"

  var_visu_attr(VP_gr_z+2,1) = "toroidal reciprocal basis vector";                  var_visu_attr(VP_gr_z+2,2) = "\\nabla\\zeta"
  var_visu_attr(VP_gr_t+2,1) = "poloidal reciprocal basis vector";                  var_visu_attr(VP_gr_t+2,2) = "\\nabla\\theta"
  var_visu_attr(VP_gr_s+2,1) = "radial reciprocal basis vector";                    var_visu_attr(VP_gr_s+2,2) = "\\nabla\\rho"

  var_visu_attr(VP_B+2,1)    = "magnetic field";                                    var_visu_attr(VP_B+2,2)    = "\mathbf{B}"

  var_visu_attr(VP_g_tt,1)   = "Metric coefficient theta/theta";                    var_visu_attr(VP_g_tt,2)   = "g_{\\theta,\\theta}"
  var_visu_attr(VP_g_tz,1)   = "Metric coefficient theta/zeta";                     var_visu_attr(VP_g_tz,2)   = "g_{\\theta,\zeta}"
  var_visu_attr(VP_g_zz,1)   = "Metric coefficient zeta/zeta";                      var_visu_attr(VP_g_zz,2)   = "g_{\\zeta,\zeta}"

  var_visu_attr(VP_F_X1,1)   = "Force residual in X1";                              var_visu_attr(VP_F_X1,2)   = "F_{X^1}"
  var_visu_attr(VP_F_X2,1)   = "Force residual in X2";                              var_visu_attr(VP_F_X2,2)   = "F_{X^2}"
  var_visu_attr(VP_F_LA,1)   = "Force residual in lambda";                          var_visu_attr(VP_F_LA,2)   = "F_{\lambda}"

  var_visu_attr(VP_dp_dr,1)  = "radial derivative of the pressure";                 var_visu_attr(VP_dp_dr,2)  = "\\frac{\partial p}{\partia\\rho}"
  var_visu_attr(VP_PRES,1)   = "pressure";                                          var_visu_attr(VP_PRES,2)   = "p"
  var_visu_attr(VP_IOTA,1)   = "rotational transform";                              var_visu_attr(VP_IOTA,2)   = "\iota"
  var_visu_attr(VP_PHI,1)    = "toroidal flux";                                     var_visu_attr(VP_PHI,2)    = "\phi"

  var_visu_attr(VP_SQRTG,1)  = "Jacobian determinant";                              var_visu_attr(VP_SQRTG,2)  = "\sqrt{g}"

  var_visu_attr(VP_LAMBDA,1) = "straight field line potential";                     var_visu_attr(VP_LAMBDA,2) = "\lambda"
  var_visu_attr(VP_X1,1)     = "first reference coordinate";                        var_visu_attr(VP_X1,2)     = "X^1"
  var_visu_attr(VP_X2,1)     = "second reference coordinate";                       var_visu_attr(VP_X2,2)     = "X^2"

  var_visu_attr(VP_s,1)      = "Logical radial coordinate on the rad-pol-tor grid"; var_visu_attr(VP_s,2)      = "\\rho"
  var_visu_attr(VP_theta,1)  = "Logical poloidal angle on the rad-pol-tor grid";    var_visu_attr(VP_theta,2)  = "\\theta"
  var_visu_attr(VP_zeta,1)   = "Logical toroidal angle on the rad-pol-tor grid";    var_visu_attr(VP_zeta,2)   = "\zeta"


  coord_attr(1,1) = "Logical radial coordinate"; coord_attr(1,2) = "\\rho"
  coord_attr(2,1) = "Logical poloidal angle"; coord_attr(2,2) = "\\theta"
  coord_attr(3,1) = "Logical toroidal angle"; coord_attr(3,2) = "\zeta"

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
!      rhopos=sgrid%sp(iElem-1)+(1.0e-06_wp+REAL(i_s-1,wp))/(2.0e-06_wp+REAL(n_s-1,wp))*sgrid%ds(iElem)
!      rhopos=MAX(1.0e-04,sgrid%sp(iElem-1)+1e-08+(REAL(i_s-1,wp))/(REAL(n_s-1,wp))*(sgrid%ds(iElem)-2*1e-8)) !for discont. data
      rhopos=MAX(1.0e-04,sgrid%sp(iElem-1)+(REAL(i_s-1,wp))/(REAL(n_s-1,wp))*(sgrid%ds(iElem)))
 
      X1_s(:)   = X1_base%s%evalDOF2D_s(rhopos,X1_base%f%modes,      0,U(0)%X1(:,:))
      dX1ds(:)  = X1_base%s%evalDOF2D_s(rhopos,X1_base%f%modes,DERIV_S,U(0)%X1(:,:))
      F_X1_s(:) = X1_base%s%evalDOF2D_s(rhopos,X1_base%f%modes,      0,F(0)%X1(:,:))
      X2_s(:)   = X2_base%s%evalDOF2D_s(rhopos,X2_base%f%modes,      0,U(0)%X2(:,:))
      dX2ds(:)  = X2_base%s%evalDOF2D_s(rhopos,X2_base%f%modes,DERIV_S,U(0)%X2(:,:))
      F_X2_s(:) = X2_base%s%evalDOF2D_s(rhopos,X2_base%f%modes,      0,F(0)%X2(:,:))
      LA_s(:)   = LA_base%s%evalDOF2D_s(rhopos,LA_base%f%modes,      0,U(0)%LA(:,:))
      F_LA_s(:) = LA_base%s%evalDOF2D_s(rhopos,LA_base%f%modes,      0,F(0)%LA(:,:))

      iota_s=iota_profile%eval_at_rho(rhopos)
      pres_s=pres_profile%eval_at_rho(rhopos)
      phiPrime_s=Phi_profile%eval_at_rho(rhopos,deriv=1)
#ifdef VISU_J_EXACT
      dX1ds_dr(:)  = X1_base%s%evalDOF2D_s(rhopos,X1_base%f%modes,DERIV_S_S,U(0)%X1(:,:))
      dX2ds_dr(:)  = X2_base%s%evalDOF2D_s(rhopos,X2_base%f%modes,DERIV_S_S,U(0)%X2(:,:))
      dLAds(:)     = LA_base%s%evalDOF2D_s(rhopos,LA_base%f%modes,DERIV_S  ,U(0)%LA(:,:))
      iota_dr=iota_profile%eval_at_rho(rhopos,deriv=1)
      phiPrime_dr=Phi_profile%eval_at_rho(rhopos,deriv=2)
#endif
      var_visu(VP_S    ,i_s,:,:,:,iElem) =rhopos
      var_visu(VP_PHI  ,i_s,:,:,:,iElem) =Phi_profile%eval_at_rho(rhopos)
      var_visu(VP_IOTA ,i_s,:,:,:,iElem) =iota_s
      var_visu(VP_PRES ,i_s,:,:,:,iElem) =pres_s
      var_visu(VP_DP_DR,i_s,:,:,:,iElem) =pres_profile%eval_at_rho(rhopos,deriv=1)
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
      X1_s_eps(:)   = X1_base%s%evalDOF2D_s(rhopos+delta_s,X1_base%f%modes,      0,U(0)%X1(:,:))
      dX1ds_eps(:)  = X1_base%s%evalDOF2D_s(rhopos+delta_s,X1_base%f%modes,DERIV_S,U(0)%X1(:,:))
      X2_s_eps(:)   = X2_base%s%evalDOF2D_s(rhopos+delta_s,X2_base%f%modes,      0,U(0)%X2(:,:))
      dX2ds_eps(:)  = X2_base%s%evalDOF2D_s(rhopos+delta_s,X2_base%f%modes,DERIV_S,U(0)%X2(:,:))
      LA_s_eps(:)   = LA_base%s%evalDOF2D_s(rhopos+delta_s,LA_base%f%modes,      0,U(0)%LA(:,:))
      iota_s_eps    = iota_profile%eval_at_rho(rhopos+delta_s)
      pres_s_eps    = pres_profile%eval_at_rho(rhopos+delta_s)
      phiPrime_s_eps= Phi_profile%eval_at_rho(rhopos+delta_s,deriv=1)
#endif
      !define theta2, which corresponds to the theta angle of a given theta_star=theta
      Itor_int = 0.
      Ipol_int = 0.

!$OMP PARALLEL DO COLLAPSE(3)     &  
!$OMP   SCHEDULE(STATIC) DEFAULT(NONE)    &
!$OMP   PRIVATE(i_m,i_n,j_s,xIP,q,q_thet,q_zeta,sqrtG,e_s,e_thet,e_zeta,theta_star, &
!$OMP           X1_visu,dX1_dr,dX1_dthet,dX1_dzeta,    &
!$OMP           X2_visu,dX2_dr,dX2_dthet,dX2_dzeta,dLA_dthet,dLA_dzeta, &
#ifdef VISU_J_FD
!$OMP           X1_eps,dX1_dr_eps,dX1_dthet_eps,dX1_dzeta_eps,    &
!$OMP           X2_eps,dX2_dr_eps,dX2_dthet_eps,dX2_dzeta_eps,dLA_dthet_eps,dLA_dzeta_eps, &
!$OMP           JCart, B_dr, B_dthet, B_dzeta, grad_Bcart, xIP_eps, delta_thet, delta_zeta,&
#endif
#ifdef VISU_J_EXACT
!$OMP           dX1_dr_dr,dX1_dr_dthet,dX1_dr_dzeta,dX1_dthet_dthet,dX1_dthet_dzeta,dX1_dzeta_dzeta, &
!$OMP           dX2_dr_dr,dX2_dr_dthet,dX2_dr_dzeta,dX2_dthet_dthet,dX2_dthet_dzeta,dX2_dzeta_dzeta, &
!$OMP                     dLA_dr_dthet,dLA_dr_dzeta,dLA_dthet_dthet,dLA_dthet_dzeta,dLA_dzeta_dzeta, &
!$OMP           dBthet_dr,dBthet_dthet,dBthet_dzeta,dBzeta_dr,dBzeta_dthet,dBzeta_dzeta, &
!$OMP           q_s,q_s_s,q_s_thet,q_s_zeta,q_thet_thet,q_thet_zeta,q_zeta_zeta, &
!$OMP           Jh,Jh_dq1,Jh_dq2,dJh_dr,dJh_dthet,dJh_dzeta, &
!$OMP           Jp              ,dJp_dr,dJp_dthet,dJp_dzeta, &
!$OMP           dsqrtg_dr,dsqrtg_dthet,dsqrtg_dzeta, &
!$OMP           g_st,g_st_dq1,g_st_dq2         ,dg_st_dthet,dg_st_dzeta, &
!$OMP           g_sz,g_sz_dq1,g_sz_dq2         ,dg_sz_dthet,dg_sz_dzeta, &
!$OMP           g_tt,g_tt_dq1,g_tt_dq2,dg_tt_dr,dg_tt_dzeta, &
!$OMP           g_tz,g_tz_dq1,g_tz_dq2,dg_tz_dr,dg_tz_dthet,dg_tz_dzeta, &
!$OMP           g_zz,g_zz_dq1,g_zz_dq2,dg_zz_dr            ,dg_zz_dthet, &
!$OMP           dBsubs_dthet,dBsubs_dzeta,dBsubthet_dr,dBsubthet_dzeta,dBsubzeta_dr,dBsubzeta_dthet, &
!$OMP           Js,Jthet,Jzeta,Jcart, &
#endif
!$OMP           Bcart, Bthet, Bzeta, grad_s, grad_thet, grad_zeta) &
!$OMP   REDUCTION(+:Itor_int,Ipol_int) &
!$OMP   SHARED(np_in,i_s,iElem,thet,zeta,SFL_theta,X1_base,X2_base,LA_base,X1_s,X2_s,LA_s,dX1ds,dX2ds,&
!$OMP          VP_LAMBDA,VP_SQRTG,VP_B,VP_F_X1,VP_F_X2,VP_F_LA, VP_Ipol,VP_Itor,VP_X1,VP_X2,&
!$OMP          VP_theta,VP_zeta,VP_g_tt,VP_g_tz,VP_g_zz,VP_gr_s,VP_gr_t,VP_gr_z,iota_s, &
#ifdef VISU_J_FD
!$OMP          X1_s_eps,X2_s_eps,LA_s_eps,dX1ds_eps,dX2ds_eps,VP_J,iota_s_eps,PhiPrime_s_eps,delta_s,&
#endif
#ifdef VISU_J_EXACT
!$OMP          dX1ds_dr,dX2ds_dr,dLAds,iota_dr,phiPrime_dr,VP_J, &
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
            
            dX1_dr     =X1_base%f%evalDOF_x(xIP,         0,dX1ds)
            dX2_dr     =X2_base%f%evalDOF_x(xIP,         0,dX2ds)
            
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

            e_s   =hmap%eval_dxdq(q,(/dX1_dr,dX2_dr,0.0_wp/))
            e_thet=hmap%eval_dxdq(q,q_thet)
            e_zeta=hmap%eval_dxdq(q,q_zeta)

           !sqrtG = hmap%eval_Jh(q)*(dX1_dr*dX2_dthet -dX2_dr*dX1_dthet) 
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
            var_visu(VP_X1   ,i_s,j_s,i_n,i_m,iElem) = X1_visu
            var_visu(VP_X2   ,i_s,j_s,i_n,i_m,iElem) = X2_visu
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
            dX1_dr_dr       = X1_base%f%evalDOF_x(xIP,         0,dX1ds_dr)
            dX2_dr_dr       = X2_base%f%evalDOF_x(xIP,         0,dX2ds_dr)

            dX1_dr_dthet    = X1_base%f%evalDOF_x(xIP,DERIV_THET,dX1ds )
            dX2_dr_dthet    = X2_base%f%evalDOF_x(xIP,DERIV_THET,dX2ds )
            dLA_dr_dthet    = LA_base%f%evalDOF_x(xIP,DERIV_THET,dLAds )

            dX1_dr_dzeta    = X1_base%f%evalDOF_x(xIP,DERIV_ZETA,dX1ds )
            dX2_dr_dzeta    = X2_base%f%evalDOF_x(xIP,DERIV_ZETA,dX2ds )
            dLA_dr_dzeta    = LA_base%f%evalDOF_x(xIP,DERIV_ZETA,dLAds )

            dX1_dthet_dthet = X1_base%f%evalDOF_x(xIP,DERIV_THET_THET,X1_s )
            dX2_dthet_dthet = X2_base%f%evalDOF_x(xIP,DERIV_THET_THET,X2_s )
            dLA_dthet_dthet = LA_base%f%evalDOF_x(xIP,DERIV_THET_THET,LA_s )

            dX1_dthet_dzeta = X1_base%f%evalDOF_x(xIP,DERIV_THET_ZETA,X1_s )
            dX2_dthet_dzeta = X2_base%f%evalDOF_x(xIP,DERIV_THET_ZETA,X2_s )
            dLA_dthet_dzeta = LA_base%f%evalDOF_x(xIP,DERIV_THET_ZETA,LA_s )

            dX1_dzeta_dzeta = X1_base%f%evalDOF_x(xIP,DERIV_ZETA_ZETA,X1_s )
            dX2_dzeta_dzeta = X2_base%f%evalDOF_x(xIP,DERIV_ZETA_ZETA,X2_s )
            dLA_dzeta_dzeta = LA_base%f%evalDOF_x(xIP,DERIV_ZETA_ZETA,LA_s )


            dBthet_dr    =  (iota_s - dLA_dzeta ) * phiPrime_dr + (iota_dr - dLA_dr_dzeta ) * phiPrime_s 
            dBthet_dthet =                                         (      - dLA_dthet_dzeta ) * phiPrime_s 
            dBthet_dzeta =                                         (      - dLA_dzeta_dzeta ) * phiPrime_s 
            dBzeta_dr    =  (1.0_wp + dLA_dthet ) * phiPrime_dr +             dLA_dr_dthet   * phiPrime_s 
            dBzeta_dthet =                                                  dLA_dthet_dthet   * phiPrime_s 
            dBzeta_dzeta =                                                  dLA_dthet_dzeta   * phiPrime_s 


            
            q_s          = (/dX1_dr         ,dX2_dr         ,0.0_wp/)
            q_s_s        = (/dX1_dr_dr      ,dX2_dr_dr      ,0.0_wp/)
            q_s_thet     = (/dX1_dr_dthet   ,dX2_dr_dthet   ,0.0_wp/)
            q_s_zeta     = (/dX1_dr_dzeta   ,dX2_dr_dzeta   ,0.0_wp/)
            q_thet_thet  = (/dX1_dthet_dthet,dX2_dthet_dthet,0.0_wp/)
            q_thet_zeta  = (/dX1_dthet_dzeta,dX2_dthet_dzeta,0.0_wp/)
            q_zeta_zeta  = (/dX1_dzeta_dzeta,dX2_dzeta_dzeta,0.0_wp/)

            Jh           = hmap%eval_Jh(q)
            Jh_dq1       = hmap%eval_Jh_dq1(q)
            Jh_dq2       = hmap%eval_Jh_dq2(q)
            dJh_dr       = q_s(   1)*Jh_dq1 + q_s(   2)*Jh_dq2
            dJh_dthet    = q_thet(1)*Jh_dq1 + q_thet(2)*Jh_dq2
            dJh_dzeta    = q_zeta(1)*Jh_dq1 + q_zeta(2)*Jh_dq2

            Jp           = q_s(     1)*q_thet(     2)-q_s(     2)*q_thet(     1)
            dJp_dr       = q_s_s(   1)*q_thet(     2)-q_s_s(   2)*q_thet(     1) &
                          +q_s(     1)*q_s_thet(   2)-q_s(     2)*q_s_thet(   1)
            dJp_dthet    = q_s_thet(1)*q_thet(     2)-q_s_thet(2)*q_thet(     1) &
                          +q_s(     1)*q_thet_thet(2)-q_s(     2)*q_thet_thet(1)
            dJp_dzeta    = q_s_zeta(1)*q_thet(     2)-q_s_zeta(2)*q_thet(     1) &
                          +q_s(     1)*q_thet_zeta(2)-q_s(     2)*q_thet_zeta(1) 
             
            sqrtg        = Jh*Jp
            dsqrtg_dr    = Jh*dJp_dr    + dJh_dr   *Jp 
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
            dg_tt_dr     = 2.0_wp*hmap%eval_gij(q_s_thet   ,q,q_thet) + q_s(   1)*g_tt_dq1 + q_s(   2)*g_tt_dq2
            dg_tt_dzeta  = 2.0_wp*hmap%eval_gij(q_thet_zeta,q,q_thet) + q_zeta(1)*g_tt_dq1 + q_zeta(2)*g_tt_dq2

            g_tz         = hmap%eval_gij(     q_thet,q,q_zeta)
            g_tz_dq1     = hmap%eval_gij_dq1( q_thet,q,q_zeta)
            g_tz_dq2     = hmap%eval_gij_dq2( q_thet,q,q_zeta)
            dg_tz_dr     = hmap%eval_gij(q_s_thet   ,q,q_zeta) + hmap%eval_gij(q_thet,q,q_s_zeta   ) + q_s(   1) * g_tz_dq1 + q_s(   2) * g_tz_dq2 
            dg_tz_dthet  = hmap%eval_gij(q_thet_thet,q,q_zeta) + hmap%eval_gij(q_thet,q,q_thet_zeta) + q_thet(1) * g_tz_dq1 + q_thet(2) * g_tz_dq2 
            dg_tz_dzeta  = hmap%eval_gij(q_thet_zeta,q,q_zeta) + hmap%eval_gij(q_thet,q,q_zeta_zeta) + q_zeta(1) * g_tz_dq1 + q_zeta(2) * g_tz_dq2 

            g_zz         = hmap%eval_gij(            q_zeta,q,q_zeta)
            g_zz_dq1     = hmap%eval_gij_dq1(        q_zeta,q,q_zeta)
            g_zz_dq2     = hmap%eval_gij_dq2(        q_zeta,q,q_zeta)
            dg_zz_dr     = 2.0_wp*hmap%eval_gij(q_s_zeta   ,q,q_zeta) + q_s(   1)*g_zz_dq1 + q_s(   2)*g_zz_dq2
            dg_zz_dthet  = 2.0_wp*hmap%eval_gij(q_thet_zeta,q,q_zeta) + q_thet(1)*g_zz_dq1 + q_thet(2)*g_zz_dq2

            dBsubs_dthet     = (-(Bthet*g_st+Bzeta*g_sz)/sqrtg*dsqrtg_dthet +Bthet*dg_st_dthet + dBthet_dthet*g_st + Bzeta*dg_sz_dthet + dBzeta_dthet*g_sz)/sqrtg
            dBsubs_dzeta     = (-(Bthet*g_st+Bzeta*g_sz)/sqrtg*dsqrtg_dzeta +Bthet*dg_st_dzeta + dBthet_dzeta*g_st + Bzeta*dg_sz_dzeta + dBzeta_dzeta*g_sz)/sqrtg

            dBsubthet_dr     = (-(Bthet*g_tt+Bzeta*g_tz)/sqrtg*dsqrtg_dr    +Bthet*dg_tt_dr    + dBthet_dr   *g_tt + Bzeta*dg_tz_dr    + dBzeta_dr   *g_tz)/sqrtg
            dBsubthet_dzeta  = (-(Bthet*g_tt+Bzeta*g_tz)/sqrtg*dsqrtg_dzeta +Bthet*dg_tt_dzeta + dBthet_dzeta*g_tt + Bzeta*dg_tz_dzeta + dBzeta_dzeta*g_tz)/sqrtg

            dBsubzeta_dr     = (-(Bthet*g_tz+Bzeta*g_zz)/sqrtg*dsqrtg_dr    +Bthet*dg_tz_dr    + dBthet_dr   *g_tz + Bzeta*dg_zz_dr    + dBzeta_dr   *g_zz)/sqrtg
            dBsubzeta_dthet  = (-(Bthet*g_tz+Bzeta*g_zz)/sqrtg*dsqrtg_dthet +Bthet*dg_tz_dthet + dBthet_dthet*g_tz + Bzeta*dg_zz_dthet + dBzeta_dthet*g_zz)/sqrtg 
           
            
            Js     = (dBsubzeta_dthet - dBsubthet_dzeta)
            Jthet  = (dBsubs_dzeta    - dBsubzeta_dr   )
            Jzeta  = (dBsubthet_dr    - dBsubs_dthet   )

            Jcart(:) = (e_s(:)*Js+ e_thet(:) * Jthet + e_zeta(:) * Jzeta)/sqrtg
            var_visu(VP_J:VP_J+2,i_s,j_s,i_n,i_m,iElem) = Jcart(:)/(2.0e-7_wp*TWOPI)  !*1/mu_0
#endif /*VISU_J_EXACT*/


#ifdef VISU_J_FD
            ! Get J components - finite difference bases
  
            
            ! Calculate ds derivative of B
            X1_eps         = X1_base%f%evalDOF_x(xIP, 0, X1_s_eps )
            X2_eps         = X2_base%f%evalDOF_x(xIP, 0, X2_s_eps )
            dX1_dr_eps     = X1_base%f%evalDOF_x(xIP, 0, dX1ds_eps)
            dX2_dr_eps     = X2_base%f%evalDOF_x(xIP, 0, dX2ds_eps)
            
            dX1_dthet_eps  = X1_base%f%evalDOF_x(xIP, DERIV_THET,X1_s_eps )
            dX2_dthet_eps  = X2_base%f%evalDOF_x(xIP, DERIV_THET,X2_s_eps )
            dLA_dthet_eps  = LA_base%f%evalDOF_x(xIP, DERIV_THET,LA_s_eps )
            
            dX1_dzeta_eps  = X1_base%f%evalDOF_x(xIP,DERIV_ZETA,X1_s_eps )
            dX2_dzeta_eps  = X2_base%f%evalDOF_x(xIP,DERIV_ZETA,X2_s_eps )
            dLA_dzeta_eps  = LA_base%f%evalDOF_x(xIP,DERIV_ZETA,LA_s_eps )

            q        = (/ X1_eps, X2_eps, xIP(2) /) !(X1,X2,zeta)
            e_s      = hmap%eval_dxdq(q,(/dX1_dr_eps   ,dX2_dr_eps   , 0.0    /)) !dxvec/ds
            e_thet   = hmap%eval_dxdq(q,(/dX1_dthet_eps,dX2_dthet_eps, 0.0    /)) !dxvec/dthet
            e_zeta   = hmap%eval_dxdq(q,(/dX1_dzeta_eps,dX2_dzeta_eps, 1.0_wp /)) !dxvec/dzeta
            sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))

            Bthet   = (iota_s_eps - dLA_dzeta_eps ) * phiPrime_s_eps   !/sqrtG
            Bzeta   = (1.0_wp  + dLA_dthet_eps ) * phiPrime_s_eps       !/sqrtG
            B_dr(:) =  (( e_thet(:) * Bthet + e_zeta(:) * Bzeta) /sqrtG - Bcart(:)) / (delta_s)      ! calculating dBx_dr, dBy_dr, dBz_dr
            
            ! Calculate dtheta derivative of B
            delta_thet = eps*SQRT(SUM(grad_thet*grad_thet))
            xIP_eps        = (/xIP(1)+delta_thet, xIP(2)/)
            X1_eps         = X1_base%f%evalDOF_x(xIP_eps, 0, X1_s )
            X2_eps         = X2_base%f%evalDOF_x(xIP_eps, 0, X2_s )
            dX1_dr_eps     = X1_base%f%evalDOF_x(xIP_eps, 0, dX1ds)
            dX2_dr_eps     = X2_base%f%evalDOF_x(xIP_eps, 0, dX2ds)
            
            dX1_dthet_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_THET, X1_s )
            dX2_dthet_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_THET, X2_s )
            dLA_dthet_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_THET, LA_s )
            
            dX1_dzeta_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X1_s )
            dX2_dzeta_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X2_s )
            dLA_dzeta_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, LA_s )

            q        = (/ X1_eps, X2_eps, xIP_eps(2) /) !(X1,X2,zeta)
            e_s      = hmap%eval_dxdq(q,(/dX1_dr_eps   ,dX2_dr_eps   , 0.0    /)) !dxvec/ds
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
            dX1_dr_eps     = X1_base%f%evalDOF_x(xIP_eps, 0, dX1ds)
            dX2_dr_eps     = X2_base%f%evalDOF_x(xIP_eps, 0, dX2ds)
            
            dX1_dthet_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_THET, X1_s )
            dX2_dthet_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_THET, X2_s )
            dLA_dthet_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_THET, LA_s )
            
            dX1_dzeta_eps  = X1_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X1_s )
            dX2_dzeta_eps  = X2_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, X2_s )
            dLA_dzeta_eps  = LA_base%f%evalDOF_x(xIP_eps, DERIV_ZETA, LA_s )

            q        = (/ X1_eps, X2_eps, xIP_eps(2) /) !(X1,X2,zeta)
            e_s      = hmap%eval_dxdq(q,(/dX1_dr_eps   ,dX2_dr_eps   , 0.0    /)) !dxvec/ds
            e_thet   = hmap%eval_dxdq(q,(/dX1_dthet_eps,dX2_dthet_eps, 0.0    /)) !dxvec/dthet
            e_zeta   = hmap%eval_dxdq(q,(/dX1_dzeta_eps,dX2_dzeta_eps, 1.0_wp /)) !dxvec/dzeta
            sqrtG    = SUM(e_s * (CROSS(e_thet,e_zeta)))

            Bthet   = (iota_s - dLA_dzeta_eps ) * phiPrime_s   !/sqrtG
            Bzeta   = (1.0_wp  + dLA_dthet_eps ) * phiPrime_s       !/sqrtG
            B_dzeta(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / (delta_zeta)    ! calculating dBx_dzeta, dBy_dzeta, dBz_dzeta

            ! Calculate B derivatives by finite difference
            grad_Bcart(1, :) = B_dr(1) * grad_s(:) + B_dthet(1) * grad_thet(:) + B_dzeta(1) * grad_zeta(:)   ! grad_BX
            grad_Bcart(2, :) = B_dr(2) * grad_s(:) + B_dthet(2) * grad_thet(:) + B_dzeta(2) * grad_zeta(:)   ! grad_BY
            grad_Bcart(3, :) = B_dr(3) * grad_s(:) + B_dthet(3) * grad_thet(:) + B_dzeta(3) * grad_zeta(:)   ! grad_BZ 

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
      tmp_nrho = (maxElem-minElem+1)*(n_s-1)+1
      ALLOCATE(coord1(tmp_nrho))
      coord1 =-42.0_wp
      coord1(1:n_s) = var_visu(VP_s,:,1,1,1,1)
      coord1(tmp_nrho-n_s+2:tmp_nrho) = var_visu(VP_s,2:,1,1,1,nElems)
      i = n_s+1
      DO iElem=2,nElems-1
        coord1(i:i+n_s-2) = var_visu(VP_s,2:n_s,1,1,1,iElem)
        i = i+n_s-1
      END DO

      tmp_ntheta = mn_IP(1)*(n_s-1)+1
      ALLOCATE(coord2(tmp_ntheta))
      coord2 =-42.0_wp
      coord2(1:n_s) = var_visu(VP_theta,1,:,1,1,1)
      coord2(tmp_ntheta-n_s+2:tmp_ntheta) = var_visu(VP_theta,1,2:,1,mn_IP(1),1)
      i = n_s+1
      DO i_m=2,mn_IP(1)-1
        coord2(i:i+n_s-2) = var_visu(VP_theta,1,2:n_s,1,i_m,1)
        i = i+n_s-1
      END DO

      coord3 = var_visu(VP_zeta,1,1,:,1,1)
      VarNames(VP_s)   = "rho_grid"
      VarNames(VP_zeta)  = "zeta_grid"
      VarNames(VP_theta) = "theta_grid"
      
      CALL WriteDataToNETCDF(3,3,nVal,(/tmp_nrho,tmp_ntheta,mn_IP(2)/),&
                          (/"rad","pol","tor"/),VarNames, &
                          tmpcoord,tmpvar, TRIM(filename),attr_values=var_visu_attr, &
                          coord1=coord1, coord2=coord2, coord3=coord3,CoordNames=(/"rho  ", "theta", "zeta "/), attr_coords=coord_attr)
      DEALLOCATE(tmpcoord,tmpvar, coord1, coord2, coord3)
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
  ALLOCATE(coord1((n_s-1)*nElems+1))
  coord1=var_visu_1d(4,:)
  CALL WriteDataToNETCDF(1,3,nVal-3,(/(n_s-1)*nElems+1/),(/"rad"/), &
       VarNames(4:nVal),var_visu_1d(1:3,:),var_visu_1d(7:3+nVal,:), TRIM(filename), coord1=coord1, CoordNames=(/"rho"/))
  DEALLOCATE(coord1)
#else
  CALL WriteDataToCSV((/CoordNames,VarNames(:)/) ,var_visu_1d,TRIM(filename)//".csv"  &
                                  ,append_in=.FALSE.,vfmt_in='E15.5')
#endif
  
  WRITE(UNIT_stdOut,'(A)') '... DONE.'
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
  REAL(wp)     ,INTENT(IN ) :: LA(1:LA_base_in%f%modes)   !< modes of Lambda (for a fixed rhopos)
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
  USE MODgvec_MHD3D_vars,     ONLY: Phi_profile, iota_profile
  USE MODgvec_fBase,          ONLY: t_fbase,sin_cos_map
  USE MODgvec_Transform_SFL,  ONLY: find_pest_angles
  USE MODgvec_SFL_Boozer,     ONLY: t_sfl_boozer,sfl_boozer_new
  USE MODgvec_output_netcdf,  ONLY: WriteDataToNETCDF
  USE MODgvec_output_vtk,     ONLY: WriteDataToVTK
  USE MODgvec_Output_vars,    ONLY: ProjectName,outputLevel
  USE MODgvec_Analyze_Vars,   ONLY: outfileType,SFLout,SFLout_nrp,SFLout_mn_pts,SFLout_mn_max,&
                                    SFLout_radialpos,SFLout_endpoint,SFLout_relambda
  USE MODgvec_sol_var_MHD3D,  ONLY: t_sol_var_mhd3d
  USE MODgvec_Globals,        ONLY: TWOPI,CROSS
  IMPLICIT NONE 
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CLASS(t_sol_var_MHD3D), INTENT(IN   ) :: Uin      !! input solution 
  INTEGER , INTENT(IN   ) :: fileID          !! added to file name before the ending
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  TYPE(t_sfl_boozer),ALLOCATABLE :: sfl_booz
  REAL(wp),ALLOCATABLE       :: coord_out(:,:,:,:),var_out(:,:,:,:),tz_pos(:,:,:,:),tz_star_pos(:,:,:)
  INTEGER                    :: i_rp,izeta,ithet,nthet_out,nzeta_out,i
  INTEGER                    :: mn_max(2),factorSFL,iVal,nfp
  REAL(wp)                   :: xp(2),sqrtG
  REAL(wp)                   :: dX1ds,dX2ds
  REAL(wp)                   :: phiPrime_int,iota_int,Itor_int,Ipol_int
  REAL(wp)                   :: X1_int,X2_int,LA_int,nu_int,dLA_dthet,dLA_dzeta
  REAL(wp)                   :: dnu_dthet,dnu_dzeta,dX1dthet,dX1dzeta,dX2dthet,dX2dzeta
  REAL(wp)                   :: dthetstar_dthet ,dthetstar_dzeta ,dzetastar_dthet ,dzetastar_dzeta,Jstar
  REAL(wp)                   :: dthet_dthetstarJ,dthet_dzetastarJ,dzeta_dthetstarJ,dzeta_dzetastarJ
  REAL(wp)                   :: Bthet,Bzeta
  REAL(wp),DIMENSION(3)      :: qvec,e_s,e_thet,e_zeta,e_thetstar,e_zetastar,Bfield
  REAL(wp),ALLOCATABLE       :: X1_s(:),dX1ds_s(:),X2_s(:),dX2ds_s(:)
  INTEGER                    :: VP_rho,VP_iota,VP_phip,VP_thetastar,VP_zetastar,VP_zeta,VP_nu,VP_lambda,VP_SQRTG,&
                                VP_SQRTGstar,VP_B,VP_modB,VP_gradrho,VP_etstar,VP_ezstar,VP_theta,VP_Itor,VP_Ipol,VP_X1,VP_X2
  INTEGER,PARAMETER          :: nVal=28
  CHARACTER(LEN=40)          :: VarNames(nval)  
  CHARACTER(LEN=10)          :: sfltype 
  CHARACTER(LEN=255)         :: filename
  INTEGER                    :: k,sflouts(2),whichSFLout
  REAL(wp)                   :: rho_pos(SFLout_nrp),iota_prof(SFLout_nrp),PhiPrime_prof(SFLout_nrp)
  REAL(wp),ALLOCATABLE        :: LA_s(:,:)
  !=================================================================================================================================
  IF(.NOT. MPIroot) RETURN
  IF(SFLout.EQ.12) THEN
     sflouts=(/1,2/)
  ELSE
     sflouts=(/SFLout,-1/)
  END IF
  !!!!! LOOP OVER WHICH SFL OUTPUT
  DO k=1,2
    whichSFLout=sflouts(k)
    IF(whichSFLout.EQ.-1) CYCLE
  
    SELECT CASE(whichSFLout)
    CASE(0) !GVEC angles(not a SFL coordinate)
      sfltype="_noSFL"
    CASE(1) !Pest
      sfltype="_pest"
    CASE(2) !Boozer
      sfltype="_boozer"
    END SELECT
    WRITE(filename,'(A,"_",I4.4,"_",I8.8,"")') & 
    TRIM(Projectname)//TRIM(sfltype),outputLevel,fileID
    WRITE(UNIT_stdOut,'(A,A,A)') 'WRITING SFL output: ',TRIM(filename),' ...'
    __PERFON("output_sfl")
    iVal=1
    VP_rho        =iVal;iVal=iVal+1; VarNames(VP_rho      )="rho"
    VP_iota       =iVal;iVal=iVal+1; VarNames(VP_iota     )="iota"
    VP_phip       =iVal;iVal=iVal+1; VarNames(VP_phip     )="phip"
    VP_Itor       =iVal;iVal=iVal+1; VarNames(VP_Itor     )="Itor"
    VP_Ipol       =iVal;iVal=iVal+1; VarNames(VP_Ipol     )="Ipol"
    VP_thetastar  =iVal;iVal=iVal+1; VarNames(VP_thetastar)="thetastar"
    VP_zetastar   =iVal;iVal=iVal+1; VarNames(VP_zetastar )="zetastar"
    VP_theta      =iVal;iVal=iVal+1; VarNames(VP_theta    )="theta"
    VP_zeta       =iVal;iVal=iVal+1; VarNames(VP_zeta     )="zeta"
    VP_nu         =iVal;iVal=iVal+1; VarNames(VP_nu       )="nu"
    VP_lambda     =iVal;iVal=iVal+1; VarNames(VP_lambda   )="lambda"
    VP_SQRTG      =iVal;iVal=iVal+1; VarNames(VP_SQRTG    )="sqrtG"
    VP_SQRTGstar  =iVal;iVal=iVal+1; VarNames(VP_SQRTGstar)="sqrtGstar"
    VP_modB       =iVal;iVal=iVal+1; VarNames(VP_modB     )="modB"
    VP_B          =iVal;iVal=iVal+3; VarNames(VP_B:VP_B+2 )=(/"BvecX","BvecY","BvecZ"/)
    VP_gradrho    =iVal;iVal=iVal+3; VarNames(VP_gradrho:VP_gradrho+2)=(/"grad_rhoX","grad_rhoY","grad_rhoZ"/)
    VP_etstar     =iVal;iVal=iVal+3; VarNames(VP_etstar:VP_etstar+2)=(/"e_thetastarX","e_thetastarY","e_thetastarZ"/)
    VP_ezstar     =iVal;iVal=iVal+3; VarNames(VP_ezstar:VP_ezstar+2)=(/"e_zetastarX","e_zetastarY","e_zetastarZ"/)
    VP_X1         =iVal;iVal=iVal+1; VarNames(VP_X1       )="X1"
    VP_X2         =iVal;iVal=iVal+1; VarNames(VP_X2       )="X2"
  
    IF(iVal.NE.Nval+1) CALL abort(__STAMP__,"nVal parameter not correctly set")
  
  
    factorSFL=4
    DO i=1,2
      IF(SFLout_mn_max(i).EQ.-1)THEN !input =-1, automatic
        mn_max(i) = factorSFL*MAXVAL((/X1_base%f%mn_max(i),X2_base%f%mn_max(i),LA_base%f%mn_max(i)/))
      ELSE 
        mn_max(i) = SFLout_mn_max(i) !user defined
      END IF
    END DO
    nfp=X1_base%f%nfp
    Nthet_out=MERGE(2*mn_max(1)+1,SFLout_mn_pts(1),SFLout_mn_pts(1).EQ.-1) !if input =-1, automatically 2*m_max+1, else user defined
    Nzeta_out=MERGE(2*mn_max(2)+1,SFLout_mn_pts(2),SFLout_mn_pts(2).EQ.-1) !if input =-1, automatically 2*n_max+1
    
    DO i_rp=1,SFLout_nrp
      !respect bounds
      rho_pos(i_rp)=MIN(MAX(1.0e-4_wp,SFLout_radialpos(i_rp)),1.0_wp)
      iota_prof(i_rp)=iota_profile%eval_at_rho(rho_pos(i_rp))
      PhiPrime_prof(i_rp)=Phi_profile%eval_at_rho(rho_pos(i_rp),deriv=1)
    END DO
    ALLOCATE(tz_star_pos(2,nthet_out,nzeta_out))
    DO ithet=1,Nthet_out
      tz_star_pos(1,ithet,:)=(TWOPI*(REAL(ithet,wp)-0.5_wp))/REAL(Nthet_out-MERGE(1,0,SFLout_endpoint),wp) 
    END DO
    DO izeta=1,Nzeta_out
      tz_star_pos(2,:,izeta)=(TWOPI*(REAL(izeta,wp)-0.5_wp))/REAL(((Nzeta_out-MERGE(1,0,SFLout_endpoint))*nfp),wp)
    END DO
    
    IF(SFLout_relambda .OR. (whichSFLout.EQ.2))THEN
      !for relambda=True, make use of the boozer transform computation
      SWRITE(UNIT_stdOut,'(A)')'recomputing lambda using boozer transform...'
      CALL sfl_boozer_new(sfl_booz,mn_max,4*mn_max+1,nfp, &  !recomputation of lambda with 4 times the number of modes
                          sin_cos_map(LA_base%f%sin_cos),hmap, &
                          SFLout_nrp,rho_pos,iota_prof,phiPrime_prof,&
                          relambda_in=SFLout_relambda)
      CALL sfl_booz%get_boozer(X1_base,X2_base,LA_base,Uin%X1,Uin%X2,Uin%LA)
    ELSE 
      ALLOCATE(LA_s(LA_base%f%modes,SFLout_nrp))
      DO i_rp=1,SFLout_nrp
        LA_s(:,i_rp)=LA_base%s%evalDOF2D_s(rho_pos(i_rp),LA_base%f%modes,0,Uin%LA)
      END DO
    END IF
    ALLOCATE(tz_pos(2,nthet_out,nzeta_out,SFLout_nrp))
    SELECT CASE(whichSFLout) !chooses which angles to use
    CASE(2) !Boozer
      CALL sfl_booz%find_angles(nthet_out*nzeta_out,tz_star_pos,tz_pos)
    CASE(1) !PEST  
      IF(SFLout_relambda)THEN
        CALL find_pest_angles(SFLout_nrp,sfl_booz%nu_fbase,sfl_booz%lambda,nthet_out*nzeta_out,tz_star_pos,tz_pos)
      ELSE 
        CALL find_pest_angles(SFLout_nrp,LA_base%f,LA_s,nthet_out*nzeta_out,tz_star_pos,tz_pos)
      END IF
    CASE(0) !no transform
      DO i_rp=1,SFLout_nrp
        tz_pos(:,:,:,i_rp)=tz_star_pos(:,:,:)
      END DO
    END SELECT
  
    ALLOCATE(coord_out(3,Nthet_out,Nzeta_out,SFLout_nrp),var_out(nVal,Nthet_out,Nzeta_out,SFLout_nrp))
    var_out=0.
  
    !use quantities given in GVEC theta and zeta:
      ALLOCATE(X1_s(X1_base%f%modes),dX1ds_s(X1_base%f%modes))
      ALLOCATE(X2_s(X2_base%f%modes),dX2ds_s(X2_base%f%modes))
      DO i_rp=1,SFLout_nrp
        Itor_int = 0.
        Ipol_int = 0.
        iota_int=iota_prof(i_rp)
        phiPrime_int=PhiPrime_prof(i_rp)
        var_out(VP_rho ,:,:,i_rp)=rho_pos(i_rp)
        var_out(VP_iota,:,:,i_rp)=iota_int
        
        !interpolate radially
        X1_s(   :) = X1_base%s%evalDOF2D_s(rho_pos(i_rp),X1_base%f%modes,       0,Uin%X1(:,:))
        dX1ds_s(:) = X1_base%s%evalDOF2D_s(rho_pos(i_rp),X1_base%f%modes, DERIV_S,Uin%X1(:,:))
      
        X2_s(   :) = X2_base%s%evalDOF2D_s(rho_pos(i_rp),X2_base%f%modes,       0,Uin%X2(:,:))
        dX2ds_s(:) = X2_base%s%evalDOF2D_s(rho_pos(i_rp),X2_base%f%modes, DERIV_S,Uin%X2(:,:))
        var_out(VP_thetastar,:,:,i_rp)=tz_star_pos(1,:,:)
        var_out(VP_zetastar ,:,:,i_rp)=tz_star_pos(2,:,:)  
        !HACK!
        tz_pos(:,:,:,i_rp)=tz_star_pos(:,:,:)
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(i_rp,iota_int,phiPrime_int,whichSFLout,SFLout_relambda,VP_theta,VP_zeta,VP_lambda,&
!$OMP              VP_nu,VP_SQRTG,VP_SQRTGstar,VP_modB,VP_B,VP_etstar,VP_ezstar,VP_gradrho,VP_X1,VP_X2) &
!$OMP PRIVATE(ithet,izeta,xp,LA_int,dLA_dthet,dLA_dzeta,nu_int,dnu_dthet,dnu_dzeta,&
!$OMP         dthetstar_dthet ,dthetstar_dzeta ,dzetastar_dthet ,dzetastar_dzeta,Jstar,&
!$OMP         dthet_dthetstarJ,dthet_dzetastarJ,dzeta_dthetstarJ,dzeta_dzetastarJ,&
!$OMP         Bthet,Bzeta,qvec,e_s,e_thet,e_zeta,e_thetstar,e_zetastar,Bfield,&
!$OMP         X1_int,X2_int,dX1ds,dX1dthet,dX1dzeta,dX2ds,dX2dthet,dX2dzeta,sqrtG) &
!$OMP REDUCTION(+:Itor_int,Ipol_int) &
!$OMP SHARED(Nzeta_out,Nthet_out,X1_base,X2_base,LA_base,hmap,sfl_booz,tz_pos,LA_s,X1_s,dX1ds_s,X2_s,dX2ds_s,coord_out,var_out)
        DO izeta=1,Nzeta_out; DO ithet=1,Nthet_out
          xp=tz_pos(:,ithet,izeta,i_rp) !=theta,zeta GVEC !!!
          IF(SFLout_relambda)THEN
            LA_int     = sfl_booz%nu_fbase%evalDOF_x(xp,         0, sfl_booz%lambda(:,i_rp))
            dLA_dthet  = sfl_booz%nu_fbase%evalDOF_x(xp,DERIV_THET, sfl_booz%lambda(:,i_rp))  
            dLA_dzeta  = sfl_booz%nu_fbase%evalDOF_x(xp,DERIV_ZETA, sfl_booz%lambda(:,i_rp))
            nu_int     = sfl_booz%nu_fbase%evalDOF_x(xp,         0, sfl_booz%nu(:,i_rp))
          ELSE
            LA_int    = LA_base%f%evalDOF_x(xp,          0, LA_s(:,i_rp) )
            dLA_dthet = LA_base%f%evalDOF_x(xp, DERIV_THET, LA_s(:,i_rp) )
            dLA_dzeta = LA_base%f%evalDOF_x(xp, DERIV_ZETA, LA_s(:,i_rp) )   
            nu_int    = 0.0_wp
          END IF
          SELECT CASE(whichSFLout)
          CASE(2)
            dnu_dthet  = sfl_booz%nu_fbase%evalDOF_x(xp,DERIV_THET, sfl_booz%nu(:,i_rp))  
            dnu_dzeta  = sfl_booz%nu_fbase%evalDOF_x(xp,DERIV_ZETA, sfl_booz%nu(:,i_rp))
            dthetstar_dthet=1.+dLA_dthet + iota_int*dnu_dthet
            dthetstar_dzeta=   dLA_dzeta + iota_int*dnu_dzeta
            dzetastar_dthet=   dnu_dthet
            dzetastar_dzeta=1.+dnu_dzeta
          CASE(1)  
            dthetstar_dthet=1.+dLA_dthet
            dthetstar_dzeta=   dLA_dzeta
            dzetastar_dthet=0.
            dzetastar_dzeta=1.
          CASE(0)!no Transform  
            dthetstar_dthet=1.
            dthetstar_dzeta=0.
            dzetastar_dthet=0.
            dzetastar_dzeta=1.
          END SELECT
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
           
          ! !transform derivative from dthet,dzeta=>dthet*,dzeta*
          ! dX1dthetstar = (dX1dthet*dthet_dthetstarJ+dX1dzeta*dzeta_dthetstarJ)/Jstar
          ! dX2dthetstar = (dX2dthet*dthet_dthetstarJ+dX2dzeta*dzeta_dthetstarJ)/Jstar
  
          ! dX1dzetastar = (dX1dthet*dthet_dzetastarJ+dX1dzeta*dzeta_dzetastarJ)/Jstar
          ! dX2dzetastar = (dX2dthet*dthet_dzetastarJ+dX2dzeta*dzeta_dzetastarJ)/Jstar
          ! IF(whichSFLout.EQ.2)THEN
          !   dnu_dthetstar=(dnu_dthet*dthet_dthetstarJ+dnu_dzeta*dzeta_dthetstarJ)/Jstar
          !   dnu_dzetastar=(dnu_dthet*dthet_dzetastarJ+dnu_dzeta*dzeta_dzetastarJ)/Jstar
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
          ! !e_s          = hmap%eval_dxdq(qvec,(/dX1ds       ,dX2ds       ,       -dnuds       /)) !dxvec/ds
          ! e_thetstar   = hmap%eval_dxdq(qvec,(/dX1dthetstar,dX2dthetstar,       -dnu_dthetstar/)) !dxvec/dthetstar
          ! e_zetastar   = hmap%eval_dxdq(qvec,(/dX1dzetastar,dX2dzetastar,1.0_wp -dnu_dzetastar/)) !dxvec/dzetastar
          ! !sqrtG        = SUM(e_s*(CROSS(e_thetstar,e_zetastar)))
          ! sqrtG        = hmap%eval_Jh(qvec)*(dX1ds*dX2dthetstar-dX1dthetstar*dX2ds)
          ! Bthetstar    = iota_int*PhiPrime_int   !/sqrtG
          ! Bzetastar    =          PhiPrime_int   !/sqrtG
          ! Bfield(:)    =  ( e_thetstar(:)*Bthetstar+e_zetastar(:)*Bzetastar) /sqrtG
          
          var_out(VP_theta    ,ithet,izeta,i_rp)=xp(1)
          var_out(VP_zeta     ,ithet,izeta,i_rp)=xp(2)
          var_out(VP_lambda   ,ithet,izeta,i_rp)=LA_int
          var_out(VP_nu       ,ithet,izeta,i_rp)=nu_int
          var_out(VP_SQRTG    ,ithet,izeta,i_rp)=sqrtG
          var_out(VP_SQRTGstar,ithet,izeta,i_rp)=sqrtG/Jstar !=sqrtGstar
          var_out(VP_B:VP_B+2 ,ithet,izeta,i_rp)=Bfield
          var_out(VP_modB     ,ithet,izeta,i_rp)=SQRT(SUM(Bfield*Bfield))
          var_out(VP_gradrho:VP_gradrho+2 ,ithet,izeta,i_rp)=CROSS(e_thet,e_zeta)/sqrtG
          var_out(VP_etstar :VP_etstar+2  ,ithet,izeta,i_rp)=e_thetstar
          var_out(VP_ezstar :VP_ezstar+2  ,ithet,izeta,i_rp)=e_zetastar
          var_out(VP_X1,ithet,izeta,i_rp)=X1_int
          var_out(VP_X2,ithet,izeta,i_rp)=X2_int
        END DO; END DO !izeta,ithet
!$OMP END PARALLEL DO
        var_out(VP_Itor ,:,:,i_rp)= Itor_int*TWOPI/(REAL((nthet_out*nzeta_out),wp)*(2.0e-7_wp*TWOPI)) !(2pi)^2/nfp /(Nt*Nz) * nfp/(2pi)
        var_out(VP_Ipol ,:,:,i_rp)= Ipol_int*TWOPI/(REAL((nthet_out*nzeta_out),wp)*(2.0e-7_wp*TWOPI))
      END DO !i_rp=1,n_rp
  
    DEALLOCATE(X1_s,dX1ds_s,X2_s,dX2ds_s,tz_pos,tz_star_pos)
    
    IF(SFLout_relambda .OR.(whichSFLout.EQ.2))THEN
      CALL sfl_booz%free(); DEALLOCATE(sfl_booz)
    ELSE 
      DEALLOCATE(LA_s)
    END IF
  
    IF((outfileType.EQ.1).OR.(outfileType.EQ.12))THEN
     CALL WriteDataToVTK(3,3,nVal,(/Nthet_out-1,Nzeta_out-1,SFLout_nrp-1/),1,VarNames, &
                        coord_out(1:3 ,1:Nthet_out,1:Nzeta_out,1:SFLout_nrp), &
                        var_out(1:nval,1:Nthet_out,1:Nzeta_out,1:SFLout_nrp),TRIM(filename)//".vtu")
    END IF
    IF((outfileType.EQ.2).OR.(outfileType.EQ.12))THEN
      CALL WriteDataToNETCDF(3,3,nVal,(/Nthet_out,Nzeta_out,SFLout_nrp/),&
                             (/"dim_theta","dim_zeta ","dim_rho  "/),VarNames, &
                             coord_out(1:3 ,1:Nthet_out,1:Nzeta_out,1:SFLout_nrp), &
                             var_out(1:nval,1:Nthet_out,1:Nzeta_out,1:SFLout_nrp), TRIM(filename))
    END IF!outfileType
    DEALLOCATE(coord_out,var_out)
    WRITE(UNIT_stdOut,'(A)') '... DONE.'
  !!! END LOOP OVER WHICH SFL OUTPUT
    __PERFOFF("output_sfl")
  END DO !k ... whichSFLout
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
  REAL(wp) :: rhopos,zeta,theta,theta0
  REAL(wp) :: UX1_s(1:X1_base%f%modes),VX1_s(1:X1_base%f%modes)
  REAL(wp) :: UX2_s(1:X2_base%f%modes),VX2_s(1:X2_base%f%modes)
  REAL(wp) :: ULA_s(1:LA_base%f%modes),VLA_s(1:LA_base%f%modes)
  REAL(wp) :: X1_visu,X2_visu,LA_visu
  REAL(wp) :: q(3),xU(3),xV(3),dist,xIP(2)
  REAL(wp),ALLOCATABLE :: theta1D(:),zeta1D(:)
  LOGICAL  :: SFL_theta=.TRUE.
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        "checkDistance should only be called by MPIroot")
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
!$OMP           UX1_s,UX2_s,ULA_s,VX1_s,VX2_s,VLA_s,rhopos,iElem,i_s) &
!$OMP   SHARED(nElems,n_s,mn_IP,theta1D,zeta1D,SFL_theta,X1_base,X2_base,LA_base,hmap,U,V,sgrid)
  DO iElem=1,nElems
    DO i_s=1,n_s
      rhopos=MAX(1.0e-06,sgrid%sp(iElem-1)+(REAL(i_s-1,wp))/(REAL(n_s,wp))*sgrid%ds(iElem)) !includes axis but not edge

      UX1_s(:) = X1_base%s%evalDOF2D_s(rhopos,X1_base%f%modes,0,U%X1(:,:))
      VX1_s(:) = X1_base%s%evalDOF2D_s(rhopos,X1_base%f%modes,0,V%X1(:,:))
      UX2_s(:) = X2_base%s%evalDOF2D_s(rhopos,X2_base%f%modes,0,U%X2(:,:))
      VX2_s(:) = X2_base%s%evalDOF2D_s(rhopos,X2_base%f%modes,0,V%X2(:,:))
      ULA_s(:) = LA_base%s%evalDOF2D_s(rhopos,LA_base%f%modes,0,U%LA(:,:))
      VLA_s(:) = LA_base%s%evalDOF2D_s(rhopos,LA_base%f%modes,0,V%LA(:,:))

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
SUBROUTINE CheckAxis(U,n_zeta,Axirhopos) 
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
  REAL(wp),INTENT(OUT)    :: Axirhopos(1:2,n_zeta)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER  :: i_n
  REAL(wp) :: zeta,UX1_s(1:X1_base%f%modes),UX2_s(1:X2_base%f%modes)
!===================================================================================================================================
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        "checkAxis should only be called by MPIroot")
  UX1_s(:) = X1_base%s%evalDOF2D_s(0.0_wp,X1_base%f%modes,0,U%X1(:,:))
  UX2_s(:) = X2_base%s%evalDOF2D_s(0.0_wp,X2_base%f%modes,0,U%X2(:,:))

  DO i_n=1,n_zeta
    zeta = TWOPI*REAL(i_n-1,wp)/REAL(n_zeta*X1_base%f%nfp,wp) !do not include periodic point 
    Axirhopos(1,i_n) = X1_base%f%evalDOF_x((/0.0_wp,zeta/),0,UX1_s(:) )
    Axirhopos(2,i_n) = X2_base%f%evalDOF_x((/0.0_wp,zeta/),0,UX2_s(:) )
  END DO !i_n
END SUBROUTINE CheckAxis

!===================================================================================================================================
!> Visualize 
!!
!===================================================================================================================================
SUBROUTINE visu_1d_modes(n_s,fileID)
! MODULES
USE MODgvec_Analyze_Vars,  ONLY: visu1D
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
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        "visu_1d_modes should only be called by MPIroot")
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
  USE MODgvec_MHD3D_Vars,    ONLY: sgrid, iota_profile, pres_profile, Phi_profile
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
    values_visu( iVar,i)=Phi_profile%eval_at_rho(s_visu(i))
  END DO !i

  iVar=iVar+1
  Varnames(iVar)='iota(Phi_norm)'
  
  DO i=1,nvisu
    values_visu(  iVar,i)=iota_profile%eval_at_rho(s_visu(i))
  END DO !i
  
  iVar=iVar+1
  Varnames(iVar)='pres(Phi_norm)'
  DO i=1,nvisu
    values_visu(  iVar,i)=pres_profile%eval_at_rho(s_visu(i))
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
  USE MODgvec_MHD3D_Vars,    ONLY: sgrid, iota_profile, pres_profile, Phi_profile
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
  IF(.NOT.MPIroot) CALL abort(__STAMP__, &
                        "writeData_MN_visu should only be called by MPIroot")
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
    values_visu(  nVal,i)=Phi_profile%eval_at_rho(s_visu(i))
  END DO !i
  
  !nVal=nVal+1
  !Varnames(   nVal)='chi'
  !values_visu(nVal,:)=0.0_wp !TODO
  
  nVal=nVal+1
  Varnames(nVal)='iota(Phi_norm)'
  
  DO i=1,nvisu
    values_visu(  nVal,i)=iota_profile%eval_at_rho(s_visu(i))
  END DO !i
  
  nVal=nVal+1
  Varnames(nVal)='pres(Phi_norm)'
  DO i=1,nvisu
    values_visu(  nVal,i)=pres_profile%eval_at_rho(s_visu(i))
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

