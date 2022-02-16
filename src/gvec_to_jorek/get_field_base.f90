!===================================================================================================================================
! Copyright (C) 2017 - 2019  Florian Hindenlang <hindenlang@gmail.com>
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
!!# Module **GET Field Basis**
!!
!! compute the 2D Fourier representation of the vector potential, magnetic field, or current density in contravariant coordinate 
!! directions.
!!
!===================================================================================================================================
MODULE MODgvec_get_field
! MODULES
USE MODgvec_Globals, ONLY:wp,abort
IMPLICIT NONE
PUBLIC

INTERFACE Get_Field_Base
  MODULE PROCEDURE Get_Field_base
END INTERFACE

CONTAINS

SUBROUTINE Get_Field_Base(mn_max,fac_nyq,field_type,vector_component,sgrid_in,field_base_out, field_representation)
! MODULES
USE MODgvec_Globals,ONLY: UNIT_stdOut,CROSS,TWOPI,PI,ProgressBar
USE MODgvec_LinAlg
USE MODgvec_base   ,ONLY: t_base,base_new
USE MODgvec_sGrid  ,ONLY: t_sgrid
USE MODgvec_fbase  ,ONLY: t_fbase,fbase_new,sin_cos_map
USE MODgvec_ReadState_vars  ,ONLY: X1_base_r,X2_base_r,LA_base_r
USE MODgvec_ReadState_vars  ,ONLY: LA_r,X1_r,X2_r 
USE MODgvec_ReadState_Vars  ,ONLY: profiles_1d,hmap_r,sbase_prof !for profiles

IMPLICIT NONE

!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER      ,INTENT(IN) :: mn_max(2)                                     !< maximum number for new variables in SFL coordinates
  INTEGER      ,INTENT(IN) :: fac_nyq                                       !< for number of integr. points  (=3...4 at least)
                                                                            !< n_IP=fac_nyq*max(mn_max_in)
  INTEGER      ,INTENT(IN) :: field_type                                    !< field to be transformed (1=B, 2=A)
  INTEGER      ,INTENT(IN) :: vector_component                              !< vector component to be transformed (1=R, 2=Z, 3=phi)
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: sgrid_in                          !< change grid for G_base_out
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base),ALLOCATABLE,INTENT(INOUT) :: field_base_out                 !< new fourier basis of function Gthet,Gzeta 
  REAL(wp)     ,ALLOCATABLE,INTENT(INOUT) :: field_representation(:,:)      !< coefficients of toroidal vector potential
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: nBase,is,iMode,modes,i_mn,mn_IP                  !< enumerators 
  INTEGER               :: mn_nyq(2),nfp,BCtype_axis(0:4)
  REAL(wp)              :: spos,dthet_dzeta
  REAL(wp)              :: Phi_int,dPhids_int,iota_int,Chi_int
  REAL(wp)              :: dPhids_int_eps,iota_int_eps
  REAL(wp)              :: theta, zeta, sqrtG
  REAL(wp)              :: xp(2),qvec(3),e_s(3), e_thet(3),e_zeta(3)                                      !< position and covariant coordinate
  REAL(wp)              :: grad_s(3), grad_thet(3),grad_zeta(3),grad_R(3),grad_Z(3)                       !< contravariant coordinates and grad R/Z vectors
  REAL(wp)              :: Acart(3), Bcart(3), Bthet, Bzeta                                               !< cartesian vector potential, magnetic field and toroidal/poloidal components
  REAL(wp)              :: Jcart(3), B_ds(3), B_dthet(3), B_dzeta(3), grad_Bcart(3, 3)                    !< cartesion current density and gradient of magnetic field components

  REAL(wp)                            :: X1_s(  1:X1_base_r%f%modes),  X1_s_eps(  1:X1_base_r%f%modes) 
  REAL(wp)                            :: dX1ds_s(1:X1_base_r%f%modes), dX1ds_s_eps(1:X1_base_r%f%modes) 
  REAL(wp)                            :: X2_s(  1:X2_base_r%f%modes),  X2_s_eps(  1:X2_base_r%f%modes) 
  REAL(wp)                            :: dX2ds_s(1:X2_base_r%f%modes), dX2ds_s_eps(1:X2_base_r%f%modes) 
  REAL(wp)                            :: LA_s(   1:LA_base_r%f%modes), LA_s_eps(   1:LA_base_r%f%modes) 
  REAL(wp),DIMENSION(:), ALLOCATABLE  :: field_representation_IP
  REAL(wp),DIMENSION(:), ALLOCATABLE  :: dLAdthet_IP,dLAdzeta_IP, dLAdthet_IP_eps,dLAdzeta_IP_eps
  REAL(wp),DIMENSION(:), ALLOCATABLE  :: X1_IP,dX1ds_IP,dX1dthet_IP,dX1dzeta_IP,X1_IP_eps,dX1ds_IP_eps,dX1dthet_IP_eps,dX1dzeta_IP_eps
  REAL(wp),DIMENSION(:), ALLOCATABLE  :: X2_IP,dX2ds_IP,dX2dthet_IP,dX2dzeta_IP,X2_IP_eps,dX2ds_IP_eps,dX2dthet_IP_eps,dX2dzeta_IP_eps
  REAL(wp),DIMENSION(:), ALLOCATABLE  :: LA_IP, LA_IP_eps
  CLASS(t_fBase),ALLOCATABLE          :: X1_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE          :: X2_fbase_nyq
  CLASS(t_fBase),ALLOCATABLE          :: LA_fbase_nyq
  
  ! Variables used in local interpolations for finite difference calculation of current density
  REAL(wp) :: X1_int,dX1ds,dX1dthet,dX1dzeta
  REAL(wp) :: X2_int,dX2ds,dX2dthet,dX2dzeta
  REAL(wp) :: dLAdthet,dLAdzeta
  
  REAL(wp) :: eps=1.0e-08                                   !< Small displacement for finite difference operations, 
                                                            !  local variables appended with _eps are used in finite different operations
  INTEGER  :: sgn
!===================================================================================================================================
  ! use maximum number of integration points from maximum mode number in both directions
  mn_nyq(1:2)=fac_nyq*MAXVAL(mn_max) 
  nfp = X1_base_r%f%nfp
  SWRITE(UNIT_StdOut,'(A,I4,3(A,2I6))')'GET FIELD, nfp=',nfp, &
                              ', mn_max_in=',LA_base_r%f%mn_max,', mn_max_out=',mn_max,', mn_int=',mn_nyq
  
  ! Initialise basis for field_representation based on existing grid representation
  CALL base_new(field_base_out,  X1_base_r%s%deg,        &
                             X1_base_r%s%continuity, &
                             sgrid_in, & !X1_base_r%s%grid,       &
                             X1_base_r%s%degGP,      &
                             mn_max,mn_nyq,nfp,      & 
                             !sin_cos_map(X1_base_r%f%sin_cos),   &
                             '_sincos_   ', &
                             X1_base_r%f%exclude_mn_zero) !exclude m=n=0
  mn_IP        = field_base_out%f%mn_IP  !total number of integration points
  modes        = field_base_out%f%modes  !number of modes in output
  nBase        = field_base_out%s%nBase  !number of radial points in output
  dthet_dzeta  = field_base_out%f%d_thet*field_base_out%f%d_zeta !integration weights

  ! Initialise bases for existing grid at higher number of integration points, based on nyquist condition
  CALL fbase_new( X1_fbase_nyq, X1_base_r%f%mn_max,  mn_nyq, &
                                X1_base_r%f%nfp, &
                    sin_cos_map(X1_base_r%f%sin_cos), &
                                X1_base_r%f%exclude_mn_zero)

  CALL fbase_new( X2_fbase_nyq, X2_base_r%f%mn_max,  mn_nyq, &
                                X2_base_r%f%nfp, &
                    sin_cos_map(X2_base_r%f%sin_cos), &
                                X2_base_r%f%exclude_mn_zero)

  CALL fbase_new(LA_fbase_nyq,  LA_base_r%f%mn_max,  mn_nyq, &
                                LA_base_r%f%nfp, &
                    sin_cos_map(LA_base_r%f%sin_cos), &
                                LA_base_r%f%exclude_mn_zero)
  
  ! Allocate arrays for storing integration point values
  ALLOCATE( X1_IP(1:mn_IP),     dX1ds_IP(1:mn_IP),     dX1dthet_IP(1:mn_IP),     dX1dzeta_IP(1:mn_IP),&
            X2_IP(1:mn_IP),     dX2ds_IP(1:mn_IP),     dX2dthet_IP(1:mn_IP),     dX2dzeta_IP(1:mn_IP),& 
            X1_IP_eps(1:mn_IP), dX1ds_IP_eps(1:mn_IP), dX1dthet_IP_eps(1:mn_IP), dX1dzeta_IP_eps(1:mn_IP),&  
            X2_IP_eps(1:mn_IP), dX2ds_IP_eps(1:mn_IP), dX2dthet_IP_eps(1:mn_IP), dX2dzeta_IP_eps(1:mn_IP))
  ALLOCATE(dLAdthet_IP(1:mn_IP), dLAdzeta_IP(1:mn_IP), dLAdthet_IP_eps(1:mn_IP), dLAdzeta_IP_eps(1:mn_IP))
  ALLOCATE(LA_IP(1:mn_IP), LA_IP_eps(1:mn_IP))
  ALLOCATE(field_representation_IP(1:mn_IP))
  ALLOCATE(field_representation(   nBase,1:modes))
  
  ! Loop over radial coordinate and evaluate modes of field_representation
  DO is=1,nBase
    ! Avoid magnetic axis and plasma boundary
    spos=MIN(MAX(1.0e-08_wp,field_base_out%s%s_IP(is)),1.0_wp-1.0e-12_wp) !interpolation points for q_in
    
    ! Evaluate grid position, derivatives and field variables at integration points and finite difference eps points
    Phi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,1))
    dPhids_int  = sbase_prof%evalDOF_s(spos, DERIV_S ,profiles_1d(:,1))
    iota_int    = sbase_prof%evalDOF_s(spos,        0,profiles_1d(:,3))
    Chi_int     = sbase_prof%evalDOF_s(spos,       0 ,profiles_1d(:,2))

    !interpolate radially
    X1_s(:)        = X1_base_r%s%evalDOF2D_s(spos    ,X1_base_r%f%modes,      0,X1_r(:,:))
    dX1ds_s(:)     = X1_base_r%s%evalDOF2D_s(spos    ,X1_base_r%f%modes,DERIV_S,X1_r(:,:))
    X2_s(:)        = X2_base_r%s%evalDOF2D_s(spos    ,X2_base_r%f%modes,      0,X2_r(:,:))
    dX2ds_s(:)     = X2_base_r%s%evalDOF2D_s(spos    ,X2_base_r%f%modes,DERIV_S,X2_r(:,:))
    
    ! Interpolate finite difference point in radial direction - direction of finite step is changed for last element to stay inside the domain
    sgn = 1
    if (is .eq. nBase) sgn=-1
    dPhids_int_eps = sbase_prof%evalDOF_s(spos+sgn*eps, DERIV_S ,profiles_1d(:,1))
    iota_int_eps   = sbase_prof%evalDOF_s(spos+sgn*eps,        0,profiles_1d(:,3))
    X1_s_eps(:)    = X1_base_r%s%evalDOF2D_s(spos+sgn*eps,X1_base_r%f%modes,      0,X1_r(:,:))
    dX1ds_s_eps(:) = X1_base_r%s%evalDOF2D_s(spos+sgn*eps,X1_base_r%f%modes,DERIV_S,X1_r(:,:))
    X2_s_eps(:)    = X2_base_r%s%evalDOF2D_s(spos+sgn*eps,X2_base_r%f%modes,      0,X2_r(:,:))
    dX2ds_s_eps(:) = X2_base_r%s%evalDOF2D_s(spos+sgn*eps,X2_base_r%f%modes,DERIV_S,X2_r(:,:))

    LA_s(:)        = LA_base_r%s%evalDOF2D_s(spos     ,LA_base_r%f%modes,      0,LA_r(:,:))
    LA_s_eps(:)    = LA_base_r%s%evalDOF2D_s(spos+sgn*eps,LA_base_r%f%modes,      0,LA_r(:,:))

    ! evaluate at integration points and finite difference eps from points
    X1_IP       = X1_fbase_nyq%evalDOF_IP(         0, X1_s(  :));     X1_IP_eps       = X1_fbase_nyq%evalDOF_IP(         0, X1_s_eps(  :))
    dX1ds_IP    = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s(:));     dX1ds_IP_eps    = X1_fbase_nyq%evalDOF_IP(         0,dX1ds_s_eps(:))
    dX1dthet_IP = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s(  :));     dX1dthet_IP_eps = X1_fbase_nyq%evalDOF_IP(DERIV_THET, X1_s_eps(  :))
    dX1dzeta_IP = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s(  :));     dX1dzeta_IP_eps = X1_fbase_nyq%evalDOF_IP(DERIV_ZETA, X1_s_eps(  :))
                                                                      
    X2_IP       = X2_fbase_nyq%evalDOF_IP(         0, X2_s(  :));     X2_IP_eps       = X2_fbase_nyq%evalDOF_IP(         0, X2_s_eps(  :))
    dX2ds_IP    = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s(:));     dX2ds_IP_eps    = X2_fbase_nyq%evalDOF_IP(         0,dX2ds_s_eps(:))
    dX2dthet_IP = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s(  :));     dX2dthet_IP_eps = X2_fbase_nyq%evalDOF_IP(DERIV_THET, X2_s_eps(  :))
    dX2dzeta_IP = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s(  :));     dX2dzeta_IP_eps = X2_fbase_nyq%evalDOF_IP(DERIV_ZETA, X2_s_eps(  :))
                                                                      
    LA_IP(:)       = LA_fbase_nyq%evalDOF_IP(         0,LA_s(:));     LA_IP_eps(:)       = LA_fbase_nyq%evalDOF_IP(         0,LA_s_eps(:))
    dLAdthet_IP(:) = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s(:));     dLAdthet_IP_eps(:) = LA_fbase_nyq%evalDOF_IP(DERIV_THET,LA_s_eps(:))
    dLAdzeta_IP(:) = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s(:));     dLAdzeta_IP_eps(:) = LA_fbase_nyq%evalDOF_IP(DERIV_ZETA,LA_s_eps(:))
    
    ! Loop over surface points and evaluate field_representation
    DO i_mn=1,mn_IP
      theta = X1_fbase_nyq%x_IP(1, i_mn)
      zeta = X1_fbase_nyq%x_IP(2, i_mn)

      ! Get the covariant basis vectors
      qvec     = (/ X1_IP(i_mn), X2_IP(i_mn), zeta /) !(X1,X2,zeta)
      e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds_IP(i_mn)   ,dX2ds_IP(i_mn)   , 0.0    /)) !dxvec/ds
      e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet_IP(i_mn),dX2dthet_IP(i_mn), 0.0    /)) !dxvec/dthet
      e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta_IP(i_mn),dX2dzeta_IP(i_mn), 1.0_wp /)) !dxvec/dzeta
      sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))
   
      ! Get contravarian basis vectors
      grad_s    = CROSS(e_thet,e_zeta) /sqrtG
      grad_thet = CROSS(e_zeta,e_s   ) /sqrtG
      grad_zeta = CROSS(e_s   ,e_thet) /sqrtG
      
      ! Get Grad R and Grad Z pol - WARNING: this implementation only works for PEST coordinates
      grad_R = dX1ds_IP(i_mn) * grad_s + dX1dthet_IP(i_mn) * grad_thet + dX1dzeta_IP(i_mn) * grad_zeta
      grad_Z = dX2ds_IP(i_mn) * grad_s + dX2dthet_IP(i_mn) * grad_thet + dX2dzeta_IP(i_mn) * grad_zeta
      
      ! Get A and X in cartesian coordinates
      Acart(:)  = (Phi_int * grad_thet(:) - (LA_IP(i_mn) * dPhids_int) * grad_s(:) - Chi_int * grad_zeta)

      ! Calculate cartesian magnetic field
      Bthet = (iota_int-dLAdzeta_IP(i_mn) )*dPhids_int   !/sqrtG
      Bzeta = (1.0_wp  +dLAdthet_IP(i_mn) )*dPhids_int   !/sqrtG
      Bcart(:) =  ( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG
        
      SELECT CASE(field_type)
      CASE(1)  ! Magnetic Field
        SELECT CASE(vector_component)
        CASE(1)
          ! Get Vertical Magnetic Field - B_X
          field_representation_IP(i_mn) = Bcart(1) 
        CASE(2)
          ! Get Vertical Magnetic Field - B_Y
          field_representation_IP(i_mn) = Bcart(2)
        CASE(3)
          ! Get Vertical Magnetic Field - B_Z
          field_representation_IP(i_mn) = Bcart(3)
        CASE(4)
          ! Get Radial Magnetic Field - B_R
          field_representation_IP(i_mn) = Bcart(1) * grad_R(1) + Bcart(2) * grad_R(2) + Bcart(3) * grad_R(3)
        CASE(5)
          ! Get Toroidal Magnetic Field - B_phi = R * (B.grad(zeta))
          field_representation_IP(i_mn) = X1_IP(i_mn) * (Bcart(1) * grad_zeta(1) + Bcart(2) * grad_zeta(2) + Bcart(3) * grad_zeta(3))
        CASE DEFAULT
          SWRITE(UNIT_StdOut,*) "Invalid vector component selected: ", vector_component
        END SELECT
      CASE(2)  ! Vector Potential
        SELECT CASE(vector_component)
        CASE(1)
          ! Get Vertical Flux - A_X
          field_representation_IP(i_mn) = Acart(1)
        CASE(2)
          ! Get Vertical Flux - A_Y
          field_representation_IP(i_mn) = Acart(2)
        CASE(3)
          ! Get Vertical Flux - A_Z
          field_representation_IP(i_mn) = Acart(3)
        CASE(4)
          ! Get Radial Flux - A_R
          field_representation_IP(i_mn) = Acart(1) * grad_R(1) + Acart(2) * grad_R(2) + Acart(3) * grad_R(3)
        CASE(5)
          ! Get Poloidal Flux - A_phi = R * (A.grad(zeta))
          field_representation_IP(i_mn) = X1_IP(i_mn) * (Acart(1) * grad_zeta(1) + Acart(2) * grad_zeta(2) + Acart(3) * grad_zeta(3))
        CASE DEFAULT
          SWRITE(UNIT_StdOut,*) "Invalid vector component selected: ", vector_component
        END SELECT
      CASE(3)  ! Current density
        ! Calculate ds derivative of B
        qvec     = (/ X1_IP_eps(i_mn), X2_IP_eps(i_mn), zeta /) !(X1,X2,zeta)
        e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds_IP_eps(i_mn)   ,dX2ds_IP_eps(i_mn)   , 0.0    /)) !dxvec/ds
        e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet_IP_eps(i_mn),dX2dthet_IP_eps(i_mn), 0.0    /)) !dxvec/dthet
        e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta_IP_eps(i_mn),dX2dzeta_IP_eps(i_mn), 1.0_wp /)) !dxvec/dzeta
        sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

        Bthet = (iota_int_eps-dLAdzeta_IP_eps(i_mn) )*dPhids_int_eps   !/sqrtG
        Bzeta = (1.0_wp  +dLAdthet_IP_eps(i_mn) )*dPhids_int_eps   !/sqrtG
        B_ds(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / (sgn*eps)      ! calculating dBx_ds, dBy_ds, dBz_ds

        ! Calculate dtheta derivative of B
        xp = (/theta+eps, zeta/)
        X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
        dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
        dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
        dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

        X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
        dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
        dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
        dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

        dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
        dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

        qvec   = (/X1_int,X2_int,zeta/)
        e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   ,dX2ds   , 0.0    /)) !dxvec/ds
        e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet,dX2dthet, 0.0    /)) !dxvec/dthet
        e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta,dX2dzeta, 1.0_wp /)) !dxvec/dzeta
        sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

        Bthet = (iota_int-dLAdzeta )*dPhids_int   !/sqrtG
        Bzeta = (1.0_wp  +dLAdthet )*dPhids_int   !/sqrtG
        B_dthet(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / eps      ! calculating dBx_dtheta, dBy_dtheta, dBz_dtheta

        ! Calculate dzeta derivative of B
        xp = (/theta, zeta+eps/)
        X1_int  =X1_base_r%f%evalDOF_x(xp,          0, X1_s  )
        dX1ds   =X1_base_r%f%evalDOF_x(xp,          0,dX1ds_s)
        dX1dthet=X1_base_r%f%evalDOF_x(xp, DERIV_THET, X1_s  )
        dX1dzeta=X1_base_r%f%evalDOF_x(xp, DERIV_ZETA, X1_s  )

        X2_int  =X2_base_r%f%evalDOF_x(xp,          0, X2_s  )
        dX2ds   =X2_base_r%f%evalDOF_x(xp,          0,dX2ds_s)
        dX2dthet=X2_base_r%f%evalDOF_x(xp, DERIV_THET, X2_s  )
        dX2dzeta=X2_base_r%f%evalDOF_x(xp, DERIV_ZETA, X2_s  )

        dLAdthet =LA_base_r%f%evalDOF_x(xp, DERIV_THET, LA_s)
        dLAdzeta =LA_base_r%f%evalDOF_x(xp, DERIV_ZETA, LA_s)

        qvec   = (/X1_int,X2_int,zeta/)
        e_s      = hmap_r%eval_dxdq(qvec,(/dX1ds   , dX2ds   , 0.0    /)) !dxvec/ds
        e_thet   = hmap_r%eval_dxdq(qvec,(/dX1dthet, dX2dthet, 0.0    /)) !dxvec/dthet
        e_zeta   = hmap_r%eval_dxdq(qvec,(/dX1dzeta, dX2dzeta, 1.0_wp /)) !dxvec/dzeta
        sqrtG    = SUM(e_s*(CROSS(e_thet,e_zeta)))

        Bthet = (iota_int-dLAdzeta )*dPhids_int   !/sqrtG
        Bzeta = (1.0_wp  +dLAdthet )*dPhids_int   !/sqrtG
        B_dzeta(:) =  (( e_thet(:)*Bthet+e_zeta(:)*Bzeta) /sqrtG - Bcart(:)) / eps      ! calculating dBx_dzeta, dBy_dzeta, dBz_dzeta

        ! Calculate B derivatives by finite difference
        grad_Bcart(1, :) = B_ds(1) * grad_s(:) + B_dthet(1) * grad_thet(:) + B_dzeta(1) * grad_zeta(:)   ! grad_BX
        grad_Bcart(2, :) = B_ds(2) * grad_s(:) + B_dthet(2) * grad_thet(:) + B_dzeta(2) * grad_zeta(:)   ! grad_BY
        grad_Bcart(3, :) = B_ds(3) * grad_s(:) + B_dthet(3) * grad_thet(:) + B_dzeta(3) * grad_zeta(:)   ! grad_BZ 

        ! Calculate current cartesian components
        Jcart(1) = grad_Bcart(3, 2) - grad_Bcart(2, 3)   ! dBZ_dY - dBY_dZ
        Jcart(2) = grad_Bcart(1, 3) - grad_Bcart(3, 1)   ! dBX_dZ - dBZ_dX
        Jcart(3) = grad_Bcart(2, 1) - grad_Bcart(1, 2)   ! dBY_dX - dBX_dY

        SELECT CASE(vector_component)
        CASE(1)
          ! Get Vertical Current Density - J_X
          field_representation_IP(i_mn) = Jcart(1)
        CASE(2)
          ! Get Vertical Current Density - J_Y
          field_representation_IP(i_mn) = Jcart(2)
        CASE(3)
          ! Get Vertical Current Density - J_Z
          field_representation_IP(i_mn) = Jcart(3)
        CASE(4)
          ! Get Radial Current Density - J_R
          field_representation_IP(i_mn) = Jcart(1) * grad_R(1) + Jcart(2) * grad_R(2) + Jcart(3) * grad_R(3)
        CASE(5)
          ! Get Toroidal Current Density - J_phi = R * (J.grad(zeta))
          field_representation_IP(i_mn) = X1_IP(i_mn) * (Jcart(1) * grad_zeta(1) + Jcart(2) * grad_zeta(2) + Jcart(3) * grad_zeta(3))
        CASE DEFAULT
          SWRITE(UNIT_StdOut,*) "Invalid vector component selected: ", vector_component
        END SELECT

      CASE DEFAULT
        SWRITE(UNIT_StdOut,*) "Invalid field type selected: ", field_type
      END SELECT

    ENDDO ! i_mn
    
    ! Convert interation points into fourier mode representation 
    CALL field_base_out%f%projectIPtoDOF(.FALSE., 1.0_wp, 0, field_representation_IP(:), field_representation(is, :))
    DO iMode=1,modes
      field_representation(is,iMode)=field_representation(is,iMode)*dthet_dzeta*field_base_out%f%snorm_base(iMode)
    END DO
  ENDDO ! is

  ! Convert radial fourier representation into radial spline
  if (vector_component .ne. 5) then
    BCtype_axis(MN_ZERO    )= BC_TYPE_DIRICHLET !=0 (should not be here!)
    BCtype_axis(M_ZERO     )= BC_TYPE_NEUMANN   ! derivative zero
    BCtype_axis(M_ODD_FIRST)= BC_TYPE_DIRICHLET !=0
    BCtype_axis(M_ODD      )= BC_TYPE_DIRICHLET !=0
    BCtype_axis(M_EVEN     )= BC_TYPE_DIRICHLET !=0
  else
    BCtype_axis(MN_ZERO    )= BC_TYPE_NEUMANN !=0 (should not be here!)
    BCtype_axis(M_ZERO     )= BC_TYPE_NEUMANN   ! derivative zero
    BCtype_axis(M_ODD_FIRST)= BC_TYPE_NEUMANN !=0
    BCtype_axis(M_ODD      )= BC_TYPE_NEUMANN !=0
    BCtype_axis(M_EVEN     )= BC_TYPE_NEUMANN !=0
  endif
  DO iMode=1, modes
    field_representation(:, iMode) = field_base_out%s%initDOF(field_representation(:, iMode))
    CALL field_base_out%s%applyBCtoDOF(field_representation(:,iMode),(/BCtype_axis(field_base_out%f%zero_odd_even(iMode)),BC_TYPE_OPEN/),(/0.,0./))
  ENDDO
  
  DEALLOCATE(X1_IP, dX1ds_IP, dX1dthet_IP, dX1dzeta_IP, X1_IP_eps, dX1ds_IP_eps, dX1dthet_IP_eps, dX1dzeta_IP_eps)
  DEALLOCATE(X2_IP, dX2ds_IP, dX2dthet_IP, dX2dzeta_IP, X2_IP_eps, dX2ds_IP_eps, dX2dthet_IP_eps, dX2dzeta_IP_eps)
  DEALLOCATE(LA_IP, dLAdthet_IP, dLAdzeta_IP, LA_IP_eps, dLAdthet_IP_eps, dLAdzeta_IP_eps)
END SUBROUTINE Get_Field_base

END MODULE MODgvec_get_field
