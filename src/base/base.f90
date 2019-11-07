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
!!# Module ** base **
!!
!! 2D Fourier base in the two angular directions: (poloidal,toroidal) ~ (m,n) ~ (theta,zeta) [0,2pi]x[0,2pi/nfp]
!!
!! explicit real fourier basis: sin(x_mn) or cos(x_mn) with x_mn=(m*theta - n*nfp*zeta)  ,
!! with mode numbers m and n 
!!
!===================================================================================================================================
MODULE MODgvec_base
! MODULES
USE MODgvec_Globals ,ONLY: wp,Unit_stdOut,abort
USE MODgvec_sBase   ,ONLY: t_sbase,sbase_new
USE MODgvec_fBase   ,ONLY: t_fbase,fbase_new
USE MODgvec_sGrid   ,ONLY: t_sgrid
IMPLICIT NONE
PUBLIC

!-----------------------------------------------------------------------------------------------------------------------------------
! TYPES 
!-----------------------------------------------------------------------------------------------------------------------------------

TYPE                 :: t_base
  !---------------------------------------------------------------------------------------------------------------------------------
  LOGICAL              :: initialized=.FALSE.      !! set to true in init, set to false in free
  !---------------------------------------------------------------------------------------------------------------------------------
  !input parameters sbase
  INTEGER              :: s_deg                      !! input parameter: degree of Spline/polynomial 
  INTEGER              :: s_degGP                    !! number of Gauss-points (degGP+1) per element >= deg
  INTEGER              :: s_continuity               !! input parameter: full spline (=deg-1) or discontinuous (=-1)
  CLASS(t_sgrid),POINTER :: s_grid                   !! pointer to grid 
  !input parameters fbase
  INTEGER              :: f_mn_max(2)   !! input parameter: maximum number of fourier modes: m_max=mn_max(1),n_max=mn_max(2)
  INTEGER              :: f_mn_nyq(2)   !! number of equidistant integration points (trapezoidal rule) in m and n
  INTEGER              :: f_nfp         !! number of field periods (toroidal repetition after 2pi/nfp)
  INTEGER              :: f_sin_cos     !! can be either only sine: _SIN_  or only cosine _COS_ or full: _SINCOS_
  !input parameters
  LOGICAL              :: f_exclude_mn_zero  !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)

  CLASS(t_sbase),ALLOCATABLE  :: s  !! container for radial basis
  CLASS(t_fbase),ALLOCATABLE  :: f  !! container for angular basis
  
  CONTAINS

  PROCEDURE :: free        => base_free
  PROCEDURE :: copy        => base_copy
  PROCEDURE :: compare     => base_compare
  PROCEDURE :: change_base => base_change_base
  PROCEDURE :: evalDOF     => base_evalDOF
  PROCEDURE :: evalDOF_x   => base_evalDOF_x

END TYPE t_base

LOGICAL  :: test_called=.FALSE.

!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> allocate and initialize the type base 
!!
!===================================================================================================================================
SUBROUTINE Base_new( sf,deg_in,continuity_in,grid_in,degGP_in, &
                     mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER       , INTENT(IN   )        :: deg_in        !! polynomial degree
  INTEGER       , INTENT(IN   )        :: continuity_in !! continuity: 
  CLASS(t_sgrid), INTENT(IN   ),TARGET :: grid_in       !! grid information
  INTEGER       , INTENT(IN   )        :: degGP_in      !! gauss quadrature points: nGP=degGP+1 
  INTEGER        , INTENT(IN   ) :: mn_max_in(2)  !! maximum mode in m and n
  INTEGER        , INTENT(IN   ) :: mn_nyq_in(2)  !! number of integration points 
  INTEGER        , INTENT(IN   ) :: nfp_in        !! number of field periods 
  CHARACTER(LEN=8),INTENT(IN   ) :: sin_cos_in    !! can be either only sine: " _sin_" only cosine: " _cos_" or full: "_sin_cos_"
  LOGICAL         ,INTENT(IN   ) :: exclude_mn_zero_in !!  =true: exclude m=n=0 mode in the basis (only important if cos is in basis)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_base), ALLOCATABLE,INTENT(INOUT)        :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  ALLOCATE(t_base :: sf)
  CALL sbase_new(sf%s,deg_in,continuity_in,grid_in,degGP_in)
  CALL fbase_new(sf%f,mn_max_in,mn_nyq_in,nfp_in,sin_cos_in,exclude_mn_zero_in)
  sf%initialized=.TRUE.

  IF(.NOT.test_called) CALL Base_test(sf)

END SUBROUTINE base_new


!===================================================================================================================================
!> finalize the type base
!!
!===================================================================================================================================
SUBROUTINE base_free( sf )
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_base), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.sf%initialized)RETURN
CALL sf%s%free()
CALL sf%f%free()

sf%initialized=.FALSE.
END SUBROUTINE base_free


!===================================================================================================================================
!> copy from input  type base to self
!!
!===================================================================================================================================
SUBROUTINE base_copy( sf , tocopy)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: tocopy
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_base), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL sf%s%copy(tocopy%s)
CALL sf%f%copy(tocopy%f)
END SUBROUTINE base_copy

!===================================================================================================================================
!> compare self and input type base
!!
!===================================================================================================================================
SUBROUTINE base_compare( sf , tocompare, is_same)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: sf !! self
  CLASS(t_base), INTENT(IN   ) :: tocompare
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  LOGICAL      , INTENT(  OUT) :: is_same
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: same_s, same_f
!===================================================================================================================================
CALL sf%s%compare(tocompare%s,is_same=same_s)
CALL sf%f%compare(tocompare%f,is_same=same_f)
is_same=same_s.AND.same_f
END SUBROUTINE base_compare

!===================================================================================================================================
!> change basis from old input base to new base, 
!!
!===================================================================================================================================
SUBROUTINE base_change_base( sf , old_base, old_data, sf_data)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: sf !! self
  CLASS(t_base), INTENT(IN   ) :: old_base
  REAL(wp)     , INTENT(IN   ) :: old_data(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)     , INTENT(  OUT) :: sf_data(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL(wp) :: tmp(old_base%s%nBase,sf%f%modes)
!===================================================================================================================================
CALL sf%f%change_base(old_base%f,1,old_data,tmp    )
CALL sf%s%change_base(old_base%s,2,tmp     ,sf_data)
END SUBROUTINE base_change_base


!===================================================================================================================================
!> evaluate all degrees of freedom at all Gauss Points (deriv=0 solution, deriv=1 first derivative d/ds)
!!
!===================================================================================================================================
SUBROUTINE base_evalDOF(sf,deriv,DOFs,y_IP_GP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: sf     !! self
  INTEGER      , INTENT(IN   ) :: deriv(2)   !! =(0,0): base, =(DERIV_S,0): ds, =(0,DERIV_THET): dtheta, =(0,DERIV_ZETA): dzeta
  REAL(wp)     , INTENT(IN   ) :: DOFs(1:sf%s%nBase,1:sf%f%modes)  !! array of all modes and all radial dofs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)     ,INTENT(OUT   ) :: y_IP_GP(sf%f%mn_IP,sf%s%nGP)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                      :: modes,nGP,nBase,mn_IP
#define PP_WHICHEVAL 110
  ! experimental optimization results, with OMP, give the ranking from slowest to fastest: 2 < 1 << 11 ~= 12 < 110  
  ! could change with the number of DOF in s and theta/zeta..., 
  ! evaluation in s does not seem to scale at all with omp
#if (PP_WHICHEVAL==1)
  !OMP loop version: first s then f
  INTEGER                      :: iMode,iGP
  REAL(wp)                     :: y_tmp(1:sf%f%modes,sf%s%nGP)

#elif (PP_WHICHEVAL==2)
  !OMP loop version: first f then s
  INTEGER                      :: i_mn,iBase
  REAL(wp)                     :: y_tmp(sf%f%mn_IP,1:sf%s%nBase)

#elif (PP_WHICHEVAL==11)
  !matrix matrix version: first s then f
  INTEGER                      :: i,j,iElem
  REAL(wp)                     :: y_tmp(sf%f%modes,1:sf%s%nGP)

#elif (PP_WHICHEVAL==110)
  !matrix matrix version: first s then f
  INTEGER                      :: i,j,iElem,iMode
  REAL(wp)                     :: y_tmp(1:sf%s%nGP,sf%f%modes)

#elif (PP_WHICHEVAL==12)
  !matrix matrix version: first f then s
  INTEGER                      :: i,j,iElem
  REAL(wp)                     :: y_tmp(sf%f%mn_IP,1:sf%s%nBase)
#elif (PP_WHICHEVAL==120)
  !matrix matrix version: first f then s
  INTEGER                      :: i,j,iElem,i_mn
  REAL(wp)                     :: y_tmp(1:sf%s%nBase,sf%f%mn_IP)
#endif /*PP_WHICHEVAL*/
!===================================================================================================================================
!WRITE(*,*)'-----------------------------'
!WRITE(*,'(A,I12)')'CHECK: #operations first s then f:', (sf%s%nGP*(sf%s%deg+1))*sf%f%modes+(sf%f%mn_IP*sf%f%modes)*sf%s%nGP
!WRITE(*,'(A,I12)')'CHECK: #operations first f then s:', (sf%f%mn_IP*sf%f%modes)*sf%s%nBase+ (sf%s%nGP*(sf%s%deg+1))*sf%f%mn_IP
!WRITE(*,'(A,F8.2)')'ratio', REAL( (sf%s%nGP*(sf%s%deg+1))*sf%f%modes+(sf%f%mn_IP*sf%f%modes)*sf%s%nGP ,wp)/ & 
!                            REAL( (sf%f%mn_IP*sf%f%modes)*sf%s%nBase+ (sf%s%nGP*(sf%s%deg+1))*sf%f%mn_IP,wp)
!WRITE(*,*)'-----------------------------'

  call perfon('BaseEval')

  nGP     = sf%s%nGP  
  modes   = sf%f%modes  
  nBase   = sf%s%nBase
  mn_IP   = sf%f%mn_IP

#if (PP_WHICHEVAL==1)
  !OMP loop version: first s then f
  call perfon('eval_s')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) & 
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iMode)  &
!$OMP   SHARED(sf,deriv,modes,y_tmp,DOFs)
  DO iMode=1,modes
    y_tmp(iMode,:)=sf%s%evalDOF_GP(deriv(1),DOFs(:,iMode))
  END DO!iMode
!$OMP END PARALLEL DO
  call perfoff('eval_s')
  call perfon('eval_f')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) & 
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iGP)  &
!$OMP   SHARED(sf,deriv,nGP,y_IP_GP,y_tmp)
  DO iGP=1,nGP
    y_IP_GP(:,iGP)=sf%f%evalDOF_IP(deriv(2),y_tmp(:,iGP))
  END DO !iGP
!$OMP END PARALLEL DO
  call perfoff('eval_f')

#elif (PP_WHICHEVAL==2)
  !OMP loop version: first f then s
  call perfon('eval_f')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) & 
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(iBase)  &
!$OMP   SHARED(sf,nBase,deriv,y_tmp,DOFs)
  DO iBase=1,nBase
    y_tmp(:,iBase)=sf%f%evalDOF_IP(deriv(2),DOFs(iBase,:))
  END DO !iBase
!$OMP END PARALLEL DO
  call perfoff('eval_f')
  call perfon('eval_s')
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) & 
!$OMP   DEFAULT(NONE)    &
!$OMP   PRIVATE(i_mn)  &
!$OMP   SHARED(sf,mn_IP,deriv,y_IP_GP,y_tmp)
  DO i_mn=1,mn_IP
    y_IP_GP(i_mn,:)=sf%s%evalDOF_GP(deriv(1),y_tmp(i_mn,:))
  END DO!iMode
!$OMP END PARALLEL DO
  call perfoff('eval_s')

#elif (PP_WHICHEVAL==11)
  ! matrix-matrix version of first s then f

  call perfon('eval_s')
  ASSOCIATE(deg=>sf%s%deg, degGP=>sf%s%degGP, nElems=>sf%s%grid%nElems)
  SELECT CASE(deriv(1))
  CASE(0)
    DO iElem=1,nElems
      j=sf%s%base_offset(iElem)
      i=(iElem-1)*(degGP+1)+1  
      !y_tmp(:,i:i+degGP)=TRANSPOSE(MATMUL(sf%s%base_GP(0:degGP,0:deg,iElem),DOFs(j:j+deg,:)))
      __MATMAT_TT(y_tmp(:,i:i+degGP),DOFs(j:j+deg,:),sf%s%base_GP(0:degGP,0:deg,iElem))
    END DO !iElem
  CASE(DERIV_S)
    DO iElem=1,nElems
      j=sf%s%base_offset(iElem)
      i=(iElem-1)*(degGP+1)+1  
      !y_tmp(:,i:i+degGP)=TRANSPOSE(MATMUL(sf%s%base_ds_GP(0:degGP,0:deg,iElem),DOFs(j:j+deg,:)))
      __MATMAT_TT(y_tmp(:,i:i+degGP),DOFs(j:j+deg,:),sf%s%base_ds_GP(0:degGP,0:deg,iElem))
    END DO !iElem
  END SELECT !deriv GP
  END ASSOCIATE
  call perfoff('eval_s')

  call perfon('eval_f')
  SELECT CASE(deriv(2))
  CASE(0)
    !Y_IP_GP = MATMUL(sf%f%base_IP(:,:),y_tmp)
    __MATMAT_NN(Y_IP_GP ,sf%f%base_IP,y_tmp)
  CASE(DERIV_THET)                        
    !Y_IP_GP = MATMUL(sf%f%base_dthet_IP(:,:),y_tmp)
    __MATMAT_NN(Y_IP_GP ,sf%f%base_dthet_IP,y_tmp)
  CASE(DERIV_ZETA)                        
    !Y_IP_GP = MATMUL(sf%f%base_dzeta_IP(:,:),y_tmp)
    __MATMAT_NN(Y_IP_GP ,sf%f%base_dzeta_IP,y_tmp)
  END SELECT !deriv IP
  call perfoff('eval_f')

#elif (PP_WHICHEVAL==110)
  ! matrix-matrix version of first s then f

  call perfon('eval_s')
  ASSOCIATE(deg=>sf%s%deg, degGP=>sf%s%degGP, nElems=>sf%s%grid%nElems)
  SELECT CASE(deriv(1))
  CASE(0)
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,iElem,j,i) 
    DO iMode=1,modes
      DO iElem=1,nElems
        j=sf%s%base_offset(iElem)
        i=(iElem-1)*(degGP+1)+1  
        y_tmp(i:i+degGP,iMode)=MATMUL(sf%s%base_GP(0:degGP,0:deg,iElem),DOFs(j:j+deg,iMode))
      END DO !iElem
    END DO!iMode
!$OMP END PARALLEL DO
  CASE(DERIV_S)
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iMode,iElem,j,i) 
    DO iMode=1,modes
      DO iElem=1,nElems
        j=sf%s%base_offset(iElem)
        i=(iElem-1)*(degGP+1)+1  
        y_tmp(i:i+degGP,iMode)=MATMUL(sf%s%base_ds_GP(0:degGP,0:deg,iElem),DOFs(j:j+deg,iMode))
      END DO !iElem
    END DO!iMode
!$OMP END PARALLEL DO
  END SELECT !deriv GP
  END ASSOCIATE
  call perfoff('eval_s')

  call perfon('eval_f')
  SELECT CASE(deriv(2))
  CASE(0)
    !Y_IP_GP = MATMUL(sf%f%base_IP(:,:),y_tmp)
    __MATMAT_NT(Y_IP_GP ,sf%f%base_IP,y_tmp)
  CASE(DERIV_THET)                        
    !Y_IP_GP = MATMUL(sf%f%base_dthet_IP(:,:),y_tmp)
    __MATMAT_NT(Y_IP_GP ,sf%f%base_dthet_IP,y_tmp)
  CASE(DERIV_ZETA)                        
    !Y_IP_GP = MATMUL(sf%f%base_dzeta_IP(:,:),y_tmp)
    __MATMAT_NT(Y_IP_GP ,sf%f%base_dzeta_IP,y_tmp)
  END SELECT !deriv IP
  call perfoff('eval_f')

#elif (PP_WHICHEVAL==12)
  ! matrix-matrix version of first f then s
  call perfon('eval_f')
  SELECT CASE(deriv(2))
  CASE(0)
    !y_tmp = MATMUL(sf%f%base_IP(:,:),TRANSPOSE(DOFs))
    __MATMAT_NT(y_tmp ,sf%f%base_IP,DOFs)
  CASE(DERIV_THET)                        
    !y_tmp = MATMUL(sf%f%base_dthet_IP(:,:),TRANSPOSE(DOFs))
    __MATMAT_NT(y_tmp,sf%f%base_dthet_IP,DOFs)
  CASE(DERIV_ZETA)                        
    !y_tmp = MATMUL(sf%f%base_dzeta_IP(:,:),TRANSPOSE(DOFs))
    __MATMAT_NT(y_tmp,sf%f%base_dzeta_IP,DOFs)
  END SELECT !deriv IP
  call perfoff('eval_f')

  call perfon('eval_s')
  ASSOCIATE(deg=>sf%s%deg, degGP=>sf%s%degGP, nElems=>sf%s%grid%nElems)
  SELECT CASE(deriv(1))
  CASE(0)
    DO iElem=1,nElems
      j=sf%s%base_offset(iElem)
      i=(iElem-1)*(degGP+1)+1  
      !y_IP_GP(:,i:i+degGP)=MATMUL(y_tmp(:,j:j+deg),TRANSPOSE(sf%s%base_GP(0:degGP,0:deg,iElem)))
      __MATMAT_NT(y_IP_GP(:,i:i+degGP),y_tmp(:,j:j+deg),sf%s%base_GP(0:degGP,0:deg,iElem))
    END DO !iElem
  CASE(DERIV_S)
    DO iElem=1,nElems
      j=sf%s%base_offset(iElem)
      i=(iElem-1)*(degGP+1)+1  
      !y_IP_GP(:,i:i+degGP)=MATMUL(y_tmp(:,j:j+deg),TRANSPOSE(sf%s%base_ds_GP(0:degGP,0:deg,iElem)))
      __MATMAT_NT(y_IP_GP(:,i:i+degGP),y_tmp(:,j:j+deg),sf%s%base_ds_GP(0:degGP,0:deg,iElem))
    END DO !iElem
  END SELECT !deriv GP
  END ASSOCIATE
  call perfoff('eval_s')

#elif (PP_WHICHEVAL==120)
  ! matrix-matrix version of first f then s
  call perfon('eval_f')
  SELECT CASE(deriv(2))
  CASE(0)
    __MATMAT_NT(y_tmp,DOFs,sf%f%base_IP)
  CASE(DERIV_THET)                        
    __MATMAT_NT(y_tmp,DOFs,sf%f%base_dthet_IP)
  CASE(DERIV_ZETA)                        
    __MATMAT_NT(y_tmp,DOFs,sf%f%base_dzeta_IP)
  END SELECT !deriv IP
  call perfoff('eval_f')

  call perfon('eval_s')
  ASSOCIATE(deg=>sf%s%deg, degGP=>sf%s%degGP, nElems=>sf%s%grid%nElems)
  SELECT CASE(deriv(1))
  CASE(0)
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i_mn,iElem,j,i) 
    DO i_mn=1,mn_IP
      DO iElem=1,nElems
        j=sf%s%base_offset(iElem)
        i=(iElem-1)*(degGP+1)+1  
        y_IP_GP(i_mn,i:i+degGP)=MATMUL(sf%s%base_GP(0:degGP,0:deg,iElem),y_tmp(j:j+deg,i_mn))
      END DO !iElem
    END DO!i_mn
!$OMP END PARALLEL DO
  CASE(DERIV_S)
!$OMP PARALLEL DO        &  
!$OMP   SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i_mn,iElem,j,i) 
    DO i_mn=1,mn_IP
      DO iElem=1,nElems
        j=sf%s%base_offset(iElem)
        i=(iElem-1)*(degGP+1)+1  
        y_IP_GP(i_mn,i:i+degGP)=MATMUL(sf%s%base_ds_GP(0:degGP,0:deg,iElem),y_tmp(j:j+deg,i_mn))
      END DO !iElem
    END DO!i_mn
!$OMP END PARALLEL DO
  END SELECT !deriv GP
  END ASSOCIATE
  call perfoff('eval_s')
#endif /*PP_WHICHEVAL*/

#undef PP_WHICHEVAL

  call perfoff('BaseEval')

END SUBROUTINE base_evalDOF


!===================================================================================================================================
!> evaluate all degrees of freedom at given s theta zeta position (deriv=0 solution, deriv=1 first derivative d/ds)
!!
!===================================================================================================================================
FUNCTION base_evalDOF_x(sf,x,deriv,DOFs) RESULT(y_IP)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CLASS(t_base), INTENT(IN   ) :: sf     !! self
  REAL(wp)     , INTENT(IN   ) :: x(3)   !! s,theta,zeta
  INTEGER      , INTENT(IN   ) :: deriv(2)   !! =(0,0): base, =(DERIV_S,0): ds, =(0,DERIV_THET): dtheta, =(0,DERIV_ZETA): dzeta
  REAL(wp)     , INTENT(IN   ) :: DOFs(1:sf%s%nBase,1:sf%f%modes)  !! array of all modes and all radial dofs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL(wp)                     :: y_IP
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                      :: iMode
  REAL(wp)                     :: y_tmp(1:sf%f%modes)
!===================================================================================================================================
  DO iMode=1,sf%f%modes
    y_tmp(iMode)=sf%s%evalDOF_s(x(1),deriv(1),DOFs(:,iMode))
  END DO!iMode
  y_IP=sf%f%evalDOF_x(x(2:3),deriv(2),y_tmp(:))

END FUNCTION base_evalDOF_x

!===================================================================================================================================
!> test base variable
!!
!===================================================================================================================================
SUBROUTINE Base_test( sf )
! MODULES
USE MODgvec_GLobals, ONLY: UNIT_StdOut,testdbg,testlevel,nfailedMsg,nTestCalled,testUnit
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CLASS(t_Base), INTENT(INOUT) :: sf !! self
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER            :: iTest,iMode,iGP
  CHARACTER(LEN=10)  :: fail
  REAL(wp),PARAMETER :: realtol=1.0E-11_wp
  REAL(wp)           :: checkreal,tmp 
  REAL(wp)           :: dofs(1:sf%s%nBase,1:sf%f%modes)
  REAL(wp)           :: g_sIP(1:sf%s%nBase)
  REAL(wp)           :: g_IP_GP(1:sf%f%mn_IP,1:sf%s%nGP)
  REAL(wp)           :: g_IP_GP_eval(1:sf%f%mn_IP,1:sf%s%nGP)
!===================================================================================================================================
  test_called=.TRUE.
  IF(testlevel.LE.0) RETURN
  IF(testdbg) THEN
     Fail=" DEBUG  !!"
  ELSE
     Fail=" FAILED !!"
  END IF
  nTestCalled=nTestCalled+1
  SWRITE(UNIT_stdOut,'(A,I4,A)')'>>>>>>>>> RUN BASE TEST ID',nTestCalled,'  >>>>>>>>>'
  ASSOCIATE(modes=>sf%f%modes,sin_range=>sf%f%sin_range,cos_range=>sf%f%cos_range, &
            deg=>sf%s%deg,nBase=>sf%s%nBase,sin_cos=>sf%f%sin_cos,nGP=>sf%s%nGP,Xmn=>sf%f%Xmn)
  IF(testlevel.GE.1)THEN

    iTest=101 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    g_IP_GP(:,:)=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +SIN(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(0.2_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +COS(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*(0.2_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode

    CALL sf%evalDOF((/0,0/),dofs,g_IP_GP_eval)
    g_IP_GP=(g_IP_GP - g_IP_GP_eval)

    checkreal=MAXVAL(ABS(g_IP_GP))
    IF(testdbg.OR.(.NOT.(checkreal .LT. realtol) )) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! BASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(4(A,E11.3))') &
      '\n =>  should be 0.0 : MAX(|g_IP_exact-g_IP(dofs)|) = ', checkreal, &
      '\n     maxval(|g_IP|)= ',MAXVAL(ABS(g_IP_GP_eval)),', minval(|g_IP|)= ',MINVAL(ABS(g_IP_GP_eval)), &
      ', avg(|g_IP|)= ',SUM(ABS(g_IP_GP_eval))/REAL(modes*nGP,wp)
    END IF !TEST

    iTest=102 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    g_IP_GP(:,:)=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(1.0_wp+0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+0.3_wp*sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +SIN(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*REAL(deg,wp)*0.3_wp*(1.0_wp+0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+0.3_wp*sf%s%s_GP(iGP))**(deg-1)
      END DO !iGP
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(1.0_wp+0.2_wp*REAL(iMode,wp)/REAL(modes,wp)+0.4_wp*sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +COS(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*REAL(deg,wp)*0.4_wp*(1.0_wp+0.2_wp*REAL(iMode,wp)/REAL(modes,wp)+0.4_wp*sf%s%s_GP(iGP))**(deg-1)
      END DO !iGP
    END DO !iMode

    CALL sf%evalDOF((/DERIV_S,0/),dofs,g_IP_GP_eval)
    g_IP_GP=(g_IP_GP - g_IP_GP_eval)

    checkreal=MAXVAL(ABS(g_IP_GP))
    IF(testdbg.OR.(.NOT.(checkreal .LT. realtol) )) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! BASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(4(A,E11.3))') &
      '\n =>  should be 0.0 : MAX(|g_IP_exact-g_IP(dofs)|) = ', checkreal, &
      '\n     maxval(|g_IP|)= ',MAXVAL(ABS(g_IP_GP_eval)),', minval(|g_IP|)= ',MINVAL(ABS(g_IP_GP_eval)), &
      ', avg(|g_IP|)= ',SUM(ABS(g_IP_GP_eval))/REAL(modes*nGP,wp)
    END IF !TEST

    iTest=103 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    g_IP_GP(:,:)=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +REAL(Xmn(1,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(0.2_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       -REAL(Xmn(1,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*(0.2_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode
    CALL sf%evalDOF((/0,DERIV_THET/),dofs,g_IP_GP_eval)
    g_IP_GP=(g_IP_GP - g_IP_GP_eval)

    checkreal=MAXVAL(ABS(g_IP_GP))
    IF(testdbg.OR.(.NOT.(checkreal .LT. realtol) )) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! BASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(4(A,E11.3))') &
      '\n =>  should be 0.0 : MAX(|g_IP_exact-g_IP(dofs)|) = ', checkreal, &
      '\n     maxval(|g_IP|)= ',MAXVAL(ABS(g_IP_GP_eval)),', minval(|g_IP|)= ',MINVAL(ABS(g_IP_GP_eval)), &
      ', avg(|g_IP|)= ',SUM(ABS(g_IP_GP_eval))/REAL(modes*nGP,wp)
    END IF !TEST

    iTest=104 ; IF(testdbg)WRITE(*,*)'iTest=',iTest
    g_IP_GP(:,:)=0.0_wp
    DO iMode=sin_range(1)+1,sin_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       -REAL(Xmn(2,iMode),wp)*COS(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode
    DO iMode=cos_range(1)+1,cos_range(2)
      tmp  = 1.0_wp/(REAL(1+Xmn(1,iMode)**2+Xmn(2,imode)**2))
      g_sIP(:)=tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_IP)**deg
      dofs(:,iMode)=sf%s%initDOF(g_sIP)
      DO iGP=1,nGP
        g_IP_GP(:,iGP)=g_IP_GP(:,iGP) &
                       +REAL(Xmn(2,iMode),wp)*SIN(REAL(Xmn(1,iMode),wp)*sf%f%x_IP(1,:)-REAL(Xmn(2,iMode),wp)*sf%f%x_IP(2,:))* &
                        tmp*(0.1_wp*REAL(iMode,wp)/REAL(modes,wp)+sf%s%s_GP(iGP))**deg
      END DO !iGP
    END DO !iMode

    CALL sf%evalDOF((/0,DERIV_ZETA/),dofs,g_IP_GP_eval)
    g_IP_GP=(g_IP_GP - g_IP_GP_eval)

    checkreal=MAXVAL(ABS(g_IP_GP))
    IF(testdbg.OR.(.NOT.(checkreal .LT. realtol) )) THEN
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(A,2(I4,A))') &
      '\n!! BASE TEST ID',nTestCalled ,': TEST ',iTest,Fail
      nfailedMsg=nfailedMsg+1 ; WRITE(testUnit,'(4(A,E11.3))') &
      '\n =>  should be 0.0 : MAX(|g_IP_exact-g_IP(dofs)|) = ', checkreal, &
      '\n     maxval(|g_IP|)= ',MAXVAL(ABS(g_IP_GP_eval)),', minval(|g_IP|)= ',MINVAL(ABS(g_IP_GP_eval)), &
      ', avg(|g_IP|)= ',SUM(ABS(g_IP_GP_eval))/REAL(modes*nGP,wp)
    END IF !TEST


  END IF !testlevel>=1
  END ASSOCIATE !sf

  test_called=.FALSE.

END SUBROUTINE Base_test



END MODULE MODgvec_base

