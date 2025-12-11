!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"

MODULE MODgvec_BiotSavart

IMPLICIT NONE
PUBLIC

CONTAINS

!===================================================================================================================================
!> Evaluate the magnetic field of a current loop discretized via line-segments.
!! The segments are defined by a list of `coil_points`, each segment being defined between points (i,i+1).
!! For a closed loop, the last point in the list is a repetition of the first point.
!! Each line segment is evaluated via the analytic compact Biot-Savart expression from
!! Hanson and Hirshman (2002) (https://doi.org/10.1063/1.1507589).
!! This implementation parallelizes over the number of evaluation positions `n_positions`.
!===================================================================================================================================
SUBROUTINE BiotSavart(n_positions, xyz, n_points, coil_points, prefactor, B)
    ! MODULES
    USE MODgvec_Globals, ONLY: TWOPI,CROSS,wp

    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    INTEGER,  INTENT(IN) :: n_positions               !! number of positions where the magnetic field is evaluated
    REAL(wp), INTENT(IN) :: xyz(3,n_positions)        !! x,y,z positions where the magnetic field is evaluated
    INTEGER,  INTENT(IN) :: n_points                  !! number of points representing the segmented coil
    REAL(wp), INTENT(IN) :: coil_points(3,n_points)   !! x,y,z positions of the segment chain (=polygon), n_segments=n_points-1
    REAL(wp), INTENT(IN) :: prefactor                 !! factor applied to the final magnetic field
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLE
    REAL(wp), INTENT(OUT):: B(3,n_positions)          !! magnetic field evaluated at xyz positions.
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER              :: iPosition,iSegment,n_segments
    REAL(wp)             :: R_i(3)                    !! vector from current position to starting point of segment
    REAL(wp)             :: R_f(3)                    !! vector from current position to end point of segment (=starting point of next segment)
    REAL(wp)             :: mod_R_i,mod_R_f
    !===================================================================================================================================

    n_segments = n_points-1

    B = 0.0_wp
    !$OMP PARALLEL DO     &
    !$OMP SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(xyz, coil_points, n_positions, n_segments) &
    !$OMP REDUCTION(+:B)
    DO iPosition=1,n_positions
        R_f = xyz(:,iPosition)-coil_points(:,1)
        mod_R_f= SQRT(SUM(R_f*R_f))
        DO iSegment=1,n_segments
            R_i = R_f
            mod_R_i = mod_R_f
            R_f = xyz(:,iPosition)-coil_points(:,iSegment+1)
            mod_R_f = SQRT(SUM(R_f*R_f))
            !! RidotRf = SUM(R_i*R_f) ;  RixRf(:) = CROSS(R_i, R_f)
            !! segment_factor = (mod_R_i+mod_R_f)/(mod_R_i*mod_R_f*(mod_R_i*mod_R_f+RidotRf))
            !! B(:,iPosition) = B(:,iPosition)+segment_factor*RixRf
            B(:,iPosition) = B(:,iPosition)+((mod_R_i+mod_R_f)/(mod_R_i*mod_R_f*(mod_R_i*mod_R_f+SUM(R_i*R_f))))*CROSS(R_i, R_f)
        END DO !iSegment
    END DO !iPosition
    !$OMP END PARALLEL DO
    B = prefactor*B
END SUBROUTINE BiotSavart

SUBROUTINE BiotSavart_VectorPotential(n_positions, xyz, n_points, n_segments, coil_points, ehat, L, prefactor, A)
    ! MODULES
    USE MODgvec_Globals, ONLY: TWOPI,CROSS,wp

    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    INTEGER,  INTENT(IN) :: n_positions               !! number of positions where the magnetic field is evaluated
    REAL(wp), INTENT(IN) :: xyz(3,n_positions)        !! x,y,z positions where the magnetic field is evaluated
    INTEGER,  INTENT(IN) :: n_points                  !! number of points representing the segmented coil
    INTEGER,  INTENT(IN) :: n_segments                !! number of coil segments
    REAL(wp), INTENT(IN) :: coil_points(3,n_points)   !! x,y,z positions of the segment chain (=polygon), n_segments=n_points-1
    REAL(wp), INTENT(IN) :: prefactor                 !! factor applied to the final magnetic field
    REAL(wp), INTENT(IN) :: ehat(3,n_segments)         !! unit vector along segment
    REAL(wp), INTENT(IN) :: L(n_segments)             !! segment length
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLE
    REAL(wp), INTENT(OUT):: A(3,n_positions)          !! magnetic field evaluated at xyz positions.
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER              :: iPosition,iSegment
    REAL(wp)             :: R_i(3)                    !! vector from current position to starting point of segment
    REAL(wp)             :: R_f(3)                    !! vector from current position to end point of segment (=starting point of next segment)
    REAL(wp)             :: norm_ehat
    REAL(wp)             :: mod_R_i,mod_R_f
    REAL(wp)             :: f_epsilon
    !===================================================================================================================================
    A = 0.0_wp
    !$OMP PARALLEL DO     &
    !$OMP SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(xyz, coil_points, n_positions, n_segments, L, ehat) &
    !$OMP REDUCTION(+:A)
    DO iPosition=1,n_positions
        R_f = xyz(:,iPosition)-coil_points(:,1)
        mod_R_f= SQRT(SUM(R_f*R_f))
        DO iSegment=1,n_segments
            R_i = R_f
            mod_R_i = mod_R_f
            R_f = xyz(:,iPosition)-coil_points(:,iSegment+1)
            mod_R_f = SQRT(SUM(R_f*R_f))
            f_epsilon = LOG( (mod_R_i + mod_R_f + L(iSegment)) / (mod_R_i + mod_R_f - L(iSegment))  )
            A(:,iPosition) = A(:,iPosition) + ehat(:,iSegment) * f_epsilon
        END DO !iSegment
    END DO !iPosition
    !$OMP END PARALLEL DO
    A = prefactor*A
END SUBROUTINE BiotSavart_VectorPotential

END MODULE MODgvec_BiotSavart
