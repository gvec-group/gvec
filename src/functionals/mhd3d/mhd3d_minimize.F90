!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.FPP"

!===================================================================================================================================
!>
!!# Module **MHD3D minimize**
!!
!! CONTAINS INITIALIZATION minimizer for MHD3D functional
!!
!===================================================================================================================================
MODULE MODgvec_MHD3D_minimize
    ! MODULES
    USE MODgvec_Globals, ONLY:wp,abort,UNIT_stdOut,fmt_sep,MPIRoot,enter_subregion,exit_subregion
    USE MODgvec_Sol_Var_MHD3D,ONLY: t_sol_var_MHD3D
    IMPLICIT NONE

    ! PRIVATE
    ! PUBLIC t_minimizer_mhd3d

    !-----------------------------------------------------------------------------------------------------------------------------------
    ! TYPES
    !-----------------------------------------------------------------------------------------------------------------------------------
    TYPE :: t_minimizer_mhd3d
        !-------------------------------------------------------------------------------------------------------------------------------
        LOGICAL   :: initialized, restart_iter, logger_is_initialized
        INTEGER   :: JacCheck, MinType
        INTEGER   :: iter,nStepDecreased,nSkip_Jac,nSkip_dw
        INTEGER   :: lastoutputIter, logiter_ramp, logscreen
        REAL(wp)  :: dt,deltaW, dW_allowed
        REAL(wp)  :: t_pseudo,Fnorm(3),Vnorm(3),Fnorm0(3),Fnorm_old(3),W_MHD3D_0

        ! Logging variables
        REAL(wp) :: min_dt_out,max_dt_out,min_dw_out,max_dw_out,sum_dW_out
        LOGICAL  :: DoCheckDistance, DoCheckAxis
        INTEGER  :: outputIter, nlogScreen, logIter, logUnit, StartTimeArray(8)

        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: dofs(:)      !! degrees of freedom at levels (k-1),(k),(k+1)
        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: force(:)     !! force
        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: temp_dofs(:) !! temporary for update

        ! acc. gradient descent variables
        REAL(wp), ALLOCATABLE  :: tau(:)
        INTEGER   :: ndamp
        REAL(wp)  :: tau_bar
        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: velocity(:)       !! 'velocity' in minimizer
        !-------------------------------------------------------------------------------------------------------------------------------
        CONTAINS
        PROCEDURE :: minimize         => MinimizeMHD3D_descent
        PROCEDURE :: reset            => MinimizeMHD3d_ResetDescent
        PROCEDURE :: StartLogging     => StartLogging_MHD3D
        PROCEDURE :: Logging          => Logging_MHD3D
        PROCEDURE :: free             => Free_minimizer
    END TYPE t_minimizer_mhd3d


    CONTAINS

    SUBROUTINE new_minimizer(&
            sf, varsize_in, dt_initial, MinType_in, dW_allowed_in,& ! Minimizer
            DoCheckDistance, DoCheckAxis,& !what to log
            outputIter, nlogScreen, logIter& !when to log
        )
        ! MODULES
        USE MODgvec_MHD3D_evalFunc, ONLY: EvalForce
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        INTEGER , INTENT(IN) :: varsize_in(:)
        REAL(wp), INTENT(IN) :: dt_initial, dW_allowed_in
        INTEGER, INTENT(IN)  :: MinType_in
        LOGICAL, INTENT(IN)  :: DoCheckDistance, DoCheckAxis
        INTEGER, INTENT(IN)  :: outputIter, nlogScreen, logIter
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        TYPE(t_minimizer_mhd3d),ALLOCATABLE, INTENT(INOUT) :: sf
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: i
        !===================================================================================================================================
        ALLOCATE(sf)

        sf%MinType = MinType_in
        sf%dt = dt_initial
        sf%dW_allowed = dW_allowed_in

        ALLOCATE(sf%dofs(-3:1))
        CALL sf%dofs(1)%init(varsize_in)
        DO i=-3,0
            CALL sf%dofs(i)%copy(sf%dofs(1))
        END DO
        ALLOCATE(sf%force(-1:0))
        DO i=-1,0
            CALL sf%force(i)%copy(sf%dofs(1))
        END DO
        ALLOCATE(sf%velocity(-1:1))
        DO i=-1,1
            CALL sf%velocity(i)%copy(sf%dofs(1))
        END DO
        ALLOCATE(sf%temp_dofs(-1:1))
        DO i=-1,1
            CALL sf%temp_dofs(i)%copy(sf%dofs(1))
        END DO

        sf%JacCheck = 1
        sf%ndamp = 10
        sf%nstepDecreased = 0
        sf%nSkip_Jac = 0
        sf%t_pseudo = 0
        sf%lastOutputIter = 0
        sf%iter = 0
        sf%Vnorm = 0.0_wp
        sf%logiter_ramp = 1
        sf%logscreen = 1
        sf%nSkip_dw = 0
        sf%JacCheck = 1

        ALLOCATE(sf%tau(sf%ndamp))

        sf%initialized = .TRUE.
        sf%restart_iter = .TRUE.
        sf%logger_is_initialized = .FALSE.

        sf%DoCheckDistance = DoCheckDistance
        sf%DoCheckAxis = DoCheckAxis

        sf%outputIter = outputIter
        sf%nlogScreen = nlogScreen
        sf%logIter = logIter

        sf%Fnorm = 10.0_wp
        sf%Fnorm0 = 10.0_wp
        sf%W_MHD3D_0 = 0.0_wp
        sf%deltaW = 0.0_wp
    END SUBROUTINE new_minimizer

    SUBROUTINE Free_minimizer(sf)
        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        INTEGER :: i

        CLOSE(sf%logUnit)
        DO i=-1,1
            IF(ALLOCATED(sf%dofs)) CALL sf%dofs(i)%free()
            IF(ALLOCATED(sf%temp_dofs)) CALL sf%temp_dofs(i)%free()
            IF(ALLOCATED(sf%velocity)) CALL sf%velocity(i)%free()
        END DO
        DO i=-1,0
            IF(ALLOCATED(sf%force)) CALL sf%force(i)%free()
        END DO

        SDEALLOCATE(sf%dofs)
        SDEALLOCATE(sf%temp_dofs)
        SDEALLOCATE(sf%force)
        SDEALLOCATE(sf%velocity)
        SDEALLOCATE(sf%tau)
    END SUBROUTINE Free_minimizer

    SUBROUTINE MinimizeMHD3d_ResetDescent(sf)
        USE MODgvec_MHD3D_EvalFunc, ONLY: EvalAux, EvalForce, EvalEnergy
        !-------------------------------------------------------------------------------------------------------------------------------
        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        !-------------------------------------------------------------------------------------------------------------------------------
        IF(.NOT.sf%logger_is_initialized)THEN
            CALL sf%StartLogging()
        END IF
        sf%JacCheck=1 !abort if detJ<0
        CALL EvalAux(sf%dofs(0), sf%JacCheck)
        sf%dofs(0)%W_MHD3D= EvalEnergy(sf%dofs(0),.FALSE.,sf%JacCheck)
        sf%W_MHD3D_0 = sf%dofs(0)%W_MHD3D
        CALL EvalForce( sf%dofs(0),.FALSE.,sf%JacCheck,sf%force(0))
        sf%Fnorm0 = SQRT(sf%force(0)%norm_2())
        sf%Fnorm = sf%Fnorm0
        sf%Fnorm_old = 1.1_wp*sf%Fnorm0
        CALL sf%dofs(-1)%set_to(sf%dofs(0)) !last state
        CALL sf%dofs(-2)%set_to(sf%dofs(0)) !state at last logging interval

        !for hirshman method
        IF(sf%MinType.EQ.10)THEN
            CALL sf%velocity(-1)%set_to(0.0_wp)
            CALL sf%velocity( 0)%set_to(0.0_wp)
            sf%tau(1:sf%ndamp)=0.15_wp/sf%dt
            sf%tau_bar = 0.075_wp
        END IF

        sf%min_dt_out=1.0e+30_wp
        sf%max_dt_out=0.0_wp
        sf%min_dW_out=1.0e+30_wp
        sf%max_dW_out=-1.0e+30_wp
        sf%sum_dW_out=0.0_wp
        sf%nSkip_dW =0
        IF(sf%restart_iter) sf%restart_iter=.FALSE.
    END SUBROUTINE MinimizeMHD3d_ResetDescent

    SUBROUTINE StartLogging_MHD3D(sf)
        USE MODgvec_Globals,     ONLY: GETFREEUNIT
        USE MODgvec_Output_Vars, ONLY: ProjectName,outputLevel
        USE MODgvec_MHD3D_visu,  ONLY: checkPos
        IMPLICIT NONE
        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        !---------------------------------------------------------------------------------------------------------------------------------
        CHARACTER(LEN=255)  :: fileString
        INTEGER             :: TimeArray(8),iLogDat
        REAL(wp)            :: Pos(2,2),W_MHD3D
        INTEGER,PARAMETER   :: nLogDat=25
        REAL(wp)            :: LogDat(1:nLogDat)
        !=================================================================================================================================
        IF(.NOT.MPIroot) RETURN
        __PERFON('log_output')
        W_MHD3D=sf%dofs(0)%W_MHD3D
        CALL DATE_AND_TIME(values=TimeArray) ! get System time
        WRITE(UNIT_stdOut,'(A,E11.4,A)')'%%%%%%%%%%  START ITERATION, dt= ',sf%dt, '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        WRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                        '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
        WRITE(UNIT_stdOut,'(A,3E21.14)') &
                '%%% dU = |Force|= ',sf%Fnorm(1:3)
        IF(sf%MinType.EQ.10) THEN
            WRITE(UNIT_stdOut,'(A,E11.4,A,3E11.4)') &
                '%%% accel.GD: tau= ',sf%tau_bar,' |vel|= ',sf%Vnorm(1:3)
        END IF

        WRITE(UNIT_stdOut,'(40(" -"))')
        !------------------------------------
        sf%StartTimeArray=TimeArray !save first time stamp

        sf%logUnit=GETFREEUNIT()
        WRITE(FileString,'("logMinimizer_",A,"_",I4.4,".csv")')TRIM(ProjectName),outputLevel
        OPEN(UNIT     = sf%logUnit       ,&
            FILE     = TRIM(FileString) ,&
            STATUS   = 'REPLACE'   ,&
            ACCESS   = 'SEQUENTIAL' )
        !header
        iLogDat=0
        WRITE(sf%logUnit,'(A)',ADVANCE="NO")'"#iterations","runtime(s)","min_dt","max_dt"'
        WRITE(sf%logUnit,'(A)',ADVANCE="NO")',"W_MHD3D","min_dW","max_dW","sum_dW"'
        WRITE(sf%logUnit,'(A)',ADVANCE="NO")',"normF_X1","normF_X2","normF_LA"'
        LogDat(ilogDat+1:iLogDat+11)=(/0.0_wp,0.0_wp,sf%dt,sf%dt,W_MHD3D,0.0_wp,0.0_wp,0.0_wp,sf%Fnorm(1:3)/)
        iLogDat=11
        IF(sf%MinType.EQ.10) THEN
            WRITE(sf%logUnit,'(A)',ADVANCE="NO")',"tau","normV_X1","normV_X2","normV_LA"'
            LogDat(ilogDat+1:iLogDat+4)=(/sf%tau_bar,sf%Vnorm(1:3)/)
            iLogDat=iLogDat+4
        END IF
        IF(sf%doCheckDistance) THEN
            WRITE(sf%logUnit,'(A)',ADVANCE="NO")',"max_Dist","avg_Dist"'
            LogDat(iLogDat+1:iLogDat+2)=(/0.0_wp,0.0_wp/)
            iLogDat=iLogDat+2
        END IF!doCheckDistance
        IF(sf%doCheckAxis) THEN
            WRITE(sf%logUnit,'(A)',ADVANCE="NO")',"X1_axis_0","X2_axis_0","X1_axis_1","X2_axis_1"'
            CALL CheckPos(sf%dofs(0),0.0_wp,2,Pos)
            LogDat(iLogDat+1:iLogDat+4)=RESHAPE(Pos,(/4/))
            iLogDat=iLogDat+4
        END IF!doCheckAxis

        WRITE(sf%logUnit,'(A)')' '
        !first data line
        WRITE(sf%logUnit,'(*(e23.15,:,","))') logDat(1:iLogDat)
        sf%logger_is_initialized = .TRUE.
        __PERFOFF('log_output')
    END SUBROUTINE StartLogging_MHD3D

    SUBROUTINE Logging_MHD3D(sf, quiet)
        USE MODgvec_MHD3D_visu, ONLY: checkDistance
        USE MODgvec_MHD3D_visu, ONLY: checkPos
        IMPLICIT NONE

        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        LOGICAL, INTENT(IN) :: quiet !! True: no screen output
        !---------------------------------------------------------------------------------------------------------------------------------
        INTEGER             :: TimeArray(8),runtime_ms,iLogDat
        REAL(wp)            :: Pos(2,2),maxDist,avgDist,W_MHD3D
        INTEGER,PARAMETER   :: nLogDat=25
        REAL(wp)            :: LogDat(1:nLogDat)
        !=================================================================================================================================
        IF(.NOT.MPIroot) RETURN
        __PERFON('log_output')
        CALL DATE_AND_TIME(values=TimeArray) ! get System time
        W_MHD3D=sf%dofs(0)%W_MHD3D
        IF(.NOT.quiet)THEN
            WRITE(UNIT_stdOut,'(80("%"))')
            WRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                            '%%% Sys date : ',timeArray(1:3),timeArray(5:7)
            WRITE(UNIT_stdOut,'(A,I8,A,2I8,A,E11.4,A,2E11.4,A,E21.14,A,3E12.4)') &
                            '%%% #ITERATIONS= ',sf%iter,', #skippedIter (Jac/dW)= ',sf%nSkip_Jac,sf%nSkip_dW, &
                    '\n%%% t_pseudo= ',sf%t_pseudo,', min/max dt= ',sf%min_dt_out,sf%max_dt_out, &
                    '\n%%% W_MHD3D= ',W_MHD3D,', min/max/sum deltaW= ' , sf%min_dW_out,sf%max_dW_out,sf%sum_dW_out
            WRITE(UNIT_stdOut,'(A,3E21.14)') &
                        '%%% dU = |Force|= ',sf%Fnorm(1:3)
            !------------------------------------
        END IF!.NOT.quiet

        iLogDat=0
        runtime_ms=MAX(0,SUM((timeArray(5:8)-sf%StartTimearray(5:8))*(/360000,6000,100,1/)))
        LogDat(ilogDat+1:iLogDat+11)=(/REAL(sf%iter,wp),REAL(runtime_ms,wp)/100.0_wp, &
                                        sf%min_dt_out,sf%max_dt_out,&
                                        W_MHD3D,sf%min_dW_out,sf%max_dW_out,sf%sum_dW_out, &
                                        sf%Fnorm(1:3)/)
        iLogDat=11
        IF(sf%MinType.EQ.10) THEN
            IF(.NOT.quiet)THEN
            WRITE(UNIT_stdOut,'(A,E11.4,A,3E11.4)') &
                    '%%% accel.GD: tau= ',sf%tau_bar,' |vel|= ',sf%Vnorm(1:3)
            END IF!.NOT.quiet
            LogDat(ilogDat+1:iLogDat+4)=(/sf%tau_bar,sf%Vnorm(1:3)/)
            iLogDat=iLogDat+4
        END IF
        IF(sf%doCheckDistance) THEN
            CALL CheckDistance(sf%dofs(0),sf%dofs(-2),maxDist,avgDist)
            CALL sf%dofs(-2)%set_to(sf%dofs(0))
            IF(.NOT.quiet)THEN
            WRITE(UNIT_stdOut,'(A,2E11.4)') &
            '               %%% Dist to last log (max/avg) : ',maxDist,avgDist
            END IF!.NOT.quiet
            LogDat(iLogDat+1:iLogDat+2)=(/maxDist,avgDist/)
            iLogDat=iLogDat+2
        END IF!doCheckDistance
        IF(sf%doCheckAxis) THEN
            CALL CheckPos(sf%dofs(0),0.0_wp,2,Pos)
            IF(.NOT.quiet)THEN
            SWRITE(UNIT_stdOut,'(2(A,2E22.14))') &
                '%%% axis position (X1,X2,zeta=0     ): ',Pos(1:2,1), &
            '\n%%% axis position (X1,X2,zeta=pi/nfp): ',Pos(1:2,2)
            END IF!.NOT.quiet
            LogDat(iLogDat+1:iLogDat+4)=RESHAPE(Pos,(/4/))
            iLogDat=iLogDat+4
        END IF !doCheckAxis

        IF(.NOT.quiet)THEN
            WRITE(UNIT_stdOut,'(40(" -"))')
        END IF!.NOT.quiet
        WRITE(sf%logUnit,'(*(e23.15,:,","))') logDat(1:iLogDat)
        __PERFOFF('log_output')
    END SUBROUTINE Logging_MHD3D

    SUBROUTINE MinimizeMHD3D_descent(sf, abstol, maxIter_in)
        USE MODgvec_MHD3D_EvalFunc, ONLY: EvalAux, EvalForce, EvalEnergy
        USE MODgvec_Analyze, ONLY:analyze
        USE MODgvec_Restart, ONLY:WriteState
        USE MODgvec_MHD3D_visu, ONLY:WriteSFLoutfile
        USE MODgvec_sol_var_MHD3D, ONLY: t_sol_var_MHD3D
        IMPLICIT NONE

        REAL(wp), INTENT(IN) :: abstol
        INTEGER, INTENT(IN)  :: maxIter_in
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: maxIter, i
        !=================================================================================================================================
        maxIter = sf%iter + maxIter_in
        DO WHILE(sf%iter.LT.maxIter)
            IF((sf%restart_iter))THEN
                CALL sf%reset()
            END IF !before first iteration or after restart Jac<0

            !COMPUTE NEW SOLUTION P(1) as a prediction
            SELECT CASE(sf%MinType)
                CASE(0) !gradient descent, previously used for minimizerType=0
                    CALL sf%temp_dofs(1)%AXBY(1.0_wp,sf%dofs(0),sf%dt,sf%force(0)) !overwrites P(1), predicts solution U(1)
                CASE(10) !hirshman method
                    !tau is damping parameter
                    sf%tau(1:sf%ndamp-1) = sf%tau(2:sf%ndamp) !save old

                    !ln(|F_n|^2/|F_{n-1}|^2), Fnorm=|F_X1|,|F_X2|,|F_LA|
                    sf%tau(sf%ndamp)  = MIN(0.15_wp,ABS(LOG(SUM(sf%Fnorm**2)/SUM(sf%Fnorm_old**2))))/sf%dt

                    sf%tau_bar = 0.5_wp*sf%dt*SUM(sf%tau)/REAL(sf%ndamp,wp)   !=1/2 * tauavg
                    CALL sf%velocity(1)%AXBY(((1.0_wp-sf%tau_bar)/(1.0_wp+sf%tau_bar)),sf%velocity(0),(sf%dt/(1.0_wp+sf%tau_bar)),sf%force(0)) !velocity V(1)
                    CALL sf%temp_dofs(1)%AXBY(1.0_wp,sf%dofs(0),sf%dt,sf%velocity(1)) !overwrites P(1), predicst solution U(1)
                    sf%Vnorm=SQRT(sf%velocity(1)%norm_2())
            END SELECT

            sf%JacCheck=2 !no abort,if detJ<0, JacCheck=-1

            sf%temp_dofs(1)%W_MHD3D=EvalEnergy(sf%temp_dofs(1),.TRUE.,sf%JacCheck)
            IF(sf%JacCheck.EQ.-1)THEN
                sf%dt=0.9_wp*sf%dt
                sf%nstepDecreased=sf%nStepDecreased+1
                sf%nSkip_Jac=sf%nSkip_Jac+1
                sf%restart_iter=.TRUE.
                CALL sf%dofs(0)%set_to(sf%dofs(-3)) !reset to initial state
                SWRITE(UNIT_stdOut,'(8X,I8,A,E11.4,A)')sf%iter,'...detJac<0, decrease stepsize to dt=',sf%dt,  ' and RESTART simulation!!!!!!!'
            ELSE
                !detJ>0
                sf%deltaW=sf%temp_dofs(1)%W_MHD3D-sf%dofs(0)%W_MHD3D!should be <=0,
                IF(sf%deltaW.LE.sf%dW_allowed*sf%W_MHD3D_0)THEN !valid step /hirshman method accept W increase!
                    IF(ALL(sf%Fnorm.LE.abstol))THEN
                        CALL sf%Logging(.FALSE.)
                        SWRITE(UNIT_stdOut,'(4x,A)')'==>Iteration finished, |force| in relative tolerance'
                        EXIT !DO LOOP
                    END IF
                    sf%iter=sf%iter+1
                    sf%t_pseudo=sf%t_pseudo+sf%dt
                    ! for simple gradient & hirshman
                    CALL sf%dofs(-1)%set_to(sf%dofs(0))
                    CALL sf%dofs(0)%set_to(sf%temp_dofs(1))
                    ! for hirshman method
                    IF(sf%MinType.EQ.10)THEN
                        CALL sf%velocity(-1)%set_to(sf%velocity(0))
                        CALL sf%velocity(0)%set_to(sf%velocity(1))
                    END IF
                    CALL EvalForce(sf%temp_dofs(1),.FALSE.,sf%JacCheck,sf%force(0)) !evalAux was already called on P(1)=U(0), so that its set false here.
                    sf%Fnorm_old=sf%Fnorm
                    sf%Fnorm=SQRT(sf%force(0)%norm_2())
                    sf%nstepDecreased=0
                    sf%min_dt_out=MIN(sf%min_dt_out,sf%dt)
                    sf%max_dt_out=MAX(sf%max_dt_out,sf%dt)
                    sf%min_dW_out=MIN(sf%min_dW_out,sf%deltaW)
                    sf%max_dW_out=MAX(sf%max_dW_out,sf%deltaW)
                    sf%sum_dW_out=sf%sum_dW_out+sf%deltaW
                    IF(MOD(sf%iter,sf%logIter_ramp).EQ.0)THEN
                        CALL sf%Logging(.NOT.((sf%logIter_ramp.GE.sf%logIter).AND.(MOD(sf%logscreen,sf%nLogScreen).EQ.0)))
                        IF(.NOT.(sf%logIter_ramp.LT.sf%logIter))THEN !only reset for logIter
                            sf%logscreen=sf%logscreen+1
                            sf%min_dt_out=1.0e+30_wp
                            sf%max_dt_out=0.0_wp
                            sf%min_dW_out=1.0e+30_wp
                            sf%max_dW_out=-1.0e+30_wp
                            sf%sum_dW_out=0.0_wp
                            sf%nSkip_dW =0
                        END IF
                        sf%logIter_ramp=MIN(sf%logIter,sf%logIter_ramp*2)
                    END IF
                ELSE !not a valid step, decrease timestep and skip P(1)
                    sf%dt=0.9_wp*sf%dt
                    sf%nstepDecreased=sf%nStepDecreased+1
                    sf%nSkip_dW=sf%nSkip_dW+1
                    sf%restart_iter=.TRUE.
                    SWRITE(UNIT_stdOut,'(8X,I8,A,E8.2,A,E8.1,A,E11.4)')sf%iter,'...deltaW=',sf%deltaW,'>',sf%dW_allowed,&
                    '*W_MHD3D_0, skip step and decrease stepsize to dt=',sf%dt
                END IF
            END IF !JacCheck
            IF(sf%nStepDecreased.GT.130) THEN ! 0.9^130 ~10^-6
                SWRITE(UNIT_stdOut,'(A,E21.11)')'Iteration stopped since timestep has been decreased by 0.9^130: ', sf%dt
                SWRITE(UNIT_stdOut,fmt_sep)
                RETURN
            END IF
            IF((MOD(sf%iter,sf%outputIter).EQ.0).AND.(sf%lastoutputIter.NE.sf%iter))THEN
                __PERFON('output')
                SWRITE(UNIT_stdOut,'(A)')'##########################  OUTPUT ##################################'
                CALL Analyze(sf%iter, sf%dofs(0), sf%force(0))
                CALL WriteState(sf%dofs(0),sf%iter)
                SWRITE(UNIT_stdOut,'(A)')'#####################################################################'
                sf%lastOutputIter=sf%iter
                __PERFOFF('output')
            END IF
        END DO !iter

        IF(sf%iter.GE.MaxIter)THEN
            SWRITE(UNIT_stdOut,'(A,E21.11)')"maximum iteration count exceeded"
        END IF
        SWRITE(UNIT_stdOut,'(A)') "... DONE."
        SWRITE(UNIT_stdOut,fmt_sep)
        IF(sf%lastoutputIter.NE.sf%iter)THEN
            CALL Analyze(MIN(sf%iter,MaxIter), sf%dofs(0), sf%force(0))
            CALL WriteState(sf%dofs(0),MIN(sf%iter,MaxIter))
        END IF
        CALL writeSFLoutfile(sf%dofs(0),MIN(sf%iter,MaxIter))
    END SUBROUTINE MinimizeMHD3D_descent

END MODULE MODgvec_MHD3D_minimize
