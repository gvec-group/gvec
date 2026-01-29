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
    TYPE, ABSTRACT :: a_minimizer_vars
        CHARACTER(LEN=40) :: MinimizerType
        LOGICAL   :: restart_iter, logger_is_initialized
        INTEGER   :: JacCheck
        INTEGER   :: iter,nStepDecreased,nSkip_Jac,nSkip_dw
        INTEGER   :: lastoutputIter, logiter_ramp, logscreen
        REAL(wp)  :: dt,deltaW, dW_allowed
        REAL(wp)  :: t_pseudo,Fnorm(3),Fnorm0(3),Fnorm_old(3),W_MHD3D_0

        ! Logging variables
        REAL(wp) :: min_dt_out,max_dt_out,min_dw_out,max_dw_out,sum_dW_out
        LOGICAL  :: DoCheckDistance, DoCheckAxis
        INTEGER  :: outputIter, nlogScreen, logIter, logUnit, StartTimeArray(8)

        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: dofs(:)      !! degrees of freedom at levels (k-1),(k),(k+1)
        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: force(:)     !! force
        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: temp_dofs(:) !! temporary for update
    END TYPE a_minimizer_vars
    TYPE, EXTENDS(a_minimizer_vars) :: t_gradient_descent_vars
    END TYPE t_gradient_descent_vars

    TYPE, EXTENDS(a_minimizer_vars) :: t_accelerated_gradient_descent_vars
        REAL(wp) :: Vnorm(3)
        REAL(wp), ALLOCATABLE  :: tau(:)
        INTEGER   :: ndamp
        REAL(wp)  :: tau_bar
        TYPE(t_sol_var_MHD3D), ALLOCATABLE :: velocity(:)       !! 'velocity' in minimizer
    END TYPE t_accelerated_gradient_descent_vars

    TYPE :: t_minimizer_mhd3d
        !-------------------------------------------------------------------------------------------------------------------------------
        LOGICAL   :: initialized
        INTEGER   :: MinType
        CLASS(a_minimizer_vars), ALLOCATABLE :: vars  !!! depend on the MinimizerType
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
        IF(MinType_in.EQ.0)THEN
            ALLOCATE(t_gradient_descent_vars :: sf%vars)
            sf%vars%MinimizerType="Gradient-Descent"
        ELSEIF(MinType_in.EQ.10)THEN
            ALLOCATE(t_accelerated_gradient_descent_vars :: sf%vars)
            sf%vars%MinimizerType="Accelerated Gradient-Descent"
        ELSE
            CALL abort(__STAMP__,&
                       "requested MinimizeType does not exist, expecting 0 or 10",intinfo=MinType_in,&
                       TypeInfo="InvalidParameterError")
        END IF
        ASSOCIATE(vars=>sf%vars)
            vars%dt = dt_initial
            vars%dW_allowed = dW_allowed_in

            ALLOCATE(vars%dofs(-3:1))
            CALL vars%dofs(1)%init(varsize_in)
            DO i=-3,0
                CALL vars%dofs(i)%copy(vars%dofs(1))
            END DO
            ALLOCATE(vars%force(-1:0))
            DO i=-1,0
                CALL vars%force(i)%copy(vars%dofs(1))
            END DO

            ALLOCATE(vars%temp_dofs(-1:1))
            DO i=-1,1
                CALL vars%temp_dofs(i)%copy(vars%dofs(1))
            END DO

            ! variables for the accelerated Gradient descent
            SELECT TYPE(vars)
            TYPE IS(t_accelerated_gradient_descent_vars)
                ALLOCATE(vars%velocity(-1:1))
                DO i=-1,1
                    CALL vars%velocity(i)%copy(vars%dofs(1))
                END DO
                vars%tau_bar = 0.075_wp
                vars%ndamp = 10
                ALLOCATE(vars%tau(vars%ndamp))
                vars%tau=0.0_wp
                vars%Vnorm = 0.0_wp
            END SELECT !Type
            vars%JacCheck = 1

            vars%nstepDecreased = 0
            vars%nSkip_Jac = 0
            vars%t_pseudo = 0
            vars%lastOutputIter = 0
            vars%iter = 0
            vars%logiter_ramp = 1
            vars%logscreen = 1
            vars%nSkip_dw = 0
            vars%JacCheck = 1



            vars%restart_iter = .TRUE.
            vars%logger_is_initialized = .FALSE.

            vars%DoCheckDistance = DoCheckDistance
            vars%DoCheckAxis = DoCheckAxis

            vars%outputIter = outputIter
            vars%nlogScreen = nlogScreen
            vars%logIter = logIter

            vars%Fnorm = 10.0_wp
            vars%Fnorm0 = 10.0_wp
            vars%W_MHD3D_0 = 0.0_wp
            vars%deltaW = 0.0_wp
        END ASSOCIATE !vars
        sf%initialized = .TRUE.
    END SUBROUTINE new_minimizer

    SUBROUTINE Free_minimizer(sf)
        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        INTEGER :: i

        ASSOCIATE(vars=>sf%vars)
            CLOSE(vars%logUnit)
            DO i=-1,1
                IF(ALLOCATED(vars%dofs)) CALL vars%dofs(i)%free()
                IF(ALLOCATED(vars%temp_dofs)) CALL vars%temp_dofs(i)%free()
            END DO
            DO i=-1,0
                IF(ALLOCATED(vars%force)) CALL vars%force(i)%free()
            END DO

            SDEALLOCATE(vars%dofs)
            SDEALLOCATE(vars%temp_dofs)
            SDEALLOCATE(vars%force)
            SELECT TYPE(vars)
            TYPE IS(t_accelerated_gradient_descent_vars)
                DO i=-1,1
                    IF(ALLOCATED(vars%velocity)) CALL vars%velocity(i)%free()
                END DO
                SDEALLOCATE(vars%velocity)
                SDEALLOCATE(vars%tau)
            END SELECT !Type
            DEALLOCATE(vars)
        END ASSOCIATE !vars
    END SUBROUTINE Free_minimizer

    SUBROUTINE MinimizeMHD3d_ResetDescent(sf)
        USE MODgvec_MHD3D_EvalFunc, ONLY: EvalAux, EvalForce, EvalEnergy
        !-------------------------------------------------------------------------------------------------------------------------------
        CLASS(t_minimizer_mhd3d), INTENT(INOUT) :: sf
        !-------------------------------------------------------------------------------------------------------------------------------
        ASSOCIATE(vars=>sf%vars)
            vars%JacCheck=1 !abort if detJ<0
            CALL EvalAux(vars%dofs(0), vars%JacCheck)
            vars%dofs(0)%W_MHD3D= EvalEnergy(vars%dofs(0),.FALSE.,vars%JacCheck)
            vars%W_MHD3D_0 = vars%dofs(0)%W_MHD3D
            CALL EvalForce( vars%dofs(0),.FALSE.,vars%JacCheck,vars%force(0))
            vars%Fnorm0 = SQRT(vars%force(0)%norm_2())
            vars%Fnorm = vars%Fnorm0
            vars%Fnorm_old = 1.1_wp*vars%Fnorm0
            CALL vars%dofs(-1)%set_to(vars%dofs(0)) !last state
            CALL vars%dofs(-2)%set_to(vars%dofs(0)) !state at last logging interval

            !for hirshman method
            SELECT TYPE(vars)
            TYPE IS(t_accelerated_gradient_descent_vars)
                CALL vars%velocity(-1)%set_to(0.0_wp)
                CALL vars%velocity( 0)%set_to(0.0_wp)
                vars%tau(1:vars%ndamp)=0.15_wp/vars%dt
                vars%tau_bar = 0.075_wp
            END SELECT

            vars%min_dt_out=1.0e+30_wp
            vars%max_dt_out=0.0_wp
            vars%min_dW_out=1.0e+30_wp
            vars%max_dW_out=-1.0e+30_wp
            vars%sum_dW_out=0.0_wp
            vars%nSkip_dW =0
            IF(vars%restart_iter) vars%restart_iter=.FALSE.
        END ASSOCIATE !vars
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
        ASSOCIATE(vars=>sf%vars)
            W_MHD3D=vars%dofs(0)%W_MHD3D
            CALL DATE_AND_TIME(values=TimeArray) ! get System time
            WRITE(UNIT_stdOut,'(A,E11.4,A)')'%%%%%  START ITERATION, dt= ',vars%dt, '  %%%%%%%%%%%%%%%%%'
            WRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                            '%% Sys date : ',timeArray(1:3),timeArray(5:7)
            WRITE(UNIT_stdOut,'(A,3E21.14)') &
                    '%% dU = |Force|= ',vars%Fnorm(1:3)
            SELECT TYPE(vars)
            TYPE IS (t_accelerated_gradient_descent_vars)
                WRITE(UNIT_stdOut,'(A,E11.4,A,3E11.4)') &
                    '%% accel.GD: tau= ',vars%tau_bar,' |vel|= ',vars%Vnorm(1:3)
            END SELECT

            WRITE(UNIT_stdOut,'(40(" -"))')
            !------------------------------------
            vars%StartTimeArray=TimeArray !save first time stamp

            vars%logUnit=GETFREEUNIT()
            WRITE(FileString,'("logMinimizer_",A,"_",I4.4,".csv")')TRIM(ProjectName),outputLevel
            OPEN(UNIT     = vars%logUnit       ,&
                FILE     = TRIM(FileString) ,&
                STATUS   = 'REPLACE'   ,&
                ACCESS   = 'SEQUENTIAL' )
            !header
            iLogDat=0
            WRITE(vars%logUnit,'(A)',ADVANCE="NO")'"#iterations","runtime(s)","min_dt","max_dt"'
            WRITE(vars%logUnit,'(A)',ADVANCE="NO")',"W_MHD3D","min_dW","max_dW","sum_dW"'
            WRITE(vars%logUnit,'(A)',ADVANCE="NO")',"normF_X1","normF_X2","normF_LA"'
            LogDat(ilogDat+1:iLogDat+11)=(/0.0_wp,0.0_wp,vars%dt,vars%dt,W_MHD3D,0.0_wp,0.0_wp,0.0_wp,vars%Fnorm(1:3)/)
            iLogDat=11
            SELECT TYPE(vars)
            TYPE IS (t_accelerated_gradient_descent_vars)
                WRITE(vars%logUnit,'(A)',ADVANCE="NO")',"tau","normV_X1","normV_X2","normV_LA"'
                LogDat(ilogDat+1:iLogDat+4)=(/vars%tau_bar,vars%Vnorm(1:3)/)
                iLogDat=iLogDat+4
            END SELECT
            IF(vars%doCheckDistance) THEN
                WRITE(vars%logUnit,'(A)',ADVANCE="NO")',"max_Dist","avg_Dist"'
                LogDat(iLogDat+1:iLogDat+2)=(/0.0_wp,0.0_wp/)
                iLogDat=iLogDat+2
            END IF!doCheckDistance
            IF(vars%doCheckAxis) THEN
                WRITE(vars%logUnit,'(A)',ADVANCE="NO")',"X1_axis_0","X2_axis_0","X1_axis_1","X2_axis_1"'
                CALL CheckPos(vars%dofs(0),0.0_wp,2,Pos)
                LogDat(iLogDat+1:iLogDat+4)=RESHAPE(Pos,(/4/))
                iLogDat=iLogDat+4
            END IF!doCheckAxis

            WRITE(vars%logUnit,'(A)')' '
            !first data line
            WRITE(vars%logUnit,'(*(e23.15,:,","))') logDat(1:iLogDat)
            vars%logger_is_initialized = .TRUE.
        END ASSOCIATE !vars
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
        ASSOCIATE(vars=>sf%vars)
            CALL DATE_AND_TIME(values=TimeArray) ! get System time
            W_MHD3D=vars%dofs(0)%W_MHD3D
            IF(.NOT.quiet)THEN
                WRITE(UNIT_stdOut,'(80("%"))')
                WRITE(UNIT_stdOut,'(A,I4.2,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
                                '%% Sys date : ',timeArray(1:3),timeArray(5:7)
                WRITE(UNIT_stdOut,'(A,I8,A,2I8,A,E11.4,A,2E11.4,A,E21.14,A,3E12.4)') &
                                '%% #ITERATIONS= ',vars%iter,', #skippedIter (Jac/dW)= ',vars%nSkip_Jac,vars%nSkip_dW, &
                        '\n%% t_pseudo= ',vars%t_pseudo,', min/max dt= ',vars%min_dt_out,vars%max_dt_out, &
                        '\n%% W_MHD3D= ',W_MHD3D,', min/max/sum deltaW= ' , vars%min_dW_out,vars%max_dW_out,vars%sum_dW_out
                WRITE(UNIT_stdOut,'(A,3E21.14)') &
                            '%% dU = |Force|= ',vars%Fnorm(1:3)
                !------------------------------------
            END IF!.NOT.quiet

            iLogDat=0
            runtime_ms=MAX(0,SUM((timeArray(5:8)-vars%StartTimearray(5:8))*(/360000,6000,100,1/)))
            LogDat(ilogDat+1:iLogDat+11)=(/REAL(vars%iter,wp),REAL(runtime_ms,wp)/100.0_wp, &
                                            vars%min_dt_out,vars%max_dt_out,&
                                            W_MHD3D,vars%min_dW_out,vars%max_dW_out,vars%sum_dW_out, &
                                            vars%Fnorm(1:3)/)
            iLogDat=11
            SELECT TYPE(vars)
            TYPE IS (t_accelerated_gradient_descent_vars)
                IF(.NOT.quiet)THEN
                WRITE(UNIT_stdOut,'(A,E11.4,A,3E11.4)') &
                        '%% accel.GD: tau= ',vars%tau_bar,' |vel|= ',vars%Vnorm(1:3)
                END IF!.NOT.quiet
                LogDat(ilogDat+1:iLogDat+4)=(/vars%tau_bar,vars%Vnorm(1:3)/)
                iLogDat=iLogDat+4
            END SELECT
            IF(vars%doCheckDistance) THEN
                CALL CheckDistance(vars%dofs(0),vars%dofs(-2),maxDist,avgDist)
                CALL vars%dofs(-2)%set_to(vars%dofs(0))
                IF(.NOT.quiet)THEN
                WRITE(UNIT_stdOut,'(A,2E11.4)') &
                '               %% Dist to last log (max/avg) : ',maxDist,avgDist
                END IF!.NOT.quiet
                LogDat(iLogDat+1:iLogDat+2)=(/maxDist,avgDist/)
                iLogDat=iLogDat+2
            END IF!doCheckDistance
            IF(vars%doCheckAxis) THEN
                CALL CheckPos(vars%dofs(0),0.0_wp,2,Pos)
                IF(.NOT.quiet)THEN
                SWRITE(UNIT_stdOut,'(2(A,2E22.14))') &
                    '%% axis position (X1,X2,zeta=0     ): ',Pos(1:2,1), &
                '\n%% axis position (X1,X2,zeta=pi/nfp): ',Pos(1:2,2)
                END IF!.NOT.quiet
                LogDat(iLogDat+1:iLogDat+4)=RESHAPE(Pos,(/4/))
                iLogDat=iLogDat+4
            END IF !doCheckAxis

            IF(.NOT.quiet)THEN
                WRITE(UNIT_stdOut,'(40(" -"))')
            END IF!.NOT.quiet
            WRITE(vars%logUnit,'(*(e23.15,:,","))') logDat(1:iLogDat)
        END ASSOCIATE !vars
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
        INTEGER :: maxIter
        !=================================================================================================================================
        ASSOCIATE(vars=>sf%vars)
            maxIter = vars%iter + maxIter_in
            DO WHILE(vars%iter.LT.maxIter)
                IF((vars%restart_iter))THEN
                    CALL sf%reset()
                END IF !before first iteration or after restart Jac<0

                IF(.NOT.vars%logger_is_initialized)THEN
                    CALL sf%StartLogging()
                END IF
                !COMPUTE NEW SOLUTION P(1) as a prediction
                SELECT TYPE(vars)
                TYPE IS(t_gradient_descent_vars) !gradient descent, previously used for minimizerType=0
                        CALL vars%temp_dofs(1)%AXBY(1.0_wp,vars%dofs(0),vars%dt,vars%force(0)) !overwrites P(1), predicts solution U(1)
                TYPE IS(t_accelerated_gradient_descent_vars) !hirshman method
                        !tau is damping parameter
                        vars%tau(1:vars%ndamp-1) = vars%tau(2:vars%ndamp) !save old

                        !ln(|F_n|^2/|F_{n-1}|^2), Fnorm=|F_X1|,|F_X2|,|F_LA|
                        vars%tau(vars%ndamp)  = MIN(0.15_wp,ABS(LOG(SUM(vars%Fnorm**2)/SUM(vars%Fnorm_old**2))))/vars%dt

                        vars%tau_bar = 0.5_wp*vars%dt*SUM(vars%tau)/REAL(vars%ndamp,wp)   !=1/2 * tauavg
                        CALL vars%velocity(1)%AXBY(((1.0_wp-vars%tau_bar)/(1.0_wp+vars%tau_bar)),vars%velocity(0),(vars%dt/(1.0_wp+vars%tau_bar)),vars%force(0)) !velocity V(1)
                        CALL vars%temp_dofs(1)%AXBY(1.0_wp,vars%dofs(0),vars%dt,vars%velocity(1)) !overwrites P(1), predicst solution U(1)
                        vars%Vnorm=SQRT(vars%velocity(1)%norm_2())
                END SELECT !Type

                vars%JacCheck=2 !no abort,if detJ<0, JacCheck=-1

                vars%temp_dofs(1)%W_MHD3D=EvalEnergy(vars%temp_dofs(1),.TRUE.,vars%JacCheck)
                IF(vars%JacCheck.EQ.-1)THEN
                    vars%dt=0.9_wp*vars%dt
                    vars%nstepDecreased=vars%nStepDecreased+1
                    vars%nSkip_Jac=vars%nSkip_Jac+1
                    vars%restart_iter=.TRUE.
                    CALL vars%dofs(0)%set_to(vars%dofs(-3)) !reset to initial state
                    SWRITE(UNIT_stdOut,'(8X,I8,A,E11.4,A)')vars%iter,'...detJac<0, decrease stepsize to dt=',vars%dt,  ' and RESTART simulation!!!!!!!'
                ELSE
                    !detJ>0
                    vars%deltaW=vars%temp_dofs(1)%W_MHD3D-vars%dofs(0)%W_MHD3D!should be <=0,
                    IF(vars%deltaW.LE.vars%dW_allowed*vars%W_MHD3D_0)THEN !valid step /hirshman method accept W increase!
                        IF(ALL(vars%Fnorm.LE.abstol))THEN
                            CALL sf%Logging(.FALSE.)
                            SWRITE(UNIT_stdOut,'(4x,A)')'==>Iteration finished, |force| in relative tolerance'
                            EXIT !DO LOOP
                        END IF
                        vars%iter=vars%iter+1
                        vars%t_pseudo=vars%t_pseudo+vars%dt
                        ! for simple gradient & hirshman
                        CALL vars%dofs(-1)%set_to(vars%dofs(0))
                        CALL vars%dofs(0)%set_to(vars%temp_dofs(1))
                        ! for hirshman method
                        SELECT TYPE(vars)
                        TYPE IS(t_accelerated_gradient_descent_vars)
                            CALL vars%velocity(-1)%set_to(vars%velocity(0))
                            CALL vars%velocity(0)%set_to(vars%velocity(1))
                        END SELECT
                        CALL EvalForce(vars%temp_dofs(1),.FALSE.,vars%JacCheck,vars%force(0)) !evalAux was already called on P(1)=U(0), so that its set false here.
                        vars%Fnorm_old=vars%Fnorm
                        vars%Fnorm=SQRT(vars%force(0)%norm_2())
                        vars%nstepDecreased=0
                        vars%min_dt_out=MIN(vars%min_dt_out,vars%dt)
                        vars%max_dt_out=MAX(vars%max_dt_out,vars%dt)
                        vars%min_dW_out=MIN(vars%min_dW_out,vars%deltaW)
                        vars%max_dW_out=MAX(vars%max_dW_out,vars%deltaW)
                        vars%sum_dW_out=vars%sum_dW_out+vars%deltaW
                        IF(MOD(vars%iter,vars%logIter_ramp).EQ.0)THEN
                            CALL sf%Logging(.NOT.((vars%logIter_ramp.GE.vars%logIter).AND.(MOD(vars%logscreen,vars%nLogScreen).EQ.0)))
                            IF(.NOT.(vars%logIter_ramp.LT.vars%logIter))THEN !only reset for logIter
                                vars%logscreen=vars%logscreen+1
                                vars%min_dt_out=1.0e+30_wp
                                vars%max_dt_out=0.0_wp
                                vars%min_dW_out=1.0e+30_wp
                                vars%max_dW_out=-1.0e+30_wp
                                vars%sum_dW_out=0.0_wp
                                vars%nSkip_dW =0
                            END IF
                            vars%logIter_ramp=MIN(vars%logIter,vars%logIter_ramp*2)
                        END IF
                    ELSE !not a valid step, decrease timestep and skip P(1)
                        vars%dt=0.9_wp*vars%dt
                        vars%nstepDecreased=vars%nStepDecreased+1
                        vars%nSkip_dW=vars%nSkip_dW+1
                        vars%restart_iter=.TRUE.
                        SWRITE(UNIT_stdOut,'(8X,I8,A,E8.2,A,E8.1,A,E11.4)')vars%iter,'...deltaW=',vars%deltaW,'>',vars%dW_allowed,&
                        '*W_MHD3D_0, skip step and decrease stepsize to dt=',vars%dt
                    END IF
                END IF !JacCheck
                IF(vars%nStepDecreased.GT.130) THEN ! 0.9^130 ~10^-6
                    SWRITE(UNIT_stdOut,'(A,E21.11)')'Iteration stopped since timestep has been decreased by 0.9^130: ', vars%dt
                    SWRITE(UNIT_stdOut,fmt_sep)
                    RETURN
                END IF
                IF((MOD(vars%iter,vars%outputIter).EQ.0).AND.(vars%lastoutputIter.NE.vars%iter))THEN
                    __PERFON('output')
                    SWRITE(UNIT_stdOut,'(A)')'##########################  OUTPUT ##################################'
                    CALL Analyze(vars%iter, vars%dofs(0), vars%force(0))
                    CALL WriteState(vars%dofs(0),vars%iter)
                    SWRITE(UNIT_stdOut,'(A)')'#####################################################################'
                    vars%lastOutputIter=vars%iter
                    __PERFOFF('output')
                END IF
            END DO !iter

            IF(vars%iter.GE.MaxIter)THEN
                SWRITE(UNIT_stdOut,'(A,E21.11)')"maximum iteration count exceeded"
            END IF
            SWRITE(UNIT_stdOut,'(A)') "... DONE."
            SWRITE(UNIT_stdOut,fmt_sep)
            IF(vars%lastoutputIter.NE.vars%iter)THEN
                CALL Analyze(MIN(vars%iter,MaxIter), vars%dofs(0), vars%force(0))
                CALL WriteState(vars%dofs(0),MIN(vars%iter,MaxIter))
            END IF
            CALL writeSFLoutfile(vars%dofs(0),MIN(vars%iter,MaxIter))
        END ASSOCIATE !vars
    END SUBROUTINE MinimizeMHD3D_descent

END MODULE MODgvec_MHD3D_minimize
