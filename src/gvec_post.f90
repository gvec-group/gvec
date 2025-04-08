!===================================================================================================================================
! Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
! License: MIT
!===================================================================================================================================
#include "defines.h"


!===================================================================================================================================
!>
!!# **GVEC** Driver program
!!
!===================================================================================================================================
PROGRAM GVEC_POST
  USE MODgvec_MPI          ,ONLY : par_Init, par_Finalize,par_bcast
  USE MODgvec_Globals
  USE MODgvec_Analyze      ,ONLY: InitAnalyze,Analyze,FinalizeAnalyze
  USE MODgvec_Output       ,ONLY: InitOutput,FinalizeOutput
  USE MODgvec_Output_vars  ,ONLY: OutputLevel
  USE MODgvec_Restart      ,ONLY: InitRestart,FinalizeRestart
  USE MODgvec_Restart      ,ONLY: RestartFromState
  USE MODgvec_Output_Vars  ,ONLY: OutputLevel,ProjectName
  USE MODgvec_ReadState_Vars,ONLY: fileID_r,outputLevel_r
  USE MODgvec_MHD3D_Vars   ,ONLY: U,F
  USE MODgvec_MHD3D_visu   ,ONLY:WriteSFLoutfile
  USE MODgvec_MHD3D_EvalFunc , ONLY: InitProfilesGP,EvalEnergy,EvalForce
  USE MODgvec_ReadInTools  ,ONLY: FillStrings,GETLOGICAL,GETINT,IgnoredStrings
  USE MODgvec_Functional
!$ USE omp_lib
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!local variables
  INTEGER                 :: iArg,nArgs
  CHARACTER(LEN=255)      :: Parameterfile
  CHARACTER(LEN=255)      :: Statefile
  INTEGER                 :: which_functional
  INTEGER                 :: JacCheck
  CLASS(t_functional),ALLOCATABLE   :: functional
  REAL(wp)                :: StartTime,EndTime
!===================================================================================================================================
  CALL par_Init()
  __PERFINIT
  __PERFON('main')
  nArgs=COMMAND_ARGUMENT_COUNT()
  IF ((nArgs.LT.2))THEN
    ! Print out error message containing valid syntax
    STOP 'ERROR - Invalid syntax. Please use: gvec_post parameter.ini [Statefiles*] '
  END IF
  CALL GET_COMMAND_ARGUMENT(1,Parameterfile)

  CALL CPU_TIME(StartTime)
!$ StartTime=OMP_GET_WTIME()

  !header
  SWRITE(Unit_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A)') "GVEC POST ! GVEC POST ! GVEC POST ! GVEC POST"
  SWRITE(Unit_stdOut,'(132("="))')
  !.only executes if compiled with OpenMP
!$ SWRITE(UNIT_stdOut,'(A,I6)')'   Number of OpenMP threads : ',OMP_GET_MAX_THREADS()
!$ SWRITE(Unit_stdOut,'(132("="))')
  !.only executes if compiled with MPI
# if MPI
  SWRITE(UNIT_stdOut,'(A,I6)')'   Number of MPI tasks : ',nRanks
  SWRITE(Unit_stdOut,fmt_sep)
  IF(nRanks.GT.1) CALL abort(__STAMP__,&
                   "GVEC post is compiled with MPI, but can only be called with 1 MPI rank." )
# endif
#include  "configuration-cmake.f90"
  SWRITE(Unit_stdOut,fmt_sep)

  CALL FillStrings(ParameterFile) !< readin parameterfile, done on MPI root + Bcast
  testdbg =.FALSE.
  testlevel=-1

  !initialization phase
  CALL InitRestart()
  CALL InitOutput()
  CALL InitAnalyze()

  which_functional=GETINT('which_functional', Proposal=1 )
  CALL InitFunctional(functional,which_functional)

  CALL IgnoredStrings()
  DO iArg=2,nArgs
    CALL GET_COMMAND_ARGUMENT(iArg,StateFile)
    SWRITE(Unit_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A,I4,A4,I4,A3,A)') 'Post-Analyze StateFile ',iArg-1,' of ',nArgs-1,' : ',TRIM(StateFile)
    SWRITE(Unit_stdOut,'(132("-"))')
    ProjectName='POST_'//TRIM(StateFile(1:INDEX(StateFile,'_State_')-1))
    CALL RestartFromState(StateFile,U(0))
    outputLevel=outputLevel_r
    JacCheck=2
    !...check this: temporarily commented for gvec_post to run with MPI version...
    CALL InitProfilesGP() !evaluate profiles once at Gauss Points (on MPIroot + BCast)
    U(0)%W_MHD3D=EvalEnergy(U(0),.TRUE.,JacCheck)
    CALL EvalForce(U(0),.FALSE.,JacCheck, F(0))
    CALL Analyze(FileID_r)
    CALL writeSFLoutfile(U(0),FileID_r)
  END DO !iArg

  CALL FinalizeFunctional(functional)

  CALL FinalizeAnalyze()
  CALL FinalizeOutput()
  CALL FinalizeRestart()

  CALL CPU_TIME(EndTime)
!$ EndTime=OMP_GET_WTIME()
  WRITE(Unit_stdOut,'(132("="))')
  WRITE(Unit_stdOut,'(A,F8.2,A)') ' GVEC POST FINISHED ! [',EndTime-StartTime,' sec ]'
  WRITE(Unit_stdOut,'(132("="))')

  __PERFOFF('main')
  __PERFOUT('main')
  CALL par_Finalize()

END PROGRAM GVEC_POST
