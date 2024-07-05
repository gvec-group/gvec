# =========================================================================
# GVEC CMake configuration
# =========================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)



# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "marconi" ON_MARCONI)
STRING(FIND "${CMAKE_HOSTNAME}" "testimony" ON_TESTIMONY)
STRING(FIND "${CMAKE_HOSTNAME}" "sharedrunner" ON_SHAREDRUNNER)
STRING(FIND "${CMAKE_HOSTNAME}" "mpcdfcirunner" ON_MPCDFCIRUNNER)
STRING(FIND "${CMAKE_HOSTNAME}" "cobra"     ON_COBRA)
STRING(FIND "${CMAKE_HOSTNAME}" "raven"     ON_RAVEN)
STRING(FIND "${CMAKE_HOSTNAME}" "draco"     ON_DRACO)
STRING(FIND "${CMAKE_HOSTNAME}" "mac_brew"  ON_MAC_BREW)
STRING(FIND "${CMAKE_HOSTNAME}" "mac_ports" ON_MAC_PORTS)
STRING(FIND "${CMAKE_HOSTNAME}" "tokp"      ON_TOKP )

IF (ON_MPCDFCIRUNNER GREATER -1)
  MESSAGE(STATUS "Compiling on MPCDF CI image shared runner")
ELSEIF (ON_SHAREDRUNNER GREATER -1)
  MESSAGE(STATUS "Compiling on MPCDF shared runner")
ELSEIF (ON_TESTIMONY GREATER -1)
  MESSAGE(STATUS "Compiling on testimony")
ELSEIF (ON_COBRA GREATER -1)
  MESSAGE(STATUS "Compiling on MPCDF cobra")
ELSEIF (ON_DRACO GREATER -1)
  MESSAGE(STATUS "Compiling on MPCDF draco")
ELSEIF (ON_RAVEN GREATER -1)
  MESSAGE(STATUS "Compiling on MPCDF raven")
ELSEIF (ON_MAC_BREW GREATER -1)
  MESSAGE(STATUS "Compiling on MacOS with homebrew!")
ELSEIF (ON_MAC_BREW_PORTS GREATER -1)
  MESSAGE(STATUS "Compiling on MacOS with macports!")
ELSEIF (ON_TOKP GREATER -1)
  MESSAGE(STATUS "Compiling on tokp")
ELSEIF (ON_MARCONI GREATER -1)
  MESSAGE(STATUS "Compiling on Marconi")
ELSE()
  MESSAGE(STATUS "Compiling on unknown hostname")
ENDIF()


# =========================================================================
# Project definitions and output paths
# =========================================================================
PROJECT(GVEC)
ENABLE_LANGUAGE(Fortran)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Set custom install dir (needs to be done after project has been defined!)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)

# =========================================================================
# Options
# =========================================================================

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

#generate lists of subfolders to be used as choices for the variables of the configuration

#SUBDIRLIST(FUNCTIONALS  ${CMAKE_CURRENT_SOURCE_DIR}/src/functionals/)
#
#SET(   FUNCTIONAL   "mhd3d"      CACHE STRING
#                        "Name of functional to be used, must be identical to subdirectory name in src/functionals folder.")
#SET_PROPERTY(CACHE FUNCTIONAL PROPERTY STRINGS ${FUNCTIONALS})


# =========================================================================
# PREPROCESSOR FLAGS
# =========================================================================

OPTION(COMPILE_GVEC "compile gvec programme for equilibrium computations (libgvec and gvec executable)" ON)
IF(COMPILE_GVEC)
  OPTION(LINK_GVEC_TO_NETCDF "Link to local serial NetCDF library, needed for VMEC input/output" ON)
ELSE()
  UNSET(LINK_GVEC_TO_NETCDF CACHE ) # if not mhd
ENDIF()
OPTION(COMPILE_GVEC_TO_GENE "compile interface gvec_to_gene library (added to libgvec)" OFF)
OPTION(COMPILE_GVEC_TO_HOPR "compile interface gvec_to_hopr library (added to libgvec)" OFF)
OPTION(COMPILE_GVEC_TO_CASTOR3D "compile converter gvec_to_castor (only executable)" OFF)
OPTION(COMPILE_GVEC_TO_JOREK "compile converter gvec_to_jorek (only executable)" OFF)

OPTION(USE_OPENMP "enable OpenMP" ON)
OPTION(USE_MPI "enable MPI" OFF)

IF(USE_MPI)
  IF(COMPILE_GVEC_TO_GENE)
    MESSAGE(SEND_ERROR " IF COMPILE_GVEC_TO_GENE=ON, set USE_MPI=OFF (MPI not needed/tested)")
  ELSEIF(COMPILE_GVEC_TO_HOPR)
    MESSAGE(SEND_ERROR " IF COMPILE_GVEC_TO_HOPR=ON, set USE_MPI=OFF (MPI not needed/tested)")
  ELSEIF(COMPILE_GVEC_TO_CASTOR3D) 
    MESSAGE(SEND_ERROR " IF COMPILE_GVEC_TO_CASTOR3D=ON, set USE_MPI=OFF (MPI not needed/tested)")
  ELSEIF(COMPILE_GVEC_TO_JOREK) 
    MESSAGE(SEND_ERROR " IF COMPILE_GVEC_TO_JOREK=ON, set USE_MPI=OFF (MPI not needed/tested)")
  ENDIF()
ENDIF()

OPTION(USE_FTIMINGS "enable profiling with ftimings" OFF)
MARK_AS_ADVANCED(USE_FTIMINGS)
SET (USE_VISU_J "off" CACHE STRING "enable visualization of current: off/ with Finite Difference / analytical")
SET_PROPERTY(CACHE USE_VISU_J PROPERTY STRINGS off FD ANA )

IF(USE_VISU_J STREQUAL "FD")
  ADD_DEFINITIONS("-DVISU_J_FD")
ELSEIF(USE_VISU_J STREQUAL "ANA")
  ADD_DEFINITIONS("-DVISU_J_EXACT")
ENDIF()
# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug / Release / Profile."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile )
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS("-DDEBUG")
ENDIF()


# =========================================================================
# MPI
# =========================================================================
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF (NOT MPI_Fortran_NO_INTERROGATE)
    FOREACH(DIR ${MPI_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    FOREACH(DIR ${MPI_Fortran_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    LIST(APPEND linkedlibs ${MPI_Fortran_LIBRARIES})
  ENDIF()
  MARK_AS_ADVANCED(MPI_LIBRARY MPI_EXTRA_LIBRARY) # don't show these ones
  MESSAGE(STATUS "Building GVEC with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
  ADD_DEFINITIONS(-DMPI=1)
  OPTION(NO_MPI_F08 "disabe MPI F2008" OFF)
  MARK_AS_ADVANCED(NO_MPI_F08)
  IF(NO_MPI_F08)
    ADD_DEFINITIONS(-DNO_MPI_F08)
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DMPI=0)
  UNSET(NO_MPI_F08 CACHE ) # if not MPI
ENDIF()


# =========================================================================
# Location of binary and filenames
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see defines.h)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")


# add basedir for location of corresponding userblock-file
#ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  #environment specific variables
  IF(ONHORNET GREATER -1)
    SET(MAIN_INSTRUCTION "-march=core-avx2")
  ELSE()
    SET(MAIN_INSTRUCTION "-march=native")
  ENDIF()
  # set Flags

  #OpenMP
  IF(USE_OPENMP)
   SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fopenmp")
  ENDIF()

  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${MAIN_INSTRUCTION} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${MAIN_INSTRUCTION} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -ffpe-trap=invalid -fbounds-check -finit-real=snan -fbacktrace -Wall -Wextra")
  # add flags only for compiling not linking!
  SET (MAIN_COMPILE_FLAGS "-xf95-cpp-input -fPIC")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  #environment specific variables
  IF((ON_TESTIMONY GREATER -1) OR (ON_SHAREDRUNNER GREATER -1) OR (ON_MPCDFCIRUNNER GREATER -1 ))
     SET(MAIN_INSTRUCTION "-xCORE-AVX-I ") #for testimony
  ELSEIF (ON_TOKP GREATER -1 )
     SET(MAIN_INSTRUCTION "") #no additional options with newest intel compiler on tokp
  ELSEIF ((ON_COBRA GREATER -1 ) OR (ON_RAVEN GREATER -1))
     SET(MAIN_INSTRUCTION "-xCORE-AVX512 -qopt-zmm-usage=high -align array64byte") #for SKL
     # add -parallel for automatic parallelization (choose OMP_NUM_THREADS=16)
  ELSEIF (ON_DRACO GREATER -1 )
     SET(MAIN_INSTRUCTION "-xCORE-AVX2 -align array64byte") #for SKL
  ELSEIF(ON_MARCONI GREATER -1)
    SET(MAIN_INSTRUCTION "-axMIC-AVX512 ") #for BW and KNL
  ELSE()
     SET(MAIN_INSTRUCTION "-xHost")
  ENDIF()
  # set Flags
  #OpenMP
  IF(USE_OPENMP)
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -qopenmp")
  ENDIF()

  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 ${MAIN_INSTRUCTION} -qopt-report=5 -diag-disable=10397 -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${MAIN_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
   # add flags only for compiling not linking!
   SET (MAIN_COMPILE_FLAGS "-fpp -assume bscc -fPIC")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -f free -s real64 -em -DCRAY")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (MAIN_COMPILE_FLAGS "-F")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -fPIC")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O3")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -g")
  # add flags only for compiling not linking!
  SET (MAIN_COMPILE_FLAGS "-Mpreprocess")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -fPIC -traceback -Wall -DNVHPC")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0  -g ") #-C option fails for LEN(*) STRINGS
  # add flags only for compiling not linking!
  SET (MAIN_COMPILE_FLAGS "-Mpreprocess")
  IF(USE_OPENMP)
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -mp=multicore") # for now not using the gpu for OpenMP!
  ELSE()
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -nomp") 
  ENDIF()


ELSE()
   MESSAGE(SEND_ERROR " Unknown compiler")
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
# use Lapack/Blas for GNU

IF(USE_OPENMP)
  SET (MKLLIB "parallel")
  SET (MKLLIBGNU "mkl_gnu_thread")
ELSE()
  SET (MKLLIB "sequential")
  SET (MKLLIBGNU "mkl_sequential")
ENDIF()


IF(ON_TESTIMONY GREATER -1)
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
    SET(MKL_DIR $ENV{MKL_HOME})
    LIST(APPEND linkedlibs -qmkl=${MKLLIB} -Wl,-rpath=${MKL_DIR}/lib/intel64 )
    MESSAGE(STATUS "Using INTEL MKL BLAS/Lapack library")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    FIND_PACKAGE(LAPACK QUIET)
    IF (LAPACK_FOUND)
      LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
      MESSAGE(STATUS "Using BLAS/Lapack library")
    ELSE()
      MESSAGE(WARNING " Lapack not found for GNU compiler!")
    ENDIF()
  ELSE()
     MESSAGE(SEND_ERROR " Unknown compiler ...")
  ENDIF()
ELSEIF((ON_SHAREDRUNNER GREATER -1) OR (ON_MPCDFCIRUNNER GREATER -1) OR (ON_COBRA GREATER -1) OR (ON_DRACO GREATER -1) OR (ON_RAVEN GREATER -1) OR (ON_TOKP GREATER -1))
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
    SET(MKL_DIR $ENV{MKL_HOME})
    LIST(APPEND linkedlibs -L${MKL_DIR}/lib/intel64 -qmkl=${MKLLIB} -Wl,-rpath=${MKL_DIR}/lib/intel64 )
    MESSAGE(STATUS "Using INTEL MKL BLAS/Lapack library with INTEL")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
    SET(MKL_DIR $ENV{MKL_HOME})
    LIST(APPEND linkedlibs -L${MKL_DIR}/lib/intel64 -lmkl_gf_lp64
		-l${MKLLIBGNU} -lmkl_core -lpthread -lm -ldl
		-Wl,-rpath=${MKL_DIR}/lib/intel64 )
    MESSAGE(STATUS "Using INTEL MKL BLAS/Lapack library with GNU")
  ELSE()
    FIND_PACKAGE(LAPACK QUIET)
    IF (LAPACK_FOUND)
      LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
      MESSAGE(STATUS "Using BLAS/Lapack library")
    ELSE()
      MESSAGE(WARNING " Lapack not found for GNU compiler!")
    ENDIF()
  ENDIF()
  MESSAGE(STATUS "Using INTEL MKL BLAS/Lapack library")
#ELSEIF(ON_MAC_PORTS GREATER -1)
#  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    SET(CMAKE_BUILD_RPATH "/opt/local/lib/libgcc")
#    LIST(APPEND linkedlibs -L/opt/local/lib -lopenblas)  #TTR: don't know why fails... needs line below to work
#    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L/opt/local/lib -lopenblas")
#    MESSAGE(STATUS "Using GNU BLAS/Lapack library")
#  ENDIF()
ELSE()
  FIND_PACKAGE(LAPACK QUIET)
  IF (LAPACK_FOUND)
    LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
    MESSAGE(STATUS "Using BLAS/Lapack library")
  ELSE()
    MESSAGE(WARNING " Lapack not found for GNU compiler!")
  ENDIF()
ENDIF()

# =========================================================================
# NETCDF (locally installed only!)
# =========================================================================
IF(LINK_GVEC_TO_NETCDF)
  IF((ON_SHAREDRUNNER GREATER -1) OR (ON_MPCDFCIRUNNER GREATER -1) OR (ON_TESTIMONY GREATER -1) OR (ON_COBRA GREATER -1) OR (ON_DRACO GREATER -1) OR (ON_RAVEN GREATER -1) OR (ON_TOKP GREATER -1))
    SET(NETCDF_DIR  $ENV{NETCDF_HOME})
  ELSEIF(ON_MAC_BREW GREATER -1)
    SET(NETCDF_DIR $ENV{HOMEBREW_PREFIX})
  ELSEIF(ON_MAC_PORTS GREATER -1)
    SET(NETCDF_DIR  /opt/local)
  ELSE()
    SET(NETCDF_DIR  /usr)
  ENDIF()
  INCLUDE_DIRECTORIES(${NETCDF_DIR}/include)
  LIST(APPEND linkedlibs -L${NETCDF_DIR}/lib -lnetcdf -lnetcdff -Wl,-rpath,${NETCDF_DIR}/lib )
  MESSAGE(STATUS "Build with NetCDF")
  ADD_DEFINITIONS(-DNETCDF=1)
ELSE()
  ADD_DEFINITIONS(-DNETCDF=0)
ENDIF()

# =========================================================================
# Ftimings
# =========================================================================
IF(USE_FTIMINGS)
  ADD_DEFINITIONS(-DPP_FTIMINGS)
  IF((ON_SHAREDRUNNER GREATER -1) OR (ON_RAVEN GREATER -1) OR (ON_COBRA GREATER -1) OR (ON_DRACO GREATER -1))
    LIST(APPEND linkedlibs -L$ENV{FTIMINGS_HOME}/lib -l$ENV{FTIMINGS_PKG}  -Wl,-rpath,$ENV{FTIMINGS_HOME}/lib ) 
    INCLUDE_DIRECTORIES($ENV{FTIMINGS_HOME}/include/$ENV{FTIMINGS_PKG})
  ENDIF()
    #TTR: don't know why above list append fails for MacOS... needs line below to work
  IF((ON_MAC_BREW GREATER -1) OR(ON_MAC_PORTS GREATER -1))
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L$ENV{FTIMINGS_HOME}/lib -l$ENV{FTIMINGS_PKG}  -Wl,-rpath,$ENV{FTIMINGS_HOME}/lib")
  ENDIF()
ENDIF()

# =========================================================================
# python 
# =========================================================================
find_package(Python COMPONENTS Interpreter REQUIRED)

# =========================================================================
# project sources & targets
# =========================================================================

function(add_exec target)
  ADD_EXECUTABLE(${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} gveclib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${MAIN_COMPILE_FLAGS})
  SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
endfunction()

# include defines.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/base/bsplines ${CMAKE_Fortran_MODULE_DIRECTORY})

FILE(GLOB_RECURSE mainF90 ./src/globals/*.f90
                          ./src/base/*.f90
                          ./src/functionals/mhd3d/*hmap*.f90
                          ./src/functionals/mhd3d/lambda_solve.f90
                          ./src/readstate/*.f90
                          ./src/analyze/analyze_vars.f90
                          ./src/io_netcdf.f90
                          ./src/output/*.f90
                          ./src/analyze/analyze_vars.f90
                          ./src/transform_sfl/*.f90
                          ./src/grid/*.f90 )

# define libs
IF(COMPILE_GVEC)
  FILE(GLOB_RECURSE gvecF90 ./src/*.f
                            ./src/functionals/*.f90
                            ./src/analyze/*.f90
                            ./src/vmec/*.f90
                            ./src/restart/*.f90 )

  IF(USE_FTIMINGS)
    LIST(APPEND mainF90 ./src/mod_timings.f90 ./src/perf2ftimings.f90)
  ENDIF()
  LIST(APPEND mainF90 ${gvecF90})
ENDIF()
IF(COMPILE_GVEC_TO_HOPR)
  FILE(GLOB tohoprF90 ./src/gvec_to_hopr/*.f90 )
  LIST(APPEND mainF90 ${tohoprF90})
ENDIF()
IF(COMPILE_GVEC_TO_GENE)
  FILE(GLOB togeneF90 ./src/gvec_to_gene/*.f90 )
  LIST(APPEND mainF90 ${togeneF90})
ENDIF()
IF(COMPILE_GVEC_TO_CASTOR3D)
  FILE(GLOB tocastorF90 ./src/gvec_to_castor3d/*.f90 )
  LIST(APPEND mainF90 ${tocastorF90})
ENDIF()
IF(COMPILE_GVEC_TO_JOREK)
  FILE(GLOB tojorekF90 ./src/gvec_to_jorek/*.f90 )
  LIST(APPEND mainF90 ${tojorekF90})
ENDIF()

ADD_LIBRARY(gveclibF90 OBJECT ${mainF90} )

OPTION(COMPILE_GVEC_AS_STATIC_LIB "compile GVEC libraray as static (OFF: compiled as shared)" ON)
MARK_AS_ADVANCED(COMPILE_GVEC_AS_STATIC_LIB)
SET_TARGET_PROPERTIES(gveclibF90 PROPERTIES COMPILE_FLAGS ${MAIN_COMPILE_FLAGS})
IF(COMPILE_GVEC_AS_STATIC_LIB)
  # generates static lib libgveclib.a
  ADD_LIBRARY(gveclib STATIC $<TARGET_OBJECTS:gveclibF90> )
ELSE()
  # generates shared lib libgveclib.so
  ADD_LIBRARY(gveclib SHARED $<TARGET_OBJECTS:gveclibF90> )
ENDIF()

  SET_TARGET_PROPERTIES(gveclib
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME "GVEC"
    DEBUG_POSTFIX "_d"
    MACOSX_RPATH ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
# Enable testing in general
include(CTest) 
#ENABLE_TESTING()

# =========================================================================
# Executable programs (simply use add_exec)
# =========================================================================
IF(COMPILE_GVEC)
  add_exec(gvec ./src/gvec.f90)
  ADD_DEPENDENCIES(gvec gveclib )
  add_exec(gvec_post ./src/gvec_post.f90)
  ADD_DEPENDENCIES(gvec_post gveclib )
  add_test(
    NAME gvec_pytest_create_rundir
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ctest_run 
  )
  add_test(
    NAME gvec_pytest_create_postdir
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ctest_post 
  )
  add_test(
    NAME gvec_pytest_dryrun
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_dryrun.txt --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --refdir=${CMAKE_BINARY_DIR}/ctest_run  --dry-run -m "run_stage or post_stage or regression_stage"
  )
  add_test(
    NAME gvec_pytest_run_shortrun_norestart
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "run_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run
  )
  add_test(
    NAME gvec_pytest_run_shortrun_restart
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_restart.txt -m "run_stage and shortrun and (restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run
  )
  add_test(
    NAME gvec_pytest_post_shortrun
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_post_shortrun.txt -m "post_stage and shortrun" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --postdir=${CMAKE_BINARY_DIR}/ctest_post
  )
  add_test(
    NAME gvec_pytest_run_debugrun
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_debugrun.txt -m "run_stage and debugrun" -k " (not highres) and (not _restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run
  )
ENDIF()
IF(COMPILE_GVEC_TO_HOPR)
  add_exec(test_gvec_to_hopr ./src/test_gvec_to_hopr.f90 )
  ADD_DEPENDENCIES(test_gvec_to_hopr gveclib )
  IF(COMPILE_GVEC)
    add_test(
      NAME gvec_to_hopr_pytest_create_convdir
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ctest_conv/to_hopr 
    )
    add_test(
      NAME gvec_to_hopr_pytest_shortrun_norestart
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --convdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_hopr"
    )
    set_tests_properties(gvec_to_hopr_pytest_shortrun_norestart PROPERTIES DEPENDS "gvec_pytest_shortrun_norestart")
  ENDIF()
ENDIF()
IF(COMPILE_GVEC_TO_GENE)
  add_exec(test_gvec_to_gene ./src/test_gvec_to_gene.f90 )
  ADD_DEPENDENCIES(test_gvec_to_gene gveclib )
  IF(COMPILE_GVEC)
    add_test(
      NAME gvec_to_gene_pytest_create_convdir
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ctest_conv/to_gene 
    )
    add_test(
      NAME gvec_to_gene_pytest_shortrun_norestart
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --convdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_gene"
    )
    set_tests_properties(gvec_to_gene_pytest_shortrun_norestart PROPERTIES DEPENDS "gvec_pytest_shortrun_norestart")
  ENDIF()
ENDIF()
IF(COMPILE_GVEC_TO_CASTOR3D)
  add_exec(convert_gvec_to_castor3d ./src/convert_gvec_to_castor3d.f90 )
  ADD_DEPENDENCIES(convert_gvec_to_castor3d gveclib )
  IF(COMPILE_GVEC)
    add_test(
      NAME gvec_to_castor3d_pytest_create_convdir
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ctest_conv/to_castor3d
    )
    add_test(
      NAME gvec_to_castor3d_pytest_shortrun_norestart
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --convdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_castor3d"
    )
    set_tests_properties(gvec_to_castor3d_pytest_shortrun_norestart PROPERTIES DEPENDS "gvec_pytest_shortrun_norestart")
  ENDIF()
ENDIF()
IF(COMPILE_GVEC_TO_JOREK)
  add_exec(convert_gvec_to_jorek ./src/convert_gvec_to_jorek.f90 )
  ADD_DEPENDENCIES(convert_gvec_to_jorek gveclib )
  IF(COMPILE_GVEC)
    add_test(
      NAME gvec_to_jorek_pytest_create_convdir
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ctest_conv/to_jorek
    )
    add_test(
      NAME gvec_to_jorek_pytest_shortrun_norestart
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --convdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_jorek"
    )
    set_tests_properties(gvec_to_jorek_pytest_shortrun_norestart PROPERTIES DEPENDS "gvec_pytest_shortrun_norestart")
  ENDIF()
ENDIF()

# =========================================================================
# Install
# =========================================================================

IF(NOT COMPILE_GVEC_TO_GENE)
  INSTALL(FILES
    ${CMAKE_Fortran_MODULE_DIRECTORY}/configuration-cmake.f90
    DESTINATION include )
ELSE()
  INCLUDE(cmake/install_gene.cmake)
ENDIF()


# =========================================================================
# USERBLOCK + Preproc_flags
# =========================================================================
# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    GET_CMAKE_PROPERTY(_vars CACHE_VARIABLES)
    STRING(REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    SET(${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# A function to get all user[A-Za-z0-9_]* defined variables with a specified string
function (getListOfVarsWith _str _varResult)
    GET_CMAKE_PROPERTY(_vars CACHE_VARIABLES)
    STRING(REGEX MATCHALL ";([A-Za-z0-9_]*${_str}[A-Za-z0-9_]*)" _matchedVars "${_vars}")
    SET(${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# Get the latest abbreviated commit hash of the working branch
 execute_process(
     COMMAND git log -n 1 --format=%H
     WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
     OUTPUT_VARIABLE GIT_HASH
     OUTPUT_STRIP_TRAILING_WHITESPACE
     )
# Get the latest abbreviated branch name
 execute_process(
     COMMAND git rev-parse --abbrev-ref HEAD
     WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
     OUTPUT_VARIABLE GIT_BRANCH
     OUTPUT_STRIP_TRAILING_WHITESPACE
     )

# A function to get all user defined variables with a specified prefix
SET(configuration ${CMAKE_Fortran_MODULE_DIRECTORY}/configuration-cmake.f90)
FILE(WRITE ${configuration} "!This file is included during compilation in gvec.f90 to show the configuration on the screen output. \n")
FILE(APPEND ${configuration} "  SWRITE(UNIT_stdOut,'(A,3(1X,A))')  \"GIT_BRANCH:\",\"'${GIT_BRANCH}'\",\",GIT_HASH:\",\"'${GIT_HASH}'\"\n")
FILE(APPEND ${configuration} "  SWRITE(UNIT_stdOut,'(A)') \" THIS IS AN OUTPUT OF THE GVEC OPTIONS USED IN CMAKE:\" \n")
GETLISTOFVARSWITH("GVEC" gvecVars)
GETLISTOFVARSSTARTINGWITH("USE_" useVars)
GETLISTOFVARSSTARTINGWITH("CMAKE_BUILD_TYPE" cbVars)
GETLISTOFVARSSTARTINGWITH("CMAKE_HOSTNAME" hnVars)
FOREACH (_var IN LISTS cbVars hnVars gvecVars useVars )
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
#  SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF")) 
#  IF(${boolian})
#    #FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
#    FILE(APPEND ${configuration} "  WRITE(UNIT_stdOut,*)  \"   CMAKE,OPTION: ${_var} ${${_var}} \" \n")
#  ELSE()
    IF(${_var})
      STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
    ENDIF()
    #FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
    FILE(APPEND ${configuration} "  SWRITE(UNIT_stdOut,'(A18,A30,A3,A40)')  \"   CMAKE-OPTION : \",\" ${_var} \",\" : \",\"${${_var}}\"\n")
#  ENDIF()
ENDFOREACH()

