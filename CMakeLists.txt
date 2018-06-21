# =========================================================================
# GVEC CMake configuration
# =========================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
# TODO: eventually require CMake 3.0 in future releases
# Advantages: Possible to split string
IF ("${CMAKE_VERSION}" VERSION_GREATER 3.0.0)
  cmake_policy (SET CMP0045 OLD)
  cmake_policy (SET CMP0046 OLD)
ENDIF()
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)



# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "marconi" ON_MARCONI)
STRING(FIND "${CMAKE_HOSTNAME}" "testimony" ON_HYDRA)
STRING(FIND "${CMAKE_HOSTNAME}" "hydra"     ON_HYDRA)
IF (ON_MARCONI GREATER -1)
  MESSAGE(STATUS "Compiling on Marconi")
ENDIF()
IF (ON_HYDRA GREATER -1)
  MESSAGE(STATUS "Compiling on hydra")
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

OPTION(BUILD_NETCDF "Link to local NetCDF library" ON)

# =========================================================================
# PREPROCESSOR FLAGS
# =========================================================================

#SET(MAIN_DISCTYPE "1" CACHE STRING "discretization type")
#SET_PROPERTY(CACHE MAIN_DISCTYPE  PROPERTY STRINGS "1" "2")
#ADD_DEFINITIONS(-DPP_DiscType=${MAIN_DISCTYPE})



# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug / Release / Profile."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS("-DDEBUG")
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
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${MAIN_INSTRUCTION} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${MAIN_INSTRUCTION} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -ffpe-trap=invalid -fbounds-check -finit-real=snan -fbacktrace -Wall -Wextra")
  # add flags only for compiling not linking!
  SET (MAIN_COMPILE_FLAGS "-xf95-cpp-input")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   #environment specific variables
   IF(ON_MARCONI GREATER -1)
      SET(MAIN_INSTRUCTION "-axMIC-AVX512 ") #for BW and KNL
   ELSE()
      SET(MAIN_INSTRUCTION "-xHost")
   ENDIF()
   # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${MAIN_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${MAIN_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
   # add flags only for compiling not linking!
   SET (MAIN_COMPILE_FLAGS "-fpp -assume bscc")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -f free -s real64 -em -DCRAY")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (MAIN_COMPILE_FLAGS "-F")

ELSE()
   MESSAGE(ERROR " Unknown compiler")
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
# use Lapack/Blas for GNU
MESSAGE(STATUS "Using BLAS/Lapack library")
IF(ON_HYDRA GREATER -1)
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
    SET(MKL_DIR $ENV{MKL_HOME})
    LIST(APPEND linkedlibs -mkl -Wl,-rpath=${MKL_DIR}/lib/intel64 )
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    FIND_PACKAGE(LAPACK QUIET)
    IF (LAPACK_FOUND)
      LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
    ELSE()
      MESSAGE(WARNING " Lapack not found for GNU compiler!")
    ENDIF()
  ELSE()
     MESSAGE(ERROR " Unknown compiler ...")
  ENDIF()
ELSE()
  FIND_PACKAGE(LAPACK QUIET)
  IF (LAPACK_FOUND)
    LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
  ELSE()
    MESSAGE(WARNING " Lapack not found for GNU compiler!")
  ENDIF()
ENDIF()

# =========================================================================
# NETCDF (locally installed only!)
# =========================================================================
IF(BUILD_NETCDF)
  IF(ON_HYDRA GREATER -1)
    SET(NETCDF_DIR  $ENV{NETCDF_HOME}) 
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

FILE(GLOB_RECURSE mainF90 ./src/*.f90 ./src/*.f )


# define libs
ADD_LIBRARY(gveclibF90 OBJECT ${mainF90} )

SET_TARGET_PROPERTIES(gveclibF90 PROPERTIES COMPILE_FLAGS ${MAIN_COMPILE_FLAGS}) 
ADD_LIBRARY(gveclib STATIC $<TARGET_OBJECTS:gveclibF90> )

# =========================================================================
# Executable programs (simply use add_exec)
# =========================================================================
add_exec(gvec ./src/gvec.f90)
ADD_DEPENDENCIES(gvec gveclib )

add_exec(test_gvec_to_hopr ./src/test_gvec_to_hopr.f90 )
ADD_DEPENDENCIES(test_gvec_to_hopr gveclib )

# =========================================================================
# Install
# =========================================================================

INSTALL(FILES    ${CMAKE_Fortran_MODULE_DIRECTORY}/configuration-cmake.f90 DESTINATION include )
