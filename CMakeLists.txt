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
STRING(FIND "${CMAKE_HOSTNAME}" "cobra"     ON_COBRA)
IF (ON_MARCONI GREATER -1)
  MESSAGE(STATUS "Compiling on Marconi")
ENDIF()
IF (ON_HYDRA GREATER -1)
  MESSAGE(STATUS "Compiling on hydra")
ENDIF()
IF (ON_COBRA GREATER -1)
  MESSAGE(STATUS "Compiling on cobra")
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

#SET(MAIN_DISCTYPE "1" CACHE STRING "discretization type")
#SET_PROPERTY(CACHE MAIN_DISCTYPE  PROPERTY STRINGS "1" "2")
#ADD_DEFINITIONS(-DPP_DiscType=${MAIN_DISCTYPE})

OPTION(COMPILE_GVEC "compile gvec programme for equilibrium computations (libgvec and gvec executable)" ON)
IF(COMPILE_GVEC)
  OPTION(LINK_GVEC_TO_NETCDF "Link to local NetCDF library, needed for VMEC input" ON)
ELSE()
  UNSET(LINK_GVEC_TO_NETCDF CACHE ) # if not mhd
ENDIF()
OPTION(COMPILE_GVEC_TO_GENE "compile interface gvec_to_gene library (added to libgvec)" ON)
OPTION(COMPILE_GVEC_TO_HOPR "compile interface gvec_to_hopr library (added to libgvec)" ON)
OPTION(COMPILE_GVEC_TO_CASTOR3D "compile converter gvec_to_castor (only executable)" ON)


# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
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
  SET (MAIN_COMPILE_FLAGS "-xf95-cpp-input -fPIC")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   #environment specific variables
   IF(ON_MARCONI GREATER -1)
      SET(MAIN_INSTRUCTION "-axMIC-AVX512 ") #for BW and KNL
   ELSEIF (ON_HYDRA GREATER -1 )
      SET(MAIN_INSTRUCTION "-xCORE-AVX-I ") #for hydra
      # add -parallel for automatic parallelization (choose OMP_NUM_THREADS=16)
   ELSE()
      SET(MAIN_INSTRUCTION "-xHost")
   ENDIF()
   # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${MAIN_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
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
ELSEIF(ON_COBRA GREATER -1)
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
    SET(MKL_DIR $ENV{MKL_HOME})
    LIST(APPEND linkedlibs -mkl -Wl,-rpath=${MKL_DIR}/lib/intel64 )
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
    SET(MKL_DIR $ENV{MKL_HOME})
    LIST(APPEND linkedlibs -L${MKL_DIR}/lib/intel64 -lmkl_gf_lp64
    		-lmkl_sequential -lmkl_core -lpthread -lm -ldl
 		-Wl,-rpath=${MKL_DIR}/lib/intel64 )
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
IF(LINK_GVEC_TO_NETCDF)
  IF(ON_HYDRA GREATER -1)
    SET(NETCDF_DIR  $ENV{NETCDF_HOME})
  ELSEIF(ON_COBRA GREATER -1)
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

FILE(GLOB_RECURSE mainF90 ./src/globals/*.f90
                          ./src/base/*.f90
                          ./src/functionals/mhd3d/*hmap*.f90
                          ./src/readstate*.f90
                          ./src/output/*.f90
                          ./src/grid/*.f90 )

# define libs
IF(COMPILE_GVEC)
  FILE(GLOB_RECURSE gvecF90 ./src/*.f
                            ./src/functionals/*.f90
                            ./src/analyze/*.f90
                            ./src/lbfgsb/*.f90
                            ./src/vmec/*.f90
                            ./src/output/*.f90
                            ./src/restart/*.f90 )
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
  FILE(GLOB tohoprF90 ./src/gvec_to_castor3d/*.f90 )
  LIST(APPEND mainF90 ${tohoprF90})
ENDIF()

ADD_LIBRARY(gveclibF90 OBJECT ${mainF90} )

SET_TARGET_PROPERTIES(gveclibF90 PROPERTIES COMPILE_FLAGS ${MAIN_COMPILE_FLAGS})
# generates shared lib libgveclib.so
ADD_LIBRARY(gveclib SHARED $<TARGET_OBJECTS:gveclibF90> )
# generates static lib libgveclib.a
#ADD_LIBRARY(gveclib STATIC $<TARGET_OBJECTS:gveclibF90> )


# =========================================================================
# Executable programs (simply use add_exec)
# =========================================================================
IF(COMPILE_GVEC)
  add_exec(gvec ./src/gvec.f90)
  ADD_DEPENDENCIES(gvec gveclib )
ENDIF()
IF(COMPILE_GVEC_TO_HOPR)
  add_exec(test_gvec_to_hopr ./src/test_gvec_to_hopr.f90 )
  ADD_DEPENDENCIES(test_gvec_to_hopr gveclib )
ENDIF()
IF(COMPILE_GVEC_TO_GENE)
  add_exec(test_gvec_to_gene ./src/test_gvec_to_gene.f90 )
  ADD_DEPENDENCIES(test_gvec_to_gene gveclib )
ENDIF()
IF(COMPILE_GVEC_TO_CASTOR3D)
  add_exec(convert_gvec_to_castor3d ./src/convert_gvec_to_castor3d.f90 )
  ADD_DEPENDENCIES(convert_gvec_to_castor3d gveclib )
ENDIF()


# =========================================================================
# Install
# =========================================================================

if(NOT COMPILE_GVEC_TO_GENE)

  INSTALL(FILES
    ${CMAKE_Fortran_MODULE_DIRECTORY}/configuration-cmake.f90
    DESTINATION include )

else()

  # Preparations to install GVEC for GENE

  message(STATUS "${PROJECT_NAME} will be installed to ${CMAKE_INSTALL_PREFIX}")

  # standartize install output directories
  include(GNUInstallDirs)

  # Offer the user the choice of overriding the installation directories
  set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH
    "Installation directory for libraries")
  set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH
    "Installation directory for executables")
  set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH
    "Installation directory for header files")
  if(WIN32 AND NOT CYGWIN)
    set(DEF_INSTALL_CMAKEDIR CMake)
  else()
    set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
  endif()
  set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH
    "Installation directory for CMake files")
  set(INSTALL_MODDIR "mod" CACHE PATH
    "Installation directory for Fortran .mod files")

  # report installation directories
  foreach(p LIB BIN INCLUDE CMAKE MOD)
    file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}DIR} _path )
    message(STATUS "Installing ${p} components to ${_path}")
    unset(_path)
  endforeach()

  target_include_directories(gveclib
    INTERFACE
    $<INSTALL_INTERFACE:${INSTALL_MODDIR}>
    )

  set_target_properties(gveclib
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME "GVEC"
    DEBUG_POSTFIX "_d"
    MACOSX_RPATH ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    )

  get_property(ARE_SHARED_LIBS_SUPPORTED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
  if(ARE_SHARED_LIBS_SUPPORTED)
    # Prepare RPATH
    file(RELATIVE_PATH _rel
      ${CMAKE_INSTALL_PREFIX}/${INSTALL_BINDIR}
      ${CMAKE_INSTALL_PREFIX})
    if(APPLE)
      set(_rpath "@loader_path/${_rel}")
    else()
      set(_rpath "\$ORIGIN/${_rel}")
    endif()
    file(TO_NATIVE_PATH "${_rpath}/${INSTALL_LIBDIR}" test_gvec_to_gene_RPATH)

    set_target_properties(test_gvec_to_gene
      PROPERTIES
      OUTPUT_NAME "test_gvec_to_gene"
      MACOSX_RPATH ON
      SKIP_BUILD_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_RPATH "${test_gvec_to_gene_RPATH}"
      INSTALL_RPATH_USE_LINK_PATH ON
      )
  endif()

  # generate gvec_to_gene_export.h
  include(GenerateExportHeader)
  generate_export_header( gveclib
    BASE_NAME gvec_to_gene)

  if(ENABLE_PYTHON_BINDINGS)
    set(_GVEC_TO_GENE_HEADERS
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/gvec_to_gene.h
      ${CMAKE_CURRENT_BINARY_DIR}/gvec_to_gene_export.h
      )

    set(_GVEC_TO_GENE_RESOURCE
      "${CMAKE_CURRENT_BINARY_DIR}/interface_file_names.cfg"
      )

    set(_GVEC_TO_GENE_PYTHON_INSTALL
      gvec_to_gene
      )

    file(RELATIVE_PATH _GVEC_TO_GENE_PYTHON_H_REL
      ${CMAKE_INSTALL_PREFIX}/${_GVEC_TO_GENE_PYTHON_INSTALL}
      ${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDEDIR}
      )

    file(RELATIVE_PATH _GVEC_TO_GENE_PYTHON_L_REL
      ${CMAKE_INSTALL_PREFIX}/${_GVEC_TO_GENE_PYTHON_INSTALL}
      ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBDIR}
      )

    file(
      GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/interface_file_names.cfg
      #INPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/interface_file_names.cfg.in
      CONTENT
      "[configuration]
    header_file_name = ${_GVEC_TO_GENE_PYTHON_H_REL}/gvec_to_gene.h
    library_file_name = ${_GVEC_TO_GENE_PYTHON_L_REL}/$<TARGET_FILE_NAME:gveclib>"
      )

    # Install python package
    install(
      FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/__init__.py
      DESTINATION
      ${_GVEC_TO_GENE_PYTHON_INSTALL}
      )

    # Add additional properties to gveclib to enable Python bindings
    set_target_properties( gveclib
      PROPERTIES
      PUBLIC_HEADER
      "${_GVEC_TO_GENE_HEADERS}"
      RESOURCE
      "${_GVEC_TO_GENE_RESOURCE}"
      )
  endif()

  # install gveclib and test_gvec_to_gene
  install(
    TARGETS
    test_gvec_to_gene
    gveclib
    EXPORT
    GVECTargets
    ARCHIVE
    DESTINATION ${INSTALL_LIBDIR}
    COMPONENT lib
    RUNTIME
    DESTINATION ${INSTALL_BINDIR}
    COMPONENT bin
    LIBRARY
    DESTINATION ${INSTALL_LIBDIR}
    COMPONENT lib
    PUBLIC_HEADER
    DESTINATION ${INSTALL_INCLUDEDIR}
    COMPONENT dev
    RESOURCE
    DESTINATION ${_GVEC_TO_GENE_PYTHON_INSTALL}
    )

  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${INSTALL_MODDIR})

  # Exporting targets for other CMake-projects

  install(
    EXPORT
    GVECTargets
    NAMESPACE
    "GVEC::"
    DESTINATION
    ${INSTALL_CMAKEDIR}
    COMPONENT
    dev
    )

  include(CMakePackageConfigHelpers)

  configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/GVECConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GVECConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CMAKEDIR}
    )

  install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/GVECConfig.cmake
    DESTINATION
    ${INSTALL_CMAKEDIR}
    )

  # turn on testing
  enable_testing()

  if(ENABLE_PYTHON_BINDINGS)
    # Testing through Python

    # require python
    find_package(PythonInterp REQUIRED)

    # define test
    add_test(
      NAME
      gvec_to_gene_python_test
      COMMAND
      ${CMAKE_COMMAND} -E env GVEC_TO_GENE_MODULE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/src
      GVEC_TO_GENE_HEADER_FILE=${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/gvec_to_gene.h
      GVEC_TO_GENE_LIBRARY_FILE=$<TARGET_FILE:gveclib>
      ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/test.py
      )
  endif()

endif()
