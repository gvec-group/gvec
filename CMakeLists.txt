# =========================================================================
# GVEC CMake configuration
# =========================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 3.22)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)

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

IF(DEFINED ENV{CONDA_PREFIX})
  MESSAGE(STATUS "Detected conda environment: $ENV{CONDA_PREFIX}")
  SET(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX} CACHE PATH "Installation prefixes to be searched for dependencies")
ENDIF()

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

# =========================================================================
# CMAKE OPTIONS
# =========================================================================

OPTION(COMPILE_GVEC "compile gvec programme for equilibrium computations (libgvec and gvec executable)" ON)
OPTION(LINK_GVEC_TO_NETCDF "Link to local serial NetCDF library, needed for VMEC input/output" ON)
OPTION(COMPILE_GVEC_TO_GENE "compile interface gvec_to_gene library (added to libgvec)" OFF)
OPTION(COMPILE_GVEC_TO_HOPR "compile interface gvec_to_hopr library (added to libgvec)" OFF)
OPTION(COMPILE_GVEC_TO_CASTOR3D "compile converter gvec_to_castor (only executable)" OFF)
OPTION(COMPILE_GVEC_TO_JOREK "compile converter gvec_to_jorek (only executable)" OFF)
OPTION(COMPILE_PYGVEC "compile python interface to gvec" OFF)
MARK_AS_ADVANCED(COMPILE_PYGVEC)
OPTION(COMPILE_SRC_PYGVEC " compile pygvec fortran sources" ON)
MARK_AS_ADVANCED(COMPILE_SRC_PYGVEC)
IF(NOT COMPILE_GVEC)
  UNSET(COMPILE_PYGVEC CACHE )
  UNSET(COMPILE_SRC_PYGVEC CACHE )
ENDIF()



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

SET (GVEC_FIX_HMAP "off" CACHE STRING "fix hmap at compile time, default is off")
SET_PROPERTY(CACHE GVEC_FIX_HMAP PROPERTY STRINGS off RZ axisNB frenet cyl knot )
IF(GVEC_FIX_HMAP STREQUAL "RZ")
  ADD_DEFINITIONS(-DPP_WHICH_HMAP=1)
ELSEIF(GVEC_FIX_HMAP STREQUAL "cyl")
  ADD_DEFINITIONS(-DPP_WHICH_HMAP=3)
ELSEIF(GVEC_FIX_HMAP STREQUAL "knot")
  ADD_DEFINITIONS(-DPP_WHICH_HMAP=10)
ELSEIF(GVEC_FIX_HMAP STREQUAL "frenet")
  ADD_DEFINITIONS(-DPP_WHICH_HMAP=20)
ELSEIF(GVEC_FIX_HMAP STREQUAL "axisNB")
  ADD_DEFINITIONS(-DPP_WHICH_HMAP=21)
ENDIF()

# =========================================================================
# Build type
# -------------------------------------------------------------------------
# make sure that the default is RELEASE
# =========================================================================

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug / Release / Profile."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
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
  MARK_AS_ADVANCED(MPI_LIBRARY MPI_EXTRA_LIBRARY)
  MESSAGE(STATUS "Building GVEC with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
  ADD_DEFINITIONS(-DMPI=1)
  OPTION(NO_MPI_F08 "disable MPI F2008" OFF)
  MARK_AS_ADVANCED(NO_MPI_F08)
  IF(NO_MPI_F08)
    ADD_DEFINITIONS(-DNO_MPI_F08)
  ENDIF()
ELSE() # if not MPI
  ADD_DEFINITIONS(-DMPI=0)
  UNSET(NO_MPI_F08 CACHE )
ENDIF()

# =========================================================================
# COMPILER FLAGS (compiler & machine specific)
# =========================================================================

IF(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 -march=native -finline-functions") #gives segfault: -fstack-arrays
  SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 -march=native -finline-functions") #gives segfault: -fstack-arrays
  SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -ffpe-trap=invalid -fbounds-check -finit-real=snan -fbacktrace -Wall -Wextra")
  SET(GVEC_COMPILE_FLAGS "-xf95-cpp-input -fPIC")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # machine specific flags
  IF(${CMAKE_HOSTNAME} MATCHES "mpcdfcirunner")
    SET(GVEC_ADD_COMPILE_FLAGS_MACHINE "-xCORE-AVX-I" CACHE STRING "additional machine dependent compile flags (for Intel)")
  ELSEIF(${CMAKE_HOSTNAME} MATCHES "raven")
    SET(GVEC_ADD_COMPILE_FLAGS_MACHINE "-xCORE-AVX512 -qopt-zmm-usage=high -align array64byte" CACHE STRING "additional machine dependent compile flags (for Intel)")
  ELSEIF(${CMAKE_HOSTNAME} MATCHES "marconi")
    SET(GVEC_ADD_COMPILE_FLAGS_MACHINE "-axMIC-AVX512" CACHE STRING "additional machine dependent compile flags (for Intel)") # for BW and KNL
  ELSE()
    SET(GVEC_ADD_COMPILE_FLAGS_MACHINE "-xHost" CACHE STRING "additional machine dependent compile flags (for Intel)")
  ENDIF()

  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 ${GVEC_ADD_COMPILE_FLAGS_MACHINE} -qopt-report=5 -diag-disable=10397 -no-prec-div")
  SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${GVEC_ADD_COMPILE_FLAGS_MACHINE} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check noarg_temp_created,noformat,nooutput_conversion,pointer -init=snan -init=arrays")
  IF(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "2021.2")
    SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all,uninit") # crash with ifx; possibly related to netCDF
  ENDIF()
  SET(GVEC_COMPILE_FLAGS "-fpp -assume bscc -fPIC")

ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -f free -s real64 -em -DCRAY")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -p . -rm")
  SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -h profile_generate -p . -rm")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -eD -rm")
  SET(GVEC_COMPILE_FLAGS "-F")

ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -fPIC")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3")
  SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O3")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -g")
  SET(GVEC_COMPILE_FLAGS "-Mpreprocess")

ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -fPIC -traceback -Wall -DNVHPC")
  SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3")
  SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0  -g ") # -C option fails for LEN(*) STRINGS
  SET(GVEC_COMPILE_FLAGS "-Mpreprocess")
  IF(USE_OPENMP)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp=multicore") # for now not using the gpu for OpenMP!
  ELSE()
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -nomp")
  ENDIF()

ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler: ${CMAKE_Fortran_COMPILER_ID}")
ENDIF()

# =========================================================================
# OpenMP
# =========================================================================

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED COMPONENTS Fortran)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
ENDIF()

# =========================================================================
# BLAS / LAPACK - MKL
# -------------------------------------------------------------------------
# rely on CMake to find BLAS / LAPACK
# on MPCDF clusters loading the module selects the correct vendor
# you can manually set BLA_VENDOR
# =========================================================================

FIND_PACKAGE(LAPACK REQUIRED)
LIST(APPEND linkedlibs LAPACK::LAPACK)

# =========================================================================
# netCDF
# -------------------------------------------------------------------------
# netCDF-fortran is not commonly available as a CMake target
# but a pkg-config is usually (ubuntu, MPCDF) available
# =========================================================================

IF(LINK_GVEC_TO_NETCDF)
  FIND_PACKAGE(PkgConfig REQUIRED)
  PKG_CHECK_MODULES(NETCDF REQUIRED IMPORTED_TARGET netcdf-fortran)
  PKG_GET_VARIABLE(NETCDF_FMODDIR netcdf-fortran fmoddir)
  INCLUDE_DIRECTORIES(${NETCDF_INCLUDEDIR} ${NETCDF_FMODDIR})
  LIST(APPEND linkedlibs PkgConfig::NETCDF)
  ADD_DEFINITIONS(-DNETCDF=1)
ENDIF()

# =========================================================================
# ftimings
# -------------------------------------------------------------------------
# Homepage/Repo: https://gitlab.mpcdf.mpg.de/loh/ftimings
# =========================================================================

IF(USE_FTIMINGS)
  FIND_PACKAGE(PkgConfig REQUIRED)
  IF(NOT DEFINED ENV{FTIMINGS_PKG})
    MESSAGE(SEND_ERROR "FTIMINGS_PKG is not set.")
  ENDIF()
  PKG_CHECK_MODULES(FTIMINGS REQUIRED IMPORTED_TARGET $ENV{FTIMINGS_PKG})
  INCLUDE_DIRECTORIES("${FTIMINGS_INCLUDEDIR}/$ENV{FTIMINGS_PKG}")
  LIST(APPEND linkedlibs PkgConfig::FTIMINGS)
  ADD_DEFINITIONS(-DPP_FTIMINGS)
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
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${GVEC_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran
                                             INSTALL_RPATH_USE_LINK_PATH TRUE) # don't strip RUNPATH during install
  IF(COMPILE_PYGVEC) # scikit-build-core = python/pip installation - install to venv/bin
    INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION ${SKBUILD_SCRIPTS_DIR})
  ELSE()
    INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  ENDIF()
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
endfunction()

# include defines.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/base/bsplines ${CMAKE_Fortran_MODULE_DIRECTORY})

FILE(GLOB_RECURSE mainF90 ./src/globals/*.f90
                          ./src/base/*.f90
                          ./src/profiles/*.f90
                          ./src/hmap/c_hmap.f90
                          ./src/hmap/hmap.f90
                          ./src/functionals/mhd3d/lambda_solve.f90
                          ./src/readstate/*.f90
                          ./src/analyze/analyze_vars.f90
                          ./src/io_netcdf.f90
                          ./src/output/*.f90
                          ./src/transform_sfl/*.f90
                          ./src/grid/*.f90 )
IF(GVEC_FIX_HMAP STREQUAL "off")
  FILE(GLOB hmapsF90 ./src/hmap/hmap_*.f90 )
ELSE()
  FILE(GLOB hmapsF90 ./src/hmap/hmap_${GVEC_FIX_HMAP}.f90)
ENDIF()
LIST(APPEND mainF90 ${hmapsF90})
# define libs
IF(COMPILE_GVEC)
  FILE(GLOB_RECURSE gvecF90 ./src/*.f
                            ./src/rungvec.f90
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
IF(COMPILE_SRC_PYGVEC)
  FILE(GLOB pygvecF90 ./src/pygvec/*.f90 )
  LIST(APPEND mainF90 ${pygvecF90})
ENDIF()

# append relative filename-macro for __FILENAME__ in stamp of abort function (see defines.h)
FOREACH(sourcefile ${mainF90})
    # inspired by https://stackoverflow.com/a/27990434
    CMAKE_PATH(NORMAL_PATH sourcefile OUTPUT_VARIABLE filename)
    CMAKE_PATH(RELATIVE_PATH filename)
    SET_PROPERTY(
        SOURCE "${sourcefile}" APPEND
        PROPERTY COMPILE_DEFINITIONS "__FILENAME__=\"${filename}\"")
ENDFOREACH()

ADD_LIBRARY(gveclibF90 OBJECT ${mainF90} )

OPTION(COMPILE_GVEC_AS_STATIC_LIB "compile GVEC library as static (OFF: compiled as shared)" ON)
MARK_AS_ADVANCED(COMPILE_GVEC_AS_STATIC_LIB)
SET_TARGET_PROPERTIES(gveclibF90 PROPERTIES COMPILE_FLAGS ${GVEC_COMPILE_FLAGS})
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
# Python bindings (pygvec), using f90wrap and f2py
# =========================================================================

IF(COMPILE_PYGVEC)
  FIND_PACKAGE(Python COMPONENTS Interpreter Development NumPy REQUIRED)

  IF(PYTHON_FOUND)
    MESSAGE(STATUS "Python version: ${Python_VERSION}")
    MESSAGE(STATUS "Python executable: ${PYTHON_EXECUTABLE}")
    MESSAGE(STATUS "Python include dir: ${Python_INCLUDE_DIRS}")
    MESSAGE(STATUS "Python libraries: ${Python_LIBRARIES}")
    MESSAGE(STATUS "Python NumPy include dir: ${Python_NumPy_INCLUDE_DIRS}")
  ELSE()
    MESSAGE(SEND_ERROR "Python not found, required for 'COMPILE_PYGVEC=ON'!")
  ENDIF()

  IF(NOT ${SKBUILD})
    MESSAGE(SEND_ERROR "pyGVEC expects scikit-build-core to be used!")
  ENDIF()

  SET(PYPKG_PATH "python/gvec")

  # --- f2py headers --- #
  SET(F2PY_INCLUDE_DIR "${Python_NumPy_INCLUDE_DIRS}/../../f2py/src")
  MESSAGE(STATUS "F2PY include dir: ${F2PY_INCLUDE_DIR}")

  # --- wrap gvec using f90wrap & f2py --- #
  # specify files to be wrapped with f90wrap
  SET(pygvec_src_files
    src/pygvec/binding.f90
    src/pygvec/state.f90
    src/pygvec/run.f90
    src/transform_sfl/sfl_boozer.f90
  )

  FOREACH(file ${pygvec_src_files})
    GET_FILENAME_COMPONENT(basename "${file}" NAME_WE)
    LIST(APPEND pygvec_wrap_files "f90wrap_${basename}.f90")
  ENDFOREACH()

  LIST(TRANSFORM pygvec_src_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
  LIST(TRANSFORM pygvec_wrap_files PREPEND "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_CUSTOM_COMMAND(
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    VERBATIM
    COMMAND "${PYTHON_EXECUTABLE}"
            -m f90wrap
            ${pygvec_src_files}
            --mod-name lib  # -> lib.py - autogenerated python wrappers
            --f90-mod-name gvec._libpygvec  # -> _libpygvec.cpython*.so - compiled python module
            --kind-map "${CMAKE_CURRENT_SOURCE_DIR}/python/kind_map.py"
            --class-names "${CMAKE_CURRENT_SOURCE_DIR}/python/class_names.py"
            --constructors
            --destructors sfl_boozer_free
    DEPENDS ${pygvec_src_files}
    OUTPUT ${pygvec_wrap_files}
    COMMENT "pyGVEC: generating the Fortran layer with f90wrap"
  )

  ADD_CUSTOM_TARGET(
    pygvec_wrap_files
    DEPENDS ${pygvec_wrap_files}
  )
  ADD_CUSTOM_COMMAND(
    DEPENDS ${pygvec_wrap_files}
    OUTPUT _libpygvecmodule.c
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    VERBATIM
    COMMAND f2py-f90wrap
            -m _libpygvec  # compiled python module
            --lower
            ${pygvec_wrap_files}
    COMMENT "pyGVEC: generating the C layer with f2py-f90wrap"
  )

  PYTHON_ADD_LIBRARY(_libpygvec MODULE WITH_SOABI
    _libpygvecmodule.c
    ${F2PY_INCLUDE_DIR}/fortranobject.c
    ${pygvec_wrap_files}
  )
  TARGET_INCLUDE_DIRECTORIES(_libpygvec PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src" ${F2PY_INCLUDE_DIR} ${Python_NumPy_INCLUDE_DIRS} ${Python_INCLUDE_DIRS})
  ADD_DEPENDENCIES(_libpygvec pygvec_wrap_files)
  TARGET_LINK_LIBRARIES(_libpygvec PRIVATE gveclib ${linkedlibs})
  SET_TARGET_PROPERTIES(_libpygvec PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE) # don't strip RUNPATH during install

  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/lib.py" DESTINATION "${SKBUILD_PROJECT_NAME}")
  INSTALL(TARGETS _libpygvec LIBRARY DESTINATION "${SKBUILD_PROJECT_NAME}")
ENDIF()

# =========================================================================
# Executable programs (simply use add_exec)
# =========================================================================
IF(COMPILE_GVEC)
  add_exec(gvec ./src/gvec.f90)
  ADD_DEPENDENCIES(gvec gveclib )
  add_exec(gvec_post ./src/gvec_post.f90)
  ADD_DEPENDENCIES(gvec_post gveclib )
  MESSAGE(STATUS "path to build directory: ${CMAKE_BINARY_DIR}")
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
    COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_dryrun.txt --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --dry-run -m "run_stage"
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
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --postdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_hopr"
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
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --postdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_gene"
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
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --postdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_castor3d"
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
      COMMAND ${Python_EXECUTABLE}  -m pytest -v -r A --log-file=${CMAKE_BINARY_DIR}/log_pytest_run_shortrun_norestart.txt -m "converter_stage and shortrun and (not restart)" --builddir=${CMAKE_BINARY_DIR} --rundir=${CMAKE_BINARY_DIR}/ctest_run --postdir=${CMAKE_BINARY_DIR}/ctest_conv -k "to_jorek"
    )
    set_tests_properties(gvec_to_jorek_pytest_shortrun_norestart PROPERTIES DEPENDS "gvec_pytest_shortrun_norestart")
  ENDIF()
ENDIF()

# =========================================================================
# Install
# =========================================================================

IF((NOT COMPILE_GVEC_TO_GENE) AND (NOT COMPILE_PYGVEC))
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
LIST(APPEND cmpVars "CMAKE_HOSTNAME")
LIST(APPEND cmpVars "CMAKE_BUILD_TYPE")
LIST(APPEND cmpVars "CMAKE_Fortran_COMPILER")
LIST(APPEND cmpVars "CMAKE_C_COMPILER")
LIST(APPEND cmpVars "CMAKE_CXX_COMPILER")
IF(COMPILE_PYGVEC)
  SET(configuration_py ${CMAKE_CURRENT_BINARY_DIR}/_compile_options.py)
  FILE(WRITE ${configuration_py} "# file generated by GVEC's CMake configuration\n")
  FILE(APPEND ${configuration_py} "# don't change, don't track in version control\n\n")
  FILE(APPEND ${configuration_py} "from pathlib import Path\n\n")
ENDIF()

FOREACH (_var IN LISTS cmpVars gvecVars useVars)
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
#  SET(boolean (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF"))
#  IF(${boolean})
#    #FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
#    FILE(APPEND ${configuration} "  WRITE(UNIT_stdOut,*)  \"   CMAKE,OPTION: ${_var} ${${_var}} \" \n")
#  ELSE()
    IF(${_var})
      STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
    ENDIF()
    #FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
  FILE(APPEND ${configuration} "  SWRITE(UNIT_stdOut,'(A18,A30,A3,A40)')  \"   CMAKE-OPTION : \",\" ${_var} \",\" : \",\"${${_var}}\"\n")
  IF(COMPILE_PYGVEC)
    STRING(TOLOWER "${${_var}}" _var_lower)
    IF(${_var_lower} MATCHES "^(on|yes|true|y)$")
      FILE(APPEND ${configuration_py} "${_var} = True\n")
    ELSEIF(${_var_lower} MATCHES "^(off|no|false|n)$")
      FILE(APPEND ${configuration_py} "${_var} = False\n")
    ELSEIF(EXISTS "${${_var}}")
      FILE(APPEND ${configuration_py} "${_var} = Path(\"${${_var}}\")\n")
    ELSE()
      FILE(APPEND ${configuration_py} "${_var} = \"${${_var}}\"\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

IF(COMPILE_PYGVEC)
  INSTALL(FILES "${configuration_py}" DESTINATION "${SKBUILD_PROJECT_NAME}")
ENDIF()
