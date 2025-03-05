  # Preparations to install GVEC for GENE

  MESSAGE(STATUS "${PROJECT_NAME} will be installed to ${CMAKE_INSTALL_PREFIX}")

  # standartize install output directories
  INCLUDE(GNUInstallDirs)

  # Offer the user the choice of overriding the installation directories
  SET(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH
    "Installation directory for libraries")
  SET(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH
    "Installation directory for executables")
  SET(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH
    "Installation directory for header files")
  IF(WIN32 AND NOT CYGWIN)
    SET(DEF_INSTALL_CMAKEDIR CMake)
  ELSE()
    SET(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
  ENDIF()
  SET(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH
    "Installation directory for CMake files")
  SET(INSTALL_MODDIR "mod" CACHE PATH
    "Installation directory for Fortran .mod files")

  # report installation directories
  FOREACH(p LIB BIN INCLUDE CMAKE MOD)
    FILE(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}DIR} _path )
    MESSAGE(STATUS "Installing ${p} components to ${_path}")
    UNSET(_path)
  ENDFOREACH()

  MARK_AS_ADVANCED(INSTALL_LIBDIR INSTALL_BINDIR INSTALL_INCLUDEDIR INSTALL_CMAKEDIR INSTALL_MODDIR)

  target_include_directories(gveclib
    INTERFACE
    $<INSTALL_INTERFACE:${INSTALL_MODDIR}>
    )

  get_property(ARE_SHARED_LIBS_SUPPORTED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
  IF(ARE_SHARED_LIBS_SUPPORTED)
    # Prepare RPATH
    FILE(RELATIVE_PATH _rel
      ${CMAKE_INSTALL_PREFIX}/${INSTALL_BINDIR}
      ${CMAKE_INSTALL_PREFIX})
    IF(APPLE)
      SET(_rpath "@loader_path/${_rel}")
    ELSE()
      SET(_rpath "\$ORIGIN/${_rel}")
    ENDIF()
    FILE(TO_NATIVE_PATH "${_rpath}/${INSTALL_LIBDIR}" test_gvec_to_gene_RPATH)

    SET_TARGET_PROPERTIES(test_gvec_to_gene
      PROPERTIES
      OUTPUT_NAME "test_gvec_to_gene"
      MACOSX_RPATH ON
      SKIP_BUILD_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_RPATH "${test_gvec_to_gene_RPATH}"
      INSTALL_RPATH_USE_LINK_PATH ON
      )
  ENDIF()

  # generate gvec_to_gene_export.h
  INCLUDE(GenerateExportHeader)
  generate_export_header( gveclib
    BASE_NAME gvec_to_gene)

  IF(ENABLE_PYTHON_BINDINGS)
    SET(_GVEC_TO_GENE_HEADERS
      ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/gvec_to_gene.h
      ${CMAKE_CURRENT_BINARY_DIR}/gvec_to_gene_export.h
      )

    SET(_GVEC_TO_GENE_RESOURCE
      "${CMAKE_CURRENT_BINARY_DIR}/interface_file_names.cfg"
      )

    SET(_GVEC_TO_GENE_PYTHON_INSTALL
      gvec_to_gene
      )

    FILE(RELATIVE_PATH _GVEC_TO_GENE_PYTHON_H_REL
      ${CMAKE_INSTALL_PREFIX}/${_GVEC_TO_GENE_PYTHON_INSTALL}
      ${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDEDIR}
      )

    FILE(RELATIVE_PATH _GVEC_TO_GENE_PYTHON_L_REL
      ${CMAKE_INSTALL_PREFIX}/${_GVEC_TO_GENE_PYTHON_INSTALL}
      ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBDIR}
      )

    FILE(
      GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/interface_file_names.cfg
      #INPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/interface_file_names.cfg.in
      CONTENT
      "[configuration]
    header_file_name = ${_GVEC_TO_GENE_PYTHON_H_REL}/gvec_to_gene.h
    library_file_name = ${_GVEC_TO_GENE_PYTHON_L_REL}/$<TARGET_FILE_NAME:gveclib>"
      )

    # Install python package
    IF(NOT COMPILE_PYGVEC)
      INSTALL(
        FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gvec_to_gene/__init__.py
        DESTINATION
        ${_GVEC_TO_GENE_PYTHON_INSTALL}
        )
    ENDIF()

    # Add additional properties to gveclib to enable Python bindings
    SET_TARGET_PROPERTIES( gveclib
      PROPERTIES
      PUBLIC_HEADER
      "${_GVEC_TO_GENE_HEADERS}"
      RESOURCE
      "${_GVEC_TO_GENE_RESOURCE}"
      )
  ENDIF()

  IF(COMPILE_PYGVEC)
    # pyGVEC installation (via scikit-build-core): only install test_gvec_to_gene for now
    INSTALL(
      TARGETS
      test_gvec_to_gene
      DESTINATION
      ${SKBUILD_SCRIPTS_DIR}
    )
  ELSE()
  # install gveclib and test_gvec_to_gene
  INSTALL(
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

  INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${INSTALL_MODDIR})

  # Exporting targets for other CMake-projects

  INSTALL(
    EXPORT
    GVECTargets
    NAMESPACE
    "GVEC::"
    DESTINATION
    ${INSTALL_CMAKEDIR}
    COMPONENT
    dev
    )

  INCLUDE(CMakePackageConfigHelpers)

  configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/GVECConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GVECConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CMAKEDIR}
    )

  INSTALL(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/GVECConfig.cmake
    DESTINATION
    ${INSTALL_CMAKEDIR}
    )

  # turn on testing
  ENABLE_TESTING()

  IF(ENABLE_PYTHON_BINDINGS)
    # Testing through Python

    # require python
    FIND_PACKAGE(PythonInterp REQUIRED)

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
  ENDIF()
  ENDIF()
