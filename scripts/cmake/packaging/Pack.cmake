if(IS_SUBPROJECT)
    include (CPack)
    return()
endif()

include(packaging/PackagingMacros)
include(packaging/ArchiveTestdata)

#### Packaging setup ####
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OGS-6 THM/C Simulator")
set(CPACK_PACKAGE_VENDOR "OpenGeoSys Community (http://www.opengeosys.org)")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "OGS-${OGS_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
# set(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/README.md")

# Package file name
if(OGS_USE_PYTHON)
    set(SUFFIX "${SUFFIX}-python-${Python3_VERSION}")
endif()
if(OGS_BUILD_GUI)
    set(SUFFIX "${SUFFIX}-de")
endif()
if(OGS_BUILD_UTILS)
    set(SUFFIX "${SUFFIX}-utils")
endif()
if(OGS_USE_MPI)
    set(SUFFIX "${SUFFIX}-mpi")
endif()

if(APPLE)
    string(REGEX MATCH "(^[0-9]*)" TMP ${CMAKE_SYSTEM_VERSION})
    math(EXPR OSX_VERSION_MINOR "${CMAKE_MATCH_1} - 4")
    set(CPACK_PACKAGE_FILE_NAME
        "ogs-${OGS_VERSION}-OSX-10.${OSX_VERSION_MINOR}-${SUFFIX}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
else()
    set(CPACK_PACKAGE_FILE_NAME "ogs-${OGS_VERSION}-${CMAKE_SYSTEM}-${SUFFIX}")
endif()

if (WIN32)
    include (packaging/PackagingWin)
endif()
if(UNIX)
    if(APPLE)
        include (packaging/PackagingMac)
    else()
        include (packaging/PackagingLinux)
    endif()
endif()

include (CPack)

cpack_add_component_group(Applications
    DISPLAY_NAME Applications
    DESCRIPTION "OpenGeoSys applications"
    EXPANDED
    BOLD_TITLE
)

cpack_add_component_group(Utilities
    DISPLAY_NAME Utilities
    DESCRIPTION "OpenGeoSys utilities"
    EXPANDED
)

cpack_add_component(ogs_extras
    DISPLAY_NAME "Extra tools"
    DESCRIPTION "Miscellaneous tools."
    GROUP Utilities
)

cpack_add_component(ogs_docs
    DISPLAY_NAME "Documentation"
    DESCRIPTION "PDF documentation."
    GROUP Utilities
)

if(OGS_USE_CONAN)
    # Install Qt platform shared libraries
    install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/platforms DESTINATION bin OPTIONAL)
endif()

if(OGS_USE_PYTHON)
    if(WIN32)
        file(GLOB PYTHON_RUNTIME_LIBS "${Python3_RUNTIME_LIBRARY_DIRS}/*.dll")
        message(STATUS "Install Python: ${PYTHON_RUNTIME_LIBS}")
        install(FILES ${PYTHON_RUNTIME_LIBS} DESTINATION bin)
    else()
        install(FILES ${Python_LIBRARIES} DESTINATION bin)
    endif()
endif()


configure_file(Documentation/README.txt.in ${PROJECT_BINARY_DIR}/README.txt)
install(FILES ${PROJECT_BINARY_DIR}/README.txt DESTINATION .)

install(FILES ${PROJECT_BINARY_DIR}/CMakeCache.txt DESTINATION ${CMAKE_INSTALL_INFODIR})
if(EXISTS ${PROJECT_BINARY_DIR}/cmake-args)
    install(FILES ${PROJECT_BINARY_DIR}/cmake-args DESTINATION ${CMAKE_INSTALL_INFODIR})
endif()

# Install dependencies via GET_RUNTIME_DEPENDENCIES. Available since CMake 3.16.
if(${CMAKE_VERSION} VERSION_LESS 3.16)
    return()
endif()
install(CODE [[
  file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES $<TARGET_FILE:ogs>
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    POST_EXCLUDE_REGEXES "/opt/local/lib/lib.*" # Disable macports zlib
  )
  file(INSTALL ${_r_deps}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    FOLLOW_SYMLINK_CHAIN
  )
  list(LENGTH _u_deps _u_length)
  if("${_u_length}" GREATER 0)
    message(WARNING "Unresolved dependencies detected!\n${_u_deps}")
  endif()
]])
