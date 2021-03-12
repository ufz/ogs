if(_IS_SUBPROJECT)
    include (CPack)
    return()
endif()

# Put ogs installs into its own component and then only install
# this component (avoids third-party installs from CPM).
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME ogs)
set(CPACK_INSTALL_CMAKE_PROJECTS
    "${PROJECT_BINARY_DIR};${PROJECT_NAME};${CMAKE_INSTALL_DEFAULT_COMPONENT_NAME};/"
)

option(OGS_INSTALL_DEPENDENCIES "Package dependencies.")
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
    list(APPEND SUFFIX_LIST "python-${Python3_VERSION}")
endif()
if(OGS_BUILD_GUI)
    list(APPEND SUFFIX_LIST "de")
endif()
if(OGS_BUILD_UTILS)
    list(APPEND SUFFIX_LIST "utils")
endif()
if(OGS_USE_MPI)
    list(APPEND SUFFIX_LIST "mpi")
endif()
string(REPLACE ";" "-" SUFFIX "${SUFFIX_LIST}")

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
if(${CMAKE_VERSION} VERSION_LESS 3.16 OR NOT OGS_INSTALL_DEPENDENCIES)
    return()
endif()
install(CODE [[
  include(GNUInstallDirs)
  if(WIN32)
    set(INSTALL_DIR ${CMAKE_INSTALL_FULL_BINDIR})
  else()
    set(INSTALL_DIR ${CMAKE_INSTALL_FULL_LIBDIR})
  endif()
  file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES $<$<TARGET_EXISTS:ogs>:$<TARGET_FILE:ogs>> $<$<TARGET_EXISTS:DataExplorer>:$<TARGET_FILE:DataExplorer>>
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    POST_EXCLUDE_REGEXES "/opt/local/lib/lib.*" # Disable macports zlib
  )
  file(INSTALL ${_r_deps}
    DESTINATION ${INSTALL_DIR}
    FOLLOW_SYMLINK_CHAIN
  )
  list(LENGTH _u_deps _u_length)
  if("${_u_length}" GREATER 0)
    message(WARNING "Unresolved dependencies detected!\n${_u_deps}")
  endif()
]])
