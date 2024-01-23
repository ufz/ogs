# Set build directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)

# cmake-lint: disable=C0103
set(Data_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Tests/Data CACHE INTERNAL "")
set(Data_BINARY_DIR ${PROJECT_BINARY_DIR}/Tests/Data CACHE INTERNAL "")

# Enable Visual Studio project folder grouping
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# RPATH setup
if(APPLE)
    set(BASEPOINT @loader_path)
else()
    set(BASEPOINT $ORIGIN)
endif()
file(RELATIVE_PATH relDir ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)

set(_py_module_rpath ${BASEPOINT}/${CMAKE_INSTALL_LIBDIR})
if(NOT OGS_BUILD_WHEEL)
    set(_py_module_rpath ${BASEPOINT}/../../..)
endif()
list(APPEND CMAKE_INSTALL_RPATH ${BASEPOINT} ${BASEPOINT}/${relDir}
     ${_py_module_rpath} # Python modules
)
list(APPEND CMAKE_BUILD_RPATH ${BASEPOINT} ${BASEPOINT}/${relDir}
     ${BASEPOINT}/${CMAKE_INSTALL_LIBDIR} # Python modules
)

# Some external dependencies always use lib instead of lib64, Fix for
# lib64-based systems, e.g. OpenSUSE:
if("${CMAKE_INSTALL_LIBDIR}" STREQUAL "lib64")
    file(RELATIVE_PATH relDirLib
         ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
         ${CMAKE_CURRENT_BINARY_DIR}/lib
    )
    list(APPEND CMAKE_INSTALL_RPATH ${BASEPOINT}/${relDirLib})
    list(APPEND CMAKE_BUILD_RPATH ${BASEPOINT}/${relDirLib})
endif()
