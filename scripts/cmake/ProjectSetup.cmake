# Set build directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
if(OGS_USE_CONAN AND MSVC)
    foreach(outputconfig ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${outputconfig} outputconfig)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${outputconfig}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        )
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${outputconfig}
            ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        )
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${outputconfig}
            ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        )
    endforeach()
endif()

# cmake-lint: disable=C0103
set(Data_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Tests/Data CACHE INTERNAL "")
set(Data_BINARY_DIR ${PROJECT_BINARY_DIR}/Tests/Data CACHE INTERNAL "")

# Enable Visual Studio project folder grouping
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# RPATH setup
if(APPLE)
    set(BASEPOINT @loader_path)
else()
    set(BASEPOINT $ORIGIN)
endif()
file(RELATIVE_PATH relDir ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_INSTALL_RPATH ${BASEPOINT} ${BASEPOINT}/${relDir})
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
