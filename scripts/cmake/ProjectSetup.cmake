# Set build directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
if(OGS_USE_CONAN AND MSVC)
    foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    endforeach(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)
endif()

set(Data_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Tests/Data CACHE INTERNAL "")
set(Data_BINARY_DIR ${PROJECT_BINARY_DIR}/Tests/Data CACHE INTERNAL "")

# Enable Visual Studio project folder grouping
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# RPATH setup
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../${CMAKE_INSTALL_LIBDIR}")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

if(NOT IS_GIT_REPO)
    return()
endif()

if(DEFINED ENV{OGS_VERSION})
    set(OGS_VERSION $ENV{OGS_VERSION})
    message(STATUS "OGS VERSION: ${OGS_VERSION} (set via environment)")
else()
    # Get version info from Git, implementation based on
    # https://github.com/tomtom-international/cpp-dependencies
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --long --dirty --always
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE DESCRIBE_RESULT
        OUTPUT_VARIABLE DESCRIBE_STDOUT
    )
    if(DESCRIBE_RESULT EQUAL 0)
        string(STRIP "${DESCRIBE_STDOUT}" DESCRIBE_STDOUT)
        message(STATUS "Git reported this project's version as '${DESCRIBE_STDOUT}'")
        if(DESCRIBE_STDOUT MATCHES "^(.*)-(dirty)$")
          set(DESCRIBE_DIRTY "${CMAKE_MATCH_2}")
          set(DESCRIBE_STDOUT "${CMAKE_MATCH_1}")
        endif()
        if(DESCRIBE_STDOUT MATCHES "^([0-9a-f]+)$")
          set(DESCRIBE_COMMIT_NAME "${CMAKE_MATCH_1}")
          set(DESCRIBE_STDOUT "")
        elseif(DESCRIBE_STDOUT MATCHES "^(.*)-g([0-9a-f]+)$")
          set(DESCRIBE_COMMIT_NAME "g${CMAKE_MATCH_2}")
          set(DESCRIBE_STDOUT "${CMAKE_MATCH_1}")
        endif()
        if(DESCRIBE_STDOUT MATCHES "^(.*)-([0-9]+)$")
          set(DESCRIBE_COMMIT_COUNT "${CMAKE_MATCH_2}")
          set(DESCRIBE_TAG "${CMAKE_MATCH_1}")
          set(DESCRIBE_STDOUT "")
        endif()

        set(OGS_VERSION ${DESCRIBE_TAG})
        if(DESCRIBE_COMMIT_COUNT GREATER 0)
          set(OGS_VERSION "${OGS_VERSION}-${DESCRIBE_COMMIT_COUNT}-${DESCRIBE_COMMIT_NAME}")
        endif()

        if(DESCRIBE_DIRTY)
          string(TIMESTAMP DESCRIBE_DIRTY_TIMESTAMP "%Y%m%d%H%M%S" UTC)
          set(OGS_VERSION "${OGS_VERSION}.dirty.${DESCRIBE_DIRTY_TIMESTAMP}")
        endif()
        message(STATUS "OGS VERSION: ${OGS_VERSION}")
    else()
        message(WARNING "Git repository contains no tags! Please run: git fetch --tags")
    endif()
endif()

# Get git commit
execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%H
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_SHA1
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_SHA1_SHORT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
