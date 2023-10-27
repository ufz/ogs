# Disallow in-source builds as the git project cluttered with generated files
# probably confuses people. source/build* is still allowed!
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(
        FATAL_ERROR
            "In-source builds are not allowed!\n"
            "Make sure to remove CMakeCache.txt and CMakeFiles/ "
            "from the source directory!"
    )
endif()

message(STATUS "Generator: ${CMAKE_GENERATOR}")
if(WIN32 AND (NOT "${CMAKE_GENERATOR}" MATCHES "Visual Studio")
   AND "$ENV{CIBUILDWHEEL}"
)
    message(FATAL_ERROR "Wheels only build in Visual Studio!")
endif()

set(_collection ${PROJECT_SOURCE_DIR}/ThirdParty/collection)
# If submodules in ThirdParty/collection are initialized and this is a Guix
# build use submodule as CPM sources.
if(EXISTS ${_collection}/ufz/vtkdiff/CMakeLists.txt AND GUIX_BUILD)
    include(${_collection}/Setup.cmake)
endif()

# Set additional CMake modules path To be replaced later. See
# https://gitlab.kitware.com/cmake/cmake/-/issues/22831
CPMAddPackage(
    NAME findmkl_cmake
    GITHUB_REPOSITORY bilke/findmkl_cmake
    GIT_TAG ee49c4f973f66bb7bfd644658d14e43459f557fa
    DOWNLOAD_ONLY YES
)
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    "${PROJECT_SOURCE_DIR}/scripts/cmake"
    "${PROJECT_SOURCE_DIR}/scripts/cmake/jedbrown"
    "${PROJECT_SOURCE_DIR}/scripts/cmake/vector-of-bool"
    "${findmkl_cmake_SOURCE_DIR}/cmake"
)

list(
    APPEND
    CMAKE_PREFIX_PATH
    $ENV{HOMEBREW_ROOT} # Homebrew package manager on Mac OS
    $ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
    ${CMAKE_LIBRARY_SEARCH_PATH}
) # CMake option, Windows

# Load additional modules
include(GNUInstallDirs)
include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
set(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

if(NOT PROJECT_IS_TOP_LEVEL)
    set(OGS_BUILD_CLI OFF CACHE BOOL "" FORCE)
endif()

if((NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   OR (NOT CMAKE_BUILD_TYPE AND MSVC AND OGS_USE_CONAN)
)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    set_property(
        CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                        "RelWithDebInfo"
    )
endif()

# Get the hostname
site_name(HOSTNAME)

# When static libraries are used in some shared libraries it is required that
# also the static libraries have position independent code.
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

# Enable Windows DLL support.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
