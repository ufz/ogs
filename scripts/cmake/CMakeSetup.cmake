# Disallow in-source builds as the git project cluttered with generated files
# probably confuses people. source/build* is still allowed!
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed!\n"
    "Make sure to remove CMakeCache.txt and CMakeFiles/ "
    "from the source directory!")
endif()

# Set additional CMake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/cmake-modules")

list(APPEND CMAKE_PREFIX_PATH
  $ENV{HOMEBREW_ROOT}             # Homebrew package manager on Mac OS
  $ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
  ${CMAKE_LIBRARY_SEARCH_PATH})   # CMake option, Windows

# Load addional modules
include(GNUInstallDirs)
include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
set(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

# Check if this project is included in another
if(NOT PROJECT_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IS_SUBPROJECT ON CACHE INTERNAL "" FORCE)
    set(OGS_BUILD_CLI OFF CACHE BOOL "" FORCE)
endif()

if((NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    OR (NOT CMAKE_BUILD_TYPE AND MSVC AND OGS_USE_CONAN))
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Get the hostname
site_name(HOSTNAME)

# Check if we are running under CI
if(DEFINED ENV{JENKINS_URL} OR DEFINED ENV{CI})
    set(IS_CI ON CACHE INTERNAL "")
endif()
