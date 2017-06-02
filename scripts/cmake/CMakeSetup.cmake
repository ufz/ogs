# Set additional CMake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/cmake-modules")

# Load addional modules
include(UseBackportedModules)
include(OptionRequires)
include(CppcheckTargets)
include(AddVtkDependency)

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
set(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

# Check if this project is included in another
if(NOT PROJECT_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IS_SUBPROJECT ON CACHE INTERNAL "" FORCE)
    set(OGS_BUILD_CLI OFF CACHE BOOL "" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Get the hostname
site_name(HOSTNAME)
