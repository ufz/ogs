# Set additional CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake")

MESSAGE(STATUS "CMAKE_MODULE_PATH is set to ${CMAKE_MODULE_PATH}")

# Load addional modules
INCLUDE(UseBackportedModules)
INCLUDE(OptionRequires)
INCLUDE(CppcheckTargets)
INCLUDE(GetCompilerInfoString)
INCLUDE(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
SET(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

# Suppress warning on setting policies
CMAKE_POLICY(SET CMP0011 OLD)

# Get the hostname
SITE_NAME(HOSTNAME)
