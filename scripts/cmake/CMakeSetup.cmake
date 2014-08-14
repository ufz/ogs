# Set additional CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake")

# Load addional modules
INCLUDE(UseBackportedModules)
INCLUDE(OptionRequires)
INCLUDE(CppcheckTargets)
INCLUDE(GetCompilerInfoString)

INCLUDE(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
SET(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

INCLUDE(GetGitRevisionDescription)
GET_GIT_HEAD_REVISION(GIT_REFSPEC GIT_SHA1)
STRING(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1_SHORT)

# Suppress warning on setting policies
CMAKE_POLICY(SET CMP0011 OLD)

# Get the hostname
SITE_NAME(HOSTNAME)
