# Set additional CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/scripts/cmake/cmake"
  "${CMAKE_SOURCE_DIR}/scripts/cmake")

# Load addional modules
INCLUDE(UseBackportedModules)
INCLUDE(OptionRequires)
INCLUDE(CppcheckTargets)
INCLUDE(GetCompilerInfoString)

# Suppress warning on setting policies
CMAKE_POLICY(SET CMP0011 OLD)

# Get the hostname
SITE_NAME(HOSTNAME)
