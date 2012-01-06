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

# Suppress warning on add_subdirectory(dir) where dir contains no CMakeLists.txt
IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.7)
	CMAKE_POLICY(SET CMP0014 OLD)
ENDIF ()
