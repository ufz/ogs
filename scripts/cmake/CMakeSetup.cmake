# Set additional CMake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/cmake-modules")

# Load addional modules
include(UseBackportedModules)
include(OptionRequires)
include(CppcheckTargets)
include(GetCompilerInfoString)
include(AddCatalystDependency)

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
set(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

# Check if this project is included in another
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(IS_SUBPROJECT ON CACHE INTERNAL "" FORCE)
	set(OGS_BUILD_CLI OFF CACHE BOOL "" FORCE)
endif()

# Get the hostname
site_name(HOSTNAME)

# Compute OS X version number
if(APPLE)
	if(CMAKE_SYSTEM_VERSION VERSION_EQUAL 12.0)
		set(OSX_VERSION 10.8 CACHE STRING "OS X version number")
		set(OSX_VERSION_NAME CACHE STRING "Mountain Lion")
	endif()
	if(CMAKE_SYSTEM_VERSION VERSION_EQUAL 13.0)
		set(OSX_VERSION 10.9 CACHE STRING "OS X version number")
		set(OSX_VERSION_NAME CACHE STRING "Mavericks")
	endif()
	if(CMAKE_SYSTEM_VERSION VERSION_EQUAL 14.0)
		set(OSX_VERSION 10.10 CACHE STRING "OS X version number")
		set(OSX_VERSION_NAME CACHE STRING "Yosemite")
	endif()
endif()
mark_as_advanced(OSX_VERSION OSX_VERSION_NAME)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
