# Set additional CMake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake")

# Load addional modules
include(UseBackportedModules)
include(OptionRequires)
include(CppcheckTargets)
include(GetCompilerInfoString)
include(AddCatalystDependency)

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
set(NUM_PROCESSORS ${NUM_PROCESSORS} CACHE STRING "Processor count")

include(GetGitRevisionDescription)
GET_GIT_HEAD_REVISION(GIT_REFSPEC GIT_SHA1)
string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1_SHORT)

# Check if this project is included in another
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(IS_SUBPROJECT ON CACHE INTERNAL "" FORCE)
	set(OGS_BUILD_CLI OFF CACHE BOOL "" FORCE)
endif()

if($ENV{CI})
	set(OGS_VERSION 6.6.6) # Dummy version for CI-environment (Travis) or subproject
	# Workaround for https://issues.jenkins-ci.org/browse/JENKINS-29329
	if(MSVC)
		set(CMAKE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} CACHE FILEPATH "workaround for jenkins cmake-plugin")
	endif()
elseif(IS_SUBPROJECT)
	set(OGS_VERSION x.x.x)
else()
	GIT_GET_TAG(GIT_DESCRIBE)
	if(GIT_DESCRIBE)
		string(REGEX MATCH ^[0-9|\\.]+ GIT_TAG ${GIT_DESCRIBE})
		set(OGS_VERSION ${GIT_TAG})

		if(GIT_DESCRIBE MATCHES ".*-.*-.*")
			# Commit is not a tag
			set(OGS_IS_RELEASE FALSE)
			string(REGEX MATCH "-([0-9]+)-" GIT_COMMITS_AFTER_TAG ${GIT_DESCRIBE})
		else()
			set(OGS_VERSION_STATUS "")
			set(OGS_IS_RELEASE TRUE)
			message(STATUS "OGS version: ${OGS_VERSION}")
		endif()
	else()
		message(WARNING "Git repository contains no tags! Please run: git fetch --tags")
	endif()
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
