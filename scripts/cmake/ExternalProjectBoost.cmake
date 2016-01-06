include(ThirdPartyLibVersions)
include(ExternalProject)

if(Boost_FOUND)
	return()
endif()

# First check for system boost
if(NOT Boost_INCLUDE_DIRS)
	if(APPLE)
		set(BOOST_ROOT $ENV{HOMEBREW_ROOT})
	endif()
	if(WIN32 AND COMPILER_IS_GCC)
		set(BOOST_INCLUDEDIR "$ENV{CMAKE_LIBRARY_SEARCH_PATH}/include/boost*")
	endif()
	if(OGS_LIB_BOOST STREQUAL "System")
		find_package(Boost 1.46.0 REQUIRED)
	elseif(OGS_LIB_BOOST STREQUAL "Default")
		find_package(Boost 1.46.0)
	endif()
	if(Boost_FOUND)
		set(Boost_FOUND TRUE CACHE BOOL "Was Boost found?" FORCE)
		set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIRS}" CACHE STRING "Boost include dirs" FORCE)
		return()
	else()
		set(Boost_INCLUDE_DIRS "")
	endif()
endif()

ExternalProject_Add(Boost
	PREFIX ${CMAKE_BINARY_DIR}/External/boost
	URL ${OGS_BOOST_URL}
	URL_MD5 ${OGS_BOOST_MD5}
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)
ExternalProject_Get_Property( Boost source_dir )

if(NOT Boost_INCLUDE_DIRS)
	set( Boost_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Boost include directories")
	# On Visual Studio Boost libs get automatically linked
	message(STATUS "Downloading Boost automatically.")
endif()
