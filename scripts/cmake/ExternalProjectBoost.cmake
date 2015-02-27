include(ThirdPartyLibVersions)
include(ExternalProject)
# Set Boost version and which libraries to compile
set(BOOST_LIBS_TO_BUILD
	# chrono
	# context
	date_time
	# exception
	filesystem
	# graph
	# graph_parallel
	# iostreams
	# locale
	# math
	# mpi
	# program_options
	# python
	# random
	# regex
	# serialization
	# signals
	system
	# test
	# thread
	# timer
	# wave
)

if(Boost_FOUND)
	return()
endif()

# Set boost toolset
if(MSVC11)
	set(BOOST_TOOLSET msvc-11.0)
elseif(MSVC12)
	set(BOOST_TOOLSET msvc-12.0)
elseif(APPLE)
	set(BOOST_TOOLSET darwin)
elseif(COMPILER_IS_CLANG)
	set(BOOST_TOOLSET clang)
elseif(COMPILER_IS_INTEL)
	# Extracts first two version numbers, e.g. 13.0 from 13.0.0.20121010
	string(REGEX MATCH ^[0-9]*.[0-9] INTEL_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
	set(BOOST_TOOLSET intel-${INTEL_COMPILER_VERSION})
elseif(COMPILER_IS_GCC)
	set(BOOST_TOOLSET gcc)
endif()

# First check for system boost
if(NOT Boost_INCLUDE_DIRS)
	if(APPLE)
		set(BOOST_ROOT $ENV{HOMEBREW_ROOT})
	endif()
	if(MSVC)
		set(Boost_USE_STATIC_LIBS ON)
		if(NOT DEFINED BOOST_LIBRARYDIR)
			set(BOOST_LIBRARYDIR "$ENV{CMAKE_LIBRARY_SEARCH_PATH}/boost/lib${BITS}-${BOOST_TOOLSET};C:/boost/lib${BITS}-${BOOST_TOOLSET};$ENV{BOOST_ROOT}/lib${BITS}-${BOOST_TOOLSET}")
			set(BOOST_INCLUDEDIR "$ENV{CMAKE_LIBRARY_SEARCH_PATH}/boost;C:/boost;$ENV{BOOST_ROOT}")
		endif()
	endif()
	if(WIN32 AND COMPILER_IS_GCC)
		set(BOOST_INCLUDEDIR "$ENV{CMAKE_LIBRARY_SEARCH_PATH}/include/boost*")
	endif()
	find_package(Boost 1.46.0 COMPONENTS ${BOOST_LIBS_TO_BUILD})
	if(Boost_FOUND)
		set(Boost_FOUND TRUE CACHE BOOL "Was Boost found?" FORCE)
		set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIRS}" CACHE STRING "Boost include dirs" FORCE)
		set(Boost_LIBRARIES "${Boost_LIBRARIES}" CACHE STRING "Boost libraries" FORCE)
		return()
	else()
		set(Boost_INCLUDE_DIRS "")
	endif()
endif()

# Check required gcc
if(COMPILER_IS_GCC AND GCC_VERSION VERSION_LESS "4.4")
	message(FATAL_ERROR "GCC version >= 4.4 is required for building Boost ${Boost_Version}!")
endif()

# Prefix with --with- for bjam (b2) build command
foreach(LIB_TO_BUILD ${BOOST_LIBS_TO_BUILD})
	set(BOOST_LIBS_TO_BUILD_CMD ${BOOST_LIBS_TO_BUILD_CMD} --with-${LIB_TO_BUILD})
endforeach()

# Prefix with boost_ for library names
foreach(LIB_TO_BUILD ${BOOST_LIBS_TO_BUILD})
	set(BOOST_LIBS_TO_BUILD_NAMES ${BOOST_LIBS_TO_BUILD_NAMES} boost_${LIB_TO_BUILD})
endforeach()

# Set update command
if(WIN32)
	set(BOOST_UPDATE_COMMAND bootstrap.bat)
else()
	set(BOOST_UPDATE_COMMAND ./bootstrap.sh)
endif()

# Set additional config options
set(BOOST_CONFIG_OPTIONS "")
if(WIN32)
	if(HAVE_64_BIT)
		set(BOOST_CONFIG_OPTIONS "architecture=x86;address-model=64")
	endif()
endif()

ExternalProject_Add(Boost
	PREFIX ${CMAKE_BINARY_DIR}/External/boost
	URL ${OGS_BOOST_URL}
	URL_MD5 ${OGS_BOOST_MD5}
	UPDATE_COMMAND "${BOOST_UPDATE_COMMAND}"
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ./b2 ${BOOST_LIBS_TO_BUILD_CMD} -j ${NUM_PROCESSORS} toolset=${BOOST_TOOLSET} link=static stage ${BOOST_CONFIG_OPTIONS}
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)
ExternalProject_Get_Property( Boost source_dir )

if(NOT Boost_INCLUDE_DIRS)
	set( Boost_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Boost include directories")
	# On Visual Studio Boost libs get automatically linked
	if(MSVC)
		set( Boost_LIBRARIES "" CACHE INTERNAL "Boost libraries")
	else()
		set( Boost_LIBRARIES ${BOOST_LIBS_TO_BUILD_NAMES} CACHE INTERNAL "Boost libraries")
	endif()
	message(STATUS "Building Boost automatically.")
endif()

link_directories( ${source_dir}/stage/lib/ )
