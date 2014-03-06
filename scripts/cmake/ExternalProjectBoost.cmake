INCLUDE(ThirdPartyLibVersions)
INCLUDE(ExternalProject)
# Set Boost version and which libraries to compile
SET(BOOST_LIBS_TO_BUILD
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
	program_options
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

IF(Boost_FOUND)
	RETURN()
ENDIF()

# First check for system boost
IF(NOT Boost_INCLUDE_DIRS)
	IF(APPLE)
		SET(BOOST_ROOT $ENV{HOMEBREW_ROOT})
	ENDIF()
	IF(WIN32)
		SET(Boost_USE_STATIC_LIBS ON)
		IF(NOT DEFINED BOOST_LIBRARYDIR)
			SET(BOOST_LIBRARYDIR "C:/boost/lib${BITS}-msvc-11.0;C:/boost/lib${BITS}-msvc-12.0")
		ENDIF()
	ENDIF()
	FIND_PACKAGE(Boost 1.46.0 COMPONENTS ${BOOST_LIBS_TO_BUILD})
	IF(Boost_FOUND)
		SET(Boost_FOUND TRUE CACHE BOOL "Was Boost found?" FORCE)
		SET(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIRS}" CACHE STRING "Boost include dirs" FORCE)
		SET(Boost_LIBRARIES "${Boost_LIBRARIES}" CACHE STRING "Boost libraries" FORCE)
		RETURN()
	ELSE()
		SET(Boost_INCLUDE_DIRS "")
	ENDIF()
ENDIF()

# Check required gcc
IF(COMPILER_IS_GCC AND GCC_VERSION VERSION_LESS "4.4")
	MESSAGE(FATAL_ERROR "GCC version >= 4.4 is required for building Boost ${Boost_Version}!")
ENDIF()

# Prefix with --with- for bjam (b2) build command
FOREACH(LIB_TO_BUILD ${BOOST_LIBS_TO_BUILD})
	SET(BOOST_LIBS_TO_BUILD_CMD ${BOOST_LIBS_TO_BUILD_CMD} --with-${LIB_TO_BUILD})
ENDFOREACH()

# Prefix with boost_ for library names
FOREACH(LIB_TO_BUILD ${BOOST_LIBS_TO_BUILD})
	SET(BOOST_LIBS_TO_BUILD_NAMES ${BOOST_LIBS_TO_BUILD_NAMES} boost_${LIB_TO_BUILD})
ENDFOREACH()

# Set boost toolset
IF(MSVC10)
	SET(BOOST_TOOLSET msvc-10.0)
ELSEIF(MSVC11)
	SET(BOOST_TOOLSET msvc-11.0)
ELSEIF(APPLE)
	SET(BOOST_TOOLSET darwin)
ELSEIF(COMPILER_IS_CLANG)
	SET(BOOST_TOOLSET clang)
ELSEIF(COMPILER_IS_INTEL)
	# Extracts first two version numbers, e.g. 13.0 from 13.0.0.20121010
	STRING(REGEX MATCH ^[0-9]*.[0-9] INTEL_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
	SET(BOOST_TOOLSET intel-${INTEL_COMPILER_VERSION})
ELSEIF(COMPILER_IS_GCC)
	SET(BOOST_TOOLSET gcc)
ENDIF()

# Set update command
IF(WIN32)
	SET(BOOST_UPDATE_COMMAND bootstrap.bat)
ELSE()
	SET(BOOST_UPDATE_COMMAND ./bootstrap.sh)
ENDIF()

# Set additional config options
SET(BOOST_CONFIG_OPTIONS "")
IF(WIN32)
	IF(HAVE_64_BIT)
		SET(BOOST_CONFIG_OPTIONS "architecture=x86;address-model=64")
	ENDIF()
ENDIF()

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

IF(NOT Boost_INCLUDE_DIRS)
	SET( Boost_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Boost include directories")
	# On Visual Studio Boost libs get automatically linked
	IF(MSVC)
		SET( Boost_LIBRARIES "" CACHE INTERNAL "Boost libraries")
	ELSE()
		SET( Boost_LIBRARIES ${BOOST_LIBS_TO_BUILD_NAMES} CACHE INTERNAL "Boost libraries")
	ENDIF()
	MESSAGE(STATUS "Building Boost automatically.")
ENDIF()

LINK_DIRECTORIES( ${source_dir}/stage/lib/ )
