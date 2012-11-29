INCLUDE(ExternalProject)
# Set Boost version and which libraries to compile
SET(Boost_Version 1.52.0)
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

# First check for system boost
SET(Boost_USE_STATIC_LIBS ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
FIND_PACKAGE(Boost COMPONENTS ${BOOST_LIBS_TO_BUILD})
IF(Boost_FOUND)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	RETURN()
ENDIF()

# Prefix with --with- for bjam (b2) build command
FOREACH(LIB_TO_BUILD ${BOOST_LIBS_TO_BUILD})
	SET(BOOST_LIBS_TO_BUILD_CMD ${BOOST_LIBS_TO_BUILD_CMD} --with-${LIB_TO_BUILD})
ENDFOREACH()

# Prefix with boost_ for library names
FOREACH(LIB_TO_BUILD ${BOOST_LIBS_TO_BUILD})
	SET(BOOST_LIBS_TO_BUILD_NAMES ${BOOST_LIBS_TO_BUILD_NAMES} boost_${LIB_TO_BUILD})
ENDFOREACH()

STRING(REPLACE "." "_" Boost_Version_Underscore ${Boost_Version})

# Set boost toolset
IF(MSVC10)
	SET(BOOST_TOOLSET msvc-10.0)
ELSEIF(MSVC11)
	SET(BOOST_TOOLSET msvc-11.0)
ELSEIF(COMPILER_IS_GCC)
	SET(BOOST_TOOLSET gcc)
ELSEIF(COMPILER_IS_INTEL)
	# Extracts first two version numbers, e.g. 13.0 from 13.0.0.20121010
	STRING(REGEX MATCH ^[0-9]*.[0-9] INTEL_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
	SET(BOOST_TOOLSET intel-${INTEL_COMPILER_VERSION})
ELSEIF(APPLE)
	SET(BOOST_TOOLSET darwin)
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

# Set archive sources
SET(BOOST_ARCHIVE_EXT "tar.bz2")
SET(BOOST_ARCHIVE_MD5 3a855e0f919107e0ca4de4d84ad3f750)

SET(BOOST_URL "http://switch.dl.sourceforge.net/project/boost/boost/${Boost_Version}/boost_${Boost_Version_Underscore}.${BOOST_ARCHIVE_EXT}")

ExternalProject_Add(Boost
	PREFIX ${CMAKE_BINARY_DIR}/External/boost
	URL ${BOOST_URL}
	URL_MD5 ${BOOST_ARCHIVE_MD5}
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
ENDIF()

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${source_dir}/stage/lib/ )
