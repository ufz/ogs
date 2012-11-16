INCLUDE(ExternalProject)
# Set Boost version and which libraries to compile
SET(Boost_Version 1.49.0)
SET(BOOST_LIBS_TO_BUILD
	--with-program_options
	--with-filesystem
	--with-system
	--with-date_time
)

STRING(REPLACE "." "_" Boost_Version_Underscore ${Boost_Version})

# Set boost toolset
IF(MSVC)
	SET(BOOST_TOOLSET msvc-10.0) # TODO: VS 2010 only
ENDIF()
IF(COMPILER_IS_GCC)
	SET(BOOST_TOOLSET gcc)
ENDIF()
IF(APPLE)
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
		SET(BOOST_CONFIG_OPTIONS "architecture=x86 address-model=64")
	ENDIF()
ENDIF()


ExternalProject_Add(Boost
	PREFIX ${CMAKE_BINARY_DIR}/External/boost
	URL http://downloads.sourceforge.net/project/boost/boost/${Boost_Version}/boost_${Boost_Version_Underscore}.zip
	URL_MD5 854dcbbff31b896c85c38247060b7713
	UPDATE_COMMAND "${BOOST_UPDATE_COMMAND}"
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ./b2 ${BOOST_LIBS_TO_BUILD} toolset=${BOOST_TOOLSET} link=static stage ${BOOST_CONFIG_OPTIONS}
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)
ExternalProject_Get_Property( Boost source_dir )

#IF(NOT Boost_INCLUDE_DIRS)
	SET( Boost_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Boost include directories")
	# TODO this has to be set manually!
	FILE( GLOB Boost_LIBRARIES "${source_dir}/stage/lib/" "${source_dir}/stage/lib/*.a")
	SET( Boost_LIBRARIES ${Boost_LIBRARIES} CACHE INTERNAL "Boost libraries")
	MESSAGE(STATUS "Boost inc: ${Boost_INCLUDE_DIRS}")
	MESSAGE(STATUS "Boost libs: ${Boost_LIBRARIES}")
#ENDIF()

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )