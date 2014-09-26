INCLUDE(ThirdPartyLibVersions)
INCLUDE(ExternalProject)

SET(CATALYST_GIT_URL https://github.com/ufz/catalyst-io.git)

FIND_PACKAGE(ParaView 4.2 COMPONENTS vtkIOXML QUIET)

IF(ParaView_FOUND)
	INCLUDE("${PARAVIEW_USE_FILE}")
	# MESSAGE("Using Catalyst in ${ParaView_FOUND}")
	RETURN()
ELSE()
	SET(ParaView_DIR ${CMAKE_BINARY_DIR}/External/catalyst/src/catalyst-build CACHE PATH "" FORCE)
ENDIF()

IF(WIN32)
	SET(CATALYST_MAKE_COMMAND
		cmake --build . --config Release --target vtkIO &&
		cmake --build . --config Debug --target vtkIO)
	SET(CATALYST_CONFIGURE_COMMAND cmake.bat)
	# MESSAGE(STATUS ${CATALYST_MAKE_COMMAND})
ELSE()
	IF($ENV{CI})
		SET(CATALYST_MAKE_COMMAND make vtkIO)
	ELSE()
		SET(CATALYST_MAKE_COMMAND make -j ${NUM_PROCESSORS} vtkIO)
	ENDIF()
	SET(CATALYST_CONFIGURE_COMMAND cmake.sh)
ENDIF()

ExternalProject_Add(Catalyst
	PREFIX ${CMAKE_BINARY_DIR}/External/catalyst
	GIT_REPOSITORY ${CATALYST_GIT_URL}
	#URL ${OGS_VTK_URL}
	#URL_MD5 ${OGS_VTK_MD5}
	CONFIGURE_COMMAND ../Catalyst/${CATALYST_CONFIGURE_COMMAND} ../Catalyst
	BUILD_COMMAND ${CATALYST_MAKE_COMMAND}
	INSTALL_COMMAND ""
)

IF(NOT ${ParaView_FOUND})
	# Rerun cmake in initial build
	ADD_CUSTOM_TARGET(VtkRescan ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} DEPENDS Catalyst)
ELSE()
	ADD_CUSTOM_TARGET(VtkRescan) # dummy target for caching
ENDIF()
