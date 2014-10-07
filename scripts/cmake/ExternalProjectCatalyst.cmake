INCLUDE(ThirdPartyLibVersions)
INCLUDE(ExternalProject)

SET(CATALYST_GIT_URL https://github.com/ufz/catalyst-io.git)

IF(NOT DEFINED ParaView_DIR AND DEFINED ENV{ParaView_DIR})
	SET(ParaView_DIR $ENV{ParaView_DIR})
ENDIF()

FIND_PACKAGE(ParaView 4.2 COMPONENTS vtkIOXML QUIET)

FIND_LIBRARY(VTKIO_LIB_FOUND vtkIOXML-pv4.2 HINTS ${ParaView_DIR}/lib PATH_SUFFIXES Release Debug)
IF(ParaView_FOUND AND VTKIO_LIB_FOUND)
	INCLUDE("${PARAVIEW_USE_FILE}")
	# MESSAGE("Using Catalyst in ${ParaView_FOUND}")
	RETURN()
ELSEIF(NOT ParaView_DIR)
	SET(ParaView_DIR ${CMAKE_BINARY_DIR}/External/catalyst/src/Catalyst-build CACHE PATH "" FORCE)
ENDIF()

SET(CATALYST_CMAKE_GENERATOR ${CMAKE_GENERATOR})
IF(WIN32)
	# Ninja tempory disabled because it builds only the Release mode.
	# FIND_PROGRAM(NINJA_TOOL_PATH ninja DOC "Ninja build tool")
	IF(NINJA_TOOL_PATH)
		SET(CATALYST_CMAKE_GENERATOR Ninja)
		SET(CATALYST_MAKE_COMMAND ninja vtkIOXML)
	ELSE()
		SET(CATALYST_MAKE_COMMAND
			cmake --build . --config Release --target vtkIOXML -- /m &&
			cmake --build . --config Debug --target vtkIOXML -- /m)
	ENDIF()
	SET(CATALYST_CONFIGURE_COMMAND cmake.bat)
ELSE()
	IF($ENV{CI})
		SET(CATALYST_MAKE_COMMAND make vtkIOXML)
	ELSE()
		SET(CATALYST_MAKE_COMMAND make -j ${NUM_PROCESSORS} vtkIOXML)
	ENDIF()
	SET(CATALYST_CONFIGURE_COMMAND cmake.sh)
ENDIF()

ExternalProject_Add(Catalyst
	PREFIX ${CMAKE_BINARY_DIR}/External/catalyst
	GIT_REPOSITORY ${CATALYST_GIT_URL}
	#URL ${OGS_VTK_URL}
	#URL_MD5 ${OGS_VTK_MD5}
	CONFIGURE_COMMAND ../Catalyst/${CATALYST_CONFIGURE_COMMAND} -G ${CATALYST_CMAKE_GENERATOR} ../Catalyst
	BUILD_COMMAND ${CATALYST_MAKE_COMMAND}
	INSTALL_COMMAND ""
)

IF(NOT ${ParaView_FOUND})
	# Rerun cmake in initial build
	ADD_CUSTOM_TARGET(VtkRescan ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} DEPENDS Catalyst)
ELSE()
	ADD_CUSTOM_TARGET(VtkRescan) # dummy target for caching
ENDIF()

