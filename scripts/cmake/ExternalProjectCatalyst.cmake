INCLUDE(ThirdPartyLibVersions)
INCLUDE(ExternalProject)

SET(CATALYST_GIT_URL https://github.com/ufz/catalyst-io.git)

IF(NOT DEFINED ParaView_DIR AND DEFINED ENV{ParaView_DIR})
	SET(ParaView_DIR $ENV{ParaView_DIR})
ENDIF()

# CLI modules
SET(PARAVIEW_MODULES vtkIOXML)

# GUI modules
IF(OGS_BUILD_GUI)
	SET(PARAVIEW_MODULES ${PARAVIEW_MODULES}
		vtkRenderingCore
		vtkRenderingOpenGL
		vtknetcdf
		vtkIOLegacy
		vtkIOImage
		vtkGUISupportQt
		vtkRenderingAnnotation
		vtkFiltersExtraction
		vtkFiltersGeometry
		vtkFiltersTexture
		vtkFiltersModeling
		vtkFiltersSources
		vtkImagingCore
		vtkInteractionWidgets
		vtkInteractionStyle
		vtkIOExport
		vtkRenderingFreeType
	)
	SET(CATALYST_GIT_URL https://github.com/ufz/catalyst-gui.git)
ENDIF()

SET(CATALYST_LIBRARIES ${PARAVIEW_MODULES} CACHE STRING "" FORCE)
IF(OGS_BUILD_GUI)
	# Replace vtknetcdf with vtkNetCDF vtkNetCDF_cxx
	LIST(REMOVE_ITEM CATALYST_LIBRARIES vtknetcdf)
	LIST(APPEND CATALYST_LIBRARIES vtkNetCDF vtkNetCDF_cxx)
ENDIF()

FIND_PACKAGE(ParaView 4.2 COMPONENTS ${PARAVIEW_MODULES} NO_MODULE QUIET)

FIND_LIBRARY(VTKIO_LIB_FOUND vtkIOXML-pv4.2 HINTS ${ParaView_DIR}/lib PATH_SUFFIXES Release Debug)
IF(ParaView_FOUND AND VTKIO_LIB_FOUND)
	FOREACH(DIR ${PARAVIEW_INCLUDE_DIRS})
		IF("${DIR}" MATCHES ".*vtknetcdf.*")
			INCLUDE_DIRECTORIES(SYSTEM ${DIR}/../cxx ${DIR}/include)
		ENDIF()
	ENDFOREACH()
	MESSAGE(STATUS "Using ParaView in ${ParaView_DIR}")
	RETURN()
ELSEIF(NOT ParaView_DIR)
	# If ParaView was not found check for VTK
	FIND_PACKAGE(VTK 6.1 COMPONENTS ${PARAVIEW_MODULES} NO_MODULE QUIET)
	IF(VTK_FOUND)
		INCLUDE( ${VTK_USE_FILE} )
		FOREACH(DIR ${VTK_INCLUDE_DIRS})
			IF("${DIR}" MATCHES ".*vtknetcdf.*")
				INCLUDE_DIRECTORIES(SYSTEM ${DIR}/../cxx ${DIR}/include)
			ELSEIF("${DIR}" MATCHES ".*vtk.*")
				INCLUDE_DIRECTORIES(SYSTEM ${DIR}/vtknetcdf/include)
			ENDIF()
		ENDFOREACH()
		INCLUDE_DIRECTORIES(SYSTEM ${VTK_DIR}/../ThirdParty/netcdf/vtknetcdf/cxx)
		MESSAGE(STATUS "Using VTK in ${VTK_DIR}")
		RETURN()
	ELSE()
		# If nothing was found build ParaView as an external project
		SET(ParaView_DIR ${CMAKE_BINARY_DIR}/External/catalyst/src/Catalyst-build CACHE PATH "" FORCE)
	ENDIF()
ENDIF()

SET(CATALYST_CMAKE_GENERATOR ${CMAKE_GENERATOR})
IF(WIN32)
	# Ninja temporary disabled because it builds only the Release mode.
	# FIND_PROGRAM(NINJA_TOOL_PATH ninja DOC "Ninja build tool")
	IF(NINJA_TOOL_PATH)
		SET(CATALYST_CMAKE_GENERATOR Ninja)
		SET(CATALYST_MAKE_COMMAND ninja ${CATALYST_LIBRARIES})
	ELSE()
		SET(CATALYST_MAKE_COMMAND
			msbuild /p:Configuration=Release /m:${NUM_PROCESSORS} ParaView.sln &&
			msbuild /p:Configuration=Debug /m:${NUM_PROCESSORS} /m ParaView.sln)
	ENDIF()
	SET(CATALYST_CONFIGURE_COMMAND cmake.bat)
ELSE()
	IF($ENV{CI})
		SET(CATALYST_MAKE_COMMAND make ${CATALYST_LIBRARIES})
	ELSE()
		SET(CATALYST_MAKE_COMMAND make -j ${NUM_PROCESSORS} ${CATALYST_LIBRARIES})
	ENDIF()
	SET(CATALYST_CONFIGURE_COMMAND cmake.sh)
ENDIF()

MESSAGE(STATUS "Building ParaView as an external project in the build directory")
IF(CMAKE_VERSION VERSION_LESS 3.0.0)
	MESSAGE(FATAL_ERROR "CMake 3.0.0 or higher is required for building VTK / ParaView!")
ENDIF()
ExternalProject_Add(Catalyst
	PREFIX ${CMAKE_BINARY_DIR}/External/catalyst
	GIT_REPOSITORY ${CATALYST_GIT_URL}
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

