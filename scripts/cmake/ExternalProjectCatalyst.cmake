include(ThirdPartyLibVersions)
include(ExternalProject)

set(CATALYST_GIT_URL https://github.com/ufz/catalyst-io.git)

if(NOT DEFINED ParaView_DIR AND DEFINED ENV{ParaView_DIR})
	set(ParaView_DIR $ENV{ParaView_DIR})
endif()

# CLI modules
set(PARAVIEW_MODULES vtkIOXML)

# GUI modules
if(OGS_BUILD_GUI)
	set(PARAVIEW_MODULES ${PARAVIEW_MODULES}
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
	set(CATALYST_GIT_URL https://github.com/ufz/catalyst-gui.git)
endif()

set(CATALYST_LIBRARIES ${PARAVIEW_MODULES} CACHE STRING "" FORCE)
if(OGS_BUILD_GUI)
	# Replace vtknetcdf with vtkNetCDF vtkNetCDF_cxx
	list(REMOVE_ITEM CATALYST_LIBRARIES vtknetcdf)
	list(APPEND CATALYST_LIBRARIES vtkNetCDF vtkNetCDF_cxx)
endif()

find_package(ParaView 4.2 COMPONENTS ${PARAVIEW_MODULES} NO_MODULE QUIET)

find_library(VTKIO_LIB_FOUND vtkIOXML-pv4.2 HINTS ${ParaView_DIR}/lib PATH_SUFFIXES Release Debug)
if(ParaView_FOUND AND VTKIO_LIB_FOUND)
	foreach(DIR ${PARAVIEW_INCLUDE_DIRS})
		if("${DIR}" MATCHES ".*vtknetcdf.*")
			include_directories(SYSTEM ${DIR}/../cxx ${DIR}/include)
		endif()
	endforeach()
	message(STATUS "Using ParaView in ${ParaView_DIR}")
	return()
elseif(NOT ParaView_DIR)
		# If nothing was found build ParaView as an external project
		set(ParaView_DIR ${CMAKE_BINARY_DIR}/External/catalyst/src/Catalyst-build CACHE PATH "" FORCE)
	endif()
endif()

set(CATALYST_CMAKE_GENERATOR ${CMAKE_GENERATOR})
if(WIN32)
	# Ninja temporary disabled because it builds only the Release mode.
	# find_program(NINJA_TOOL_PATH ninja DOC "Ninja build tool")
	if(NINJA_TOOL_PATH)
		set(CATALYST_CMAKE_GENERATOR Ninja)
		set(CATALYST_MAKE_COMMAND ninja ${CATALYST_LIBRARIES})
	else()
		set(CATALYST_MAKE_COMMAND
			msbuild /p:Configuration=Release /m:${NUM_PROCESSORS} ParaView.sln &&
			msbuild /p:Configuration=Debug /m:${NUM_PROCESSORS} /m ParaView.sln)
	endif()
	set(CATALYST_CONFIGURE_COMMAND cmake.bat)
else()
	if($ENV{CI})
		set(CATALYST_MAKE_COMMAND make ${CATALYST_LIBRARIES})
	else()
		set(CATALYST_MAKE_COMMAND make -j ${NUM_PROCESSORS} ${CATALYST_LIBRARIES})
	endif()
	set(CATALYST_CONFIGURE_COMMAND cmake.sh)
endif()

message(STATUS "Building ParaView as an external project in the build directory")

ExternalProject_Add(Catalyst
	PREFIX ${CMAKE_BINARY_DIR}/External/catalyst
	GIT_REPOSITORY ${CATALYST_GIT_URL}
	CONFIGURE_COMMAND ../Catalyst/${CATALYST_CONFIGURE_COMMAND} -G ${CATALYST_CMAKE_GENERATOR}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ../Catalyst
	BUILD_COMMAND ${CATALYST_MAKE_COMMAND}
	INSTALL_COMMAND ""
)

if(NOT ${ParaView_FOUND})
	# Rerun cmake in initial build
	add_custom_target(VtkRescan ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} DEPENDS Catalyst)
else()
	add_custom_target(VtkRescan) # dummy target for caching
endif()

