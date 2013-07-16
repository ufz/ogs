INCLUDE(ExternalProject)
# Set Boost version and which libraries to compile
SET(OGS_VTK_VERSION 6.0.0)
SET(OGS_VTK_REQUIRED_LIBS
	vtkIOXML
)
IF(OGS_BUILD_GUI)
	SET(OGS_VTK_REQUIRED_LIBS ${OGS_VTK_REQUIRED_LIBS}
		vtkRenderingCore
		vtkGUISupportQt
		vtkInteractionWidgets
		vtkFiltersTexture
		vtkIONetCDF
		vtkIOLegacy
		vtkIOExport
	)
ENDIF()

IF(VTK_FOUND)
	INCLUDE( ${VTK_USE_FILE} )
	RETURN()
ENDIF()

IF(NOT DEFINED VTK_DIR)
	SET(VTK_DIR ${CMAKE_BINARY_DIR}/External/vtk/src/Vtk-build)
ENDIF()
FIND_PACKAGE(VTK ${OGS_VTK_VERSION} COMPONENTS ${OGS_VTK_REQUIRED_LIBS} NO_MODULE)
IF(VTK_FOUND)
	INCLUDE( ${VTK_USE_FILE} )
	RETURN()
ENDIF()

# Set archive sources
SET(VTK_ARCHIVE_MD5 72ede4812c90bdc55172702f0cad02bb)
SET(VTK_URL "http://www.vtk.org/files/release/6.0/vtk-${OGS_VTK_VERSION}.tar.gz")

IF(OGS_BUILD_GUI)
	SET(OGS_VTK_CMAKE_ARGS "-DVTK_Group_Qt:BOOL=ON")
ENDIF()

ExternalProject_Add(VTK
	PREFIX ${CMAKE_BINARY_DIR}/External/vtk
	URL /Users/bilke/tmp/vtk-6.0.0.tar.gz #${VTK_URL}
	URL_MD5 ${VTK_ARCHIVE_MD5}
	CMAKE_ARGS
		-DBUILD_TESTING:BOOL=OFF
		-DCMAKE_BUILD_TYPE:STRING=Release
		${OGS_VTK_CMAKE_ARGS}
	BUILD_COMMAND make -j ${NUM_PROCESSORS} ${OGS_VTK_REQUIRED_LIBS}
	INSTALL_COMMAND ""
)


IF(NOT ${VTK_FOUND})
	# Rerun cmake in initial build
	ADD_CUSTOM_TARGET(VtkRescan ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR} DEPENDS VTK)
ELSE()
	ADD_CUSTOM_TARGET(VtkRescan) # dummy target for caching
ENDIF()
