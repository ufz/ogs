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

# First check for system boost
IF(NOT Boost_INCLUDE_DIRS)
	FIND_PACKAGE(VTK ${OGS_VTK_VERSION} REQUIRED COMPONENTS ${OGS_VTK_REQUIRED_LIBS} NO_MODULE)
	IF(VTK_FOUND)
		INCLUDE( ${VTK_USE_FILE} )
		RETURN()
	ENDIF()
ENDIF()

# Set archive sources
SET(VTK_ARCHIVE_MD5 72ede4812c90bdc55172702f0cad02bb)
SET(VTK_URL "http://www.vtk.org/files/release/6.0/vtk-${OGS_VTK_VERSION}.tar.gz")

IF(OGS_BUILD_GUI)
	SET(OGS_VTK_CMAKE_ARGS "-DVTK_Group_Qt:BOOL=ON")
ENDIF()

ExternalProject_Add(Vtk
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

set(VTK_DIR ${CMAKE_BINARY_DIR}/External/vtk/src/Vtk-build)
#ExternalProject_Get_Property( Vtk CMAKE_BUILD_TYPE )
#MESSAGE("${VTK_DEFINITIONS}")
#INCLUDE( ${CMAKE_BINARY_DIR}/External/vtk/src/Vtk/CMake/UseVTK.cmake )
#SET(VTK_USE_FILE "${CMAKE_BINARY_DIR}/External/vtk/src/Vtk/CMake/UseVTK.cmake")
