# Source files
SET( SOURCES
	MeshFromRasterDialog.cpp
	QVtkDataSetMapper.cpp
	VisPrefsDialog.cpp
	VisualizationWidget.cpp
	VtkAddFilterDialog.cpp
	VtkAlgorithmPropertyCheckbox.cpp
	VtkAlgorithmPropertyLineEdit.cpp
	VtkAlgorithmPropertyVectorEdit.cpp
	VtkAppendArrayFilter.cpp
	VtkBGImageSource.cpp
	VtkColorByHeightFilter.cpp
	VtkColorLookupTable.cpp
	VtkCompositeFilter.cpp
	VtkCompositeColorByHeightFilter.cpp
	VtkCompositeColormapToImageFilter.cpp
	VtkCompositeContourFilter.cpp
	VtkCompositeElementSelectionFilter.cpp
	VtkCompositeGeoObjectFilter.cpp
	VtkCompositeImageToCylindersFilter.cpp
	VtkCompositeLineToTubeFilter.cpp
	VtkCompositeNodeSelectionFilter.cpp
	VtkCompositePointToGlyphFilter.cpp
	VtkCompositeTextureOnSurfaceFilter.cpp
	VtkCompositeThresholdFilter.cpp
	VtkFilterFactory.cpp
	VtkGeoImageSource.cpp
	VtkImageDataToLinePolyDataFilter.cpp
	VtkMeshConverter.cpp
	VtkMeshSource.cpp
	VtkPolylinesSource.cpp
	VtkPointsSource.cpp
	VtkRaster.cpp
	VtkStationSource.cpp
	VtkSurfacesSource.cpp
	VtkTextureOnSurfaceFilter.cpp
	VtkTrackedCamera.cpp
	VtkVisHelper.cpp
	VtkVisImageItem.cpp
	VtkVisPipeline.cpp
	VtkVisPipelineItem.cpp
	VtkVisPipelineView.cpp
	VtkVisPointSetItem.cpp
	VtkVisTabWidget.cpp
	VtkAlgorithmProperties.cpp
)

# Moc headers
SET( MOC_HEADERS
	MeshFromRasterDialog.h
	QVtkDataSetMapper.h
	VisPrefsDialog.h
	VisualizationWidget.h
	VtkAddFilterDialog.h
	VtkAlgorithmProperties.h
	VtkAlgorithmPropertyLineEdit.h
	VtkAlgorithmPropertyCheckbox.h
	VtkAlgorithmPropertyVectorEdit.h
	VtkTrackedCamera.h
	VtkVisPipeline.h
	VtkVisPipelineView.h
	VtkVisTabWidget.h
)

# Header files
SET( HEADERS
	VtkAppendArrayFilter.h
	VtkBGImageSource.h
	VtkColorByHeightFilter.h
	VtkColorLookupTable.h
	VtkCompositeFilter.h
	VtkCompositeColorByHeightFilter.h
	VtkCompositeColormapToImageFilter.h
	VtkCompositeContourFilter.h
	VtkCompositeElementSelectionFilter.h
	VtkCompositeGeoObjectFilter.h
	VtkCompositeImageToCylindersFilter.h
	VtkCompositeLineToTubeFilter.h
	VtkCompositeNodeSelectionFilter.h
	VtkCompositePointToGlyphFilter.h
	VtkCompositeTextureOnSurfaceFilter.h
	VtkCompositeThresholdFilter.h
	VtkFilterFactory.h
	VtkGeoImageSource.h
	VtkImageDataToLinePolyDataFilter.h
	VtkMeshConverter.h
	VtkMeshSource.h
	VtkPolylinesSource.h
	VtkPointsSource.h
	VtkRaster.h
	VtkStationSource.h
	VtkSurfacesSource.h
	VtkTextureOnSurfaceFilter.h
	VtkVisHelper.h
	VtkVisImageItem.h
	VtkVisPipelineItem.h
	VtkVisPointSetItem.h
)

# UI files
SET ( UIS
	MeshFromRaster.ui
	VtkAddFilterDialogBase.ui
	VtkVisTabWidgetBase.ui
	VisualizationWidgetBase.ui
	VisPrefs.ui
)

# Put moc files in a project folder
SOURCE_GROUP("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
SOURCE_GROUP("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Put filter in a folder
SOURCE_GROUP("Filter Header Files" REGULAR_EXPRESSION "[.]*Filter.h")
SOURCE_GROUP("Filter Source Files" REGULAR_EXPRESSION "[.]*Filter.cpp")

# Put vtk source files in a folder
SOURCE_GROUP("Source Header Files" REGULAR_EXPRESSION "[.]*Source.h")
SOURCE_GROUP("Source Source Files" REGULAR_EXPRESSION "[.]*Source.cpp")

# Run Qts user interface compiler uic on .ui files
QT4_WRAP_UI( UI_HEADERS ${UIS} )

# Run Qts meta object compiler moc on header files
QT4_WRAP_CPP( MOC_SOURCES ${MOC_HEADERS} )

INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/../VtkAct
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/../../BaseLib
	${CMAKE_CURRENT_BINARY_DIR}/../../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../GeoLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../FileIO
	${CMAKE_CURRENT_SOURCE_DIR}/../../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../Base
	${CMAKE_CURRENT_SOURCE_DIR}/../DataView
	${CMAKE_CURRENT_BINARY_DIR}/../DataView
)

# Create the library
ADD_LIBRARY( VtkVis STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
)

IF(TARGET VtkRescan)
	ADD_DEPENDENCIES(VtkVis VtkRescan)
ENDIF()

IF (GEOTIFF_FOUND)
	INCLUDE_DIRECTORIES( ${GEOTIFF_INCLUDE_DIRS} )
	TARGET_LINK_LIBRARIES( VtkVis ${GEOTIFF_LIBRARIES} )
ENDIF () # GEOTIFF_FOUND

TARGET_LINK_LIBRARIES( VtkVis ${QT_LIBRARIES} )

ADD_DEPENDENCIES(VtkVis QtDataView)

IF(OGS_USE_VRPN)
	INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Vrpn)
ENDIF() # OGS_USE_VRPN

IF (VTKOSGCONVERTER_FOUND)
	USE_OPENSG(VtkVis)
	INCLUDE_DIRECTORIES( ${VTKOSGCONVERTER_INCLUDE_DIRS} )
ENDIF ()

set_property(TARGET VtkVis PROPERTY FOLDER "DataExplorer")
