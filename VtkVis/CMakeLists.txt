# Source files
SET( SOURCES
	GridAdapter.cpp
	QVtkDataSetMapper.cpp
	VtkAddFilterDialog.cpp
	VtkAlgorithmPropertyCheckbox.cpp
	VtkAlgorithmPropertyLineEdit.cpp
	VtkAlgorithmPropertyVectorEdit.cpp
	VtkApplyColorTableFilter.cpp
	VtkBGImageSource.cpp
	VtkColorByHeightFilter.cpp
	VtkColorLookupTable.cpp
	VtkCompositeFilter.cpp
	VtkCompositeColorByHeightFilter.cpp
	VtkCompositeColormapToImageFilter.cpp
	VtkCompositeImageToCylindersFilter.cpp
	VtkCompositeLineToTubeFilter.cpp
	VtkCompositePointToGlyphFilter.cpp
	VtkCompositeSelectionFilter.cpp	
	VtkCompositeTextureOnSurfaceFilter.cpp
	VtkCompositeThresholdFilter.cpp
	VtkConditionSource.cpp
	VtkFilterFactory.cpp
	VtkGeoImageSource.cpp
	VtkImageDataToLinePolyDataFilter.cpp
	VtkMeshSource.cpp
	VtkPolylinesSource.cpp
	VtkPointsSource.cpp
	VtkSelectionFilter.cpp
	VtkStationSource.cpp
	VtkSurfacesSource.cpp
	VtkTextureOnSurfaceFilter.cpp
	VtkTrackedCamera.cpp
	VtkVisPipeline.cpp
	VtkVisPipelineItem.cpp
	VtkVisPipelineView.cpp
	VtkVisTabWidget.cpp
	VtkVisHelper.cpp
	VisualizationWidget.cpp
)

# Moc headers
SET( MOC_HEADERS
	QVtkDataSetMapper.h
	VtkAddFilterDialog.h
	VtkAlgorithmProperties.h
	VtkAlgorithmPropertyLineEdit.h
	VtkAlgorithmPropertyCheckbox.h
	VtkAlgorithmPropertyVectorEdit.h
	VtkTrackedCamera.h
	VtkVisPipeline.h
	VtkVisPipelineView.h
	VtkVisTabWidget.h
	VisualizationWidget.h
)

# Header files
SET( HEADERS
	GridAdapter.h
	VtkApplyColorTableFilter.h
	VtkBGImageSource.h
	VtkColorByHeightFilter.h
	VtkColorLookupTable.h
	VtkCompositeFilter.h
	VtkCompositeColorByHeightFilter.h
	VtkCompositeColormapToImageFilter.h
	VtkCompositeImageToCylindersFilter.h
	VtkCompositeLineToTubeFilter.h
	VtkCompositePointToGlyphFilter.h
	VtkCompositeSelectionFilter.h
	VtkCompositeTextureOnSurfaceFilter.h
	VtkCompositeThresholdFilter.h
	VtkConditionSource.h
	VtkFilterFactory.h
	VtkGeoImageSource.h
	VtkImageDataToLinePolyDataFilter.h
	VtkMeshSource.h
	VtkPolylinesSource.h
	VtkPointsSource.h
	VtkSelectionFilter.h
	VtkStationSource.h
	VtkSurfacesSource.h
	VtkTextureOnSurfaceFilter.h
	VtkVisHelper.h
	VtkVisPipelineItem.h
)

# UI files
SET ( UIS
	VtkAddFilterDialogBase.ui
	VtkVisTabWidgetBase.ui
	VisualizationWidgetBase.ui
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
	.
	../VtkAct
	${CMAKE_BINARY_DIR}/Qt/VtkVis
	../../Base
	../../FEM
	../../GEO
	../../MathLib
	../../MSH
	../Base
	../DataView
)

# Create the library
ADD_LIBRARY( VtkVis STATIC
	${SOURCES}
	${HEADERS}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
)

TARGET_LINK_LIBRARIES( VtkVis ${QT_LIBRARIES} )

IF (OGS_USE_OPENSG)
	USE_OPENSG(VtkVis)
	INCLUDE_DIRECTORIES( ../OpenSG )
ENDIF (OGS_USE_OPENSG)
