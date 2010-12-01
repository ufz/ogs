# Source files
SET( SOURCES
	GridAdapter.cpp
	VtkAddFilterDialog.cpp
	VtkApplyColorTableFilter.cpp
	VtkColorByHeightFilter.cpp
	VtkColorLookupTable.cpp
	VtkMeshSource.cpp
	VtkTextureOnSurfaceFilter.cpp
	VtkPolylinesSource.cpp
	VtkPointsSource.cpp
	VtkStationSource.cpp
	VtkSurfacesSource.cpp
	VtkVisPipeline.cpp
	VtkVisPipelineItem.cpp
	VtkVisPipelineView.cpp
	VtkVisTabWidget.cpp
	VtkVisHelper.cpp
	VtkTrackedCamera.cpp
	VtkGeoImageSource.cpp
	VtkBGImageSource.cpp
	VtkImageDataToLinePolyDataFilter.cpp
	VtkAlgorithmPropertyLineEdit.cpp
	VtkCompositeFilter.cpp
	VtkCompositeImageToCylindersFilter.cpp
	VtkAlgorithmPropertyCheckbox.cpp
	VtkFilterFactory.cpp
	VtkCompositeLineToTubeFilter.cpp
	VtkCompositePointToGlyphFilter.cpp
	VtkCompositeColorByHeightFilter.cpp
	VtkCompositeColormapToImageFilter.cpp
	VtkCompositeTextureOnSurfaceFilter.cpp
	VtkCompositeThresholdFilter.cpp
	VtkAlgorithmPropertyVectorEdit.cpp
	QVtkDataSetMapper.cpp
)

# Moc headers
SET( MOC_HEADERS
	VtkAddFilterDialog.h
	VtkVisPipeline.h
	VtkVisPipelineView.h
	VtkVisTabWidget.h
	VtkTrackedCamera.h
	VtkAlgorithmPropertyLineEdit.h
	VtkAlgorithmPropertyCheckbox.h
	VtkAlgorithmPropertyVectorEdit.h
	VtkAlgorithmProperties.h
	QVtkDataSetMapper.h
)

# Header files
SET( HEADERS
	GridAdapter.h
	#OGSFilterInfo.h
	VtkApplyColorTableFilter.h
	VtkColorByHeightFilter.h
	VtkColorLookupTable.h
	VtkMeshSource.h
	#VtkOGSFilter.h
	VtkTextureOnSurfaceFilter.h
	VtkPolylinesSource.h
	VtkPointsSource.h
	VtkStationSource.h
	VtkSurfacesSource.h
	VtkVisHelper.h
	VtkGeoImageSource.h
	VtkBGImageSource.h
	VtkImageDataToLinePolyDataFilter.h
	VtkCompositeFilter.h
	VtkCompositeImageToCylindersFilter.h
	VtkFilterFactory.h
	VtkCompositeLineToTubeFilter.h
	VtkCompositePointToGlyphFilter.h
	VtkCompositeColorByHeightFilter.h
	VtkCompositeColormapToImageFilter.h
	VtkCompositeTextureOnSurfaceFilter.h
	VtkCompositeThresholdFilter.h
	VtkVisPipelineItem.h
)

# UI files
SET ( UIS
	VtkVisTabWidgetBase.ui
	VtkAddFilterDialogBase.ui
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
	${CMAKE_BINARY_DIR}/Qt/VtkVis
	../../Base
	../../FEM
	../../GEO
	../../MathLib
	../../MSH
	../Base
	../DataView
	../StationView
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
