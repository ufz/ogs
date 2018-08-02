set(SOURCES
    MeshFromRasterDialog.cpp
    QVtkDataSetMapper.cpp
    VisPrefsDialog.cpp
    VisualizationWidget.cpp
    VtkAddFilterDialog.cpp
    VtkAlgorithmProperties.cpp
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
    VtkCompositeImageToPointCloudFilter.cpp
    VtkCompositeImageToSurfacePointsFilter.cpp
    VtkCompositeLineToTubeFilter.cpp
    VtkCompositeNodeSelectionFilter.cpp
    VtkCompositePointToGlyphFilter.cpp
    VtkCompositeTextureOnSurfaceFilter.cpp
    VtkCompositeThresholdFilter.cpp
    VtkConsoleOutputWindow.cpp
    VtkCustomInteractorStyle.cpp
    VtkFilterFactory.cpp
    VtkGeoImageSource.cpp
    VtkImageDataToLinePolyDataFilter.cpp
    VtkImageDataToPointCloudFilter.cpp
    VtkImageDataToSurfacePointsFilter.cpp
    VtkPickCallback.cpp
    VtkPolylinesSource.cpp
    VtkPointsSource.cpp
    VtkRaster.cpp
    VtkStationSource.cpp
    VtkSurfacesSource.cpp
    VtkTextureOnSurfaceFilter.cpp
    VtkVisHelper.cpp
    VtkVisImageItem.cpp
    VtkVisPipeline.cpp
    VtkVisPipelineItem.cpp
    VtkVisPipelineView.cpp
    VtkVisPointSetItem.cpp
    VtkVisTabWidget.cpp
)

set(HEADERS
    MeshFromRasterDialog.h
    QVtkDataSetMapper.h
    VisPrefsDialog.h
    VisualizationWidget.h
    VtkAddFilterDialog.h
    VtkAlgorithmProperties.h
    VtkAlgorithmPropertyLineEdit.h
    VtkAlgorithmPropertyCheckbox.h
    VtkAlgorithmPropertyVectorEdit.h
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
    VtkCompositeImageToPointCloudFilter.h
    VtkCompositeImageToSurfacePointsFilter.h
    VtkCompositeLineToTubeFilter.h
    VtkCompositeNodeSelectionFilter.h
    VtkCompositePointToGlyphFilter.h
    VtkCompositeTextureOnSurfaceFilter.h
    VtkCompositeThresholdFilter.h
    VtkConsoleOutputWindow.h
    VtkCustomInteractorStyle.h
    VtkFilterFactory.h
    VtkGeoImageSource.h
    VtkImageDataToLinePolyDataFilter.h
    VtkImageDataToPointCloudFilter.h
    VtkImageDataToSurfacePointsFilter.h
    VtkPickCallback.h
    VtkPolylinesSource.h
    VtkPointsSource.h
    VtkRaster.h
    VtkStationSource.h
    VtkSurfacesSource.h
    VtkTextureOnSurfaceFilter.h
    VtkVisHelper.h
    VtkVisImageItem.h
    VtkVisPipeline.h
    VtkVisPipelineItem.h
    VtkVisPipelineView.h
    VtkVisPointSetItem.h
    VtkVisTabWidget.h
)

# Visual Studio folder
source_group("Filter Header Files" REGULAR_EXPRESSION "[.]*Filter.h")
source_group("Filter Source Files" REGULAR_EXPRESSION "[.]*Filter.cpp")
source_group("Source Header Files" REGULAR_EXPRESSION "[.]*Source.h")
source_group("Source Source Files" REGULAR_EXPRESSION "[.]*Source.cpp")
file(GLOB_RECURSE UIS *.ui)
source_group("UI Files" FILES ${UIS})

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/..)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SOURCE_DIR_REL}/Applications/FileIO
    ${SOURCE_DIR_REL}/BaseLib
    ${SOURCE_DIR_REL}/GeoLib
    ${SOURCE_DIR_REL}/MathLib
    ${SOURCE_DIR_REL}/MeshLib
    ${GUI_SOURCE_DIR_REL}/Base
    ${GUI_SOURCE_DIR_REL}/DataView
    ${GUI_SOURCE_DIR_REL}/VtkModules/Qt
)

add_library(VtkVis
    ${SOURCES}
    ${HEADERS}
    ${UIS}
)

if(GEOTIFF_FOUND)
    include_directories(${GEOTIFF_INCLUDE_DIRS})
    target_link_libraries(VtkVis INTERFACE ${GEOTIFF_LIBRARIES})
endif() # GEOTIFF_FOUND

target_link_libraries(VtkVis
    PUBLIC BaseLib GeoLib MeshLib DataHolderLib QtBase vtkIOImage
    PRIVATE MathLib ApplicationsFileIO Qt5::Gui logog NetCdfDialogLib
)
set_property(TARGET VtkVis PROPERTY FOLDER "DataExplorer")

if(OGS_USE_PCH)
    cotire(VtkVis)
endif()
