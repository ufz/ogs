# ---- VtkVisFilter ----
set(FILTER_SOURCES
    VtkAlgorithmProperties.cpp
    VtkAppendArrayFilter.cpp
    VtkBGImageSource.cpp
    VtkColorByHeightFilter.cpp
    VtkColorLookupTable.cpp
    VtkCompositeColorByHeightFilter.cpp
    VtkCompositeColormapToImageFilter.cpp
    VtkCompositeContourFilter.cpp
    VtkCompositeElementSelectionFilter.cpp
    VtkCompositeFilter.cpp
    VtkCompositeGeoObjectFilter.cpp
    VtkCompositeImageToCylindersFilter.cpp
    VtkCompositeImageToPointCloudFilter.cpp
    VtkCompositeImageToSurfacePointsFilter.cpp
    VtkCompositeLineToTubeFilter.cpp
    VtkCompositeNodeSelectionFilter.cpp
    VtkCompositePointToGlyphFilter.cpp
    VtkCompositeTextureOnSurfaceFilter.cpp
    VtkCompositeThresholdFilter.cpp
    VtkGeoImageSource.cpp
    VtkImageDataToLinePolyDataFilter.cpp
    VtkImageDataToPointCloudFilter.cpp
    VtkImageDataToSurfacePointsFilter.cpp
    VtkPointsSource.cpp
    VtkPolylinesSource.cpp
    VtkRaster.cpp
    VtkStationSource.cpp
    VtkSurfacesSource.cpp
    VtkTextureOnSurfaceFilter.cpp
)
if(OGS_USE_NETCDF)
    list(APPEND FILTER_SOURCES NetCdfConfigureDialog.cpp NetCdfConfigure.ui)
    set(UIS NetCdfConfigure.ui)
endif()

ogs_add_library(VtkVisFilter ${FILTER_SOURCES})
target_include_directories(
    VtkVisFilter PUBLIC ${VTK_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}
)
target_link_libraries(VtkVisFilter PRIVATE GeoLib Qt5::Core Qt5::Widgets ${VTK_LIBRARIES})

if(OGS_USE_NETCDF)
    target_link_libraries(
        VtkVisFilter PUBLIC ${NETCDF_LIBRARIES_CXX} ${NETCDF_LIBRARIES_C}
                            ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES}
    )
    if(OGS_USE_CONAN AND MSVC)
        target_link_libraries(VtkVisFilter PUBLIC ${CONAN_LIBS})
    endif()
    # Workaround for system installed VTK (tested on arch)
    if(NOT OGS_USE_CONAN)
        target_include_directories(
            VtkVisFilter SYSTEM PUBLIC ${NETCDF_INCLUDES_C}
                                       ${NETCDF_INCLUDES_CXX}
        )
    endif()
    add_autogen_include(VtkVisFilter)
endif()

# ---- VtkVis ----
set(SOURCES
    MeshFromRasterDialog.cpp
    QVtkDataSetMapper.cpp
    VisPrefsDialog.cpp
    VisualizationWidget.cpp
    VtkAddFilterDialog.cpp
    VtkAlgorithmPropertyCheckbox.cpp
    VtkAlgorithmPropertyLineEdit.cpp
    VtkAlgorithmPropertyVectorEdit.cpp
    VtkConsoleOutputWindow.cpp
    VtkCustomInteractorStyle.cpp
    VtkFilterFactory.cpp
    VtkPickCallback.cpp
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

list(
    APPEND
    UIS
    MeshFromRaster.ui
    VisPrefs.ui
    VisualizationWidgetBase.ui
    VtkAddFilterDialogBase.ui
    VtkVisTabWidgetBase.ui
)

# Visual Studio folder
source_group("Filter Header Files" REGULAR_EXPRESSION "[.]*Filter.h")
source_group("Filter Source Files" REGULAR_EXPRESSION "[.]*Filter.cpp")
source_group("Source Header Files" REGULAR_EXPRESSION "[.]*Source.h")
source_group("Source Source Files" REGULAR_EXPRESSION "[.]*Source.cpp")

source_group("UI Files" FILES ${UIS})

# TODO: Workaround for #3172.
set(CMAKE_AUTOUIC NO)
qt5_wrap_ui(genUiHeaders ${UIS})

ogs_add_library(VtkVis ${SOURCES} ${HEADERS} ${UIS} ${genUiHeaders})
target_include_directories(VtkVis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(GEOTIFF_FOUND)
    target_compile_definitions(VtkVis PRIVATE GEOTIFF_FOUND)
    include_directories(${GEOTIFF_INCLUDE_DIRS})
    target_link_libraries(VtkVis PUBLIC ${GEOTIFF_LIBRARIES})
endif() # GEOTIFF_FOUND

target_link_libraries(
    VtkVis PUBLIC BaseLib GeoLib MeshLib DataHolderLib QtBase VtkVisFilter
    PRIVATE MathLib ApplicationsFileIO Qt5::Gui
)

set_property(TARGET VtkVis PROPERTY FOLDER "DataExplorer")

add_autogen_include(VtkVis)
