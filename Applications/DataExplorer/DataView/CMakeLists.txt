set(SOURCES
    AddFaultsToVoxelGridDialog.cpp
    AddLayerToMeshDialog.cpp
    ColorTableModel.cpp
    ColorTableView.cpp
    CondFromRasterDialog.cpp
    CreateStructuredGridDialog.cpp
    DataExplorerSettingsDialog.cpp
    DirectConditionGenerator.cpp
    ElementTreeModel.cpp
    ElementTreeView.cpp
    FemConditionModel.cpp
    FemConditionView.cpp
    GEOModels.cpp
    GeoOnMeshMappingDialog.cpp
    GeoTabWidget.cpp
    GeoTreeModel.cpp
    GeoTreeView.cpp
    GMSHPrefsDialog.cpp
    Layers2GridDialog.cpp
    LicenseDialog.cpp
    LinearEditDialog.cpp
    LineEditDialog.cpp
    MergeGeometriesDialog.cpp
    MeshAnalysisDialog.cpp
    MeshElementRemovalDialog.cpp
    MeshItem.cpp
    MeshLayerEditDialog.cpp
    MeshMapping2DDialog.cpp
    MeshModel.cpp
    MeshQualitySelectionDialog.cpp
    MeshTabWidget.cpp
    MeshValueEditDialog.cpp
    MeshView.cpp
    ModellingTabWidget.cpp
    ModelTreeItem.cpp
    ProcessModel.cpp
    ProcessView.cpp
    RasterDataToMeshDialog.cpp
    SaveMeshDialog.cpp
    SelectMeshDialog.cpp
    SetNameDialog.cpp
    SHPImportDialog.cpp
    StationTabWidget.cpp
    StationTreeModel.cpp
    StationTreeView.cpp
    SurfaceExtractionDialog.cpp
    TranslateDataDialog.cpp
    Vtu2GridDialog.cpp
)

set(HEADERS
    AddFaultsToVoxelGridDialog.h
    AddLayerToMeshDialog.h
    BaseItem.h
    ColorTableModel.h
    ColorTableView.h
    CondFromRasterDialog.h
    CondItem.h
    CreateStructuredGridDialog.h
    DataExplorerSettingsDialog.h
    DirectConditionGenerator.h
    ElementTreeModel.h
    ElementTreeView.h
    FemConditionModel.h
    FemConditionView.h
    GEOModels.h
    GeoObjectListItem.h
    GeoOnMeshMappingDialog.h
    GeoTabWidget.h
    GeoTreeItem.h
    GeoTreeModel.h
    GeoTreeView.h
    GMSHPrefsDialog.h
    Layers2GridDialog.h
    LicenseDialog.h
    LinearEditDialog.h
    LineEditDialog.h
    MergeGeometriesDialog.h
    MeshAnalysisDialog.h
    MeshElementRemovalDialog.h
    MeshItem.h
    MeshLayerEditDialog.h
    MeshMapping2DDialog.h
    MeshModel.h
    MeshQualitySelectionDialog.h
    MeshTabWidget.h
    MeshValueEditDialog.h
    MeshView.h
    ModellingTabWidget.h
    ModelTreeItem.h
    ProcessModel.h
    ProcessVarItem.h
    ProcessView.h
    RasterDataToMeshDialog.h
    SaveMeshDialog.h
    SelectMeshDialog.h
    SetNameDialog.h
    SHPImportDialog.h
    StationTabWidget.h
    StationTreeModel.h
    StationTreeView.h
    SurfaceExtractionDialog.h
    Vtu2GridDialog.h
)

# Visual Studio folder
source_group("Dialog Header Files" REGULAR_EXPRESSION "[.]*Dialog.h")
source_group("Dialog Source Files" REGULAR_EXPRESSION "[.]*Dialog.cpp")
source_group(
    "Data Model Header Files"
    REGULAR_EXPRESSION "[.]*Item.h|[.]*Model.h|[.]*View.h"
)
source_group(
    "Data Model Source Files"
    REGULAR_EXPRESSION "[.]*Item.cpp|[.]*Model.cpp|[.]*View.cpp"
)
file(GLOB UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/..)

if(GEOTIFF_FOUND)
    include_directories(${GEOTIFF_INCLUDE_DIRS})
endif() # GEOTIFF_FOUND

# TODO: Workaround for #3172.
set(CMAKE_AUTOUIC NO)
qt5_wrap_ui(genUiHeaders ${UIS})

ogs_add_library(QtDataView ${SOURCES} ${HEADERS} ${UIS} ${genUiHeaders})
target_include_directories(
    QtDataView PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)
target_link_libraries(
    QtDataView
    PUBLIC GeoLib MeshGeoToolsLib MeshLib QtBase VtkVis
    PRIVATE BaseLib
            MathLib
            ApplicationsFileIO
            QtDiagramView
            QtStratView
            DataHolderLib
            Qt5::Core
            Qt5::Gui
)

if(GEOTIFF_FOUND)
    target_compile_definitions(QtDataView PRIVATE GEOTIFF_FOUND)
    target_link_libraries(QtDataView PRIVATE ${GEOTIFF_LIBRARIES})
endif() # GEOTIFF_FOUND

add_autogen_include(QtDataView)

set_property(TARGET QtDataView PROPERTY FOLDER "DataExplorer")
