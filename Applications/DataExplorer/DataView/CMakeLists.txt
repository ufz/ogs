set(SOURCES
    AddLayerToMeshDialog.cpp
    ColorTableModel.cpp
    ColorTableView.cpp
    CondFromRasterDialog.cpp
    CreateStructuredGridDialog.cpp
    DataExplorerSettingsDialog.cpp
    DirectConditionGenerator.cpp
    ElementTreeModel.cpp
    ElementTreeView.cpp
    GEOModels.cpp
    GeoOnMeshMappingDialog.cpp
    GeoTabWidget.cpp
    GeoTreeModel.cpp
    GeoTreeView.cpp
    GMSHPrefsDialog.cpp
    LicenseDialog.cpp
    LinearEditDialog.cpp
    LineEditDialog.cpp
    MergeGeometriesDialog.cpp
    MeshAnalysisDialog.cpp
    MeshElementRemovalDialog.cpp
    MeshMapping2DDialog.cpp
    MeshQualitySelectionDialog.cpp
    MeshValueEditDialog.cpp
    ModellingTabWidget.cpp
    ModelTreeItem.cpp
    MeshLayerEditDialog.cpp
    MshItem.cpp
    MshModel.cpp
    MshTabWidget.cpp
    MshView.cpp
    NetCdfConfigureDialog.cpp
    SaveMeshDialog.cpp
    SelectMeshDialog.cpp
    SetNameDialog.cpp
    SHPImportDialog.cpp
    StationTabWidget.cpp
    StationTreeModel.cpp
    StationTreeView.cpp
    SurfaceExtractionDialog.cpp
)

set(HEADERS
    AddLayerToMeshDialog.h
    BaseItem.h
    ColorTableModel.h
    ColorTableView.h
    CondFromRasterDialog.h
    CondItem.h
    CondObjectListItem.h
    CreateStructuredGridDialog.h
    DataExplorerSettingsDialog.h
    DirectConditionGenerator.h
    ElementTreeModel.h
    ElementTreeView.h
    GEOModels.h
    GeoObjectListItem.h
    GeoOnMeshMappingDialog.h
    GeoTabWidget.h
    GeoTreeItem.h
    GeoTreeModel.h
    GeoTreeView.h
    GMSHPrefsDialog.h
    LicenseDialog.h
    LinearEditDialog.h
    LineEditDialog.h
    MergeGeometriesDialog.h
    MeshAnalysisDialog.h
    MeshElementRemovalDialog.h
    MeshMapping2DDialog.h
    MeshQualitySelectionDialog.h
    MeshValueEditDialog.h
    ModellingTabWidget.h
    ModelTreeItem.h
    MeshLayerEditDialog.h
    MshItem.h
    MshModel.h
    MshTabWidget.h
    MshView.h
    NetCdfConfigureDialog.h
    SaveMeshDialog.h
    SelectMeshDialog.h
    SetNameDialog.h
    SHPImportDialog.h
    StationTabWidget.h
    StationTreeModel.h
    StationTreeView.h
    SurfaceExtractionDialog.h
)

# Visual Studio folder
source_group("Dialog Header Files" REGULAR_EXPRESSION "[.]*Dialog.h")
source_group("Dialog Source Files" REGULAR_EXPRESSION "[.]*Dialog.cpp")
source_group("Data Model Header Files" REGULAR_EXPRESSION "[.]*Item.h|[.]*Model.h|[.]*View.h")
source_group("Data Model Source Files" REGULAR_EXPRESSION "[.]*Item.cpp|[.]*Model.cpp|[.]*View.cpp")
file(GLOB UIS *.ui)
source_group("UI Files" FILES ${UIS})

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/..)
include_directories(
    ${SOURCE_DIR_REL}/Applications/FileIO
    ${SOURCE_DIR_REL}/BaseLib
    ${SOURCE_DIR_REL}/MathLib
    ${SOURCE_DIR_REL}/GeoLib
    ${SOURCE_DIR_REL}/MeshGeoToolsLib
    ${SOURCE_DIR_REL}/MeshLib
    ${GUI_SOURCE_DIR_REL}/Base
    ${GUI_SOURCE_DIR_REL}/VtkVis

    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/DiagramView
    ${CMAKE_CURRENT_SOURCE_DIR}/StratView

    # Qt generated file includes
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/DiagramView
    ${CMAKE_CURRENT_BINARY_DIR}/StratView

    # Workaround for CMake 3.8
    ${CMAKE_CURRENT_BINARY_DIR}/DiagramView/QtDiagramView_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/StratView/QtStratView_autogen/include
)

if(GEOTIFF_FOUND)
    include_directories(${GEOTIFF_INCLUDE_DIRS})
endif() # GEOTIFF_FOUND

add_library(QtDataView
    ${SOURCES}
    ${HEADERS}
    ${UIS}
)

target_link_libraries(QtDataView
    PUBLIC GeoLib MeshGeoToolsLib MeshLib QtBase VtkVis
    PRIVATE BaseLib MathLib ApplicationsFileIO QtDiagramView QtStratView
        DataHolderLib logog Qt5::Core Qt5::Gui
)

ADD_VTK_DEPENDENCY(QtDataView)

if(GEOTIFF_FOUND)
    target_link_libraries(QtDataView PRIVATE ${GEOTIFF_LIBRARIES})
endif () # GEOTIFF_FOUND

set_property(TARGET QtDataView PROPERTY FOLDER "DataExplorer")

if(OGS_USE_PCH)
    cotire(QtDataView)
endif()
