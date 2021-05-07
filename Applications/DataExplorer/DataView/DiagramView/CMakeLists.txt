set(SOURCES
    DetailWindow.cpp
    DiagramList.cpp
    DiagramPrefsDialog.cpp
    DiagramScene.cpp
    DiagramView.cpp
    QArrow.cpp
    QGraphicsGrid.cpp
)

set(HEADERS
    DiagramList.h
    GetDateTime.h
    QArrow.h
    QGraphicsGrid.h
    DiagramScene.h
    DetailWindow.h
    DiagramPrefsDialog.h
    DiagramView.h
)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
file(GLOB_RECURSE UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

ogs_add_library(QtDiagramView ${SOURCES} ${HEADERS} ${UIS})
target_include_directories(QtDiagramView PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
    QtDiagramView PRIVATE BaseLib GeoLib QtBase DataHolderLib Qt5::Gui
                          spdlog::spdlog
)

add_autogen_include(QtDiagramView)

set_property(TARGET QtDiagramView PROPERTY FOLDER "DataExplorer")
