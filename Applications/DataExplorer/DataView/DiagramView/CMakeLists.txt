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
    getDateTime.h
    QArrow.h
    QGraphicsGrid.h
    DiagramScene.h
    DetailWindow.h
    DiagramPrefsDialog.h
    DiagramView.h
)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(
    ${SOURCE_DIR_REL}/BaseLib
    ${SOURCE_DIR_REL}/GeoLib
    ${SOURCE_DIR_REL}/MathLib
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GUI_SOURCE_DIR_REL}/Base
    ${GUI_SOURCE_DIR_REL}/DataView
)

file(GLOB_RECURSE UIS *.ui)
source_group("UI Files" FILES ${UIS})

add_library(QtDiagramView
    ${SOURCES}
    ${HEADERS}
    ${UIS}
)

target_link_libraries(QtDiagramView
    PRIVATE BaseLib GeoLib QtBase DataHolderLib Qt5::Gui logog
)

set_property(TARGET QtDiagramView PROPERTY FOLDER "DataExplorer")

if(OGS_USE_PCH)
    cotire(QtDiagramView)
endif()
