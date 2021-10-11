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

file(GLOB UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

# TODO: Workaround for #3172.
set(CMAKE_AUTOUIC NO)
qt5_wrap_ui(genUiHeaders ${UIS})

ogs_add_library(QtDiagramView ${SOURCES} ${HEADERS} ${UIS} ${genUiHeaders})
target_include_directories(
    QtDiagramView PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(
    QtDiagramView PRIVATE BaseLib GeoLib QtBase DataHolderLib Qt5::Gui
)

add_autogen_include(QtDiagramView)

set_property(TARGET QtDiagramView PROPERTY FOLDER "DataExplorer")
