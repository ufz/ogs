set(SOURCES StratBar.cpp StratScene.cpp StratView.cpp StratWindow.cpp)

set(HEADERS StratBar.h StratScene.h StratView.h StratWindow.h)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../../../..)
set(GUI_SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
file(GLOB UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

ogs_add_library(QtStratView ${SOURCES} ${HEADERS} ${UIS})
target_include_directories(QtStratView PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
    QtStratView PRIVATE BaseLib GeoLib DataHolderLib QtBase Qt5::Gui
)
add_autogen_include(QtStratView)

set_property(TARGET QtStratView PROPERTY FOLDER "DataExplorer")
