set(SOURCES StratBar.cpp StratScene.cpp StratView.cpp StratWindow.cpp)

set(HEADERS StratBar.h StratScene.h StratView.h StratWindow.h)

file(GLOB UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

# TODO: Workaround for #3172.
set(CMAKE_AUTOUIC NO)
qt5_wrap_ui(genUiHeaders ${UIS})

ogs_add_library(QtStratView ${SOURCES} ${HEADERS} ${UIS})
target_include_directories(
    QtStratView PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(
    QtStratView PRIVATE BaseLib GeoLib DataHolderLib QtBase Qt5::Gui
)
add_autogen_include(QtStratView)

set_property(TARGET QtStratView PROPERTY FOLDER "DataExplorer")
