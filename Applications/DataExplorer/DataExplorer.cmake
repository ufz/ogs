# Source files
set(SOURCES mainwindow.cpp mainwindow.h
            ${CMAKE_CURRENT_SOURCE_DIR}/Img/icons.qrc
)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Put moc files in a project folder
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")
file(GLOB UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

# Application icon
set(APP_ICON ${SOURCE_DIR_REL}/scripts/packaging/ogs-de-icon.icns)

# Create the executable
ogs_add_executable(
    DataExplorer main.cpp ${SOURCES} ${UIS} ${APP_ICON} exe-icon.rc
)
target_compile_definitions(
    DataExplorer PUBLIC $<$<BOOL:${VTK_ADDED}>:VTK_VIA_CPM>
)

target_link_libraries(
    DataExplorer
    BaseLib
    GeoLib
    GitInfoLib
    MeshLib
    ApplicationsFileIO
    DataHolderLib
    OGSFileConverterLib
    QtBase
    QtDataView
    QtDiagramView
    VtkVis
    Threads::Threads
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    Qt5::Xml
    Qt5::Network
    spdlog::spdlog
    ${VTK_LIBRARIES}
)

if(UNIX AND NOT APPLE)
    target_link_libraries(DataExplorer Qt5::X11Extras)
endif()

if(GEOTIFF_FOUND)
    target_link_libraries(DataExplorer ${GEOTIFF_LIBRARIES})
endif()

if(MSVC)
    # Set linker flags
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG
        "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:MSVCRT /IGNORE:4099"
    )
    target_link_libraries(DataExplorer winmm)
endif()

set_property(TARGET DataExplorer PROPERTY FOLDER "DataExplorer")

# ---- Installation ----
install(TARGETS DataExplorer RUNTIME DESTINATION bin)

cpack_add_component(
    ogs_gui
    DISPLAY_NAME "OGS Data Explorer"
    DESCRIPTION "The graphical user interface for OpenGeoSys."
    GROUP Applications
)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "DataExplorer"
                              "OGS Data Explorer" PARENT_SCOPE
)
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/DataExplorer.exe"
                          "Data Explorer" PARENT_SCOPE
)
