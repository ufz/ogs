# Source files
set(SOURCES
    mainwindow.cpp
    mainwindow.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Img/icons.qrc
)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/Base
    ${CMAKE_CURRENT_SOURCE_DIR}/DataView
    ${CMAKE_CURRENT_SOURCE_DIR}/DataView/StratView
    ${CMAKE_CURRENT_SOURCE_DIR}/DataView/DiagramView
    ${CMAKE_CURRENT_SOURCE_DIR}/VtkVis
)

# Put moc files in a project folder
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")
file(GLOB UIS CONFIGURE_DEPENDS *.ui)
source_group("UI Files" FILES ${UIS})

# Application icon
set(APP_ICON ${SOURCE_DIR_REL}/scripts/packaging/ogs-de-icon.icns)

# Create the executable
add_executable(DataExplorer
    main.cpp
    ${SOURCES}
    ${UIS}
    ${APP_ICON}
    exe-icon.rc
)

target_link_libraries(DataExplorer
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
    logog
    ${VTK_LIBRARIES}
)


if(OGS_USE_NETCDF)
    add_definitions(-DOGS_USE_NETCDF)
    add_subdirectory(NetCdfDialog)
    target_link_libraries(DataExplorer NetCdfDialogLib)
endif()

if(NOT APPLE AND OGS_USE_CONAN)
    if(UNIX)
        target_link_libraries(DataExplorer Qt5::X11Extras)
    endif()
endif()

# Workaround for Windows conan tiff-package
if(OGS_USE_CONAN AND WIN32)
    find_package(ZLIB REQUIRED)
    target_link_libraries(DataExplorer ${ZLIB_LIBRARIES})
endif()

if(GEOTIFF_FOUND)
    target_link_libraries(DataExplorer ${GEOTIFF_LIBRARIES} )
endif()

if(MSVC)
    # Set linker flags
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:MSVCRT /IGNORE:4099")
    target_link_libraries(DataExplorer winmm)
endif()

if(VTKFBXCONVERTER_FOUND)
    target_link_libraries(DataExplorer ${VTKFBXCONVERTER_LIBRARIES})
endif()

set_property(TARGET DataExplorer PROPERTY FOLDER "DataExplorer")

if(OGS_USE_PCH)
    cotire(DataExplorer)
endif()

# ---- Installation ----
install(TARGETS DataExplorer RUNTIME DESTINATION bin COMPONENT ogs_gui)

cpack_add_component(ogs_gui
    DISPLAY_NAME "OGS Data Explorer"
    DESCRIPTION "The graphical user interface for OpenGeoSys."
    GROUP Applications
)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "DataExplorer" "OGS Data Explorer" PARENT_SCOPE)
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/DataExplorer.exe" "Data Explorer" PARENT_SCOPE)

