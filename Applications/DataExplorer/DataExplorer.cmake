# Source files
set(SOURCES
    mainwindow.cpp
    mainwindow.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Img/icons.qrc
)

set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(
    ${SOURCE_DIR_REL}/Applications/Utils/OGSFileConverter
    ${SOURCE_DIR_REL}/Applications/FileIO
    ${SOURCE_DIR_REL}/BaseLib
    ${SOURCE_DIR_REL}/MathLib
    ${SOURCE_DIR_REL}/GeoLib
    ${SOURCE_DIR_REL}/FileIO
    ${SOURCE_DIR_REL}/MeshLib
    ${SOURCE_DIR_REL}/MeshLibGEOTOOLS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Base
    ${CMAKE_CURRENT_SOURCE_DIR}/DataView
    ${CMAKE_CURRENT_SOURCE_DIR}/DataView/StratView
    ${CMAKE_CURRENT_SOURCE_DIR}/DataView/DiagramView
    ${CMAKE_CURRENT_SOURCE_DIR}/VtkVis
    ${CMAKE_CURRENT_SOURCE_DIR}/VtkAct

    # Qt generated file includes
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/../Utils/OGSFileConverter
    ${CMAKE_CURRENT_BINARY_DIR}/DataView
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/DiagramView
    ${CMAKE_CURRENT_BINARY_DIR}/VtkVis

    # Workaround for CMake 3.8
    ${CMAKE_CURRENT_BINARY_DIR}/../Utils/OGSFileConverter/OGSFileConverterLib_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/DiagramView/QtDiagramView_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/QtDataView_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/VtkVis/VtkVis_autogen/include
)

# Put moc files in a project folder
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")
file(GLOB UIS *.ui)
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
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    Qt5::Xml
    Qt5::Network
    ApplicationsLib
    ApplicationsFileIO
    DataHolderLib
    QtDataView
    QtDiagramView
    QtStratView
    VtkVis
    VtkAct
    Threads::Threads
    OGSFileConverterLib
    vtkGUISupportQt
    logog
)

# Workaround for Windows conan tiff-package
if(USE_CONAN AND WIN32)
    find_package(ZLIB REQUIRED)
    target_link_libraries(DataExplorer ${ZLIB_LIBRARIES})
endif()

if(CMAKE_CROSSCOMPILING)
    target_link_libraries(DataExplorer
        ${QT_XML_DEPS_LIBRARIES}
        ${QT_GUI_DEPS_LIBRARIES}
        ${QT_NETWORK_DEPS_LIBRARIES}
    )
endif()

if(VTK_NETCDF_FOUND)
    target_link_libraries(DataExplorer vtkNetCDF vtkNetCDF_cxx )
endif()

if(GEOTIFF_FOUND)
    target_link_libraries(DataExplorer ${GEOTIFF_LIBRARIES} )
endif()

add_dependencies (DataExplorer VtkVis)

if(MSVC)
    # Set linker flags
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:MSVCRT /IGNORE:4099")
    target_link_libraries(DataExplorer winmm)
endif()

if(VTKFBXCONVERTER_FOUND)
    target_link_libraries(DataExplorer ${VTKFBXCONVERTER_LIBRARIES})
endif()

ADD_VTK_DEPENDENCY(DataExplorer)

set_property(TARGET DataExplorer PROPERTY FOLDER "DataExplorer")

cotire(DataExplorer)

####################
### Installation ###
####################
install(TARGETS DataExplorer RUNTIME DESTINATION bin COMPONENT ogs_gui)

cpack_add_component(ogs_gui
    DISPLAY_NAME "OGS Data Explorer"
    DESCRIPTION "The graphical user interface for OpenGeoSys."
    GROUP Applications
)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "DataExplorer" "OGS Data Explorer" PARENT_SCOPE)
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/DataExplorer.exe" "Data Explorer" PARENT_SCOPE)
