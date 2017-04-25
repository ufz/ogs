set(SOURCES
    FileListDialog.cpp
    OGSFileConverter.cpp
)

include_directories(
    ${PROJECT_SOURCE_DIR}/BaseLib
    ${PROJECT_SOURCE_DIR}/FileIO
    ${PROJECT_SOURCE_DIR}/GeoLib
    ${PROJECT_SOURCE_DIR}/MeshLib
    ${PROJECT_SOURCE_DIR}/Applications/DataExplorer/Base
    ${PROJECT_SOURCE_DIR}/Applications/DataExplorer/DataView
    ${PROJECT_SOURCE_DIR}/Applications/DataExplorer/VtkVis

    # Qt generated file includes
    ${CMAKE_CURRENT_BINARY_DIR}

    # Workaround for CMake 3.8
    ${CMAKE_CURRENT_BINARY_DIR}/OGSFileConverterLib_autogen/include
)

file(GLOB_RECURSE UIS *.ui)
source_group("UI Files" FILES ${UIS})

add_library(OGSFileConverterLib
    ${SOURCES}
    ${HEADERS}
    ${UIS}
)
target_link_libraries(OGSFileConverterLib
    QtBase
    MeshLib)

add_executable(OGSFileConverter main.cpp)

target_link_libraries(OGSFileConverter
    ApplicationsLib
    OGSFileConverterLib
    QtBase
    GeoLib)
ADD_VTK_DEPENDENCY(OGSFileConverter)

set_property(TARGET OGSFileConverter PROPERTY FOLDER "Utilities")

if(OGS_USE_PCH)
    cotire(OGSFileConverter)
endif()

####################
### Installation ###
####################

install(TARGETS OGSFileConverter RUNTIME DESTINATION bin COMPONENT ogs_cli)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "OGSFileConverter" "File Converter")
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/OGSFileConverter.exe" "File Converter" PARENT_SCOPE)
