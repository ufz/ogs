cmake_minimum_required(VERSION 2.8.6)

project(OGSFileConverter)

cmake_policy(SET CMP0020 NEW)

set(CMAKE_AUTOMOC TRUE)
set(CMAKE_AUTOUIC TRUE)
set(CMAKE_AUTORCC TRUE)

set(SOURCES
    FileListDialog.cpp
    OGSFileConverter.cpp
)

find_package(Qt5 QUIET REQUIRED Gui Widgets Xml XmlPatterns)

include_directories(
    ${CMAKE_BINARY_DIR}/Applications/Utils/OGSFileConverter
    ${CMAKE_SOURCE_DIR}/BaseLib
    ${CMAKE_SOURCE_DIR}/FileIO
    ${CMAKE_SOURCE_DIR}/GeoLib
    ${CMAKE_SOURCE_DIR}/MeshLib
    ${CMAKE_SOURCE_DIR}/Applications/DataExplorer/Base
    ${CMAKE_SOURCE_DIR}/Applications/DataExplorer/DataView
    ${CMAKE_SOURCE_DIR}/Applications/DataExplorer/VtkVis
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

set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE TYPE INTERNAL FORCE)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

set_property(TARGET OGSFileConverter PROPERTY FOLDER "Utilities")

####################
### Installation ###
####################

install(TARGETS OGSFileConverter RUNTIME DESTINATION bin COMPONENT ogs_cli)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "OGSFileConverter" "File Converter")
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/OGSFileConverter.exe" "File Converter" PARENT_SCOPE)
