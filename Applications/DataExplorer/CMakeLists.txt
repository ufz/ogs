set(BUILD_SHARED_LIBS_SAVED "${BUILD_SHARED_LIBS}")
# Always built Data Explorer statically as there are cyclic dependencies
# between QtDataView and VtkVis
set(BUILD_SHARED_LIBS OFF)

if(MSVC)
    add_definitions(-DVTKGUISUPPORTQT_STATIC_DEFINE)
endif()

if(VTKFBXCONVERTER_FOUND)
    add_definitions(-DVTKFBXCONVERTER_FOUND)
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}

    # Qt generated file includes
    ${CMAKE_CURRENT_BINARY_DIR}/Base
    ${CMAKE_CURRENT_BINARY_DIR}/DataView
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/DiagramView
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/StratView
    ${CMAKE_CURRENT_BINARY_DIR}/NetCdfDialog

    # Workaround for CMake 3.8
    ${CMAKE_CURRENT_BINARY_DIR}/Base/QtBase_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/QtDataView_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/DiagramView/QtDiagramView_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/DataView/StratView/QtStratView_autogen/include
    ${CMAKE_CURRENT_BINARY_DIR}/NetCdfDialog/NetCdfDialogLib_autogen/include
)

set(CMAKE_AUTOMOC TRUE)
set(CMAKE_AUTOUIC TRUE)
set(CMAKE_AUTORCC TRUE)

# Add subprojects
add_subdirectory(Base)
add_subdirectory(DataView/StratView)
add_subdirectory(DataView)
add_subdirectory(DataView/DiagramView)
add_subdirectory(NetCdfDialog)
add_subdirectory(VtkVis)
include(DataExplorer.cmake)

set(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS_SAVED}")
