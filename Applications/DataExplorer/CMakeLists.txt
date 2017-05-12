set(BUILD_SHARED_LIBS_SAVED "${BUILD_SHARED_LIBS}")
# Always built Data Explorer statically as there are cyclic dependencies
# between QtDataView and VtkVis
set(BUILD_SHARED_LIBS OFF)

if(MSVC)
    add_definitions(-DVTKGUISUPPORTQT_STATIC_DEFINE)
endif()

if(VTKFBXCONVERTER_FOUND)
    add_definitions(-DVTKFBXCONVERTER_FOUND)
    include_directories(${VTKFBXCONVERTER_INCLUDE_DIRS})
endif()

# Add subprojects
add_subdirectory(Base)
add_subdirectory(DataView/StratView)
add_subdirectory(DataView)
add_subdirectory(DataView/DiagramView)
add_subdirectory(VtkVis)
include(DataExplorer.cmake)

set(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS_SAVED}")
