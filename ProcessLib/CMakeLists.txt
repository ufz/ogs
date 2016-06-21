include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Source files
APPEND_SOURCE_FILES(SOURCES)

add_subdirectory(GroundwaterFlow)
APPEND_SOURCE_FILES(SOURCES GroundwaterFlow)

add_subdirectory(TES)
APPEND_SOURCE_FILES(SOURCES TES)

APPEND_SOURCE_FILES(SOURCES Utils)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/../MathLib
    ${CMAKE_CURRENT_BINARY_DIR}/../NumLib
    ${CMAKE_CURRENT_BINARY_DIR}/../MaterialsLib
)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib
    MaterialsLib
    MeshGeoToolsLib
    NumLib # for shape matrices
    ${VTK_LIBRARIES}
)

ADD_VTK_DEPENDENCY(ProcessLib)

if(TARGET Eigen)
    add_dependencies(ProcessLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(ProcessLib Boost)
endif()
