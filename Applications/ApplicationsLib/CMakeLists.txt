# Source files
GET_SOURCE_FILES(SOURCES_APPLICATIONSLIB)
set(LIB_SOURCES ${SOURCES_APPLICATIONSLIB})

if(MSVC)
    include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}/../../NumLib
        ${CMAKE_CURRENT_BINARY_DIR}/../../MathLib
        ${CMAKE_CURRENT_BINARY_DIR}/../../MaterialsLib
        ${CMAKE_CURRENT_BINARY_DIR}/../../MeshLib
    )
endif()

# Library
add_library(ApplicationsLib ${LIB_SOURCES})

target_link_libraries(ApplicationsLib
    MeshGeoToolsLib
    NumLib
    ProcessLib
    ${VTK_LIBRARIES}
)
ADD_VTK_DEPENDENCY(ApplicationsLib)

if(TARGET Eigen)
    add_dependencies(ApplicationsLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(ApplicationsLib Boost)
endif()
