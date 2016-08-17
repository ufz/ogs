add_executable(ComputeNodeAreasFromSurfaceMesh ComputeNodeAreasFromSurfaceMesh.cpp)
set_target_properties(ComputeNodeAreasFromSurfaceMesh PROPERTIES FOLDER Utilities)
target_link_libraries(ComputeNodeAreasFromSurfaceMesh
    MeshLib
    ${OGS_VTK_REQUIRED_LIBS}
)

add_executable(scaleProperty scaleProperty.cpp)
set_target_properties(scaleProperty PROPERTIES FOLDER Utilities)
target_link_libraries(scaleProperty
    MeshLib
    ${OGS_VTK_REQUIRED_LIBS}
)
