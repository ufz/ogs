ogs_add_executable(
    mergeMeshToBulkMesh
    MergeMeshToBulkMesh.cpp
)
target_link_libraries(mergeMeshToBulkMesh GitInfoLib tclap MeshToolsLib NumLib)
install(TARGETS mergeMeshToBulkMesh RUNTIME DESTINATION bin)
