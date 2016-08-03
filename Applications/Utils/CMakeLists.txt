add_subdirectory(FileConverter)
add_subdirectory(GeoTools)
add_subdirectory(MeshGeoTools)
add_subdirectory(MeshEdit)
add_subdirectory(ModelPreparation)
add_subdirectory(SimpleMeshCreation)

if(OGS_BUILD_GUI)
    add_subdirectory(OGSFileConverter)
else()
    # When OGS_BUILD_GUI=ON, vtk libraries are linked to every object in
    # building, which lets the compilation of mpmetis (in the project
    # of PartitionMesh) fail. Therefore, the following subdirectory is
    # excluded when OGS_BUILD_GUI=ON.
    add_subdirectory(ModelPreparation/PartitionMesh)
endif() # OGS_BUILD_GUI

if(OGS_BUILD_SWMM)
    add_subdirectory(SWMMConverter)
endif()
