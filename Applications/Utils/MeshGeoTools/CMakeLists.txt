add_executable( ComputeSurfaceNodeIDsInPolygonalRegion ComputeSurfaceNodeIDsInPolygonalRegion.cpp )
set_target_properties( ComputeSurfaceNodeIDsInPolygonalRegion PROPERTIES FOLDER Utilities)
target_link_libraries( ComputeSurfaceNodeIDsInPolygonalRegion
    MeshLib
    GeoLib
    MathLib
    ApplicationsFileIO
    ${OGS_VTK_REQUIRED_LIBS}
)

add_executable(constructMeshesFromGeometry ConstructMeshesFromGeometry.cpp )
set_target_properties(constructMeshesFromGeometry PROPERTIES FOLDER Utilities)
target_link_libraries(constructMeshesFromGeometry
    MeshGeoToolsLib
    ApplicationsFileIO
    ${OGS_VTK_REQUIRED_LIBS}
)
