add_executable( ComputeSurfaceNodeIDsInPolygonalRegion ComputeSurfaceNodeIDsInPolygonalRegion.cpp )
set_target_properties( ComputeSurfaceNodeIDsInPolygonalRegion PROPERTIES FOLDER Utilities)
target_link_libraries( ComputeSurfaceNodeIDsInPolygonalRegion
    MeshLib
    GeoLib
    MathLib
    ApplicationsFileIO
    ${OGS_VTK_REQUIRED_LIBS}
)

