ADD_EXECUTABLE( ComputeSurfaceNodeIDsInPolygonalRegion ComputeSurfaceNodeIDsInPolygonalRegion.cpp )
SET_TARGET_PROPERTIES( ComputeSurfaceNodeIDsInPolygonalRegion PROPERTIES FOLDER Utilities)
TARGET_LINK_LIBRARIES( ComputeSurfaceNodeIDsInPolygonalRegion
    MeshLib
    GeoLib
    MathLib
    ApplicationsFileIO
    ${OGS_VTK_REQUIRED_LIBS}
)

