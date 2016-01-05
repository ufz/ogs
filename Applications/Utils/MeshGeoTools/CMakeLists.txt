ADD_EXECUTABLE( ComputeSurfaceNodeIDsInPolygonalRegion ComputeSurfaceNodeIDsInPolygonalRegion.cpp )
SET_TARGET_PROPERTIES( ComputeSurfaceNodeIDsInPolygonalRegion PROPERTIES FOLDER Utils)
TARGET_LINK_LIBRARIES( ComputeSurfaceNodeIDsInPolygonalRegion
	MeshLib
	FileIO
	GeoLib
	MathLib
	${OGS_VTK_REQUIRED_LIBS}
)

