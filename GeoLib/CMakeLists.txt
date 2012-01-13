# Source files
SET ( GeoLib_Files
	AxisAlignedBoundingBox.cpp
	BruteForceClosestPair.cpp
	Color.cpp
	GEOObjects.cpp
	GeoType.cpp
	Point.cpp
	PointVec.cpp
	Polygon.cpp
	Polyline.cpp
#	ProjectData.cpp # uses CFEMesh class
#	SimplePolygonTree.cpp # uses FileIO/MeshIO
	Station.cpp
	Surface.cpp
	Triangle.cpp
	Raster.cpp
	AxisAlignedBoundingBox.h
	BruteForceClosestPair.h
	ClosestPair.h
	Color.h
	GeoObject.h
	GEOObjects.h
	GeoType.h
	Point.h
	PointVec.h
	PointWithID.h
	Polygon.h
	Polyline.h
	PolylineVec.h
#	ProjectData.h # uses CFEMesh class
	PropertyBounds.h
	QuadTree.h
#	SimplePolygonTree.h # uses FileIO/MeshIO
	Station.h
	Surface.h
	SurfaceVec.h
	TemplatePoint.h
	TemplateVec.h
	Triangle.h
	Raster.h
)
SOURCE_GROUP(GeoLib FILES ${GeoLib_Files})

# Create the library
ADD_LIBRARY(GeoLib STATIC ${GeoLib_Files})


include_directories(
	.
	../Base
	../MathLib
)


target_link_libraries (
	GeoLib
	Base
	MathLib
)

