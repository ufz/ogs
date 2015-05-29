# Source files
set(LIB_SOURCES ProjectData.cpp)

# Library
add_library(ApplicationsLib STATIC ${LIB_SOURCES})

target_link_libraries(ApplicationsLib
	GeoLib
	MeshGeoToolsLib
	ProcessLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(ApplicationsLib)
