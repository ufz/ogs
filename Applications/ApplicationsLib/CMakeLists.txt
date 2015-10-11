# Source files
set(LIB_SOURCES ProjectData.cpp)

# Library
add_library(ApplicationsLib STATIC ${LIB_SOURCES})

target_link_libraries(ApplicationsLib INTERFACE
	GeoLib
	MeshGeoToolsLib
	NumLib
	ProcessLib
	${CATALYST_LIBRARIES}
	logog
)
ADD_CATALYST_DEPENDENCY(ApplicationsLib)

if(TARGET Eigen)
	add_dependencies(ApplicationsLib Eigen)
endif()
