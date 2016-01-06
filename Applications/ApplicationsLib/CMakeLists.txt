# Source files
set(LIB_SOURCES ProjectData.cpp)

# Library
add_library(ApplicationsLib STATIC ${LIB_SOURCES})

target_link_libraries(ApplicationsLib INTERFACE
	GeoLib
	MeshGeoToolsLib
	NumLib
	ProcessLib
	${VTK_LIBRARIES}
	logog
)
ADD_VTK_DEPENDENCY(ApplicationsLib)

if(TARGET Eigen)
	add_dependencies(ApplicationsLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(ApplicationsLib Boost)
endif()
