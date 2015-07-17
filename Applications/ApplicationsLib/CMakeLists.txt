# Source files
set(LIB_SOURCES ProjectData.cpp)

# Library
add_library(ApplicationsLib STATIC ${LIB_SOURCES})

target_link_libraries(ApplicationsLib INTERFACE
	GeoLib
	MeshGeoToolsLib
	ProcessLib
	${CATALYST_LIBRARIES}
	logog
)
ADD_CATALYST_DEPENDENCY(ApplicationsLib)

if(OGS_BUILD_GUI)
	target_link_libraries(ApplicationsLib PUBLIC Qt4::QtCore)
endif()
