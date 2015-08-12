if(WIN32)
	set(ADDITIONAL_LIBS Winmm.lib)
endif()

# Create the executable
add_executable(MeshRead
	MeshRead.cpp
	${SOURCES}
	${HEADERS}
)

target_link_libraries(MeshRead
	MeshLib
	FileIO
	MathLib
	BaseLib
	GeoLib
	logog
	zlib
	${ADDITIONAL_LIBS}
	${BOOST_LIBRARIES}
)

# Create MeshSearchTest executable
add_executable(MeshSearchTest
	MeshSearchTest.cpp
	${SOURCES}
	${HEADERS}
)

target_link_libraries(MeshSearchTest
	MeshLib
	FileIO
	MathLib
	BaseLib
	GeoLib
	logog
	${ADDITIONAL_LIBS}
	${BOOST_LIBRARIES}
)
