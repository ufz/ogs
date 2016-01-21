include_directories(
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
)

if(Qt5XmlPatterns_FOUND)
	add_executable(TriangulatePolyline TriangulatePolyline.cpp)
	target_link_libraries(TriangulatePolyline FileIO)
	set_target_properties(TriangulatePolyline PROPERTIES FOLDER Utilities)
	install(TARGETS
		TriangulatePolyline
		RUNTIME DESTINATION bin
		COMPONENT Utilities
	)

	add_executable(MoveGeometry MoveGeometry.cpp)
	target_link_libraries(MoveGeometry FileIO)
	set_target_properties(MoveGeometry PROPERTIES FOLDER Utilities)
	install(TARGETS
		MoveGeometry
		RUNTIME DESTINATION bin
		COMPONENT Utilities
	)
endif()
