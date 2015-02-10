include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
)

if(QT4_FOUND)
	add_executable( TriangulatePolyline TriangulatePolyline.cpp )
	target_link_libraries( TriangulatePolyline
		BaseLib
		FileIO
		MathLib
		${QT_LIBRARIES}
	)
	set_target_properties(TriangulatePolyline PROPERTIES FOLDER Utilities)
endif() # QT4_FOUND
