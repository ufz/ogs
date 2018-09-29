if(Qt5XmlPatterns_FOUND)
    add_executable(TriangulatePolyline TriangulatePolyline.cpp)
    target_link_libraries(TriangulatePolyline GeoLib ApplicationsFileIO)
    set_target_properties(TriangulatePolyline PROPERTIES FOLDER Utilities)
    install(TARGETS
        TriangulatePolyline
        RUNTIME DESTINATION bin
        COMPONENT Utilities
    )

    add_executable(MoveGeometry MoveGeometry.cpp)
    target_link_libraries(MoveGeometry GeoLib)
    set_target_properties(MoveGeometry PROPERTIES FOLDER Utilities)
    install(TARGETS
        MoveGeometry
        RUNTIME DESTINATION bin
        COMPONENT Utilities
    )
endif()
