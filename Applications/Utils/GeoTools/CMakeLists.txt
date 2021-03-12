if(NOT OGS_BUILD_GUI)
    return()
endif()

set(TOOLS MoveGeometry TriangulatePolyline)
foreach(tool ${TOOLS})
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(
        ${tool} GeoLib GitInfoLib ApplicationsFileIO tclap Qt5::Core
    )
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
