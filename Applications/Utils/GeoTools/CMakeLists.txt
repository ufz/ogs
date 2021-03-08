if(NOT OGS_BUILD_GUI)
    return()
endif()

set(TOOLS MoveGeometry TriangulatePolyline)
foreach(TOOL ${TOOLS})
    add_executable(${TOOL} ${TOOL}.cpp)
    target_link_libraries(${TOOL}
                          GeoLib
                          GitInfoLib
                          ApplicationsFileIO
                          tclap
                          Qt5::Core)
endforeach()
set_target_properties(${TOOLS} PROPERTIES FOLDER Utilities)
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
