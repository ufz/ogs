set(TOOLS MoveGeometry)
foreach(tool ${TOOLS})
    ogs_add_executable(${tool} ${tool}.cpp)
    target_link_libraries(
        ${tool} GeoLib GitInfoLib ApplicationsFileIO tclap
    )
endforeach()
install(TARGETS ${TOOLS} RUNTIME DESTINATION bin)
