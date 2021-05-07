ogs_add_executable(SWMMConverter SWMMConverter.cpp)
target_link_libraries(
    SWMMConverter
    ApplicationsFileIO
    GeoLib
    GitInfoLib
    MeshLib
    SwmmInterface
    tclap
)
install(TARGETS SWMMConverter RUNTIME DESTINATION bin)
