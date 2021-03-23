add_executable(SWMMConverter SWMMConverter.cpp)
target_link_libraries(
    SWMMConverter
    ApplicationsFileIO
    GeoLib
    GitInfoLib
    MeshLib
    SwmmInterface
    tclap
)
include_directories(${PROJECT_SOURCE_DIR}/Applications/FileIO)
install(TARGETS SWMMConverter RUNTIME DESTINATION bin)
