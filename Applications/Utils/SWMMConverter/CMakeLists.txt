add_executable(SWMMConverter SWMMConverter.cpp)
target_link_libraries(SWMMConverter
                      ApplicationsFileIO
                      GeoLib
                      GitInfoLib
                      MeshLib
                      SwmmInterface
                      tclap
)
include_directories(${PROJECT_SOURCE_DIR}/Applications/FileIO)
set_target_properties(SWMMConverter PROPERTIES FOLDER Utilities)
install(TARGETS SWMMConverter RUNTIME DESTINATION bin)
