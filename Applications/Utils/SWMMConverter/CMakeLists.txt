add_executable(SWMMConverter
    SWMMConverter.cpp
)

target_link_libraries(SWMMConverter
    ApplicationsFileIO
    GeoLib
    MeshLib
    SwmmInterface
)

include_directories( ${CMAKE_SOURCE_DIR}/Applications/FileIO )

ADD_VTK_DEPENDENCY(SWMMConverter)

set_target_properties(SWMMConverter PROPERTIES FOLDER Utilities)
set_target_properties(SWMMConverter PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcmt")
