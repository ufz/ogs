add_executable(SWMMConverter
    SWMMConverter.cpp
    SwmmInterface.h
    SwmmInterface.cpp
)

target_link_libraries(SWMMConverter
    ApplicationsFileIO
    GeoLib
    MeshLib
    swmm5interface
)

ADD_VTK_DEPENDENCY(SWMMConverter)

set_target_properties(SWMMConverter PROPERTIES FOLDER Utilities)
