add_executable(SWMMInterface SWMMInterface.cpp)

target_link_libraries(SWMMInterface
    MeshLib
    swmm5interface
)

ADD_VTK_DEPENDENCY(SWMMInterface)

set_target_properties(SWMMInterface PROPERTIES FOLDER Utilities)
