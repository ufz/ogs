
add_executable(postLIE postLIE.cpp)
target_link_libraries(postLIE MeshLib ProcessLib)
ADD_VTK_DEPENDENCY(postLIE)
set_target_properties(postLIE PROPERTIES FOLDER Utilities)

