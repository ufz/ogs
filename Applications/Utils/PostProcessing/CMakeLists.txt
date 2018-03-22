
add_executable(postLIE postLIE.cpp)
target_link_libraries(postLIE MeshLib LIE)
set_target_properties(postLIE PROPERTIES FOLDER Utilities)

####################
### Installation ###
####################
install(TARGETS
    postLIE
    RUNTIME DESTINATION bin
    COMPONENT Utilities
)
