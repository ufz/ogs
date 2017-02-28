
add_executable(ogs
    ogs.cpp
)

target_link_libraries(ogs
    ApplicationsLib
    MeshGeoToolsLib
    ProcessLib
    NumLib
)

ADD_VTK_DEPENDENCY(ogs)

if(OGS_USE_PETSC)
    target_link_libraries(ogs ${PETSC_LIBRARIES})
endif()

if(OGS_INSITU)
    target_link_libraries(ogs InSituLib)
endif()

####################
### Tests ##########
####################
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

####################
### Installation ###
####################
install(TARGETS ogs RUNTIME DESTINATION bin COMPONENT ogs_cli)

set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ogs" "OGS Simulator")
cpack_add_component(ogs_cli
    DISPLAY_NAME "OGS THMC Simulator"
    DESCRIPTION "The command line interface for OpenGeoSys."
    GROUP Applications
)
