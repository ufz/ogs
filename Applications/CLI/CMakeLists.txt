
add_executable(ogs
    ogs.cpp
)

target_link_libraries(ogs
    ApplicationsLib
    BaseLib
    FileIO
    MeshLib
    ${Boost_LIBRARIES}
)
if(NOT OGS_DONT_USE_VTK)
    target_link_libraries(ogs
        InSituLib
        ${CATALYST_LIBRARIES}
    )
endif()

ADD_CATALYST_DEPENDENCY(ogs)

if(OGS_USE_PETSC)
    target_link_libraries(ogs ${PETSC_LIBRARIES})
endif()

if(NOT OGS_USE_MPI)
    include(Tests.cmake)
endif()

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

include(packaging/InstallDependencies)
InstallDependencies(ogs ogs_cli)
