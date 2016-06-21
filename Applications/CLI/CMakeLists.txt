
include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}/../../BaseLib
    ${CMAKE_CURRENT_BINARY_DIR}/../../MathLib
    ${CMAKE_CURRENT_BINARY_DIR}/../../NumLib
    ${CMAKE_CURRENT_BINARY_DIR}/../../MaterialsLib
)

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

include(Tests.cmake)

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
