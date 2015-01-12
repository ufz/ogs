
ADD_EXECUTABLE(ogs
    ogs.cpp
)

TARGET_LINK_LIBRARIES(ogs
    ApplicationsLib
    BaseLib
    FileIO
    MeshLib
    InSituLib
    ${CATALYST_LIBRARIES}
    ${Boost_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(ogs)

IF(OGS_USE_PETSC)
    TARGET_LINK_LIBRARIES(ogs ${PETSC_LIBRARIES})
ENDIF()

INCLUDE(Tests.cmake)

####################
### Installation ###
####################
INSTALL (TARGETS ogs RUNTIME DESTINATION bin COMPONENT ogs_cli)

set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ogs" "OGS Simulator")
cpack_add_component(ogs_cli
	DISPLAY_NAME "OGS THMC Simulator"
	DESCRIPTION "The OGS simulation as a command line executable."
	GROUP Applications
)
