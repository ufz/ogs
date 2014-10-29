
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
