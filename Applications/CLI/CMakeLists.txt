
ADD_EXECUTABLE(ogs
    ogs.cpp
)

TARGET_LINK_LIBRARIES(ogs
    ApplicationsLib
    BaseLib
    FileIO
    MeshLib
    ${Boost_LIBRARIES}
)

IF(OGS_USE_PETSC)
    TARGET_LINK_LIBRARIES(ogs ${PETSC_LIBRARIES})
ENDIF()

INCLUDE(Tests.cmake)
