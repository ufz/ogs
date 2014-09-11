
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

IF(OGS_BUILD_GUI)
    TARGET_LINK_LIBRARIES(ogs QtDataView ${VTK_LIBRARIES} VtkVis)
ENDIF()
