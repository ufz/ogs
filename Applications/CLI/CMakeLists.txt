
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
