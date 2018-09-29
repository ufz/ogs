# Create the executable
add_executable(MeshRead
    MeshRead.cpp
    ${SOURCES}
    ${HEADERS}
)

target_link_libraries(MeshRead
    MeshLib
    MathLib
    BaseLib
    GeoLib
    ${BOOST_LIBRARIES}
)

# Create MeshSearchTest executable
add_executable(MeshSearchTest
    MeshSearchTest.cpp
    ${SOURCES}
    ${HEADERS}
)

target_link_libraries(MeshSearchTest
    MeshLib
    MathLib
    BaseLib
    GeoLib
    ${BOOST_LIBRARIES}
)
