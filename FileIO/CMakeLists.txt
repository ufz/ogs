GET_SOURCE_FILES(SOURCES_FileIO)
set(SOURCES ${SOURCES_FileIO})

if(MSVC)
    include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}/../MathLib
        ${CMAKE_CURRENT_BINARY_DIR}/../MeshLib
        )
endif()

# Create the library
add_library(FileIO ${SOURCES})
target_link_libraries(FileIO GeoLib MeshLib)

ADD_VTK_DEPENDENCY(FileIO)

if(TARGET Eigen)
    add_dependencies(FileIO Eigen)
endif()

if(TARGET Boost)
    add_dependencies(FileIO Boost)
endif()

