# Source files
GET_SOURCE_FILES(SOURCES_APPLICATIONSLIB)
set(LIB_SOURCES ${SOURCES_APPLICATIONSLIB})

# Library
add_library(ApplicationsLib ${LIB_SOURCES})

target_link_libraries(ApplicationsLib
    PUBLIC BaseLib GeoLib Processes logog
    PRIVATE MathLib MeshLib
)

if(TARGET Eigen)
    add_dependencies(ApplicationsLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(ApplicationsLib Boost)
endif()

if(OGS_USE_PCH)
    cotire(ApplicationsLib)
endif()
