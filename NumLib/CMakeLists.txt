GET_SOURCE_FILES(SOURCES)

APPEND_SOURCE_FILES(SOURCES Assembler)
APPEND_SOURCE_FILES(SOURCES DOF)
APPEND_SOURCE_FILES(SOURCES Fem)
APPEND_SOURCE_FILES(SOURCES Fem/CoordinatesMapping)
APPEND_SOURCE_FILES(SOURCES Fem/FiniteElement)
APPEND_SOURCE_FILES(SOURCES Fem/Integration)
APPEND_SOURCE_FILES(SOURCES Fem/ShapeFunction)
APPEND_SOURCE_FILES(SOURCES TimeStepping)
APPEND_SOURCE_FILES(SOURCES TimeStepping/Algorithms)
APPEND_SOURCE_FILES(SOURCES Function)
APPEND_SOURCE_FILES(SOURCES Distribution)
APPEND_SOURCE_FILES(SOURCES ODESolver)
APPEND_SOURCE_FILES(SOURCES Extrapolation)

# Create the library
add_library(NumLib ${SOURCES})

set_target_properties(NumLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(NumLib
    PUBLIC BaseLib GeoLib MathLib MeshLib logog
    PRIVATE MeshGeoToolsLib
)

if(TARGET Eigen)
    add_dependencies(NumLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(NumLib Boost)
endif()

if(OGS_USE_PCH)
    cotire(NumLib)
endif()
