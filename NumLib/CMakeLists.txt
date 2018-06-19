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
include(GenerateExportHeader)
generate_export_header(NumLib)
target_include_directories(NumLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(NumLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(NumLib
    PUBLIC BaseLib GeoLib MathLib MeshLib logog
    PRIVATE MeshGeoToolsLib
)

if(OGS_USE_PCH)
    cotire(NumLib)
endif()
