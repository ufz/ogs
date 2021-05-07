get_source_files(SOURCES)

append_source_files(SOURCES Assembler)
append_source_files(SOURCES DOF)
append_source_files(SOURCES Fem)
append_source_files(SOURCES Fem/CoordinatesMapping)
append_source_files(SOURCES Fem/FiniteElement)
append_source_files(SOURCES Fem/Integration)
append_source_files(SOURCES Fem/ShapeFunction)
append_source_files(SOURCES TimeStepping)
append_source_files(SOURCES TimeStepping/Algorithms)
append_source_files(SOURCES Function)
append_source_files(SOURCES ODESolver)
append_source_files(SOURCES Extrapolation)

# Create the library
ogs_add_library(NumLib ${SOURCES})
if(OGS_USE_PETSC)
    target_compile_definitions(NumLib PUBLIC USE_PETSC)
endif()
set_target_properties(NumLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(
    NumLib PUBLIC BaseLib GeoLib MathLib MeshLib PRIVATE MeshGeoToolsLib
                                                         spdlog::spdlog
)
