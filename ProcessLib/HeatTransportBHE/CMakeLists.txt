append_source_files(SOURCES)
append_source_files(SOURCES BHE)
append_source_files(SOURCES BoundaryConditions)
append_source_files(SOURCES LocalAssemblers)

ogs_add_library(HeatTransportBHE ${SOURCES})

target_link_libraries(HeatTransportBHE PUBLIC ProcessLib)
if (OGS_USE_PYTHON)
    target_link_libraries(HeatTransportBHE PRIVATE pybind11::pybind11)
endif()

target_link_libraries(HeatTransportBHE PUBLIC ProcessLib PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
