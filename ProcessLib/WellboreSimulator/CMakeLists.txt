append_source_and_header_files(_sources _headers _private DIRECTORIES .)

ogs_add_library(WellboreSimulator ${_sources} PUBLIC_HEADERS ${_headers})
target_link_libraries(WellboreSimulator PUBLIC ProcessLib PRIVATE ParameterLib)

if(OGS_BUILD_TESTING)
    include(Tests.cmake)
endif()
