APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BHE)
APPEND_SOURCE_FILES(SOURCES LocalAssemblers)

add_library(HeatTransportBHE ${SOURCES})
target_link_libraries(HeatTransportBHE PUBLIC ProcessLib)

# include(Tests.cmake)
