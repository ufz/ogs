# Source files grouped by a directory
get_source_files(SOURCES)
append_source_files(SOURCES PhreeqcIOData)
append_source_files(SOURCES PhreeqcKernelData)
append_source_files(SOURCES Common)

# Create the library
add_library(ChemistryLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ChemistryLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

include(GenerateExportHeader)
generate_export_header(ChemistryLib)

target_link_libraries(ChemistryLib PUBLIC iphreeqc PRIVATE NumLib)

set_target_properties(ChemistryLib PROPERTIES CXX_STANDARD 14)

if(OGS_USE_PCH)
    cotire(ChemistryLib)
endif()
