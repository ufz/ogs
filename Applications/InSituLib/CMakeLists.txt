# Source files
append_source_and_header_files(_sources _headers _private DIRECTORIES .)

# Library
ogs_add_library(InSituLib ${_sources} PUBLIC_HEADERS ${_headers})
if(BUILD_SHARED_LIBS)
    install(TARGETS InSituLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(
    InSituLib
    PUBLIC BaseLib
    INTERFACE VTK::PythonUsed
    PRIVATE MeshLib ParaView::PythonCatalyst VTK::CommonDataModel
)
