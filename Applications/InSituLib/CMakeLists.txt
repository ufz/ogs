# Source files
get_source_files(SOURCES)

# Library
ogs_add_library(InSituLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS InSituLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(InSituLib
    INTERFACE
      VTK::PythonUsed
    PRIVATE
      ParaView::PythonCatalyst
      VTK::CommonDataModel)

target_link_libraries(InSituLib PUBLIC BaseLib PRIVATE MeshLib)
