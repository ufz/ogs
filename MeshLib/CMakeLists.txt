get_source_files(SOURCES)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
    list(REMOVE_ITEM SOURCES NodePartitionedMesh.h)
endif()

append_source_files(SOURCES MeshEditing)
append_source_files(SOURCES MeshGenerators)
append_source_files(SOURCES MeshSearch)
append_source_files(SOURCES Elements)
append_source_files(SOURCES IO)
append_source_files(SOURCES IO/Legacy)
append_source_files(SOURCES IO/VtkIO)
append_source_files(SOURCES MeshQuality)
append_source_files(SOURCES Vtk)

# It could be used for other MPI based DDC approach in future.
if(OGS_USE_PETSC)
    append_source_files(SOURCES IO/MPI_IO)
endif()

# Create the library
add_library(MeshLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS MeshLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
if(OGS_USE_PCH)
    cotire(MeshLib)
endif()

target_link_libraries(MeshLib
                      PUBLIC BaseLib
                             GeoLib
                             MathLib
                             logog
                             ${VTK_LIBRARIES})
