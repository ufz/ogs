# Create OpenGeoSys python module
# https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-cmake
pybind11_add_module(
    mesh MODULE ogs_mesh_python_module.cpp OGSMesh.cpp
)

# lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
get_target_property(_link_options mesh LINK_OPTIONS)
if(_link_options)
    list(REMOVE_ITEM _link_options -fuse-ld=lld)
    set_target_properties(mesh PROPERTIES LINK_OPTIONS "${_link_options}")
endif()

target_link_libraries(
    mesh PRIVATE MeshLib BaseLib CMakeInfoLib GitInfoLib tclap
)
target_include_directories(mesh  PRIVATE ../../../MeshLib)

# Install into Python module root dir (enables 'import ogs.simulator')
install(TARGETS mesh LIBRARY DESTINATION ogs)
