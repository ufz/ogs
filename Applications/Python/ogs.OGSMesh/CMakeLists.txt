# Helper library which provides a better view on the mesh data. This allows to
# be used by several modules, e.g. this and the simulator module.
add_library(_python_mesh STATIC OGSMesh.cpp)
target_link_libraries(
    _python_mesh PUBLIC MeshLib BaseLib GitInfoLib pybind11::pybind11
)
target_include_directories(_python_mesh PRIVATE ${PROJECT_SOURCE_DIR})

# The actual Python module should only contain the binding code.
pybind11_add_module(OGSMesh MODULE ogs_mesh_python_module.cpp)

ogs_remove_lld_linker_option(OGSMesh)

target_link_libraries(OGSMesh PRIVATE _python_mesh)
target_include_directories(OGSMesh PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS OGSMesh LIBRARY DESTINATION ${_py_install_location})
