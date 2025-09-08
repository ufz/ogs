# Helper library which provides a better view on the mesh data. This allows to
# be used by several modules, e.g. this and the simulator module.
add_library(_python_mesh STATIC OGSMesh.cpp)
target_link_libraries(
    _python_mesh PUBLIC MeshLib BaseLib GitInfoLib pybind11::pybind11
)
target_include_directories(_python_mesh PRIVATE ${PROJECT_SOURCE_DIR})

# The actual Python module should only contain the binding code.
pybind11_add_module(mesh MODULE ogs_mesh_python_module.cpp)

ogs_remove_lld_linker_option(mesh)

target_link_libraries(mesh PRIVATE _python_mesh)
target_include_directories(mesh PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS mesh LIBRARY DESTINATION ${_py_install_location})
