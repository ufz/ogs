# Helper library which provides a better view on the mesh data. This allows to
# be used by several modules, e.g. this and the simulator module.
add_library(_python_mesh STATIC OGSMesh.cpp)
target_link_libraries(
    _python_mesh PUBLIC MeshLib BaseLib GitInfoLib pybind11::pybind11
)

# The actual Python module should only contain the binding code.
pybind11_add_module(mesh MODULE ogs_mesh_python_module.cpp)

# lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
get_target_property(_link_options mesh LINK_OPTIONS)
if(_link_options)
    list(REMOVE_ITEM _link_options -fuse-ld=lld)
    set_target_properties(mesh PROPERTIES LINK_OPTIONS "${_link_options}")
endif()

target_link_libraries(mesh PRIVATE _python_mesh)

install(TARGETS mesh LIBRARY DESTINATION ${_py_install_location})
