# Create OpenGeoSys python module
# https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-cmake
pybind11_add_module(
    simulator MODULE ogs_python_module.cpp
    ../../CLI/CommandLineArgumentParser.cpp
)

# lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
get_target_property(_link_options simulator LINK_OPTIONS)
if(_link_options)
    list(REMOVE_ITEM _link_options -fuse-ld=lld)
    set_target_properties(simulator PROPERTIES LINK_OPTIONS "${_link_options}")
endif()

target_link_libraries(
    simulator PRIVATE ApplicationsLib BaseLib CMakeInfoLib GitInfoLib tclap
                      _python_mesh
)
target_include_directories(simulator PRIVATE ../../CLI)

install(TARGETS simulator LIBRARY DESTINATION ${_py_install_location})
