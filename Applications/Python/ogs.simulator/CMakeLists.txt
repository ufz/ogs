# Create OpenGeoSys python module
# https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-cmake
pybind11_add_module(
    simulator MODULE ogs_python_module.cpp
    ../../CLI/CommandLineArgumentParser.cpp
)

ogs_remove_lld_linker_option(simulator)

target_link_libraries(
    simulator PRIVATE ApplicationsLib BaseLib CMakeInfoLib GitInfoLib tclap
                      _python_mesh
)
target_include_directories(simulator PRIVATE ../../CLI ${PROJECT_SOURCE_DIR})

install(TARGETS simulator LIBRARY DESTINATION ${_py_install_location})
