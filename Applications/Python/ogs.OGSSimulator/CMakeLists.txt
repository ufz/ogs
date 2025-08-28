# Create OpenGeoSys python module
# https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-cmake
pybind11_add_module(
    OGSSimulator MODULE ogs_python_module.cpp
    ../../CLI/CommandLineArgumentParser.cpp
)

ogs_remove_lld_linker_option(OGSSimulator)

target_link_libraries(
    OGSSimulator PRIVATE ApplicationsLib BaseLib CMakeInfoLib GitInfoLib tclap
                      _python_mesh
)
target_include_directories(OGSSimulator PRIVATE ../../CLI ${PROJECT_SOURCE_DIR})

install(TARGETS OGSSimulator LIBRARY DESTINATION ${_py_install_location})
