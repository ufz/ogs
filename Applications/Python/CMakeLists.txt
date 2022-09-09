# Create OpenGeoSys python module
# https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-cmake
pybind11_add_module(
    _cli MODULE ogs_python_module.cpp ../CLI/CommandLineArgumentParser.cpp
)

# lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
get_target_property(_link_options _cli LINK_OPTIONS)
if(_link_options)
    list(REMOVE_ITEM _link_options -fuse-ld=lld)
    set_target_properties(_cli PROPERTIES LINK_OPTIONS "${_link_options}")
endif()

target_link_libraries(_cli PRIVATE ApplicationsLib BaseLib CMakeInfoLib GitInfoLib tclap)
target_include_directories(_cli PRIVATE ../CLI)

install(TARGETS _cli LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
