pybind11_add_module(callbacks MODULE ogs_callbacks_module.cpp)

# lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
get_target_property(_link_options callbacks LINK_OPTIONS)
if(_link_options)
    list(REMOVE_ITEM _link_options -fuse-ld=lld)
    set_target_properties(callbacks PROPERTIES LINK_OPTIONS "${_link_options}")
endif()

target_link_libraries(
    callbacks PRIVATE ProcessLibBoundaryConditionAndSourceTermPythonModule
)

install(TARGETS callbacks LIBRARY DESTINATION ${_py_install_location})
