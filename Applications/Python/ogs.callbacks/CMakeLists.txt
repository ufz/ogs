pybind11_add_module(callbacks MODULE ogs_callbacks_module.cpp)

ogs_remove_lld_linker_option(callbacks)

target_link_libraries(
    callbacks PRIVATE ProcessLibBoundaryConditionAndSourceTermPythonModule
)
target_include_directories(callbacks PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS callbacks LIBRARY DESTINATION ${_py_install_location})
