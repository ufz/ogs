pybind11_add_module(mpl MODULE ogs_mpl.cpp)

ogs_remove_lld_linker_option(mpl)

target_link_libraries(
    mpl PRIVATE MaterialLib BaseLib range-v3
)
target_include_directories(mpl PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS mpl LIBRARY DESTINATION ${_py_install_location})
