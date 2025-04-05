pybind11_add_module(mpl MODULE ogs_mpl.cpp)

# lld linker strips out PyInit_OpenGeoSys symbol. Use standard linker.
get_target_property(_link_options mpl LINK_OPTIONS)
if(_link_options)
    list(REMOVE_ITEM _link_options -fuse-ld=lld)
    set_target_properties(mpl PROPERTIES LINK_OPTIONS "${_link_options}")
endif()

target_link_libraries(
    mpl PRIVATE MaterialLib BaseLib range-v3
)

install(TARGETS mpl LIBRARY DESTINATION ${_py_install_location})
