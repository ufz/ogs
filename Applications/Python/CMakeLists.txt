# wheel: Install into Python module root dir (enables 'import ogs.simulator')
set(_py_install_location ogs)
if(NOT OGS_BUILD_WHEEL)
    set(_py_install_location
        "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/ogs"
    )
endif()
add_subdirectory(ogs)
add_subdirectory(ogs.simulator)
add_subdirectory(ogs.mesh)
add_subdirectory(ogs.callbacks)
