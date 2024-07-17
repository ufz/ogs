set(_py_site_packages_dir
    "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages"
)
if(WIN32)
    set(_py_site_packages_dir "Lib/site-packages")
endif()

# wheel: Install into Python module root dir (enables 'import ogs.simulator')
set(_py_install_location ogs)
if(NOT OGS_BUILD_WHEEL)
    set(_py_install_location "{_py_site_packages_dir}/ogs")
endif()

add_subdirectory(ogs)
add_subdirectory(ogs.simulator)
add_subdirectory(ogs.mesh)
add_subdirectory(ogs.callbacks)

if(OGS_USE_PIP)
    set_target_properties(
        simulator mesh callbacks
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                   ${LOCAL_VIRTUALENV_DIR}/${_py_site_packages_dir}/ogs
    )
    file(
        COPY ogs/.
        DESTINATION ${LOCAL_VIRTUALENV_DIR}/${_py_site_packages_dir}/ogs
        PATTERN "__pycache__" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
    )
endif()
