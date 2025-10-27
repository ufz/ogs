set(_config_file ${PROJECT_BINARY_DIR}/site-packages/ogs/config.py)
file(CONFIGURE OUTPUT "${_config_file}" CONTENT "
OGS_USE_PETSC = '@OGS_USE_PETSC@'
OGS_VERSION = '@OGS_VERSION@'
OGS_USE_MKL = '@OGS_USE_MKL@'
OGS_COVERAGE = '@OGS_COVERAGE@'
OGS_GUIX_BUILD = '@GUIX_BUILD@'
")

install(FILES ${_config_file} DESTINATION ${_py_install_location})

if(OGS_BUILD_WHEEL)
    return()
endif()

install(
    DIRECTORY .
    DESTINATION ${_py_install_location}
    PATTERN "__pycache__" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
)
