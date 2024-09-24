# cmake-lint: disable=E1126
set(_config_content
    "OGS_USE_PETSC = '@OGS_USE_PETSC@'" "OGS_VERSION = '@OGS_VERSION@'"
    "OGS_USE_MKL = '@OGS_USE_MKL@'"
)
string(REPLACE ";" "\n" _config_content "${_config_content}")
set(_config_file ${PROJECT_BINARY_DIR}/site-packages/ogs/config.py)
file(CONFIGURE OUTPUT ${_config_file} CONTENT "${_config_content}")

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
