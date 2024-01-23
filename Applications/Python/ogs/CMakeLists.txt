if(OGS_BUILD_WHEEL)
    return()
endif()

install(
    DIRECTORY .
    DESTINATION ${_py_install_location}
    PATTERN "__pycache__" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
)
