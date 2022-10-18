if(OGS_BUILD_PYTHON_MODULE)
    add_subdirectory(ogs.simulator)
endif()

if(OGS_USE_PYTHON)
    add_subdirectory(ogs.callbacks)
endif()
