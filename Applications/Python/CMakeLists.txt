# wheel: Install into Python module root dir (enables 'import ogs.simulator')
set(_py_install_location ogs)
set(_py_build_location ogs)
if(NOT OGS_BUILD_WHEEL)
    set(_py_install_location
        "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/ogs"
    )
    set(_py_build_location "${PROJECT_BINARY_DIR}/site-packages/ogs")
endif()

add_subdirectory(ogs)
add_subdirectory(ogs.simulator)
add_subdirectory(ogs.mesh)
add_subdirectory(ogs.callbacks)

if(OGS_USE_PIP)
    set_target_properties(
        simulator mesh callbacks PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                            ${_py_build_location}
    )
    file(
        COPY ogs/.
        DESTINATION ${_py_build_location}
        PATTERN "__pycache__" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
    )
    # Fails with libGitInfoLib.so: undefined symbol: __asan_report_load8
    if(NOT ENABLE_ASAN)
        add_test(
            NAME pytest
            # TODO: currently failing with:
            #
            # Can't add new modules after the interpreter has been initialized
            COMMAND pytest -c ${PROJECT_SOURCE_DIR}/pyproject.toml
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
        set_tests_properties(
            pytest PROPERTIES LABELS "default;python" COST 10
        )

        add_dependencies(
            ctest
            GMSH2OGS
            Layers2Grid
            AddFaultToVoxelGrid
            ExtractBoundary
            vtkdiff
            generateStructuredMesh
        )
    endif()

endif()
