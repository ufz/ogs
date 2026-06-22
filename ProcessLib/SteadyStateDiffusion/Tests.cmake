# CUBE 1x1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3)
    if (${mesh_size} STREQUAL "1e2")
        AddTest(
            NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}.xml --write-prj
            TESTER vtkdiff
            REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu ${mesh_size}_cube_1x1x1_hex_${mesh_size}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_0_t_0.000000.vtu ${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_0_t_0.000000.vtu pressure pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_0_t_0.000000.vtu ${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_0_t_0.000000.vtu v v 5e-13 3e-13
            meshes/${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_1_t_1.000000.vtu pressure pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_hex_${mesh_size}_0_1_ts_1_t_1.000000.vtu v v 5e-13 3e-13
            meshes/${mesh_size}_cube_1x1x1_geometry_back_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_back_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_back_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_back_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_back_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_back_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_front_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_front_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_front_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_front_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_back_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_back_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_back_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_back_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_front_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_front_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_front_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_front_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_back_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_back_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_back_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_back_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
        )
    else (${mesh_size} STREQUAL "1e2")
        AddTest(
            NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}.xml --write-prj
            TESTER vtkdiff
            REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu ${mesh_size}_cube_1x1x1_hex_${mesh_size}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_back_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_back_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_back_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_front_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_left_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_horizontal_polyline_right_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_back_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_back_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_back_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_back_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_front_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_front_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_front_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_front_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_left_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_left_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_back_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_back_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_back_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_back_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_front_bottom_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_front_bottom_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_front_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_front_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_right_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_right_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_top_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_top_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_back_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_back_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_left_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_back_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_back_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
            meshes/${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_front_ts_1_t_1.000000.vtu ${mesh_size}_cube_1x1x1_geometry_vertical_polyline_right_front_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
        )
    endif()

    if(TEST ogs-SteadyStateDiffusion_cube_1x1x1_${mesh_size} AND DIFF_TOOL_PATH)
        set(_processed_path
            Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}_processed.prj
        )
        add_test(
            NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}_prj_diff
            COMMAND ${DIFF_TOOL_PATH} ${Data_SOURCE_DIR}/${_processed_path}
                    ${Data_BINARY_DIR}/${_processed_path}
        )
        set_tests_properties(
            SteadyStateDiffusion_cube_1x1x1_${mesh_size}_prj_diff
            PROPERTIES LABELS "default" DEPENDS
                       ogs-SteadyStateDiffusion_cube_1x1x1_${mesh_size}
        )
    endif()

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}_Newton
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        # `-m .` just for testing input mesh dir parameter:
        EXECUTABLE_ARGS -m . cube_${mesh_size}_newton.prj
        TESTER vtkdiff
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_newton_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure 1e-1 1e-1
    )
endforeach()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
        set(RUNTIME 10)
        if("${mesh_size}" STREQUAL "1e6")
            set(RUNTIME 75)
        endif()
        OgsTest(
            PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}.prj
            RUNTIME ${RUNTIME}
            NAME_SUFFIX ${mesh_size}
        )

        OgsTest(
            PROJECTFILE
                Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}_neumann.prj
            RUNTIME ${RUNTIME}
            NAME_SUFFIX ${mesh_size}_Neumann
        )
    endforeach()
endif()

# CG TriangularMatrix test
foreach(matrix lower upper lowerupper)
    set(RUNTIME 3)
    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_1e0_${matrix}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        RUNTIME ${RUNTIME}
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e0_${matrix}.xml --write-prj
        TESTER vtkdiff
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
        DIFF_DATA
        cube_1x1x1_hex_1e0.vtu 1e0_cube_1x1x1_hex_1e0_${matrix}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
    )
endforeach()



# Test FixedTimeStepping and fixed output times
if (NOT OGS_USE_LIS)
    AddTest(
        NAME SteadyStateDiffusion_square_1x1_quad_1e1_FixedTimeStepping_FixedOutputTimes
        PATH Elliptic/square_1x1_SteadyStateDiffusion/FixedTimeSteppingFixedOutputTimes
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e1-fixed_timestepping-fixed_output_times.prj -m ../
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        GLOB square_1e1_ts_*.vtu pressure pressure 1e-15 1e-15
    )
endif()

# Quadratic hex element.
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e0_quadratic_hex.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    # SQUARE 1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}.xml
            NAME_SUFFIX ${mesh_size}
        )

        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}_neumann.xml
            NAME_SUFFIX ${mesh_size}_Neumann
        )
    endforeach()

    foreach(mesh_size 1e5 1e6)
        set(RUNTIME 30)
        if("${mesh_size}" STREQUAL "1e6")
            set(RUNTIME 90)
        endif()
        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}_neumann.xml
            RUNTIME ${RUNTIME}
            NAME_SUFFIX ${mesh_size}_Neumann
        )
    endforeach()

    OgsTest(
        PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
        PATCH_FILES square_1e5.xml
        RUNTIME 2
        NAME_SUFFIX 1e5
    )

    # The largest test is less accurate
    OgsTest(
        PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
        PATCH_FILES square_1e6.xml
        RUNTIME 15
        NAME_SUFFIX 1e6
    )

    # LINE 1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e1)
        OgsTest(
            PROJECTFILE Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}.prj
            NAME_SUFFIX ${mesh_size}
        )

        OgsTest(
            PROJECTFILE Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}_neumann.prj
            NAME_SUFFIX ${mesh_size}_Neumann
        )
        OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/drainage_excavation.prj)
        OgsTest(
            PROJECTFILE
                Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}_robin_right_picard.prj
            NAME_SUFFIX ${mesh_size}_Robin_Right_Picard
        )

        OgsTest(
            PROJECTFILE
                Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}_robin_left_picard.prj
            NAME_SUFFIX ${mesh_size}_Robin_Left_Picard
        )

        OgsTest(
            PROJECTFILE
                Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}_time_dep_dirichlet.prj
            NAME_SUFFIX ${mesh_size}_Time_Dep_Dirichlet
        )

        OgsTest(
            PROJECTFILE
                Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}_time_dep_neumann.prj
            NAME_SUFFIX ${mesh_size}_Time_Dep_Neumann
        )
    endforeach()
endif()

# Some Neumann BC tests
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e3_top_neumann.prj)
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e3_bottom_neumann.prj)
    # Some Neumann BC tests -- Newton
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e3_top_neumann_newton.prj)
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e3_bottom_neumann_newton.prj)
    # test SurfaceFlux
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e3_calculatesurfaceflux.prj)
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e3_neumann_calculatesurfaceflux.prj)
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_2e3_prism_surfaceflux_left_right.prj)
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_2e3_prism_surfaceflux_front_back.prj)
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_2e3_prism_surfaceflux_top_bottom.prj)
    OgsTest(PROJECTFILE Elliptic/wedge_1x1x1_SteadyStateDiffusion/wedge_1e3_prism_surfaceflux_diagonal.prj)
endif()

# SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC, 1e2
# The results were compared to an analytical solution (method of manufactured
# solutions). The vtkdiff comparison is against the numerical solution.
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e2_axi.prj)
endif()
# WEDGE 1x1 GROUNDWATER FLOW TEST -- same setup as above test but in cartesian coordinates
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/wedge_1e2_axi_ang_0.02.prj)
    # SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC, 1e4
    # The results were compared to an analytical solution (method of manufactured
    # solutions). The vtkdiff comparison is against the numerical solution.
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e4_axi.prj)
    # WEDGE 1x1 GROUNDWATER FLOW TEST -- same setup as above test but in cartesian coordinates
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/wedge_1e4_axi_ang_0.02.prj)
endif()

# Serial XDMF output
if(TARGET xdmfdiff)
    OgsTest(
        PROJECTFILE
            Elliptic/cube_1x1x1_SteadyStateDiffusion/xdmf/cube_1e4_anisotropic.prj
    )
endif()

# MPI groundwater flow tests
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE EllipticPETSc/quad_20x10_GroundWaterFlow.prj
            WRAPPER mpirun -np 3
    )
endif()
if(TEST ogs-EllipticPETSc/quad_20x10_GroundWaterFlow-mpi
   AND XMLSTARLET_TOOL_PATH AND BASH_TOOL_PATH
)
    get_test_property(
        ogs-EllipticPETSc/quad_20x10_GroundWaterFlow-mpi
        WORKING_DIRECTORY
        _groundwater_flow_2d_output_dir
    )
    set(_groundwater_flow_2d_pvtu
        "${_groundwater_flow_2d_output_dir}/quad_20x10_GroundWaterFlow_result_ts_0_t_0_000000.pvtu"
    )

    # Just checks if there are 3 <Piece>-elements in the pvtu
    add_test(
        NAME ParallelFEM_GroundWaterFlow2D_pvtu
        COMMAND
            ${BASH_TOOL_PATH} -c
            "if [[ $(${XMLSTARLET_TOOL_PATH} sel -t -v 'count(/VTKFile/PUnstructuredGrid/Piece)' \"${_groundwater_flow_2d_pvtu}\") == '3' ]] ; then exit 0; else cat \"${_groundwater_flow_2d_pvtu}\"; exit 1; fi"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    set_tests_properties(
        ParallelFEM_GroundWaterFlow2D_pvtu
        PROPERTIES LABELS "default" DEPENDS
                   ogs-EllipticPETSc/quad_20x10_GroundWaterFlow-mpi
    )
endif()
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE EllipticPETSc/cube_1e3.prj WRAPPER mpirun -np 3)
    OgsTest(PROJECTFILE EllipticPETSc/cube_1e3_neumann.prj WRAPPER mpirun -np 3)
endif()

if(OGS_USE_MPI AND TARGET xdmfdiff)
    OgsTest(
        PROJECTFILE EllipticPETSc/cube_1e3_XDMF_np3.prj
        WRAPPER mpirun -np 3
    )
    OgsTest(
        PROJECTFILE EllipticPETSc/cube_1e3_XDMF_np3_2files.prj
        WRAPPER mpirun -np 3
    )
    OgsTest(
        PROJECTFILE EllipticPETSc/cube_1e3_XDMF_np3_3files.prj
        WRAPPER mpirun -np 3
    )
    OgsTest(
        PROJECTFILE EllipticPETSc/cube_1e3_XDMF_np2.prj
        WRAPPER mpirun -np 2
    )
endif()

if(OGS_USE_MPI)
    OgsTest(PROJECTFILE EllipticPETSc/square_1e1_neumann.prj WRAPPER mpirun -np 2)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE "Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e4_anisotropic.prj")
endif()

if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE
            EllipticPETSc/cube_1x1x1_SteadyStateDiffusion/2/cube_1e4_anisotropic_mpi.xml
        WRAPPER mpirun -np 2
    )
    # Single core
    # CUBE 1x1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3)
        if("${mesh_size}" STREQUAL "1e0")
            OgsTest(
                PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e0.prj
                WRAPPER mpirun -np 1
            )
        else()
            OgsTest(
                PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}_processed.prj
                NAME_SUFFIX ${mesh_size}
                WRAPPER mpirun -np 1
            )
        endif()
        OgsTest(
            PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}_neumann.prj
            WRAPPER mpirun -np 1
        )
    endforeach()

    # TODO: Parallel LARGE tests not tested!
    foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
        set(RUNTIME 10)
        if("${mesh_size}" STREQUAL "1e5")
            set(RUNTIME 55)
        endif()
        if("${mesh_size}" STREQUAL "1e6")
            set(RUNTIME 430)
        endif()
        OgsTest(
            PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}.prj
            RUNTIME ${RUNTIME}
            WRAPPER mpirun -np 1
        )
        OgsTest(
            PROJECTFILE
                Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}_neumann.prj
            RUNTIME ${RUNTIME}
            WRAPPER mpirun -np 1
        )
    endforeach()

    # SQUARE 1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}.xml
            WRAPPER mpirun -np 1
            NAME_SUFFIX ${mesh_size}
        )

        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}_neumann.xml
            WRAPPER mpirun -np 1
            NAME_SUFFIX ${mesh_size}_Neumann
        )
    endforeach()

    foreach(mesh_size 1e5 1e6)
        set(RUNTIME 65)
        if("${mesh_size}" STREQUAL "1e6")
            set(RUNTIME 450)
        endif()
        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}.xml
            RUNTIME ${RUNTIME}
            WRAPPER mpirun -np 1
            NAME_SUFFIX ${mesh_size}
        )

        OgsTest(
            PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
            PATCH_FILES square_${mesh_size}_neumann.xml
            RUNTIME ${RUNTIME}
            WRAPPER mpirun -np 1
            NAME_SUFFIX ${mesh_size}_Neumann
        )
    endforeach()

    # LINE 1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e1)
        OgsTest(
            PROJECTFILE Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}.prj
            WRAPPER mpirun -np 1
            NAME_SUFFIX ${mesh_size}
        )

        OgsTest(
            PROJECTFILE Elliptic/line_1_SteadyStateDiffusion/line_${mesh_size}_neumann.prj
            WRAPPER mpirun -np 1
            NAME_SUFFIX ${mesh_size}_Neumann
        )
    endforeach()
endif()
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Elliptic/nonuniform_bc_SteadyStateDiffusion/inhomogeneous_permeability.prj)
    OgsTest(PROJECTFILE Elliptic/nonuniform_bc_SteadyStateDiffusion/neumann_nonuniform.prj)
    OgsTest(PROJECTFILE Elliptic/nonuniform_bc_SteadyStateDiffusion/dirichlet_nonuniform.prj)
    # tests for nodal source term implementation
    # For the special setup with a dirac source term at position (xi, eta) the
    # analytical solution in 2 dimensions is valid:
    # u(x,y) = ln(sqrt((x-xi)^2+(y-eta)^2))/(2 * Pi)
    OgsTest(PROJECTFILE Elliptic/circle_radius_1/circle_1e1_axi.prj)
    OgsTest(PROJECTFILE Elliptic/circle_radius_1/circle_1e2_axi.prj)
    OgsTest(PROJECTFILE Elliptic/circle_radius_1/circle_1e3_axi.prj)
    OgsTest(PROJECTFILE Elliptic/circle_radius_1/circle_1e4_axi.prj)
    OgsTest(PROJECTFILE Elliptic/circle_radius_1/circle_1e5_axi.prj)
    OgsTest(PROJECTFILE Elliptic/circle_radius_1/circle_1e6_axi.prj RUNTIME 4)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e6_with_nodal_sources.prj RUNTIME 30)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e2_volumetricsourceterm.prj)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e3_volumetricsourceterm.prj)
    OgsTest(PROJECTFILE Elliptic/quarter_disc/quarter_disc_neumann.prj)
    OgsTest(PROJECTFILE Elliptic/quarter_circle/quarter_circle_neumann.prj)
    # the analytical solution is: sin(2*Pi*x-Pi/2)*sin(2*Pi*y-Pi/2)
    # the source term in the data array was set to: -2*(2*Pi)^2 * sin(2*Pi*x-Pi/2)*sin(2*Pi*y-Pi/2)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e3_volumetricsourcetermdataarray.prj)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion_Python/square_1e3_laplace_eq.prj)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion_Python/square_1e3_poisson_sin_x_sin_y.prj)
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion_Python/square_1e5_poisson_sin_x_sin_y.prj)
endif()

if(OGS_USE_EIGEN_UNSUPPORTED AND NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Elliptic/square_1x1_SteadyStateDiffusion/square_1e2_GMRES.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e2_3d_submesh_output.xml)
endif()

if(TARGET xdmfdiff AND OGS_USE_EIGEN_UNSUPPORTED AND NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE
            Elliptic/square_1x1_SteadyStateDiffusion/square_1e2_GMRES_GML_output_xdmf-hdf5.prj
    )
endif()

if(OGS_USE_MPI)
    NotebookTest(NOTEBOOKFILE Notebooks/SimplePETSc.py RUNTIME 10)
    OgsTest(PROJECTFILE EllipticPETSc/cube_1x1x1_SteadyStateDiffusion/gml_output/3/cube_hex_27.prj WRAPPER mpirun -np 3)
else()
    NotebookTest(NOTEBOOKFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/ssd-cube.py RUNTIME 6)
endif()
