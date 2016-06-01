
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

if(NOT OGS_USE_MPI)
    # CUBE 1x1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3)
        AddTest(
            NAME GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA cube_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
            DATA cube_${mesh_size}.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
        )

        AddTest(
            NAME GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA cube_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure
            DATA cube_${mesh_size}_neumann.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
        )
    endforeach()

    foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
        AddTest(
            NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA cube_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
            DATA cube_${mesh_size}.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
        )

        AddTest(
            NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA cube_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure
            DATA cube_${mesh_size}_neumann.prj cube_1x1x1_hex_${mesh_size}.vtu cube_1x1x1.gml
        )
    endforeach()

    # SQUARE 1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
        AddTest(
            NAME GroundWaterFlowProcess_square_1x1_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA square_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
            DATA square_${mesh_size}.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
        )

        AddTest(
            NAME GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure
            DATA square_${mesh_size}_neumann.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
        )
    endforeach()

    foreach(mesh_size 1e5 1e6)
        AddTest(
            NAME LARGE_GroundWaterFlowProcess_square_1x1_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA square_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
            DATA square_${mesh_size}.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
        )

        AddTest(
            NAME LARGE_GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure
            DATA square_${mesh_size}_neumann.prj square_1x1_quad_${mesh_size}.vtu square_1x1.gml
        )
    endforeach()

    # LINE 1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e1)
        AddTest(
            NAME GroundWaterFlowProcess_line_1_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}.prj
            WRAPPER time
            TESTER vtkdiff
            ABSTOL 1e-16 RELTOL 1e-16
            DIFF_DATA line_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
            DATA line_${mesh_size}.prj line_1_line_${mesh_size}.vtu line_1.gml
        )

        AddTest(
                    NAME GroundWaterFlowProcess_line_1_Neumann_${mesh_size}
                    PATH Elliptic/line_1_GroundWaterFlow
                    EXECUTABLE ogs
                    EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
                    WRAPPER time
                    TESTER vtkdiff
                    ABSTOL 1e-16 RELTOL 1e-16
                    DIFF_DATA line_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure
                    DATA line_${mesh_size}_neumann.prj line_1_line_${mesh_size}.vtu line_1.gml
                )
        endforeach()


    AddTest(
        NAME TES_zeolite_discharge_small
        PATH Parabolic/TES/1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-small.prj
        WRAPPER time
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu pressure pressure
        tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu temperature temperature
        tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu v_mass_frac v_mass_frac
        tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu solid_density solid_density
        DATA line_0.1.gml line_0.1_37.msh tes-1D-zeolite-discharge.prj
    )

    AddTest(
        NAME LARGE_TES_zeolite_discharge
        PATH Parabolic/TES/1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-large.prj
        WRAPPER time
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu pressure pressure
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu temperature temperature
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu v_mass_frac v_mass_frac
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu solid_density solid_density
        DATA line_0.1.gml line_0.1_100.msh tes-1D-zeolite-discharge.prj
    )

else()
    # MPI groundwater flow tests
    AddTest(
        NAME ParallelFEM_GroundWaterFlow2D
        PATH EllipticPETSc
        EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER diff
        DIFF_DATA
            quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1.000000_0.vtu
            quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1.000000_1.vtu
            quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1.000000_2.vtu
    )

    AddTest(
        NAME ParallelFEM_GroundWaterFlow3D_DirichletBC
        PATH EllipticPETSc
        EXECUTABLE_ARGS cube_1e3.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER diff
        DIFF_DATA
            cube_1e3_pcs_0_ts_1_t_1.000000_0.vtu
            cube_1e3_pcs_0_ts_1_t_1.000000_1.vtu
            cube_1e3_pcs_0_ts_1_t_1.000000_2.vtu
    )

    AddTest(
        NAME ParallelFEM_GroundWaterFlow3D_NeumannBC
        PATH EllipticPETSc
        EXECUTABLE_ARGS cube_1e3_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER diff
        DIFF_DATA
            cube_1e3_neumann_pcs_0_ts_1_t_1.000000_0.vtu
            cube_1e3_neumann_pcs_0_ts_1_t_1.000000_1.vtu
            cube_1e3_neumann_pcs_0_ts_1_t_1.000000_2.vtu
    )
endif()
