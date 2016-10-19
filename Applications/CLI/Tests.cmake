
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
            TESTER vtkdiff
            ABSTOL 1e-15 RELTOL 1e-15
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_cube_1x1x1_${mesh_size}_Newton
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}_newton.prj
            TESTER vtkdiff
            ABSTOL 1e-15 RELTOL 1e-15
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_newton_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
            TESTER vtkdiff
            ABSTOL 1e-1 RELTOL 1e-1
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure
        )
    endforeach()

    foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
        AddTest(
            NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}.prj
            TESTER vtkdiff
            ABSTOL 1e-13 RELTOL 1e-13
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
            TESTER vtkdiff
            ABSTOL 1e-2 RELTOL 1e-2
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure
        )
    endforeach()

    # SQUARE 1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
        AddTest(
            NAME GroundWaterFlowProcess_square_1x1_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}.prj
            TESTER vtkdiff
            ABSTOL 1e-13 RELTOL 1e-13
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
            TESTER vtkdiff
            ABSTOL 1e-1 RELTOL 1e-1
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure
        )
    endforeach()

    foreach(mesh_size 1e5 1e6)
        AddTest(
            NAME LARGE_GroundWaterFlowProcess_square_1x1_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}.prj
            TESTER vtkdiff
            ABSTOL 1e-12 RELTOL 1e-16
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME LARGE_GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
            TESTER vtkdiff
            ABSTOL 1e-02 RELTOL 1e-02
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure
        )
    endforeach()

    # LINE 1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e1)
        AddTest(
            NAME GroundWaterFlowProcess_line_1_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}.prj
            TESTER vtkdiff
            ABSTOL 1e-15 RELTOL 1e-15
            DIFF_DATA
            line_1_line_${mesh_size}.vtu line_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_line_1_Neumann_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
            TESTER vtkdiff
            ABSTOL 1e-14 RELTOL 1e-14
            DIFF_DATA
            line_1_line_${mesh_size}.vtu line_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_line_1_Robin_Right_Picard_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}_robin_right_picard.prj
            TESTER vtkdiff
            ABSTOL 4e-14 RELTOL 2e-14
            DIFF_DATA
            line_1_line_${mesh_size}.vtu line_${mesh_size}_robin_right_picard_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_line_1_Robin_Left_Picard_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}_robin_left_picard.prj
            TESTER vtkdiff
            ABSTOL 1e-14 RELTOL 1e-14
            DIFF_DATA
            line_1_line_${mesh_size}.vtu line_${mesh_size}_robin_left_picard_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_line_1_Time_Dep_Dirichlet_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}_time_dep_dirichlet.prj
            TESTER vtkdiff
            ABSTOL 1e-14 RELTOL 1e-14
            DIFF_DATA
            line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_pcs_0_ts_1_t_1.000000.vtu t_1s pressure
            line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_pcs_0_ts_5_t_5.000000.vtu t_5s pressure
            line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_pcs_0_ts_10_t_10.000000.vtu t_10s pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_line_1_Time_Dep_Neumann_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE ogs
            EXECUTABLE_ARGS line_${mesh_size}_time_dep_neumann.prj
            TESTER vtkdiff
            ABSTOL 1e-14 RELTOL 1e-14
            DIFF_DATA
            line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_pcs_0_ts_1_t_1.000000.vtu t_1s pressure
            line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_pcs_0_ts_5_t_5.000000.vtu t_5s pressure
            line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_pcs_0_ts_10_t_10.000000.vtu t_10s pressure
        )
    endforeach()

    # Some Neumann BC tests
    AddTest(
        NAME GroundWaterFlowProcess_cube_top
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3_top_neumann.prj
        TESTER vtkdiff
        ABSTOL 1e-14 RELTOL 1e-14
        DIFF_DATA
        cube_1e3_top_neumann.vtu cube_1e3_top_neumann_pcs_0_ts_1_t_1.000000.vtu pressure pressure
    )
    AddTest(
        NAME GroundWaterFlowProcess_cube_bottom
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3_bottom_neumann.prj
        TESTER vtkdiff
        ABSTOL 1e-14 RELTOL 1e-14
        DIFF_DATA
        cube_1e3_bottom_neumann.vtu cube_1e3_bottom_neumann_pcs_0_ts_1_t_1.000000.vtu pressure pressure
    )
    # Some Neumann BC tests -- Newton
    AddTest(
        NAME GroundWaterFlowProcess_cube_top_Newton
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3_top_neumann_newton.prj
        TESTER vtkdiff
        ABSTOL 1e-14 RELTOL 1e-14
        DIFF_DATA
        cube_1e3_top_neumann.vtu cube_1e3_top_neumann_newton_pcs_0_ts_1_t_1.000000.vtu pressure pressure
    )
    AddTest(
        NAME GroundWaterFlowProcess_cube_bottom_Newton
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3_bottom_neumann_newton.prj
        TESTER vtkdiff
        ABSTOL 1e-14 RELTOL 1e-14
        DIFF_DATA
        cube_1e3_bottom_neumann.vtu cube_1e3_bottom_neumann_newton_pcs_0_ts_1_t_1.000000.vtu pressure pressure
    )

    # test CalculateSurfaceFlux
    AddTest(
        NAME GroundWaterFlowProcess_cube_1x1x1_1e3_dirichlet_calculatesurfaceflux
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3_calculatesurfaceflux.prj
        TESTER vtkdiff
        ABSTOL 1e-15 RELTOL 1e-15
        DIFF_DATA
        cube_1x1x1_hex_1e3_complete_surface.vtu cube_1x1x1_hex_1e3_complete_surface_left_right_dirichlet_surfaceflux.vtu surfaceflux_left_right_dirichlet_reference surfaceflux
    )
    AddTest(
        NAME GroundWaterFlowProcess_cube_1x1x1_1e3_neumann_calculatesurfaceflux
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3_neumann_calculatesurfaceflux.prj
        TESTER vtkdiff
        ABSTOL 1e-15 RELTOL 1e-15
        DIFF_DATA
        cube_1x1x1_hex_1e3_complete_surface.vtu cube_1x1x1_hex_1e3_complete_surface_neumann_surfaceflux.vtu surfaceflux_neumann_reference surfaceflux
    )
    # TES tests
    AddTest(
        NAME TES_zeolite_discharge_small
        PATH Parabolic/TES/1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-small.prj
        TESTER vtkdiff
        ABSTOL 1e-7 RELTOL 5e-9
        DIFF_DATA
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu pressure pressure
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu temperature temperature
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu vapour_partial_pressure vapour_partial_pressure
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu solid_density solid_density
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu solid_density fct_solid_density
    )

    AddTest(
        NAME LARGE_TES_zeolite_discharge
        PATH Parabolic/TES/1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-large.prj
        TESTER vtkdiff
        ABSTOL 1e-3 RELTOL 1e-4
        DIFF_DATA
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu pressure pressure
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu temperature temperature
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu vapour_partial_pressure vapour_partial_pressure
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu solid_density solid_density
        tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu solid_density fct_solid_density
    )

    AddTest(
        NAME LARGE_TES_zeolite_discharge_Newton
        PATH Parabolic/TES/1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-small-newton.prj
        TESTER vtkdiff
        ABSTOL 1.5e-3 RELTOL 1.5e-3
        DIFF_DATA
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu pressure pressure
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu temperature temperature
        # tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu vapour_partial_pressure vapour_partial_pressure
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu solid_density solid_density
    )

    AddTest(
         NAME 1D_HeatConduction_dirichlet
         PATH Parabolic/T/1D_dirichlet
         EXECUTABLE ogs
         EXECUTABLE_ARGS line_60_heat.prj
         TESTER vtkdiff
         ABSTOL 1e-5 RELTOL 1e-5
         DIFF_DATA
         temperature_analytical.vtu line_60_heat_pcs_0_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature
         temperature_analytical.vtu line_60_heat_pcs_0_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature
    )

    AddTest(
         NAME 1D_HeatConduction_neumann
         PATH Parabolic/T/1D_neumann
         EXECUTABLE ogs
         EXECUTABLE_ARGS line_60_heat.prj
         TESTER vtkdiff
         ABSTOL 1e-4 RELTOL 1e-4
         DIFF_DATA
         temperature_analytical.vtu line_60_heat_pcs_0_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature
         temperature_analytical.vtu line_60_heat_pcs_0_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature
    )
    AddTest(
         NAME 2D_RichardsFlow_h_us_quad
         PATH Parabolic/Richards
         EXECUTABLE ogs
         EXECUTABLE_ARGS RichardsFlow_2d.prj
         TESTER vtkdiff
         ABSTOL 1e-1 RELTOL 1e-1
         DIFF_DATA
         h_us_quad_1000.vtu richards_pcs_0_ts_100_t_100.000000.vtu PRESSURE1 pressure
    )
    # Mechanics; Small deformations, linear (SDL)
    AddTest(
        NAME Mechanics_SDL_square_1e0_displacementBC
        PATH Mechanics/Linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e0.prj
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        square_1e0_expected_pcs_0_ts_4_t_1.000000.vtu square_1e0_pcs_0_ts_4_t_1.000000.vtu displacement displacement
    )
    AddTest(
        NAME Mechanics_SDL_square_1e2_tractionBC
        PATH Mechanics/Linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e2.prj
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        square_1e2_expected_pcs_0_ts_4_t_1.000000.vtu square_1e2_pcs_0_ts_4_t_1.000000.vtu displacement displacement
    )
    AddTest(
        NAME LARGE_Mechanics_SDL_disc_with_hole
        PATH Mechanics/Linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS disc_with_hole.prj
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        disc_with_hole_expected_pcs_0_ts_4_t_1.000000.vtu disc_with_hole_pcs_0_ts_4_t_1.000000.vtu displacement displacement
    )
    AddTest(
        NAME LARGE_Mechanics_SDL_square_1e5_tractionBC
        PATH Mechanics/Linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e5.prj
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        square_1e5_expected_pcs_0_ts_4_t_1.000000.vtu square_1e5_pcs_0_ts_4_t_1.000000.vtu displacement displacement
    )
    AddTest(
        NAME Mechanics_SDL_square_1e2_quad8_traction_topBC
        PATH Mechanics/Linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e2_quad8_traction_top.prj
        WRAPPER time
        TESTER vtkdiff
        ABSTOL 5e-16 RELTOL 1e-15
        DIFF_DATA
        expected_square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu displacement displacement
    )

    # Mechanics; Small deformations, Burgers (SDB)
    AddTest(
        NAME Mechanics_SDB_cube_1e0_tractionBC
        PATH Mechanics/Burgers
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e0.prj
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        cube_1e0_expected_pcs_0_ts_1_t_0.000100.vtu cube_1e0_pcs_0_ts_1_t_0.000100.vtu displacement displacement
        cube_1e0_expected_pcs_0_ts_101_t_1.000000.vtu cube_1e0_pcs_0_ts_101_t_1.000000.vtu displacement displacement
    )
    AddTest(
        NAME LARGE_Mechanics_SDB_cube_1e3_tractionBC
        PATH Mechanics/Burgers
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_1e3.prj
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        cube_1e3_expected_pcs_0_ts_1_t_0.000100.vtu cube_1e3_pcs_0_ts_1_t_0.000100.vtu displacement displacement
        cube_1e3_expected_pcs_0_ts_101_t_1.000000.vtu cube_1e3_pcs_0_ts_101_t_1.000000.vtu displacement displacement
    )

    # SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC
    # test results are compared to 3D simulation on a wedge-shaped domain
    AddTest(
        NAME GroundWaterFlowProcess_square_1x1_1e2_axi
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e2_axi.prj
        TESTER vtkdiff
        ABSTOL 1.6e-5 RELTOL 1e-5
        DIFF_DATA
        wedge-1e2-ang-0.02-surface.vtu square_1e2_axi_pcs_0_ts_1_t_1.000000.vtu temperature temperature
    )
    # # WEDGE 1x1 GROUNDWATER FLOW TEST -- computes reference results for the above test
    # AddTest(
    #     NAME GroundWaterFlowProcess_wedge_1e2_ang_0.02
    #     PATH Elliptic/square_1x1_GroundWaterFlow
    #     EXECUTABLE ogs
    #     EXECUTABLE_ARGS wedge_1e2_axi_ang_0.02.prj
    # )

    # SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC
    # test results are compared to 3D simulation on a wedge-shaped domain
    AddTest(
        NAME GroundWaterFlowProcess_square_1x1_1e4_axi_ang_0.02
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e4_axi.prj
        TESTER vtkdiff
        ABSTOL 1.6e-5 RELTOL 1e-5
        DIFF_DATA
        wedge-1e4-ang-0.02-surface.vtu square_1e4_axi_pcs_0_ts_1_t_1.000000.vtu temperature temperature
    )
    # # WEDGE 1x1 GROUNDWATER FLOW TEST -- computes reference results for the above test
    # AddTest(
    #     NAME GroundWaterFlowProcess_wedge_1e4_ang_0.02
    #     PATH Elliptic/square_1x1_GroundWaterFlow
    #     EXECUTABLE ogs
    #     EXECUTABLE_ARGS wedge_1e4_axi_ang_0.02.prj
    # )

    # SMALL DEFORMATION TEST -- AXIALLY SYMMETRIC
    AddTest(
        NAME LARGE_SmallDeformation_ring_plane_strain_axi
        PATH Mechanics/Linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS ring_plane_strain.prj
        TESTER vtkdiff
        ABSTOL 6e-4 RELTOL 1e-4
        DIFF_DATA
        ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu u displacement
        ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu sigma_rr sigma_rr
        ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu sigma_ff sigma_ff
    )

    # SQUARE 1x1 HEAT CONDUCTION TEST -- AXIALLY SYMMETRIC
    # test results are compared to 3D simulation on a wedge-shaped domain
    AddTest(
         NAME 2D_HeatConduction_axi
         PATH Parabolic/T/2D_axially_symmetric
         EXECUTABLE ogs
         EXECUTABLE_ARGS square_1e2_axi.prj
         TESTER vtkdiff
         ABSTOL 1.7e-5 RELTOL 1e-5
         DIFF_DATA
         wedge_1e2_axi_ang_0.02_t_2s_extracted_surface.vtu square_1e2_axi_pcs_0_ts_2_t_2.000000.vtu temperature temperature
         wedge_1e2_axi_ang_0.02_t_2s_extracted_surface.vtu square_1e2_axi_pcs_0_ts_2_t_2.000000.vtu heat_flux_x heat_flux_x
    )
    # # WEDGE 1x1 HEATCONDUCTION TEST -- computes reference results for the above test
    # AddTest(
    #      NAME 2D_HeatConduction_wedge
    #      PATH Parabolic/T/2D_axially_symmetric
    #      EXECUTABLE ogs
    #      EXECUTABLE_ARGS wedge_1e2_axi_ang_0.02.prj
    # )

    # SQUARE 1x1 TES TEST -- AXIALLY SYMMETRIC
    # test results are compared to 3D simulation on a wedge-shaped domain
    AddTest(
        NAME LARGE_TES_inert_axi
        PATH Parabolic/TES/2D
        EXECUTABLE ogs
        EXECUTABLE_ARGS tes-inert-axi.prj
        TESTER vtkdiff
        # Note: Since the temperature and pressure only vary by a factor of ~ 1.e-6 in x-direction
        # the relative tolerance has to be much smaller than 1.e-6
        ABSTOL 1e-12 RELTOL 2e-9
        DIFF_DATA
        inert-wedge-extracted-surface-t-1s.vtu tes_inert_axi_pcs_0_ts_4_t_1.000000.vtu pressure pressure
        inert-wedge-extracted-surface-t-1s.vtu tes_inert_axi_pcs_0_ts_4_t_1.000000.vtu temperature temperature
        inert-wedge-extracted-surface-t-1s.vtu tes_inert_axi_pcs_0_ts_4_t_1.000000.vtu v_mass_frac v_mass_frac
    )
    # # WEDGE 1x1 TES TEST -- computes reference results for the above test
    # AddTest(
    #     NAME TES_inert_wedge
    #     PATH Parabolic/TES/2D
    #     EXECUTABLE ogs
    #     EXECUTABLE_ARGS tes-inert-wedge.prj
    # )

else()
    # MPI groundwater flow tests
    AddTest(
        NAME ParallelFEM_GroundWaterFlow2D
        PATH EllipticPETSc
        EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER vtkdiff
        ABSTOL 1e-15 RELTOL 1e-14
        DIFF_DATA
        quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_0.vtu quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_0.vtu pressure pressure
        quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_1.vtu quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_1.vtu pressure pressure
        quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_2.vtu quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_2.vtu pressure pressure
    )

    AddTest(
        NAME ParallelFEM_GroundWaterFlow3D_DirichletBC
        PATH EllipticPETSc
        EXECUTABLE_ARGS cube_1e3.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER vtkdiff
        ABSTOL 1e-14 RELTOL 1e-6
        DIFF_DATA
        cube_1e3_pcs_0_ts_1_t_1_000000_0.vtu cube_1e3_pcs_0_ts_1_t_1_000000_0.vtu pressure pressure
        cube_1e3_pcs_0_ts_1_t_1_000000_1.vtu cube_1e3_pcs_0_ts_1_t_1_000000_1.vtu pressure pressure
        cube_1e3_pcs_0_ts_1_t_1_000000_2.vtu cube_1e3_pcs_0_ts_1_t_1_000000_2.vtu pressure pressure
    )

    AddTest(
        NAME ParallelFEM_GroundWaterFlow3D_NeumannBC
        PATH EllipticPETSc
        EXECUTABLE_ARGS cube_1e3_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 3
        TESTER vtkdiff
        ABSTOL 1e-14 RELTOL 1e-14
        DIFF_DATA
        cube_1e3_neumann_pcs_0_ts_1_t_1_000000_0.vtu cube_1e3_neumann_pcs_0_ts_1_t_1_000000_0.vtu pressure pressure
        cube_1e3_neumann_pcs_0_ts_1_t_1_000000_1.vtu cube_1e3_neumann_pcs_0_ts_1_t_1_000000_1.vtu pressure pressure
        cube_1e3_neumann_pcs_0_ts_1_t_1_000000_2.vtu cube_1e3_neumann_pcs_0_ts_1_t_1_000000_2.vtu pressure pressure
    )

    # Single core
    # CUBE 1x1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3)
        AddTest(
            NAME GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE_ARGS cube_${mesh_size}.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-15 RELTOL 1e-15
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-1 RELTOL 1e-1
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure
        )
    endforeach()


    foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
        AddTest(
            NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE_ARGS cube_${mesh_size}.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-7 RELTOL 1e-7
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
            PATH Elliptic/cube_1x1x1_GroundWaterFlow
            EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-2 RELTOL 1e-2
            DIFF_DATA
            cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure
        )
    endforeach()

    # SQUARE 1x1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
        AddTest(
            NAME GroundWaterFlowProcess_square_1x1_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE_ARGS square_${mesh_size}.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-13 RELTOL 1e-13
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-1 RELTOL 1e-1
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure
        )
    endforeach()

    foreach(mesh_size 1e5 1e6)
        AddTest(
            NAME LARGE_GroundWaterFlowProcess_square_1x1_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE_ARGS square_${mesh_size}.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-7 RELTOL 1e-7
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME LARGE_GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
            PATH Elliptic/square_1x1_GroundWaterFlow
            EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-02 RELTOL 1e-02
            DIFF_DATA
            square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure
        )
    endforeach()

    # LINE 1 GROUNDWATER FLOW TESTS
    foreach(mesh_size 1e1)
        AddTest(
            NAME GroundWaterFlowProcess_line_1_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE_ARGS line_${mesh_size}.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-15 RELTOL 1e-15
            DIFF_DATA
            line_1_line_${mesh_size}.vtu line_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure
        )

        AddTest(
            NAME GroundWaterFlowProcess_line_1_Neumann_${mesh_size}
            PATH Elliptic/line_1_GroundWaterFlow
            EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
            WRAPPER mpirun
            WRAPPER_ARGS -np 1
            TESTER vtkdiff
            ABSTOL 1e-14 RELTOL 1e-14
            DIFF_DATA
            line_1_line_${mesh_size}.vtu line_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_N1_right pressure
        )
    endforeach()

    AddTest(
        NAME TES_zeolite_discharge_small
        PATH Parabolic/TES/1D
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-small.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        ABSTOL 1e-7 RELTOL 5e-9
        DIFF_DATA
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0_100000_0.vtu pressure pressure
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0_100000_0.vtu temperature temperature
        tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0_100000_0.vtu v_mass_frac v_mass_frac
#        tes_zeolite_discharge_small_ts_19_t_0.100000.vtu solid_density solid_density
    )

    AddTest(
        NAME LARGE_TES_zeolite_discharge
        PATH Parabolic/TES/1D
        EXECUTABLE_ARGS tes-1D-zeolite-discharge-large.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        ABSTOL 1e-16 RELTOL 1e-16
        DIFF_DATA
        tes_zeolite_discharge_large_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000_0.vtu pressure pressure
        tes_zeolite_discharge_large_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000_0.vtu temperature temperature
        tes_zeolite_discharge_large_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000_0.vtu v_mass_frac v_mass_frac
#        tes_zeolite_discharge_large_ts_28_t_1_0.vtu solid_density solid_density
    )
endif()
