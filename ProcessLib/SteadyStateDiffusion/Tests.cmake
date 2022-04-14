# CUBE 1x1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3)
    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}.xml --write-prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    if(TEST ogs-SteadyStateDiffusion_cube_1x1x1_${mesh_size} AND DIFF_TOOL_PATH)
        set(_processed_path Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_${mesh_size}_processed.prj)
        add_test(NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}_prj_diff
            COMMAND ${DIFF_TOOL_PATH} ${Data_SOURCE_DIR}/${_processed_path} ${Data_BINARY_DIR}/${_processed_path})
        set_tests_properties(SteadyStateDiffusion_cube_1x1x1_${mesh_size}_prj_diff
            PROPERTIES LABELS "default" DEPENDS SteadyStateDiffusion_cube_1x1x1_${mesh_size})
    endif()

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}_Newton
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        # `-m .` just for testing input mesh dir parameter:
        EXECUTABLE_ARGS -m . cube_${mesh_size}_newton.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_newton_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure 1e-1 1e-1
    )
endforeach()

foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
    set(RUNTIME 10)
    if("${mesh_size}" STREQUAL "1e6")
        set(RUNTIME 75)
    endif()
    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        RUNTIME ${RUNTIME}
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
    )

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        RUNTIME ${RUNTIME}
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    )
endforeach()

# Quadratic hex element.
AddTest(
    NAME SteadyStateDiffusion_cube_1x1x1_1e0_QuadraticHex
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_quadratic_hex.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex20_1e0.vtu cube_1e0_quadratic_hex_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
)

# SQUARE 1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
    AddTest(
        NAME SteadyStateDiffusion_square_1x1_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS -p square_${mesh_size}.xml square_1e0.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
    )

    AddTest(
        NAME SteadyStateDiffusion_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS -p square_${mesh_size}.xml -p square_neumann.xml square_1e0.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure 1e-1 1e-1
    )
endforeach()

foreach(mesh_size 1e5 1e6)
    set(RUNTIME 30)
    if("${mesh_size}" STREQUAL "1e6")
        set(RUNTIME 90)
    endif()
    AddTest(
        NAME SteadyStateDiffusion_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        RUNTIME ${RUNTIME}
        EXECUTABLE ogs
        EXECUTABLE_ARGS -p square_${mesh_size}.xml -p square_neumann.xml square_1e0.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure 1e-02 1e-02
    )
endforeach()

AddTest(
    NAME SteadyStateDiffusion_square_1x1_1e5
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    RUNTIME 2
    EXECUTABLE ogs
    EXECUTABLE_ARGS -p square_1e5.xml square_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1x1_quad_1e5.vtu square_1e5_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-16
)

# The largest test is less accurate
AddTest(
    NAME SteadyStateDiffusion_square_1x1_1e6
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    RUNTIME 75
    EXECUTABLE ogs
    EXECUTABLE_ARGS -p square_1e6.xml square_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1x1_quad_1e6.vtu square_1e6_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 3e-12 1e-16
)

# LINE 1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e1)
    AddTest(
        NAME SteadyStateDiffusion_line_1_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME SteadyStateDiffusion_line_1_Neumann_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_neumann_ts_1_t_1.000000.vtu D1_left_N1_right pressure 1e-14 1e-14
    )
    if (NOT OGS_USE_MPI)
        OgsTest(PROJECTFILE Elliptic/cube_1x1x1_SteadyStateDiffusion/drainage_excavation.prj)
    endif()
    AddTest(
        NAME SteadyStateDiffusion_line_1_Robin_Right_Picard_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_robin_right_picard.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_robin_right_picard_ts_1_t_1.000000.vtu D1_left_N1_right pressure 4e-14 2e-14
    )

    AddTest(
        NAME SteadyStateDiffusion_line_1_Robin_Left_Picard_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_robin_left_picard.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_robin_left_picard_ts_1_t_1.000000.vtu D1_left_N1_right pressure 1e-14 1e-14
    )

    AddTest(
        NAME SteadyStateDiffusion_line_1_Time_Dep_Dirichlet_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_time_dep_dirichlet.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_ts_1_t_1.000000.vtu t_1s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_ts_5_t_5.000000.vtu t_5s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_ts_10_t_10.000000.vtu t_10s pressure 1e-14 1e-14
    )

    AddTest(
        NAME SteadyStateDiffusion_line_1_Time_Dep_Neumann_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_time_dep_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_ts_1_t_1.000000.vtu t_1s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_ts_5_t_5.000000.vtu t_5s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_ts_10_t_10.000000.vtu t_10s pressure 1e-14 1e-14
    )
endforeach()

# Some Neumann BC tests
AddTest(
    NAME SteadyStateDiffusion_cube_top
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_top_neumann.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_top_neumann.vtu cube_1e3_top_neumann_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)
AddTest(
    NAME SteadyStateDiffusion_cube_bottom
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_bottom_neumann.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_bottom_neumann.vtu cube_1e3_bottom_neumann_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)
# Some Neumann BC tests -- Newton
AddTest(
    NAME SteadyStateDiffusion_cube_top_Newton
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_top_neumann_newton.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_top_neumann.vtu cube_1e3_top_neumann_newton_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)
AddTest(
    NAME SteadyStateDiffusion_cube_bottom_Newton
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_bottom_neumann_newton.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_bottom_neumann.vtu cube_1e3_bottom_neumann_newton_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)

# test SurfaceFlux
AddTest(
    NAME SteadyStateDiffusion_cube_1x1x1_1e3_dirichlet_calculatesurfaceflux
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_calculatesurfaceflux.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_complete_surface_left_right_dirichlet_specific_flux_t_1.000000.vtu cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_1_t_1.000000.vtu specific_flux specific_flux 5e-15 5e-15
)
AddTest(
    NAME SteadyStateDiffusion_cube_1x1x1_1e3_neumann_calculatesurfaceflux
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_neumann_calculatesurfaceflux.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_complete_surface_neumann_specific_flux_t_1.000000.vtu cube_1e3_neumann_balance_cube_1x1x1_hex_1e3_complete_surface_ts_1_t_1.000000.vtu specific_flux specific_flux 2e-14 2e-14
)
AddTest(
    NAME SteadyStateDiffusion_cube_1x1x1_2e3_prism_surfaceflux_left_right
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_2e3_prism_surfaceflux_left_right.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_prism_2e3_complete_surface_left_right_dirichlet_specific_flux_t_1.000000.vtu cube_2e3_surface_flux_left_right_cube_1x1x1_prism_2e3_complete_surface_ts_1_t_1.000000.vtu specific_flux specific_flux 1e-14 1e-14
)

AddTest(
    NAME SteadyStateDiffusion_cube_1x1x1_2e3_prism_surfaceflux_front_back
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_2e3_prism_surfaceflux_front_back.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_prism_2e3_complete_surface_front_back_dirichlet_specific_flux_t_1.000000.vtu cube_2e3_surface_flux_front_back_cube_1x1x1_prism_2e3_complete_surface_ts_1_t_1.000000.vtu specific_flux specific_flux  1e-14 1e-14
)

AddTest(
    NAME SteadyStateDiffusion_cube_1x1x1_2e3_prism_surfaceflux_top_bottom
    PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_2e3_prism_surfaceflux_top_bottom.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_prism_2e3_complete_surface_top_bottom_dirichlet_specific_flux_t_1.000000.vtu cube_2e3_surface_flux_top_bottom_cube_1x1x1_prism_2e3_complete_surface_ts_1_t_1.000000.vtu specific_flux specific_flux 1e-14 1e-14
)

AddTest(
    NAME SteadyStateDiffusion_wedge_1x1x1_1e3_prism_surfaceflux
    PATH Elliptic/wedge_1x1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS wedge_1e3_prism_surfaceflux_diagonal.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    wedge_1x1x1_1e3_prism_complete_surface_specific_flux_t_1.000000.vtu wedge_1e3_surfaceflux_wedge_1x1x1_1e3_prism_complete_surface_ts_1_t_1.000000.vtu specific_flux specific_flux 2e-14 0
)

# SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC
# test results are compared to 3D simulation on a wedge-shaped domain
AddTest(
    NAME SteadyStateDiffusion_square_1x1_1e2_axi
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_axi.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    wedge-1e2-ang-0.02-surface.vtu square_1e2_axi_ts_1_t_1.000000.vtu temperature temperature 1.6e-5 1e-5
)
# # WEDGE 1x1 GROUNDWATER FLOW TEST -- computes reference results for the above test
# AddTest(
#     NAME SteadyStateDiffusion_wedge_1e2_ang_0.02
#     PATH Elliptic/square_1x1_SteadyStateDiffusion
#     EXECUTABLE ogs
#     EXECUTABLE_ARGS wedge_1e2_axi_ang_0.02.prj
# )

# SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC
# test results are compared to 3D simulation on a wedge-shaped domain
AddTest(
    NAME SteadyStateDiffusion_square_1x1_1e4_axi_ang_0.02
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e4_axi.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    wedge-1e4-ang-0.02-surface.vtu square_1e4_axi_ts_1_t_1.000000.vtu temperature temperature 1.6e-5 1e-5
)
# # WEDGE 1x1 GROUNDWATER FLOW TEST -- computes reference results for the above test
# AddTest(
#     NAME SteadyStateDiffusion_wedge_1e4_ang_0.02
#     PATH Elliptic/square_1x1_SteadyStateDiffusion
#     EXECUTABLE ogs
#     EXECUTABLE_ARGS wedge_1e4_axi_ang_0.02.prj
# )

# MPI groundwater flow tests
AddTest(
    NAME ParallelFEM_GroundWaterFlow2D
    PATH EllipticPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    quad_20x10_GroundWaterFlow_result_ts_1_t_1_000000_0.vtu quad_20x10_GroundWaterFlow_result_ts_1_t_1_000000_0.vtu pressure pressure 2e-15 1e-16
    quad_20x10_GroundWaterFlow_result_ts_1_t_1_000000_1.vtu quad_20x10_GroundWaterFlow_result_ts_1_t_1_000000_1.vtu pressure pressure 2e-15 1e-16
    quad_20x10_GroundWaterFlow_result_ts_1_t_1_000000_2.vtu quad_20x10_GroundWaterFlow_result_ts_1_t_1_000000_2.vtu pressure pressure 2e-15 1e-16
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_DirichletBC
    PATH EllipticPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_ts_1_t_1_000000_0.vtu cube_1e3_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 2e-15 1e-16
    cube_1e3_ts_1_t_1_000000_1.vtu cube_1e3_ts_1_t_1_000000_1.vtu Linear_1_to_minus1 pressure 2e-15 1e-16
    cube_1e3_ts_1_t_1_000000_2.vtu cube_1e3_ts_1_t_1_000000_2.vtu Linear_1_to_minus1 pressure 2e-15 1e-16
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_NeumannBC
    PATH EllipticPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_neumann.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_neumann_ts_1_t_1_000000_0.vtu cube_1e3_neumann_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    cube_1e3_neumann_ts_1_t_1_000000_1.vtu cube_1e3_neumann_ts_1_t_1_000000_1.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    cube_1e3_neumann_ts_1_t_1_000000_2.vtu cube_1e3_neumann_ts_1_t_1_000000_2.vtu D1_left_front_N1_right pressure 1e-2 1e-2
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_NeumannBC_XDMF_np3_1file
    PATH EllipticPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_XDMF_np3.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER xdmfdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_np3.xdmf cube_1e3_np3_cube_1x1x1_hex_1e3.xdmf pressure pressure 1e-3 1e-3
    cube_1e3_np3.xdmf cube_1e3_np3_cube_1x1x1_hex_1e3.xdmf v v 1e-3 1e-3
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_NeumannBC_XDMF_np3_2files
    PATH EllipticPETSc/XDMF_NP3_2
    EXECUTABLE ogs
    EXECUTABLE_ARGS ../cube_1e3_XDMF_np3_2files.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER xdmfdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_np3_2files_0.xdmf cube_1e3_np3_cube_1x1x1_hex_1e3.xdmf pressure pressure 1e-3 1e-3
    cube_1e3_np3_2files_0.xdmf cube_1e3_np3_cube_1x1x1_hex_1e3.xdmf v v 1e-3 1e-3
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_NeumannBC_XDMF_np3_3files
    PATH EllipticPETSc/XDMF_NP3_3
    EXECUTABLE ogs
    EXECUTABLE_ARGS ../cube_1e3_XDMF_np3_3files.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER xdmfdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_np3_3files_0.xdmf cube_1e3_np3_cube_1x1x1_hex_1e3.xdmf pressure pressure 1e-3 1e-3
    cube_1e3_np3_3files_0.xdmf cube_1e3_np3_cube_1x1x1_hex_1e3.xdmf v v 1e-3 1e-3
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_NeumannBC_XDMF_np2
    PATH EllipticPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_XDMF_np2.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER xdmfdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_np2.xdmf cube_1e3_np2_cube_1x1x1_hex_1e3.xdmf pressure pressure 1e-3 1e-3
    cube_1e3_np2.xdmf cube_1e3_np2_cube_1x1x1_hex_1e3.xdmf v v 1e-3 1e-3
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow2D_NeumannBC
    PATH EllipticPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e1_neumann.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    square_1e1_neumann_ts_1_t_1_000000_0.vtu square_1e1_neumann_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure 1e-2 0
    square_1e1_neumann_ts_1_t_1_000000_1.vtu square_1e1_neumann_ts_1_t_1_000000_1.vtu D1_left_bottom_N1_right pressure 1e-2 0
    square_1e1_neumann_ts_1_t_1_000000_0.vtu square_1e1_neumann_ts_1_t_1_000000_0.vtu pressure pressure 1e-14 0
    square_1e1_neumann_ts_1_t_1_000000_1.vtu square_1e1_neumann_ts_1_t_1_000000_1.vtu pressure pressure 1e-14 0
)

AddTest(
    NAME ParallelFEM_SteadyStateDiffusion_cube_2
    PATH EllipticPETSc/cube_1x1x1_SteadyStateDiffusion/2
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e4_anisotropic.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e4_anisotropic_ts_1_t_1_000000_0.vtu cube_1e4_anisotropic_ts_1_t_1_000000_0.vtu pressure pressure 1e-14 0
    cube_1e4_anisotropic_ts_1_t_1_000000_1.vtu cube_1e4_anisotropic_ts_1_t_1_000000_1.vtu pressure pressure 1e-14 0
)
#OgsTest(PROJECTFILE "EllipticPETSc/cube_1x1x1_SteadyStateDiffusion/2/cube_1e4_anisotropic.prj")

# Single core
# CUBE 1x1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3)
    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}.xml
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure 1e-1 1e-1
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
    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        RUNTIME ${RUNTIME}
        EXECUTABLE_ARGS cube_${mesh_size}.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-7 1e-7
    )

    AddTest(
        NAME SteadyStateDiffusion_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        RUNTIME ${RUNTIME}
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    )
endforeach()

# SQUARE 1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
    AddTest(
        NAME SteadyStateDiffusion_square_1x1_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS -p square_${mesh_size}.xml square_1e0.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
    )

    AddTest(
        NAME SteadyStateDiffusion_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS -p square_${mesh_size}.xml -p square_neumann.xml square_1e0.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure 1e-1 1e-1
    )
endforeach()

foreach(mesh_size 1e5 1e6)
    set(RUNTIME 65)
    if("${mesh_size}" STREQUAL "1e6")
        set(RUNTIME 450)
    endif()
    AddTest(
        NAME SteadyStateDiffusion_square_1x1_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        RUNTIME ${RUNTIME}
        EXECUTABLE_ARGS -p square_${mesh_size}.xml square_1e0.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-7 1e-7
    )

    AddTest(
        NAME SteadyStateDiffusion_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_SteadyStateDiffusion
        EXECUTABLE ogs
        RUNTIME ${RUNTIME}
        EXECUTABLE_ARGS -p square_${mesh_size}.xml -p square_neumann.xml square_1e0.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure 1e-02 1e-02
    )
endforeach()

# LINE 1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e1)
    AddTest(
        NAME SteadyStateDiffusion_line_1_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME SteadyStateDiffusion_line_1_Neumann_${mesh_size}
        PATH Elliptic/line_1_SteadyStateDiffusion
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_neumann_ts_1_t_1_000000_0.vtu D1_left_N1_right pressure 1e-14 1e-14
    )
endforeach()

AddTest(
    NAME SteadyStateDiffusion_Inhomogeneous_permeability
    PATH Elliptic/nonuniform_bc_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS inhomogeneous_permeability.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    inhomogeneous_permeability.vtu inhomogeneous_permeability_ts_1_t_1.000000.vtu mass_flux_ref mass_flux 4e-2 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_Neumann_nonuniform_cosY
    PATH Elliptic/nonuniform_bc_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS neumann_nonuniform.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_neumann_nonuniform_ts_1_t_1.000000.vtu neumann_nonuniform_ts_1_t_1.000000.vtu pressure pressure 2e-14 0
    expected_neumann_nonuniform_ts_1_t_1.000000.vtu neumann_nonuniform_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-12 0
)

AddTest(
    NAME SteadyStateDiffusion_Dirichlet_nonuniform_linearY
    PATH Elliptic/nonuniform_bc_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS dirichlet_nonuniform.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_dirichlet_nonuniform_ts_1_t_1.000000.vtu dirichlet_nonuniform_ts_1_t_1.000000.vtu pressure pressure 1e-14 0
)

# tests for nodal source term implementation
# For the special setup with a dirac source term at position (xi, eta) the
# analytical solution in 2 dimensions is valid:
# u(x,y) = ln(sqrt((x-xi)^2+(y-eta)^2))/(2 * Pi)
AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_circle_1e1
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e1_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e1_expected.vtu circle_1e1_axi_ts_1_t_1.000000.vtu analytical_solution pressure 0.7 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_circle_1e2
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e2_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e2_expected.vtu circle_1e2_axi_ts_1_t_1.000000.vtu analytical_solution pressure 1.1 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_circle_1e3
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e3_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e3_expected.vtu circle_1e3_axi_ts_1_t_1.000000.vtu analytical_solution pressure 1.6 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_circle_1e4
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e4_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e4_expected.vtu circle_1e4_axi_ts_1_t_1.000000.vtu analytical_solution pressure 1.8 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_circle_1e5
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e5_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e5_expected.vtu circle_1e5_axi_ts_1_t_1.000000.vtu analytical_solution pressure 2.15 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_circle_1e6
    PATH Elliptic/circle_radius_1
    RUNTIME 20
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e6_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e6_expected.vtu circle_1e6_axi_ts_1_t_1.000000.vtu analytical_solution pressure 2.52 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_NodalSourceTerm_square_1e6
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    RUNTIME 200
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e6_with_nodal_sources.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e6_nodal_sources_expected.vtu square_1e6_with_nodal_sources_ts_1_t_1.000000.vtu analytical_solution pressure 1.41 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_VolumetricSourceTerm_square_1e2
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_volumetricsourceterm.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e2_volumetricsourceterm_analytical_solution.vtu square_1e2_volumetricsourceterm_ts_1_t_1.000000.vtu analytical_solution pressure 1e-14 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_VolumetricSourceTerm_square_1e3
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e3_volumetricsourceterm.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e3_volumetricsourceterm_analytical_solution.vtu square_1e3_volumetricsourceterm_ts_1_t_1.000000.vtu analytical_solution pressure 1e-10 1e-11
)

AddTest(
    NAME SteadyStateDiffusion_NeumannBC_Along_Line_in_3D_domain
    PATH Elliptic/quarter_disc
    EXECUTABLE ogs
    EXECUTABLE_ARGS quarter_disc_neumann.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    quarter_disc_r_1.vtu neumann_along_line_ts_1_t_1.000000.vtu analytical_solution pressure 6e-5 0
)

AddTest(
    NAME SteadyStateDiffusion_NeumannBC_In_Center_Point_2D_domain
    PATH Elliptic/quarter_circle
    EXECUTABLE ogs
    EXECUTABLE_ARGS quarter_circle_neumann.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    quarter_circle_r_1.vtu neumann_in_center_ts_1_t_1.000000.vtu pressure_nodal_source_term pressure 1e-14 1e-14
)

AddTest(
    NAME SteadyStateDiffusion_VolumetricSourceTerm_sin_x_sin_y_square_1e3
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e3_volumetricsourcetermdataarray.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    # the analytical solution is: sin(2*Pi*x-Pi/2)*sin(2*Pi*y-Pi/2)
    # the source term in the data array was set to: -2*(2*Pi)^2 * sin(2*Pi*x-Pi/2)*sin(2*Pi*y-Pi/2)
    DIFF_DATA
    square_1x1_quad_1e3_volumetricsourcetermdataarray.vtu square_1e3_volumetricsourcetermdataarray_ts_1_t_1.000000.vtu analytical_solution pressure 2e-2 1e-16
)

AddTest(
    NAME PythonBCSteadyStateDiffusionLaplaceEqDirichletNeumann
    PATH Elliptic/square_1x1_SteadyStateDiffusion_Python
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e3_laplace_eq.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_PYTHON AND NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    python_laplace_eq_ref.vtu square_1e3_neumann_ts_1_t_1.000000.vtu pressure_expected pressure 4e-4 1e-16
)

AddTest(
    NAME PythonSourceTermPoissonSinAXSinBYDirichlet_square_1e3
    PATH Elliptic/square_1x1_SteadyStateDiffusion_Python
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e3_poisson_sin_x_sin_y.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_PYTHON AND NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e3.vtu square_1e3_volumetricsourceterm_ts_1_t_1.000000.vtu analytical_solution pressure 0.7e-2 1e-16
)

AddTest(
    NAME PythonSourceTermPoissonSinAXSinBYDirichlet_square_1e5
    PATH Elliptic/square_1x1_SteadyStateDiffusion_Python
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e5_poisson_sin_x_sin_y.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_PYTHON AND NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e5.vtu square_1e5_volumetricsourceterm_ts_1_t_1.000000.vtu analytical_solution pressure 0.75e-4 1e-16
)

AddTest(
    NAME SteadyStateDiffusion_square_1x1_1e2_GMRES
    PATH Elliptic/square_1x1_SteadyStateDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_GMRES.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e2.vtu square_1e2_GMRES_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-14 1e-14
)

if(OGS_USE_MPI)
    NotebookTest(NOTEBOOKFILE Notebooks/SimplePETSc.ipynb RUNTIME 10)
else()
    OgsTest(PROJECTFILE "Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e4_anisotropic.prj")
endif() # OGS_USE_MPI
