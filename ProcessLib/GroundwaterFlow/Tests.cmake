# CUBE 1x1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3)
    AddTest(
        NAME GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME GroundWaterFlowProcess_cube_1x1x1_${mesh_size}_Newton
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_newton.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_newton_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure 1e-1 1e-1
    )
endforeach()

foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
    if("${mesh_size}" STREQUAL "1e6")
        set(benchmark_type LARGE_)
    endif()
    AddTest(
        NAME ${benchmark_type}GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
    )

    AddTest(
        NAME ${benchmark_type}GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    )
endforeach()

# Quadratic hex element.
AddTest(
    NAME GroundWaterFlowProcess_cube_1x1x1_1e0_QuadraticHex
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_quadratic_hex.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex20_1e0.vtu cube_1e0_quadratic_hex_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
)

# SQUARE 1x1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e0 1e1 1e2 1e3 1e4)
    AddTest(
        NAME GroundWaterFlowProcess_square_1x1_${mesh_size}
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_${mesh_size}.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
        VIS square_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu
    )

    AddTest(
        NAME GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure 1e-1 1e-1
        VIS square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu
    )
endforeach()

foreach(mesh_size 1e5 1e6)
    if("${mesh_size}" STREQUAL "1e6")
        set(benchmark_type LARGE_)
    endif()
    AddTest(
        NAME ${benchmark_type}GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_bottom_N1_right pressure 1e-02 1e-02
    )
endforeach()

AddTest(
    NAME LARGE_GroundWaterFlowProcess_square_1x1_1e5
    PATH Elliptic/square_1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e5.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1x1_quad_1e5.vtu square_1e5_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-13 1e-16
)

# The largest test is less accurate
AddTest(
    NAME LARGE_GroundWaterFlowProcess_square_1x1_1e6
    PATH Elliptic/square_1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e6.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1x1_quad_1e6.vtu square_1e6_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 3e-12 1e-16
)

# LINE 1 GROUNDWATER FLOW TESTS
foreach(mesh_size 1e1)
    AddTest(
        NAME GroundWaterFlowProcess_line_1_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME GroundWaterFlowProcess_line_1_Neumann_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_neumann_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure 1e-14 1e-14
    )

    AddTest(
        NAME GroundWaterFlowProcess_line_1_Robin_Right_Picard_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_robin_right_picard.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_robin_right_picard_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure 4e-14 2e-14
    )

    AddTest(
        NAME GroundWaterFlowProcess_line_1_Robin_Left_Picard_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_robin_left_picard.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_robin_left_picard_pcs_0_ts_1_t_1.000000.vtu D1_left_N1_right pressure 1e-14 1e-14
    )

    AddTest(
        NAME GroundWaterFlowProcess_line_1_Time_Dep_Dirichlet_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_time_dep_dirichlet.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_pcs_0_ts_1_t_1.000000.vtu t_1s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_pcs_0_ts_5_t_5.000000.vtu t_5s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_dirichlet_pcs_0_ts_10_t_10.000000.vtu t_10s pressure 1e-14 1e-14
    )

    AddTest(
        NAME GroundWaterFlowProcess_line_1_Time_Dep_Neumann_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_${mesh_size}_time_dep_neumann.prj
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_pcs_0_ts_1_t_1.000000.vtu t_1s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_pcs_0_ts_5_t_5.000000.vtu t_5s pressure 1e-14 1e-14
        line_1_time_dep_dirichlet.vtu line_${mesh_size}_time_dep_neumann_pcs_0_ts_10_t_10.000000.vtu t_10s pressure 1e-14 1e-14
    )
endforeach()

# Some Neumann BC tests
AddTest(
    NAME GroundWaterFlowProcess_cube_top
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_top_neumann.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_top_neumann.vtu cube_1e3_top_neumann_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)
AddTest(
    NAME GroundWaterFlowProcess_cube_bottom
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_bottom_neumann.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_bottom_neumann.vtu cube_1e3_bottom_neumann_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)
# Some Neumann BC tests -- Newton
AddTest(
    NAME GroundWaterFlowProcess_cube_top_Newton
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_top_neumann_newton.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_top_neumann.vtu cube_1e3_top_neumann_newton_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)
AddTest(
    NAME GroundWaterFlowProcess_cube_bottom_Newton
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_bottom_neumann_newton.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_bottom_neumann.vtu cube_1e3_bottom_neumann_newton_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-14 1e-14
)

# test CalculateSurfaceFlux
AddTest(
    NAME GroundWaterFlowProcess_cube_1x1x1_1e3_dirichlet_calculatesurfaceflux
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_calculatesurfaceflux.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_complete_surface.vtu cube_1x1x1_hex_1e3_complete_surface_left_right_dirichlet_surfaceflux.vtu surfaceflux_left_right_dirichlet_reference surfaceflux 1e-15 1e-15
)
AddTest(
    NAME GroundWaterFlowProcess_cube_1x1x1_1e3_neumann_calculatesurfaceflux
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_neumann_calculatesurfaceflux.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_complete_surface.vtu cube_1x1x1_hex_1e3_complete_surface_neumann_surfaceflux.vtu surfaceflux_neumann_reference surfaceflux 1e-15 1e-15
)
AddTest(
    NAME GroundWaterFlowProcess_cube_1x1x1_2e3_prism_surfaceflux_left_right
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_2e3_prism_surfaceflux_left_right.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_prism_2e3_complete_surface.vtu cube_1x1x1_prism_2e3_complete_surface_left_right_dirichlet_surfaceflux.vtu surfaceflux_left_right_reference surfaceflux_left_right 1e-15 1e-15
)

AddTest(
    NAME GroundWaterFlowProcess_cube_1x1x1_2e3_prism_surfaceflux_front_back
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_2e3_prism_surfaceflux_front_back.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_prism_2e3_complete_surface.vtu cube_1x1x1_prism_2e3_complete_surface_front_back_dirichlet_surfaceflux.vtu surfaceflux_front_back_reference surfaceflux_front_back 1e-15 1e-15
)

AddTest(
    NAME GroundWaterFlowProcess_cube_1x1x1_2e3_prism_surfaceflux_top_bottom
    PATH Elliptic/cube_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_2e3_prism_surfaceflux_top_bottom.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_prism_2e3_complete_surface.vtu cube_1x1x1_prism_2e3_complete_surface_top_bottom_dirichlet_surfaceflux.vtu surfaceflux_top_bottom_reference surfaceflux_top_bottom 1e-15 1e-15
)

AddTest(
    NAME GroundWaterFlowProcess_wedge_1x1x1_1e3_prism_surfaceflux
    PATH Elliptic/wedge_1x1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS wedge_1e3_prism_surfaceflux_diagonal.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    wedge_1x1x1_1e3_prism_complete_surface.vtu wedge_1x1x1_1e3_prism_complete_surface_surfaceflux.vtu surfaceflux_reference surfaceflux 1e-15 1e-15
)

# SQUARE 1x1 GROUNDWATER FLOW TEST -- AXIALLY SYMMETRIC
# test results are compared to 3D simulation on a wedge-shaped domain
AddTest(
    NAME GroundWaterFlowProcess_square_1x1_1e2_axi
    PATH Elliptic/square_1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_axi.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    wedge-1e2-ang-0.02-surface.vtu square_1e2_axi_pcs_0_ts_1_t_1.000000.vtu temperature temperature 1.6e-5 1e-5
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
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    wedge-1e4-ang-0.02-surface.vtu square_1e4_axi_pcs_0_ts_1_t_1.000000.vtu temperature temperature 1.6e-5 1e-5
)
# # WEDGE 1x1 GROUNDWATER FLOW TEST -- computes reference results for the above test
# AddTest(
#     NAME GroundWaterFlowProcess_wedge_1e4_ang_0.02
#     PATH Elliptic/square_1x1_GroundWaterFlow
#     EXECUTABLE ogs
#     EXECUTABLE_ARGS wedge_1e4_axi_ang_0.02.prj
# )

# MPI groundwater flow tests
AddTest(
    NAME ParallelFEM_GroundWaterFlow2D
    PATH EllipticPETSc
    EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_0.vtu quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_0.vtu pressure pressure 2e-15 1e-16
    quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_1.vtu quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_1.vtu pressure pressure 2e-15 1e-16
    quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_2.vtu quad_20x10_GroundWaterFlow_result_pcs_0_ts_1_t_1_000000_2.vtu pressure pressure 2e-15 1e-16
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_DirichletBC
    PATH EllipticPETSc
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_pcs_0_ts_1_t_1_000000_0.vtu cube_1e3_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 2e-15 1e-16
    cube_1e3_pcs_0_ts_1_t_1_000000_1.vtu cube_1e3_pcs_0_ts_1_t_1_000000_1.vtu Linear_1_to_minus1 pressure 2e-15 1e-16
    cube_1e3_pcs_0_ts_1_t_1_000000_2.vtu cube_1e3_pcs_0_ts_1_t_1_000000_2.vtu Linear_1_to_minus1 pressure 2e-15 1e-16
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow3D_NeumannBC
    PATH EllipticPETSc
    EXECUTABLE_ARGS cube_1e3_neumann.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    cube_1e3_neumann_pcs_0_ts_1_t_1_000000_0.vtu cube_1e3_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    cube_1e3_neumann_pcs_0_ts_1_t_1_000000_1.vtu cube_1e3_neumann_pcs_0_ts_1_t_1_000000_1.vtu D1_left_front_N1_right pressure 1e-2 1e-2
    cube_1e3_neumann_pcs_0_ts_1_t_1_000000_2.vtu cube_1e3_neumann_pcs_0_ts_1_t_1_000000_2.vtu D1_left_front_N1_right pressure 1e-2 1e-2
)

AddTest(
    NAME ParallelFEM_GroundWaterFlow2D_NeumannBC
    PATH EllipticPETSc
    EXECUTABLE_ARGS square_1e1_neumann.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    square_1e1_neumann_pcs_0_ts_1_t_1_000000_0.vtu square_1e1_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure 1e-2 0
    square_1e1_neumann_pcs_0_ts_1_t_1_000000_1.vtu square_1e1_neumann_pcs_0_ts_1_t_1_000000_1.vtu D1_left_bottom_N1_right pressure 1e-2 0
    square_1e1_neumann_pcs_0_ts_1_t_1_000000_0.vtu square_1e1_neumann_pcs_0_ts_1_t_1_000000_0.vtu pressure pressure 1e-14 0
    square_1e1_neumann_pcs_0_ts_1_t_1_000000_1.vtu square_1e1_neumann_pcs_0_ts_1_t_1_000000_1.vtu pressure pressure 1e-14 0
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
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure 1e-1 1e-1
    )
endforeach()

# TODO: Parallel LARGE tests not tested!
foreach(mesh_size 1e4 2e4 3e4 4e4 5e4 1e5 1e6)
    AddTest(
        NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE_ARGS cube_${mesh_size}.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-7 1e-7
    )

    AddTest(
        NAME LARGE_GroundWaterFlowProcess_cube_1x1x1_Neumann_${mesh_size}
        PATH Elliptic/cube_1x1x1_GroundWaterFlow
        EXECUTABLE_ARGS cube_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        cube_1x1x1_hex_${mesh_size}.vtu cube_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_front_N1_right pressure 1e-2 1e-2
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
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-13 1e-13
    )

    AddTest(
        NAME GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure 1e-1 1e-1
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
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-7 1e-7
    )

    AddTest(
        NAME LARGE_GroundWaterFlowProcess_square_1x1_Neumann_${mesh_size}
        PATH Elliptic/square_1x1_GroundWaterFlow
        EXECUTABLE_ARGS square_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        square_1x1_quad_${mesh_size}.vtu square_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_bottom_N1_right pressure 1e-02 1e-02
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
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_pcs_0_ts_1_t_1_000000_0.vtu Linear_1_to_minus1 pressure 1e-15 1e-15
    )

    AddTest(
        NAME GroundWaterFlowProcess_line_1_Neumann_${mesh_size}
        PATH Elliptic/line_1_GroundWaterFlow
        EXECUTABLE_ARGS line_${mesh_size}_neumann.prj
        WRAPPER mpirun
        WRAPPER_ARGS -np 1
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        line_1_line_${mesh_size}.vtu line_${mesh_size}_neumann_pcs_0_ts_1_t_1_000000_0.vtu D1_left_N1_right pressure 1e-14 1e-14
    )
endforeach()

AddTest(
    NAME GroundWaterFlowProcess_Inhomogeneous_permeability
    PATH Elliptic/nonuniform_bc_Groundwaterflow
    EXECUTABLE ogs
    EXECUTABLE_ARGS inhomogeneous_permeability.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    inhomogeneous_permeability.vtu inhomogeneous_permeability_pcs_0_ts_1_t_1.000000.vtu mass_flux_ref mass_flux 4e-2 1e-16
)

AddTest(
    NAME GroundWaterFlowProcess_Neumann_nonuniform_cosY
    PATH Elliptic/nonuniform_bc_Groundwaterflow
    EXECUTABLE ogs
    EXECUTABLE_ARGS neumann_nonuniform.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_neumann_nonuniform_pcs_0_ts_1_t_1.000000.vtu neumann_nonuniform_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-14 0
    expected_neumann_nonuniform_pcs_0_ts_1_t_1.000000.vtu neumann_nonuniform_pcs_0_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-12 0
)

AddTest(
    NAME GroundWaterFlowProcess_Dirichlet_nonuniform_linearY
    PATH Elliptic/nonuniform_bc_Groundwaterflow
    EXECUTABLE ogs
    EXECUTABLE_ARGS dirichlet_nonuniform.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_dirichlet_nonuniform_pcs_0_ts_1_t_1.000000.vtu dirichlet_nonuniform_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-14 0
)

# tests for nodal source term implementation
# For the special setup with a dirac source term at position (xi, eta) the
# analytical solution in 2 dimensions is valid:
# u(x,y) = ln(sqrt((x-xi)^2+(y-eta)^2))/(2 * Pi)
AddTest(
    NAME GroundWaterFlowProcess_NodalSourceTerm_circle_1e1
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e1_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e1_expected.vtu circle_1e1_axi_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 0.7 1e-16
    VIS circle_1e1_axi_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME GroundWaterFlowProcess_NodalSourceTerm_circle_1e2
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e2_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e2_expected.vtu circle_1e2_axi_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 1.1 1e-16
    VIS circle_1e2_axi_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME GroundWaterFlowProcess_NodalSourceTerm_circle_1e3
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e3_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e3_expected.vtu circle_1e3_axi_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 1.6 1e-16
    VIS circle_1e3_axi_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME GroundWaterFlowProcess_NodalSourceTerm_circle_1e4
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e4_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e4_expected.vtu circle_1e4_axi_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 1.8 1e-16
    VIS circle_1e4_axi_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME GroundWaterFlowProcess_NodalSourceTerm_circle_1e5
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e5_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e5_expected.vtu circle_1e5_axi_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 2.15 1e-16
    VIS circle_1e5_axi_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME LARGE_GroundWaterFlowProcess_NodalSourceTerm_circle_1e6
    PATH Elliptic/circle_radius_1
    EXECUTABLE ogs
    EXECUTABLE_ARGS circle_1e6_axi.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    line_1_lines_1e6_expected.vtu circle_1e6_axi_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 2.52 1e-16
    VIS circle_1e6_axi_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME LARGE_GroundWaterFlowProcess_NodalSourceTerm_square_1e6
    PATH Elliptic/square_1x1_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e6_with_nodal_sources.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_1x1_quad_1e6_nodal_sources_expected.vtu square_1e6_with_nodal_sources_pcs_0_ts_1_t_1.000000.vtu analytical_solution pressure 1.41 1e-16
    VIS square_1e6__nodal_sources_expected_pcs_0_ts_1_t_1.000000.vtu
)

