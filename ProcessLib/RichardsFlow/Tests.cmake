if (NOT OGS_USE_MPI)
    # Comparison test for richards mechanics w/o deformations.
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_richardsflow.prj)
    OgsTest(PROJECTFILE Parabolic/Richards/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE Parabolic/Richards/infiltration_1D/1d-infiltration-RM.prj)
endif()

AddTest(
        NAME 2D_RichardsFlow_h_us_quad_ogs5
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_compare_ogs5.prj
        TESTER vtkdiff
        DIFF_DATA
        h_us_quad_1000.vtu richards_pcs_0_ts_100_t_100.000000.vtu PRESSURE1 pressure 1e-1 1e-1
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME 2D_RichardsFlow_h_us_quad_small
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_small.prj
        TESTER vtkdiff
        DIFF_DATA
        ref_t_1600.000000.vtu richards_pcs_0_ts_1100_t_1600.000000.vtu pressure pressure 1e-8 1e-8
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME LARGE_2D_RichardsFlow_h_us_quad
        PATH Parabolic/Richards
        RUNTIME 90
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_large.prj
        TESTER vtkdiff
        DIFF_DATA
        ref_t_20000.000000.vtu richards_pcs_0_ts_18200_t_20000.000000.vtu pressure pressure 1e-8 1e-8
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small_PID_adaptive_dt
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small_PID_adaptive_dt.prj
    TESTER vtkdiff
    DIFF_DATA
    ref_t_1600.000000.vtu richards_pcs_0_ts_803_t_1600.000000.vtu pressure pressure 1e-8 1e-3
# The following three comparisons are used just to check whether the output is
# made at the fixed times of 50, 100 and 500, which are given in the project
# file of RichardsFlow_2d_small_adaptive_dt.prj
    richards_pcs_0_ts_28_spec_t_50.000000.vtu richards_pcs_0_ts_28_t_50.000000.vtu pressure pressure 1e-10 1e-10
    richards_pcs_0_ts_53_spec_t_100.000000.vtu richards_pcs_0_ts_53_t_100.000000.vtu pressure pressure 1e-10 1e-10
    richards_pcs_0_ts_253_spec_t_500.000000.vtu richards_pcs_0_ts_253_t_500.000000.vtu pressure pressure 1e-10 1e-10
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small_iteration_adaptive_dt.prj
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 51
    # No vtkdiff comparison here, because of the different file names for
    # different machines, which again is due to the adaptive time stepping
    # scheme. When the output file format can be specified in the project files,
    # e.g. in the form Richards_%t where %t is the current time, the output will
    # no longer be ambiguous.
)

#PETSc/MPI
#AddTest(
#    NAME 2D_RichardsFlow_h_us_quad_small_PID_adaptive_dt
#    PATH Parabolic/Richards
#    EXECUTABLE_ARGS RichardsFlow_2d_small_PID_adaptive_dt.prj
#    WRAPPER mpirun
#    WRAPPER_ARGS -np 1
#    TESTER vtkdiff
#    REQUIREMENTS OGS_USE_MPI
#    RUNTIME 220
#    DIFF_DATA
#    ref_t_1600.000000.vtu richards_pcs_0_ts_803_t_1600_000000_0.vtu pressure pressure 1e-8 1e-3
#)

# Comparison test for richards mechanics w/o deformations.
AddTest(
    NAME Parallel_RichardsMechanics_RichardsFlow_2d_richardsflow
    PATH RichardsMechanics
    EXECUTABLE_ARGS RichardsFlow_2d_richardsflow.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900.000000.vtu RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900_000000_0.vtu pressure pressure 5e-8 1e-10
    RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900.000000.vtu RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900_000000_0.vtu saturation saturation 1e-10 1e-11
)

AddTest(
    NAME 1d-infiltration
    PATH Parabolic/Richards/infiltration_1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS 1d-infiltration.prj
    RUNTIME 3
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB 1d-infiltration_t_*.vtu pressure pressure 1e-15 0
    GLOB 1d-infiltration_t_*.vtu darcy_velocity darcy_velocity 1e-15 0
    GLOB 1d-infiltration_t_*.vtu saturation saturation 1e-15 0
)
