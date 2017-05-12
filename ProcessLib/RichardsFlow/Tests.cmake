AddTest(
        NAME 2D_RichardsFlow_h_us_quad_ogs5
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_compare_ogs5.prj
        TESTER vtkdiff
        ABSTOL 1e-1 RELTOL 1e-1
        DIFF_DATA
        h_us_quad_1000.vtu richards_pcs_0_ts_100_t_100.000000.vtu PRESSURE1 pressure
    REQUIREMENTS NOT OGS_USE_MPI
)
AddTest(
        NAME 2D_RichardsFlow_h_us_quad_small
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_small.prj
        TESTER vtkdiff
        ABSTOL 1e-8 RELTOL 1e-8
        DIFF_DATA
        ref_t_1600.000000.vtu richards_pcs_0_ts_1100_t_1600.000000.vtu pressure pressure
    REQUIREMENTS NOT OGS_USE_MPI
)
AddTest(
        NAME LARGE_2D_RichardsFlow_h_us_quad
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_large.prj
        TESTER vtkdiff
        ABSTOL 1e-8 RELTOL 1e-8
        DIFF_DATA
        ref_t_20000.000000.vtu richards_pcs_0_ts_18200_t_20000.000000.vtu pressure pressure
    REQUIREMENTS NOT OGS_USE_MPI
)
