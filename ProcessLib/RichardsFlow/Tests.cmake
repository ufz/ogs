AddTest(
        NAME LARGE_2D_RichardsFlow_h_us_quad
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d.prj
        TESTER vtkdiff
        ABSTOL 1e-1 RELTOL 1e-1
        DIFF_DATA
        h_us_quad_1000.vtu richards_pcs_0_ts_100_t_100.000000.vtu PRESSURE1 pressure
    REQUIREMENTS NOT OGS_USE_MPI
)
