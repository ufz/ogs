AddTest(
    NAME RichardsComponentTransport_1D_Padilla_NaCl1
    PATH Parabolic/RichardsComponentTransport/Padilla/Padilla_NaCl1
    EXECUTABLE ogs
    EXECUTABLE_ARGS Padilla_NaCl1.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB Padilla_NaCl1_pcs_0_ts_*.vtu concentration concentration 1e-12 0
    GLOB Padilla_NaCl1_pcs_0_ts_*.vtu pressure pressure 1e-10 0
)

AddTest(
    NAME RichardsComponentTransport_1D_Padilla_NaCl1_quadratic
    PATH Parabolic/RichardsComponentTransport/Padilla/Padilla_NaCl1
    EXECUTABLE ogs
    EXECUTABLE_ARGS Padilla_NaCl1_quadratic.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB Padilla_NaCl1_quadratic_pcs_0_ts_*.vtu concentration concentration 5e-8 0
    GLOB Padilla_NaCl1_quadratic_pcs_0_ts_*.vtu pressure pressure 5e-7 0
)

AddTest(
    NAME RichardsComponentTransport_1D_Padilla_NaCl6
    PATH Parabolic/RichardsComponentTransport/Padilla/Padilla_NaCl6
    EXECUTABLE ogs
    EXECUTABLE_ARGS Padilla_NaCl6.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB Padilla_NaCl6_pcs_0_ts_*.vtu concentration concentration 1e-13 0
    GLOB Padilla_NaCl6_pcs_0_ts_*.vtu pressure pressure 1e-10 0
)

