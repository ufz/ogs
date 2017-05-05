# LIE; HydroMechanics
AddTest(
    NAME LIE_HM_single_fracture
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-12 RELTOL 1e-12
    DIFF_DATA
    expected_single_fracture_pcs_0_ts_10_t_100.000000.vtu single_fracture_pcs_0_ts_10_t_100.000000.vtu pressure pressure
    expected_single_fracture_pcs_0_ts_10_t_100.000000.vtu single_fracture_pcs_0_ts_10_t_100.000000.vtu displacement displacement
    expected_single_fracture_pcs_0_ts_10_t_100.000000.vtu single_fracture_pcs_0_ts_10_t_100.000000.vtu displacement_jump1 displacement_jump1
)

AddTest(
    NAME LIE_HM_TaskB
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TaskB.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-12 RELTOL 1e-12
    DIFF_DATA
    expected_TaskB_pcs_0_ts_4_t_18.000000.vtu TaskB_pcs_0_ts_4_t_18.000000.vtu pressure pressure
    expected_TaskB_pcs_0_ts_4_t_18.000000.vtu TaskB_pcs_0_ts_4_t_18.000000.vtu displacement displacement
    expected_TaskB_pcs_0_ts_4_t_18.000000.vtu TaskB_pcs_0_ts_4_t_18.000000.vtu displacement_jump1 displacement_jump1
)
