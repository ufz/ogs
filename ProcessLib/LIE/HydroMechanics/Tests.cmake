# LIE; HydroMechanics
AddTest(
    NAME LIE_HM_single_fracture
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_pcs_0_ts_10_t_100.000000.vtu single_fracture_pcs_0_ts_10_t_100.000000.vtu pressure pressure 1e-12 1e-12
    expected_single_fracture_pcs_0_ts_10_t_100.000000.vtu single_fracture_pcs_0_ts_10_t_100.000000.vtu displacement displacement 1e-12 1e-12
    expected_single_fracture_pcs_0_ts_10_t_100.000000.vtu single_fracture_pcs_0_ts_10_t_100.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
)

AddTest(
    NAME LARGE_LIE_HM_single_fracture_3D
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_3D_pcs_0_ts_10_t_100.000000.vtu single_fracture_3D_pcs_0_ts_10_t_100.000000.vtu pressure pressure 1e-12 1e-12
    expected_single_fracture_3D_pcs_0_ts_10_t_100.000000.vtu single_fracture_3D_pcs_0_ts_10_t_100.000000.vtu displacement displacement 1e-12 1e-12
    expected_single_fracture_3D_pcs_0_ts_10_t_100.000000.vtu single_fracture_3D_pcs_0_ts_10_t_100.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
)

AddTest(
    NAME LIE_HM_TaskB
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TaskB.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_TaskB_pcs_0_ts_4_t_18.000000.vtu TaskB_pcs_0_ts_4_t_18.000000.vtu pressure pressure 1e-12 1e-12
    expected_TaskB_pcs_0_ts_4_t_18.000000.vtu TaskB_pcs_0_ts_4_t_18.000000.vtu displacement displacement 1e-12 1e-12
    expected_TaskB_pcs_0_ts_4_t_18.000000.vtu TaskB_pcs_0_ts_4_t_18.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
)
