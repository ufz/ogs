AddTest(
    NAME ThermoMechanicalPhaseField_3D_non-isothermal
    PATH ThermoMechanicalPhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu displacement displacement 1e-6 1e-6
    expected_cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu temperature temperature 1e-6 1e-6
    analytical_cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu phasefield phasefield 1e-5 1e-5
    mechanical_cube_1e0_pcs_1_ts_1_t_100.000000.vtu cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu phasefield phasefield 1e-5 1e-5
   )

AddTest(
    NAME ThermoMechanicalPhaseField_3D_beam_LARGE
    PATH ThermoMechanicalPhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_pcs_2_ts_1_t_10.000000.vtu displacement displacement 1e-5 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_pcs_2_ts_1_t_10.000000.vtu phasefield phasefield 1e-6 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_pcs_2_ts_1_t_10.000000.vtu temperature temperature 1e-6 0
   )
