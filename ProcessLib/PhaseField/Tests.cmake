AddTest(
    NAME PhaseField_2D_StaticHydraulicFracture
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS 2D_static.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_2D_StaticCrack_ts_1_t_1.000000.vtu 2D_StaticCrack_ts_1_t_1.000000.vtu displacement displacement 1e-15 0
    expected_2D_StaticCrack_ts_1_t_1.000000.vtu 2D_StaticCrack_ts_1_t_1.000000.vtu phasefield phasefield 5e-15 0
   )

AddTest(
    NAME PhaseField_3D_beam
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d_stag_1pcs.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 18
    DIFF_DATA
    expected_beam3_stag1pcsAT2_ts_10_t_1.000000.vtu beam3_stag1pcsAT2_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
    expected_beam3_stag1pcsAT2_ts_10_t_1.000000.vtu beam3_stag1pcsAT2_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
   )
