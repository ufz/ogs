AddTest(
    NAME PhaseField_3D_Unconfined_Compression
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_cube_1e0_pcs_0_ts_10000_t_1.000000.vtu cube_1e0_pcs_0_ts_10000_t_1.000000.vtu displacement displacement 1e-6 1e-6
    expected_cube_1e0_pcs_0_ts_10000_t_1.000000.vtu cube_1e0_pcs_0_ts_10000_t_1.000000.vtu phasefield phasefield 1e-6 1e-6
   )

AddTest(
    NAME PhaseField_2D_StaticHydraulicFracture
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS 2D_static.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu 2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu displacement displacement 1e-15 1e-15
    expected_2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu 2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu phasefield phasefield 1e-15 1e-15
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
    expected_beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
    expected_beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
   )


AddTest(
    NAME LARGE_PhaseField_Staggered_square_line
    PATH PhaseField
    RUNTIME 70
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_line_h_400.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu displacement displacement 1e-16 0
    square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu phasefield phasefield 5e-15 0
    square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu history_field history_field 1e-16 0
    square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu epsilon epsilon 1e-16 0
    square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu square_line_h_400_pcs_0_ts_2000_t_0.200000.vtu sigma sigma 1e-14 0
   )

AddTest(
    NAME LARGE_PhaseField_Staggered_square_shear
    PATH PhaseField
    RUNTIME 200
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_shear_h_400.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu displacement displacement 1e-16 0
    square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu phasefield phasefield 1e-15 0
    square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu sigma sigma 1e-10 0
    square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu square_shear_h_400_pcs_0_ts_4000_t_400.000000.vtu epsilon epsilon 1e-15 0
   )
