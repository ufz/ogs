# HydroMechanics; Small deformations, linear poroelastic (HML)
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad8_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-15
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu velocity_x velocity_x
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu velocity_x velocity_x
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu velocity_x velocity_x
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu velocity_x velocity_x
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu velocity_y velocity_y
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu velocity_y velocity_y
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu velocity_y velocity_y
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu velocity_y velocity_y
)
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad9_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_quad9.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-15
    DIFF_DATA
    expected_square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu displacement displacement
    expected_square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu pressure pressure
)
AddTest(
    NAME HydroMechanics_HML_square_1e2_tri6_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_tri.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-15
    DIFF_DATA
    expected_square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu displacement displacement
    expected_square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu pressure pressure
)
AddTest(
    NAME HydroMechanics_HML_cube_1e3_hex20_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-15
    DIFF_DATA
    expected_cube_1e3_pcs_0_ts_1_t_5.000000.vtu cube_1e3_pcs_0_ts_1_t_5.000000.vtu displacement displacement
    expected_cube_1e3_pcs_0_ts_1_t_5.000000.vtu cube_1e3_pcs_0_ts_1_t_5.000000.vtu pressure pressure
)

# HydroMechanics; Small deformation, linear poroelastic (unconfined compression early) The drainage process is ongoing and the displacement behaviour is related to water pressure and solid properties.
AddTest(
    NAME HydroMechanics_HML_square_1e2_unconfined_compression_early
    PATH HydroMechanics/Linear/Unconfined_Compression_early
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_UC_early.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-3 RELTOL 1e-3
    DIFF_DATA
    UnconfinedCompressionAnalytical_1s.vtu square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu displacement_ana displacement
)
# HydroMechanics; Small deformation, linear poroelastic (unconfined compression late) the drainage process is finished and the displacement of the porous media is only a result of solid properties.
AddTest(
    NAME HydroMechanics_HML_square_1e2_unconfined_compression_late
    PATH HydroMechanics/Linear/Unconfined_Compression_late
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_UC_late.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-3 RELTOL 1e-3
    DIFF_DATA
    UnconfinedCompressionAnalytical_1000s.vtu square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu displacement_ana displacement
)
