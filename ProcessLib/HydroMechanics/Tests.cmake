# HydroMechanics; Small deformations, linear poroelastic (HML)

### With monolithic scheme
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad8_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu pressure pressure 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu pressure pressure 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu pressure pressure 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_1_t_5.000000.vtu square_1e2_pcs_0_ts_1_t_5.000000.vtu velocity velocity 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_20_t_100.000000.vtu square_1e2_pcs_0_ts_20_t_100.000000.vtu velocity velocity 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_120_t_1000.000000.vtu square_1e2_pcs_0_ts_120_t_1000.000000.vtu velocity velocity 1e-15 1e-15
    expected_square_1e2_pcs_0_ts_420_t_4000.000000.vtu square_1e2_pcs_0_ts_420_t_4000.000000.vtu velocity velocity 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad9_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_quad9.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_HML_square_1e2_tri6_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_tri.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_HML_cube_1e3_hex20_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e3_pcs_0_ts_1_t_5.000000.vtu cube_1e3_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_cube_1e3_pcs_0_ts_1_t_5.000000.vtu cube_1e3_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
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
    DIFF_DATA
    expected_square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu displacement displacement 1e-11 1e-16
    expected_square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu pressure pressure 1e-11 1e-16
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
    DIFF_DATA
    expected_square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu displacement displacement 1e-13 1e-16
    expected_square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu pressure pressure 1e-13 1e-16
)

## Test as the reference of InjectionProduction1D
AddTest(
    NAME MonolithicInjectionProduction1D
    PATH HydroMechanics/StaggeredScheme/InjectionProduction1D/RerefenceSolutionByMonolithicScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS InjectionProduction1DMono.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu displacement displacement 1e-11 1e-11
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu pressure pressure 1e-11 1e-11
)

### With staggered scheme
AddTest(
    NAME StaggeredInjectionProduction1D
    PATH HydroMechanics/StaggeredScheme/InjectionProduction1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS InjectionProduction1D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu displacement displacement 1e-11 1e-11
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu pressure pressure 1e-11 1e-11
)
