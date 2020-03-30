# ThermoHydroMechanics; Small deformation, linear poroelastic, homogeneous
AddTest(
    NAME ThermoHydroMechanics_square_1e0
    PATH ThermoHydroMechanics/Linear/Square_sealed_homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu epsilon epsilon 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu sigma sigma 1e-8 1e-8
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, sealed, bimaterial
AddTest(
    NAME ThermoHydroMechanics_square_1e2_sealed_bimaterial
    PATH ThermoHydroMechanics/Linear/Beam_sealed_bimaterial
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu epsilon epsilon 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu sigma sigma 1e-8 1e-8
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, unsealed, bimaterial
AddTest(
    NAME ThermoHydroMechanics_square_1e2_unsealed_bimaterial
    PATH ThermoHydroMechanics/Linear/Beam_unsealed_bimaterial
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu epsilon epsilon 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu sigma sigma 1e-8 1e-8
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation
AddTest(
    NAME LARGE_ThermoHydroMechanics_square_1e2_point_heat_injection
    PATH ThermoHydroMechanics/Linear/Point_injection
    RUNTIME 400
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_10_t_50000.000000.vtu square_1e0_pcs_0_ts_10_t_50000.000000.vtu displacement displacement 1e-5 1e-5
    expected_square_1e0_pcs_0_ts_10_t_50000.000000.vtu square_1e0_pcs_0_ts_10_t_50000.000000.vtu pressure pressure 1e-5 1e-5
    expected_square_1e0_pcs_0_ts_10_t_50000.000000.vtu square_1e0_pcs_0_ts_10_t_50000.000000.vtu temperature temperature 1e-5 1e-5
    expected_square_1e0_pcs_0_ts_10_t_50000.000000.vtu square_1e0_pcs_0_ts_10_t_50000.000000.vtu epsilon epsilon 1e-5 1e-5
    expected_square_1e0_pcs_0_ts_10_t_50000.000000.vtu square_1e0_pcs_0_ts_10_t_50000.000000.vtu sigma sigma 1e-5 1e-5
)
# ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation, linear elements for displacement
AddTest(
    NAME LARGE_ThermoHydroMechanics_square_1e2_point_heat_injection_lin
    PATH ThermoHydroMechanics/Linear/Point_injection
    RUNTIME 400
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_lin.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu displacement displacement 1e-5 1e-5
    expected_square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu pressure pressure 1e-5 1e-5
    expected_square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu temperature temperature 1e-5 1e-5
    expected_square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu epsilon epsilon 1e-5 1e-5
    expected_square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu square_1e0_lin_pcs_0_ts_10_t_50000.000000.vtu sigma sigma 1e-5 1e-5
)
