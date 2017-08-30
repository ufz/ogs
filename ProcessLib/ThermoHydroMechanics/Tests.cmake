# ThermoHydroMechanics; Small deformation, linear poroelastic, homogeneous
AddTest(
    NAME ThermoHydroMechanics_square_1e0
    PATH ThermoHydroMechanics/Linear/Square_sealed_homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu displacement displacement
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu pressure pressure
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu temperature temperature
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu epsilon_xx epsilon_xx
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu epsilon_yy epsilon_yy
)

AddTest(
    NAME ThermoHydroMechanics_square_1e0_analyitcal
    PATH ThermoHydroMechanics/Linear/Square_sealed_homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-3 RELTOL 1e-3
    DIFF_DATA
    analytical.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu pressure_ana pressure
    analytical.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu sigma_xx_ana sigma_xx
    analytical.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu sigma_yy_ana sigma_yy
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, sealed, bimaterial
AddTest(
    NAME ThermoHydroMechanics_square_1e2_sealed_bimaterial
    PATH ThermoHydroMechanics/Linear/Beam_sealed_bimaterial
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu sigma_yy sigma_yy
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, unsealed, bimaterial
AddTest(
    NAME ThermoHydroMechanics_square_1e2_unsealed_bimaterial
    PATH ThermoHydroMechanics/Linear/Beam_unsealed_bimaterial
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu sigma_yy sigma_yy
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation
AddTest(
    NAME ThermoHydroMechanics_square_1e2_point_heat_injection
    PATH ThermoHydroMechanics/Linear/Point_injection
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-6 RELTOL 1e-6
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu displacement displacement
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu pressure pressure
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu temperature temperature
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu sigma_xx sigma_xx
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu sigma_yy sigma_yy
)
