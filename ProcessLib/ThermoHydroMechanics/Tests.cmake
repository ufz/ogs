# ThermoHydroMechanics; Small deformation, linear poroelastic, homogeneous
AddTest(
    NAME ThermoHydroMechanics_square_1e0
    PATH ThermoHydroMechanics/Linear/Square_Homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_2_t_200.000000.vtu square_1e0_pcs_0_ts_2_t_200.000000.vtu displacement displacement
    expected_square_1e0_pcs_0_ts_4_t_400.000000.vtu square_1e0_pcs_0_ts_4_t_400.000000.vtu displacement displacement
    expected_square_1e0_pcs_0_ts_6_t_600.000000.vtu square_1e0_pcs_0_ts_6_t_600.000000.vtu displacement displacement
    expected_square_1e0_pcs_0_ts_8_t_800.000000.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu displacement displacement
    expected_square_1e0_pcs_0_ts_2_t_200.000000.vtu square_1e0_pcs_0_ts_2_t_200.000000.vtu pressure pressure
    expected_square_1e0_pcs_0_ts_4_t_400.000000.vtu square_1e0_pcs_0_ts_4_t_400.000000.vtu pressure pressure
    expected_square_1e0_pcs_0_ts_6_t_600.000000.vtu square_1e0_pcs_0_ts_6_t_600.000000.vtu pressure pressure
    expected_square_1e0_pcs_0_ts_8_t_800.000000.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu pressure pressure
    expected_square_1e0_pcs_0_ts_2_t_200.000000.vtu square_1e0_pcs_0_ts_2_t_200.000000.vtu temperature temperature
    expected_square_1e0_pcs_0_ts_4_t_400.000000.vtu square_1e0_pcs_0_ts_4_t_400.000000.vtu temperature temperature
    expected_square_1e0_pcs_0_ts_6_t_600.000000.vtu square_1e0_pcs_0_ts_6_t_600.000000.vtu temperature temperature
    expected_square_1e0_pcs_0_ts_8_t_800.000000.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu temperature temperature
    expected_square_1e0_pcs_0_ts_2_t_200.000000.vtu square_1e0_pcs_0_ts_2_t_200.000000.vtu epsilon_xx epsilon_xx
    expected_square_1e0_pcs_0_ts_4_t_400.000000.vtu square_1e0_pcs_0_ts_4_t_400.000000.vtu epsilon_xx epsilon_xx
    expected_square_1e0_pcs_0_ts_6_t_600.000000.vtu square_1e0_pcs_0_ts_6_t_600.000000.vtu epsilon_xx epsilon_xx
    expected_square_1e0_pcs_0_ts_8_t_800.000000.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu epsilon_xx epsilon_xx
    expected_square_1e0_pcs_0_ts_2_t_200.000000.vtu square_1e0_pcs_0_ts_2_t_200.000000.vtu epsilon_yy epsilon_yy
    expected_square_1e0_pcs_0_ts_4_t_400.000000.vtu square_1e0_pcs_0_ts_4_t_400.000000.vtu epsilon_yy epsilon_yy
    expected_square_1e0_pcs_0_ts_6_t_600.000000.vtu square_1e0_pcs_0_ts_6_t_600.000000.vtu epsilon_yy epsilon_yy
    expected_square_1e0_pcs_0_ts_8_t_800.000000.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu epsilon_yy epsilon_yy
)

AddTest(
    NAME ThermoHydroMechanics_square_1e0_analyitcal
    PATH ThermoHydroMechanics/Linear/Square_Homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-3 RELTOL 1e-3
    DIFF_DATA
    analytical.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu pressure_ana pressure
    analytical.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu sigma_xx_ana sigma_xx
    analytical.vtu square_1e0_pcs_0_ts_8_t_800.000000.vtu sigma_yy_ana sigma_yy
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, sealed, heterogeneous
AddTest(
    NAME ThermoHydroMechanics_square_1e2_heterogeneous
    PATH ThermoHydroMechanics/Linear/Beam_sealed_Heterogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_6_t_300.000000.vtu square_1e2_pcs_0_ts_6_t_300.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_12_t_600.000000.vtu square_1e2_pcs_0_ts_12_t_600.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_18_t_900.000000.vtu square_1e2_pcs_0_ts_18_t_900.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_24_t_1200.000000.vtu square_1e2_pcs_0_ts_24_t_1200.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_30_t_1500.000000.vtu square_1e2_pcs_0_ts_30_t_1500.000000.vtu displacement displacement
    expected_square_1e2_pcs_0_ts_6_t_300.000000.vtu square_1e2_pcs_0_ts_6_t_300.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_12_t_600.000000.vtu square_1e2_pcs_0_ts_12_t_600.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_18_t_900.000000.vtu square_1e2_pcs_0_ts_18_t_900.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_24_t_1200.000000.vtu square_1e2_pcs_0_ts_24_t_1200.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_30_t_1500.000000.vtu square_1e2_pcs_0_ts_30_t_1500.000000.vtu pressure pressure
    expected_square_1e2_pcs_0_ts_6_t_300.000000.vtu square_1e2_pcs_0_ts_6_t_300.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_12_t_600.000000.vtu square_1e2_pcs_0_ts_12_t_600.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_18_t_900.000000.vtu square_1e2_pcs_0_ts_18_t_900.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_24_t_1200.000000.vtu square_1e2_pcs_0_ts_24_t_1200.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_30_t_1500.000000.vtu square_1e2_pcs_0_ts_30_t_1500.000000.vtu temperature temperature
    expected_square_1e2_pcs_0_ts_6_t_300.000000.vtu square_1e2_pcs_0_ts_6_t_300.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_12_t_600.000000.vtu square_1e2_pcs_0_ts_12_t_600.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_18_t_900.000000.vtu square_1e2_pcs_0_ts_18_t_900.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_24_t_1200.000000.vtu square_1e2_pcs_0_ts_24_t_1200.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_30_t_1500.000000.vtu square_1e2_pcs_0_ts_30_t_1500.000000.vtu sigma_xx sigma_xx
    expected_square_1e2_pcs_0_ts_6_t_300.000000.vtu square_1e2_pcs_0_ts_6_t_300.000000.vtu sigma_yy sigma_yy
    expected_square_1e2_pcs_0_ts_12_t_600.000000.vtu square_1e2_pcs_0_ts_12_t_600.000000.vtu sigma_yy sigma_yy
    expected_square_1e2_pcs_0_ts_18_t_900.000000.vtu square_1e2_pcs_0_ts_18_t_900.000000.vtu sigma_yy sigma_yy
    expected_square_1e2_pcs_0_ts_24_t_1200.000000.vtu square_1e2_pcs_0_ts_24_t_1200.000000.vtu sigma_yy sigma_yy
    expected_square_1e2_pcs_0_ts_30_t_1500.000000.vtu square_1e2_pcs_0_ts_30_t_1500.000000.vtu sigma_yy sigma_yy
)
