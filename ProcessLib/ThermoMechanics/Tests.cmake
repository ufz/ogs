AddTest(
    NAME 3D_ThermoElastic_Stress_Analysis
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-10 RELTOL 1e-12
    DIFF_DATA
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu sigma_xx sigma_xx
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu sigma_yy sigma_yy
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu sigma_zz sigma_zz
    expected_cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu displacement displacement
    expected_cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_80_t_72000.000000.vtu temperature temperature
)
