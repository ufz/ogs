AddTest(
    NAME 3D_ThermoElastic_Stress_Analysis
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma_xx sigma_xx 1e-10 1e-12
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma_yy sigma_yy 1e-10 1e-12
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma_zz sigma_zz 1e-10 1e-12
    expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu displacement displacement 1e-10 1e-12
    expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu temperature temperature 1e-10 1e-12
)

AddTest(
    NAME 2D_ThermoElastic_IGLU_Plane_Strain
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_quarter_plane_strain.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 1e-9 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 2e-6 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu sigma_xx sigma_xx 5e-6 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu sigma_yy sigma_yy 5e-6 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu sigma_zz sigma_zz 5e-6 1e-15
)

AddTest(
    NAME 2D_ThermoElastic_IGLU_Axisymmetric_Plane_Strain
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_axisymmetric_plane_strain.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 1e-10 1e-15
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 1e-10 1e-15
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu sigma_xx sigma_xx 1e-10 1e-15
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu sigma_yy sigma_yy 1e-10 1e-15
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu sigma_zz sigma_zz 1e-10 1e-15
)

AddTest(
    NAME LARGE_2D_ThermoElastic_IGLU_Plane_Strain_Quadratic_Mesh
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_quarter_plane_strain_quad.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 5e-10 1e-15
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 2e-6 1e-15
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu sigma_xx sigma_xx 5e-6 0
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu sigma_yy sigma_yy 6e-6 0
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu sigma_zz sigma_zz 5e-6 0
)

AddTest(
    NAME 2D_ThermoElastic_IGLU_Axisymmetric_Plane_Strain_Quadratic_Mesh
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_axisymmetric_plane_strain_quad.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 5e-10 1e-15
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 5e-10 1e-15
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu sigma_xx sigma_xx 5e-10 1e-15
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu sigma_yy sigma_yy 5e-10 1e-15
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu sigma_zz sigma_zz 5e-10 1e-15
)

AddTest(
    NAME ThermoMechanics_CreepBGRa_SimpleAxisymmetricCreep
    PATH ThermoMechanics/CreepBGRa/SimpleAxisymmetricCreep
    EXECUTABLE ogs
    EXECUTABLE_ARGS SimpleAxisymmetricCreep.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu  SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu displacement displacement 1e-16 1e-10
    expected_SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu  SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu sigma_yy sigma_yy 1e-16 1e-10
)

AddTest(
    NAME ThermoMechanics_CreepBGRa_SimpleAxisymmetricCreepWithAnalyticSolution
    PATH ThermoMechanics/CreepBGRa/SimpleAxisymmetricCreep
    EXECUTABLE ogs
    EXECUTABLE_ARGS SimpleAxisymmetricCreepWithAnalyticSolution.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    SimpleAxisymmetricCreepWithAnalyticSolution.vtu SimpleAxisymmetricCreepWithAnalyticalSolution_pcs_0_ts_1000_t_100.000000.vtu   analytic_strain_yy  epsilon_yy 1e-10 1e-4
)
