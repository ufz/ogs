if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1D1bt/m2_1D1bt.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1D2bt/m2_1D2bt.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Dcreep/m2_1Dcreep.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Dlozenge/m2_1Dlozenge.prj RUNTIME 23)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Dlozengebt/m2_1Dlozengebt.prj RUNTIME 82)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_1Drelax/m2_1Drelax.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_2Dload/m2_2Dload.prj)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_2Dload/m2_2Dload_ym45.prj RUNTIME 43)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_2Dloadbt/m2_2Dloadbt.prj RUNTIME 64)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_3Dload/m2_3Dload.prj RUNTIME 24)
    OgsTest(PROJECTFILE ThermoMechanics/CreepBGRa/Verification/m2_3Dloadbt/m2_3Dloadbt.prj RUNTIME 67)
endif()

AddTest(
    NAME ThermoMechanics_3D_ThermoElastic_Stress_Analysis
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    stress_analytical.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma sigma 1e-5 1e-12
    expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu displacement displacement 1e-10 1e-12
    expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu temperature temperature 1e-10 1e-12
    expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma sigma 1e-6 1e-12
    expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu epsilon epsilon 1e-16 0
)

AddTest(
    NAME ThermoMechanics_2D_ThermoElastic_IGLU_Plane_Strain
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_quarter_plane_strain.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 39
    DIFF_DATA
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 1e-9 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 2e-6 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu sigma sigma 5e-6 1e-15
    expected_tm_q_pcs_0_ts_20_t_20000.000000.vtu tm_q_pcs_0_ts_20_t_20000.000000.vtu epsilon epsilon 5e-6 1e-15
)

AddTest(
    NAME ThermoMechanics_2D_ThermoElastic_IGLU_Axisymmetric_Plane_Strain
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_axisymmetric_plane_strain.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 1e-9 1e-15
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 1e-10 1e-8
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu sigma sigma 1e-6 0
    expected_tm_a_pcs_0_ts_20_t_20000.000000.vtu tm_a_pcs_0_ts_20_t_20000.000000.vtu epsilon epsilon 1e-10 0
)

AddTest(
    NAME ThermoMechanics_LARGE_2D_ThermoElastic_IGLU_Plane_Strain_Quadratic_Mesh
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_quarter_plane_strain_quad.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 5e-10 1e-15
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 2e-6 1e-15
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu sigma sigma 5e-6 0
    expected_tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu tm_q_quad_pcs_0_ts_20_t_20000.000000.vtu epsilon epsilon 6e-6 0
)

AddTest(
    NAME ThermoMechanics_2D_ThermoElastic_IGLU_Axisymmetric_Plane_Strain_Quadratic_Mesh
    PATH ThermoMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS iglu_axisymmetric_plane_strain_quad.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu displacement displacement 5e-10 1e-15
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu temperature temperature 5e-10 1e-7
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu sigma sigma 1e-6 0
    expected_tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu tm_a_quad_pcs_0_ts_20_t_20000.000000.vtu epsilon epsilon 5e-10 1e-8
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
    expected_SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu displacement displacement 1e-14 1e-10
    expected_SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu sigma sigma 1e-7 0
    expected_SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu SimpleAxisymmetricCreep_pcs_0_ts_370_t_360.000000.vtu epsilon epsilon 1e-12 0
)

AddTest(
    NAME ThermoMechanics_CreepBGRa_SimpleAxisymmetricCreepWithAnalyticSolution
    PATH ThermoMechanics/CreepBGRa/SimpleAxisymmetricCreep
    EXECUTABLE ogs
    EXECUTABLE_ARGS SimpleAxisymmetricCreepWithAnalyticSolution.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 22
    DIFF_DATA
    SimpleAxisymmetricCreepWithAnalyticSolution.vtu SimpleAxisymmetricCreepWithAnalyticalSolution_pcs_0_ts_1000_t_100.000000.vtu analytic_strain epsilon 1e-7 0
)

AddTest(
    NAME ThermoMechanics_CreepAfterExcavation
    PATH ThermoMechanics/CreepBGRa/CreepAfterExcavation
    EXECUTABLE ogs
    EXECUTABLE_ARGS CreepAfterExcavation.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 24
    DIFF_DATA
    ExpectedCreepAfterExcavation_pcs_0_ts_61_t_4320000.000000.vtu CreepAfterExcavation_pcs_0_ts_61_t_4320000.000000.vtu sigma sigma 5e-6 0
    ExpectedCreepAfterExcavation_pcs_0_ts_61_t_4320000.000000.vtu CreepAfterExcavation_pcs_0_ts_61_t_4320000.000000.vtu epsilon epsilon 1e-15 0
    ExpectedCreepAfterExcavation_pcs_0_ts_61_t_4320000.000000.vtu CreepAfterExcavation_pcs_0_ts_61_t_4320000.000000.vtu displacement displacement 1e-16 1e-9
)

# Basic test that MFront models work for TM.
# Linear elastic, no internal state variables, but external temperature.
AddTest(
    NAME ThermoMechanics_confined_thermal_expansion_mfront
    PATH ThermoMechanics/LinearMFront
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_lin.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MFRONT AND NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_1_t_1.000000.vtu   sigma_1  sigma 1e-16 0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_2_t_2.000000.vtu   sigma_2  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_3_t_3.000000.vtu   sigma_3  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_4_t_4.000000.vtu   sigma_4  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_5_t_5.000000.vtu   sigma_5  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_6_t_6.000000.vtu   sigma_6  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_7_t_7.000000.vtu   sigma_7  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_8_t_8.000000.vtu   sigma_8  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_9_t_9.000000.vtu   sigma_9  sigma 1e-8  0
    cthex_ref.vtu  cube_1e0_lin_pcs_0_ts_10_t_10.000000.vtu sigma_10 sigma 1e-8  0
)

# Test of a creep law.
AddTest(
    NAME ThermoMechanics_BDT_mfront
    PATH ThermoMechanics/BDT
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_bdt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MFRONT AND NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_51_t_300.000000.vtu     epsilon_300   epsilon  1e-8   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_51_t_300.000000.vtu     sigma_300     sigma    1e-3   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_151_t_900.000000.vtu    epsilon_900   epsilon  1e-9   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_151_t_900.000000.vtu    sigma_900     sigma    1e-3   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_251_t_1500.000000.vtu   epsilon_1500  epsilon  1e-8   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_251_t_1500.000000.vtu   sigma_1500    sigma    1e-3   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_501_t_3000.000000.vtu   epsilon_3000  epsilon  1e-10  0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_501_t_3000.000000.vtu   sigma_3000    sigma    1e-3   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_1001_t_6000.000000.vtu  epsilon_6000  epsilon  1e-7   0
    bdt_ref.vtu  cube_1e0_bdt_pcs_0_ts_1001_t_6000.000000.vtu  sigma_6000    sigma    1e-3   0
)
