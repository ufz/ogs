if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/LinearMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/FullySaturatedFlowMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/RichardsFlow2D/RichardsFlow_2d_small.prj RUNTIME 10)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/A2/A2.prj RUNTIME 18)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/OrthotropicSwelling/orthotropic_swelling_xy.xml)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/OrthotropicSwelling/orthotropic_swelling_xyz.xml)
endif()

AddTest(
    NAME ThermoRichardsMechanics_anisotropic_thermal_expansion_vector
    PATH ThermoRichardsMechanics/anisotropic_thermal_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS aniso_expansion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu sigma sigma 4e-11 0
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu displacement displacement 1e-14 1e-14
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu saturation saturation 1e-12 1e-12
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu pressure pressure 1e-10 1e-10
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu epsilon epsilon 1e-14 1e-14
)

AddTest(
    NAME ThermoRichardsMechanics_anisotropic_thermal_expansion_expansivity_matrix
    PATH ThermoRichardsMechanics/anisotropic_thermal_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS aniso_expansion_expansivity_matrix.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_ts_1_t_1000000.000000.vtu sigma sigma 4e-11 0
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_ts_1_t_1000000.000000.vtu displacement displacement 1e-14 1e-14
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_ts_1_t_1000000.000000.vtu saturation saturation 1e-14 1e-14
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_ts_1_t_1000000.000000.vtu pressure pressure 1e-10 1e-10
    expected_anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_ts_1_t_1000000.000000.vtu epsilon epsilon 1e-12 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_anisotropic_thermal_expansion_expansivity_matrix_z90
    PATH ThermoRichardsMechanics/anisotropic_thermal_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS aniso_expansion_expansivity_matrix_z90.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_anisotropic_thermal_expansion_z90_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_z90_ts_1_t_1000000.000000.vtu sigma sigma 4e-11 0
    expected_anisotropic_thermal_expansion_z90_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_z90_ts_1_t_1000000.000000.vtu displacement displacement 1e-14 1e-14
    expected_anisotropic_thermal_expansion_z90_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_z90_ts_1_t_1000000.000000.vtu saturation saturation 1e-12 1e-12
    expected_anisotropic_thermal_expansion_z90_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_z90_ts_1_t_1000000.000000.vtu pressure pressure 1e-10 1e-10
    expected_anisotropic_thermal_expansion_z90_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_expansivity_matrix_z90_ts_1_t_1000000.000000.vtu epsilon epsilon 1e-12 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_anisotropic_thermal_expansion_x45
    PATH ThermoRichardsMechanics/anisotropic_thermal_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS aniso_expansion_x45.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_anisotropic_thermal_expansion_x45_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_x45_ts_1_t_1000000.000000.vtu epsilon epsilon 1e-12 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_anisotropic_thermal_expansion_y45
    PATH ThermoRichardsMechanics/anisotropic_thermal_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS aniso_expansion_y45.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_anisotropic_thermal_expansion_y45_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_y45_ts_1_t_1000000.000000.vtu epsilon epsilon 1e-12 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_anisotropic_thermal_expansion_z45
    PATH ThermoRichardsMechanics/anisotropic_thermal_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS aniso_expansion_z45.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_anisotropic_thermal_expansion_z45_ts_1_t_1000000.000000.vtu anisotropic_thermal_expansion_z45_ts_1_t_1000000.000000.vtu epsilon epsilon 1e-12 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_liakopoulosHM
    PATH ThermoRichardsMechanics/LiakopoulosHM
    EXECUTABLE ogs
    EXECUTABLE_ARGS liakopoulos.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    GLOB liakopoulos_t_*.vtu sigma sigma 1e-9 1e-12
    GLOB liakopoulos_t_*.vtu displacement displacement 1e-10 1e-12
    GLOB liakopoulos_t_*.vtu saturation saturation 1e-10 1e-12
    GLOB liakopoulos_left_t_*.vtu sigma sigma 1e-9 1e-12
    GLOB liakopoulos_left_t_*.vtu displacement displacement 1e-10 1e-12
    GLOB liakopoulos_left_t_*.vtu saturation saturation 1e-10 1e-12
    GLOB liakopoulos_right_t_*.vtu sigma sigma 1e-9 1e-12
    GLOB liakopoulos_right_t_*.vtu displacement displacement 1e-10 1e-12
    GLOB liakopoulos_right_t_*.vtu saturation saturation 1e-10 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_liakopoulosHM_restart
    PATH ThermoRichardsMechanics/LiakopoulosHM
    EXECUTABLE ogs
    EXECUTABLE_ARGS liakopoulos_restart.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    liakopoulos_t_600.vtu liakopoulos_restart_liakopoulos_t_600_t_600.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_600.vtu liakopoulos_restart_liakopoulos_t_600_t_600.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_600.vtu liakopoulos_restart_liakopoulos_t_600_t_600.vtu saturation saturation 1e-10 1e-12
    liakopoulos_t_7200.vtu liakopoulos_restart_liakopoulos_t_600_t_7200.vtu sigma sigma 1.5 0
    liakopoulos_t_7200.vtu liakopoulos_restart_liakopoulos_t_600_t_7200.vtu displacement displacement 3.4e-7 0
    liakopoulos_t_7200.vtu liakopoulos_restart_liakopoulos_t_600_t_7200.vtu saturation saturation 4e-5 0
)

AddTest(
    NAME ThermoRichardsMechanics_LiakopoulosMixedElementsPETSc
    PATH ThermoRichardsMechanics/LiakopoulosPETSc
    EXECUTABLE ogs
    EXECUTABLE_ARGS liakopoulos_mixElem_mumps.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    # Run on envinf only because PETSc MUMPS solver is used
    REQUIREMENTS OGS_USE_PETSC
    LABELS "petsc-mumps"
    RUNTIME 2
    DIFF_DATA
    GLOB liakopoulosBulk_mixElem_t_*.vtu sigma sigma 1e-9 1e-12
    GLOB liakopoulosBulk_mixElem_t_*.vtu displacement displacement 1e-10 1e-12
    GLOB liakopoulosBulk_mixElem_t_*.vtu saturation saturation 1e-10 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_3D_ThermoElastic_Stress_Analysis
    PATH ThermoRichardsMechanics/Simple3DThermoMechanicsFromTM
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    stress_analytical.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu sigma sigma 1e-5 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu displacement displacement 1e-10 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu temperature temperature 1e-10 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu sigma sigma 1e-6 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu epsilon epsilon 1e-16 0
)

AddTest(
    NAME ThermoRichardsMechanics_HeatTransportInStationaryFlow
    PATH ThermoRichardsMechanics/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS HeatTransportInStationaryFlow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature 1e-6 1e-10
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure 1e-10 1e-10
)

AddTest(
    NAME ThermoRichardsMechanics_point_heat_injection
    PATH ThermoRichardsMechanics/PointHeatSource
    RUNTIME 25
    EXECUTABLE ogs
    EXECUTABLE_ARGS point_heat_source_2D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu displacement displacement 1e-6 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu pressure pressure 2.9e4 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu temperature temperature 0.2 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu epsilon epsilon 5e-6 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu sigma sigma 1.4e4 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu HeatFlowRate HeatFlowRate 4e-11 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu MassFlowRate MassFlowRate 2e-8 0
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu PointHeatSource_ts_10_t_50000.000000.vtu NodalForces NodalForces 3.6e2 0
    # submesh residuum output quarter_002_2nd_r_gt_2
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu displacement displacement 1e-15 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu pressure pressure 8e-8 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu temperature temperature 1.2e-13 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu epsilon epsilon 1e-15 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu sigma sigma 7e-8 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu HeatFlowRate HeatFlowRate 4e-11 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu MassFlowRate MassFlowRate 1e-15 0
    PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_gt_2_ts_10_t_50000.000000.vtu NodalForces NodalForces 6e-8 0
    # submesh residuum output quarter_002_2nd_r_lt_2
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu displacement displacement 1e-15 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu pressure pressure 5e-7 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu temperature temperature 3e-11 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu epsilon epsilon 1e-15 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu sigma sigma 1.1e-5 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu HeatFlowRate HeatFlowRate 6.2e-12 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu MassFlowRate MassFlowRate 1e-15 0
    PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu PointHeatSource_quarter_002_2nd_r_lt_2_ts_10_t_50000.000000.vtu NodalForces NodalForces 2.1e-8 0
)

AddTest(
    NAME ThermoRichardsMechanics_TaskCDECOVALEX2023
    PATH ThermoRichardsMechanics/TaskCDECOVALEX2023
    EXECUTABLE ogs
    EXECUTABLE_ARGS Decovalex-0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-0_ts_10_t_864000.000000.vtu sigma sigma 1e-9 1e-8
    Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-0_ts_10_t_864000.000000.vtu displacement displacement 1e-10 1e-12
    Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-0_ts_10_t_864000.000000.vtu saturation saturation 1e-10 1e-12
    Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-0_ts_10_t_864000.000000.vtu temperature temperature 1e-10 1e-12
    Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-0_ts_10_t_864000.000000.vtu velocity velocity 1e-10 1e-12
    Decovalex-0_ts_10_t_864000.000000.vtu Decovalex-0_ts_10_t_864000.000000.vtu liquid_density liquid_density 1e-10 1e-11
)

AddTest(
    NAME ThermoRichardsMechanics_CTF1
    PATH ThermoRichardsMechanics/CTF1
    EXECUTABLE ogs
    EXECUTABLE_ARGS CTF1.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    CTF1_14.000000.vtu CTF1_14.000000.vtu sigma sigma 1e-9 1e-8
    CTF1_14.000000.vtu CTF1_14.000000.vtu displacement displacement 1e-10 1e-12
    CTF1_14.000000.vtu CTF1_14.000000.vtu saturation saturation 1e-10 1e-12
    CTF1_14.000000.vtu CTF1_14.000000.vtu temperature temperature 1e-10 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_CreepBGRa_SimpleAxisymmetricCreep
    PATH ThermoRichardsMechanics/SimpleAxisymmetricCreep
    EXECUTABLE ogs
    EXECUTABLE_ARGS SimpleAxisymmetricCreep.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 6
    DIFF_DATA
    expected_SimpleAxisymmetricCreep_ts_370_t_360.000000.vtu SimpleAxisymmetricCreep_ts_370_t_360.000000.vtu displacement displacement 1e-13 1e-10
    expected_SimpleAxisymmetricCreep_ts_370_t_360.000000.vtu SimpleAxisymmetricCreep_ts_370_t_360.000000.vtu sigma sigma 2e-7 0
    expected_SimpleAxisymmetricCreep_ts_370_t_360.000000.vtu SimpleAxisymmetricCreep_ts_370_t_360.000000.vtu epsilon epsilon 1e-11 0
)

#PETSc
AddTest(
    NAME ParallelFEM_ThermoRichardsMechanics_3D_ThermoElastic_Stress_Analysis
    PATH ThermoRichardsMechanics/Simple3DThermoMechanicsFromTM
    RUNTIME 40
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    LABELS "petsc-mumps"
    DIFF_DATA
    cube_1e3_tm_ts_17_t_72000_000000_0.vtu cube_1e3_tm_ts_17_t_72000_000000_0.vtu velocity velocity 1e-10 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_1.vtu cube_1e3_tm_ts_17_t_72000_000000_1.vtu velocity velocity 1e-10 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_0.vtu cube_1e3_tm_ts_17_t_72000_000000_0.vtu displacement displacement 1e-10 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_1.vtu cube_1e3_tm_ts_17_t_72000_000000_1.vtu displacement displacement 1e-10 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_0.vtu cube_1e3_tm_ts_17_t_72000_000000_0.vtu sigma sigma 1e-7 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_1.vtu cube_1e3_tm_ts_17_t_72000_000000_1.vtu sigma sigma 1e-7 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_0.vtu cube_1e3_tm_ts_17_t_72000_000000_0.vtu epsilon epsilon 1e-10 1e-9
    cube_1e3_tm_ts_17_t_72000_000000_1.vtu cube_1e3_tm_ts_17_t_72000_000000_1.vtu epsilon epsilon 1e-10 1e-9
)

AddTest(
    NAME ParallelFEM_ThermoRichardsMechanics_FullySaturatedFlowMechanics
    PATH ThermoRichardsMechanics/FullySaturatedFlowMechanics/PETSc
    RUNTIME 10
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_fully_saturated_petsc.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    flow_fully_saturated_ts_2_t_2_000000_0.vtu flow_fully_saturated_ts_2_t_2_000000_0.vtu pressure pressure 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_0.vtu flow_fully_saturated_ts_2_t_2_000000_0.vtu velocity velocity 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_1.vtu flow_fully_saturated_ts_2_t_2_000000_1.vtu velocity velocity 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_0.vtu flow_fully_saturated_ts_2_t_2_000000_0.vtu displacement displacement 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_1.vtu flow_fully_saturated_ts_2_t_2_000000_1.vtu displacement displacement 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_0.vtu flow_fully_saturated_ts_2_t_2_000000_0.vtu sigma sigma 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_1.vtu flow_fully_saturated_ts_2_t_2_000000_1.vtu sigma sigma 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_0.vtu flow_fully_saturated_ts_2_t_2_000000_0.vtu epsilon epsilon 1e-10 1e-9
    flow_fully_saturated_ts_2_t_2_000000_1.vtu flow_fully_saturated_ts_2_t_2_000000_1.vtu epsilon epsilon 1e-10 1e-9
)

AddTest(
    NAME ParallelFEM_ThermoRichardsMechanics_point_heat_injection
    PATH ThermoRichardsMechanics/PointHeatSource
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS -p point_heat_source_2D_non_submesh_r_output.xml point_heat_source_2D.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    LABELS "petsc-mumps"
    DIFF_DATA
    PointHeatSource_ts_10_t_50000_000000_0.vtu PointHeatSource_ts_10_t_50000_000000_0.vtu displacement displacement 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_0.vtu PointHeatSource_ts_10_t_50000_000000_0.vtu pressure pressure 1e-10 1.0e-6
    PointHeatSource_ts_10_t_50000_000000_0.vtu PointHeatSource_ts_10_t_50000_000000_0.vtu temperature temperature 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_0.vtu PointHeatSource_ts_10_t_50000_000000_0.vtu epsilon epsilon 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_0.vtu PointHeatSource_ts_10_t_50000_000000_0.vtu sigma sigma 1e-10 1.0e-6
#
    PointHeatSource_ts_10_t_50000_000000_1.vtu PointHeatSource_ts_10_t_50000_000000_1.vtu displacement displacement 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_1.vtu PointHeatSource_ts_10_t_50000_000000_1.vtu pressure pressure 1e-10 1.0e-6
    PointHeatSource_ts_10_t_50000_000000_1.vtu PointHeatSource_ts_10_t_50000_000000_1.vtu temperature temperature 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_1.vtu PointHeatSource_ts_10_t_50000_000000_1.vtu epsilon epsilon 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_1.vtu PointHeatSource_ts_10_t_50000_000000_1.vtu sigma sigma 1e-10 1.0e-6
#
    PointHeatSource_ts_10_t_50000_000000_2.vtu PointHeatSource_ts_10_t_50000_000000_2.vtu displacement displacement 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_2.vtu PointHeatSource_ts_10_t_50000_000000_2.vtu pressure pressure 1e-10 1.0e-6
    PointHeatSource_ts_10_t_50000_000000_2.vtu PointHeatSource_ts_10_t_50000_000000_2.vtu temperature temperature 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_2.vtu PointHeatSource_ts_10_t_50000_000000_2.vtu epsilon epsilon 1e-10 1.0e-9
    PointHeatSource_ts_10_t_50000_000000_2.vtu PointHeatSource_ts_10_t_50000_000000_2.vtu sigma sigma 1e-10 1.0e-6
)

AddTest(
    NAME ParallelFEM_ThermoRichardsMechanics_point_heat_injection_gml
    PATH ThermoRichardsMechanics/PointHeatSource
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS point_heat_source_2D_gml.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    LABELS "petsc-mumps"
    DIFF_DATA
    PointHeatSource_gml_ts_10_t_50000_000000_0.vtu PointHeatSource_gml_ts_10_t_50000_000000_0.vtu displacement displacement 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_0.vtu PointHeatSource_gml_ts_10_t_50000_000000_0.vtu pressure pressure 1e-10 1.0e-6
    PointHeatSource_gml_ts_10_t_50000_000000_0.vtu PointHeatSource_gml_ts_10_t_50000_000000_0.vtu temperature temperature 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_0.vtu PointHeatSource_gml_ts_10_t_50000_000000_0.vtu epsilon epsilon 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_0.vtu PointHeatSource_gml_ts_10_t_50000_000000_0.vtu sigma sigma 1e-10 1.0e-6
#
    PointHeatSource_gml_ts_10_t_50000_000000_1.vtu PointHeatSource_gml_ts_10_t_50000_000000_1.vtu displacement displacement 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_1.vtu PointHeatSource_gml_ts_10_t_50000_000000_1.vtu pressure pressure 1e-10 1.0e-6
    PointHeatSource_gml_ts_10_t_50000_000000_1.vtu PointHeatSource_gml_ts_10_t_50000_000000_1.vtu temperature temperature 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_1.vtu PointHeatSource_gml_ts_10_t_50000_000000_1.vtu epsilon epsilon 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_1.vtu PointHeatSource_gml_ts_10_t_50000_000000_1.vtu sigma sigma 1e-10 1.0e-6
#
    PointHeatSource_gml_ts_10_t_50000_000000_2.vtu PointHeatSource_gml_ts_10_t_50000_000000_2.vtu displacement displacement 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_2.vtu PointHeatSource_gml_ts_10_t_50000_000000_2.vtu pressure pressure 1e-10 1.0e-6
    PointHeatSource_gml_ts_10_t_50000_000000_2.vtu PointHeatSource_gml_ts_10_t_50000_000000_2.vtu temperature temperature 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_2.vtu PointHeatSource_gml_ts_10_t_50000_000000_2.vtu epsilon epsilon 1e-10 1.0e-9
    PointHeatSource_gml_ts_10_t_50000_000000_2.vtu PointHeatSource_gml_ts_10_t_50000_000000_2.vtu sigma sigma 1e-10 1.0e-6
)

AddTest(
    NAME ParallelFEM_ThermoRichardsMechanics_TaskCDECOVALEX2023
    PATH ThermoRichardsMechanics/TaskCDECOVALEX2023
    RUNTIME 10
    EXECUTABLE ogs
    EXECUTABLE_ARGS Decovalex-0.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    LABELS "petsc-mumps"
    DIFF_DATA
#
    Decovalex-0_ts_10_t_864000_000000_0.vtu Decovalex-0_ts_10_t_864000_000000_0.vtu displacement displacement 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_0.vtu Decovalex-0_ts_10_t_864000_000000_0.vtu pressure pressure 1e-10 1.0e-6
    Decovalex-0_ts_10_t_864000_000000_0.vtu Decovalex-0_ts_10_t_864000_000000_0.vtu saturation saturation 1e-10 1.0e-6
    Decovalex-0_ts_10_t_864000_000000_0.vtu Decovalex-0_ts_10_t_864000_000000_0.vtu temperature temperature 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_0.vtu Decovalex-0_ts_10_t_864000_000000_0.vtu epsilon epsilon 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_0.vtu Decovalex-0_ts_10_t_864000_000000_0.vtu sigma sigma 1e-10 1.0e-6
#
    Decovalex-0_ts_10_t_864000_000000_1.vtu Decovalex-0_ts_10_t_864000_000000_1.vtu displacement displacement 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_1.vtu Decovalex-0_ts_10_t_864000_000000_1.vtu pressure pressure 1e-10 1.0e-6
    Decovalex-0_ts_10_t_864000_000000_1.vtu Decovalex-0_ts_10_t_864000_000000_1.vtu saturation saturation 1e-10 1.0e-6
    Decovalex-0_ts_10_t_864000_000000_1.vtu Decovalex-0_ts_10_t_864000_000000_1.vtu temperature temperature 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_1.vtu Decovalex-0_ts_10_t_864000_000000_1.vtu epsilon epsilon 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_1.vtu Decovalex-0_ts_10_t_864000_000000_1.vtu sigma sigma 1e-10 1.0e-6
#
    Decovalex-0_ts_10_t_864000_000000_2.vtu Decovalex-0_ts_10_t_864000_000000_2.vtu displacement displacement 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_2.vtu Decovalex-0_ts_10_t_864000_000000_2.vtu pressure pressure 1e-10 1.0e-6
    Decovalex-0_ts_10_t_864000_000000_2.vtu Decovalex-0_ts_10_t_864000_000000_2.vtu saturation saturation 1e-10 1.0e-6
    Decovalex-0_ts_10_t_864000_000000_2.vtu Decovalex-0_ts_10_t_864000_000000_2.vtu temperature temperature 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_2.vtu Decovalex-0_ts_10_t_864000_000000_2.vtu epsilon epsilon 1e-10 1.0e-9
    Decovalex-0_ts_10_t_864000_000000_2.vtu Decovalex-0_ts_10_t_864000_000000_2.vtu sigma sigma 1e-10 1.0e-6
)
# ThermoRichardsMechanics; thermo_osmosis and thermo_filtration effects, linear poroelastic, column consolidation
AddTest(
    NAME ThermoRichardsMechanics_thermo_osmosis_filtration_effects_Column
    PATH ThermoRichardsMechanics/ThermoOsmosis
    RUNTIME 15
    EXECUTABLE ogs
    EXECUTABLE_ARGS Column.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_Column_ts_68_t_7200000.000000.vtu Column_ts_68_t_7200000.000000.vtu displacement displacement 1e-5 1e-5
    expected_Column_ts_68_t_7200000.000000.vtu Column_ts_68_t_7200000.000000.vtu pressure pressure 1e-5 1e-5
    expected_Column_ts_68_t_7200000.000000.vtu Column_ts_68_t_7200000.000000.vtu temperature temperature 1e-5 1e-5
    expected_Column_ts_68_t_7200000.000000.vtu Column_ts_68_t_7200000.000000.vtu epsilon epsilon 1e-5 1e-5
    expected_Column_ts_68_t_7200000.000000.vtu Column_ts_68_t_7200000.000000.vtu sigma sigma 1e-5 1e-5
)
# ThermoRichardsMechanics; test for removing body force from displacement equation
AddTest(
    NAME ThermoRichardsMechanics_dont_apply_body_force_for_deformation
    PATH ThermoRichardsMechanics/BodyForce
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS square.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_ts_10_t_864000.000000.vtu displacement displacement 1e-10 1e-10
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_ts_10_t_864000.000000.vtu pressure pressure 5e-2 1e-8
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_ts_10_t_864000.000000.vtu temperature temperature 1e-10 1e-10
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_ts_10_t_864000.000000.vtu epsilon epsilon 1e-10 1e-10
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_ts_10_t_864000.000000.vtu sigma sigma 5e-2 1e-8
)

AddTest(
    NAME ThermoRichardsMechanics_total_initial_stress_dont_apply_body_force_for_deformation
    PATH ThermoRichardsMechanics/BodyForce
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_total_stress_test.xml
    WRAPPER time
    TESTER vtkdiff
    DIFF_DATA
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_total_stess0_test_ts_10_t_864000.000000.vtu displacement displacement 1e-10 1e-10
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_total_stess0_test_ts_10_t_864000.000000.vtu pressure pressure 5e-2 1e-8
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_total_stess0_test_ts_10_t_864000.000000.vtu temperature temperature 1e-10 1e-10
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_total_stess0_test_ts_10_t_864000.000000.vtu epsilon epsilon 1e-10 1e-10
    expected_dont_apply_body_force_for_deformation_ts_0_t_0.000000.vtu dont_apply_body_force_for_deformation_total_stess0_test_ts_10_t_864000.000000.vtu sigma sigma 5e-2 1e-8
)

if(OGS_USE_MFRONT)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MultiMaterialEhlers/square_1e1_2_matIDs.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MultiMaterialEhlers/square_1e1_2_matIDs_restart.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/ThermoPoroElasticity/uniaxial_isothermal_drainage_imbibition_basic_mfront_model_ctest.xml)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/ThermoPoroElasticity/uniaxial_isothermal_drainage_imbibition_extended_mfront_model_ctest.xml)

    if (NOT OGS_USE_MPI)
        OgsTest(PROJECTFILE ThermoRichardsMechanics/MFront/A2/A2.xml RUNTIME 18)
        AddTest(
            NAME ThermoRichardsMechanics_A2_effective_initial_stress
            PATH ThermoRichardsMechanics/MFront/A2
            RUNTIME 1
            EXECUTABLE ogs
            EXECUTABLE_ARGS A2_effective_stress0.xml
            WRAPPER time
            TESTER vtkdiff
            DIFF_DATA
            A2_ts_76_t_2764800.000000.vtu A2_effective_stess0_test_ts_76_t_2764800.000000.vtu displacement displacement 1e-10 1e-10
            A2_ts_76_t_2764800.000000.vtu A2_effective_stess0_test_ts_76_t_2764800.000000.vtu pressure pressure 1e-9 1e-8
            A2_ts_76_t_2764800.000000.vtu A2_effective_stess0_test_ts_76_t_2764800.000000.vtu temperature temperature 1e-10 1e-10
            A2_ts_76_t_2764800.000000.vtu A2_effective_stess0_test_ts_76_t_2764800.000000.vtu epsilon epsilon 1e-10 1e-10
            mfront_A2_ts_76_t_2764800.000000.vtu A2_effective_stess0_test_ts_76_t_2764800.000000.vtu sigma_total sigma_total 6e-8 1e-8
        )
    endif()
endif()

if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/Mockup2D/mockup.prj RUNTIME 60)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/Mockup2D/mockup_restart.xml RUNTIME 30)
endif()
