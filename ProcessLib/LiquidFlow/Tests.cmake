# Liquid flow
AddTest(
    NAME LiquidFlow_LineDirichletNeumannBC
    PATH Parabolic/LiquidFlow/LineDirichletNeumannBC
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_dirichlet_neumannBC.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    sat1D.vtu sat_1D_ts_1_t_1.000000.vtu AnalyticPressure pressure 1e-8 1e-8
    sat1D.vtu sat_1D_ts_1_t_1.000000.vtu AnalyticVec v 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_PressureBCatCornerOfAnisotropicSquare
    PATH Parabolic/LiquidFlow/PressureBCatCornerOfAnisotropicSquare
    EXECUTABLE ogs
    EXECUTABLE_ARGS pressureBC_at_corner_of_anisotropic_square.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    mesh2D.vtu sat_2D_lflow_ts_1_t_1.000000.vtu OGS5_Results pressure 1e-8 1e-8
)

if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/DrainageExcavation/drainage_LiquidFlow.prj)
endif()

AddTest(
    NAME LiquidFlow_GravityDriven
    PATH Parabolic/LiquidFlow/GravityDriven
    EXECUTABLE ogs
    EXECUTABLE_ARGS gravity_driven.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    mesh2D.vtu gravity_driven_ts_1_t_1.000000.vtu AnalyticPressure pressure 1e-8 1e-8
    mesh2D.vtu gravity_driven_ts_1_t_1.000000.vtu v_ref v 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_AxisymTheis
    PATH Parabolic/LiquidFlow/AxiSymTheis
    EXECUTABLE ogs
    EXECUTABLE_ARGS axisym_theis.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    axisym_theis.vtu liquid_pcs_ts_30_t_1728.000000.vtu OGS5_pressure pressure 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_BuildupTest
    PATH Parabolic/LiquidFlow/BuildupTest
    EXECUTABLE ogs
    EXECUTABLE_ARGS buildup_test.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    standard_solution_buildup_test_ts_107_t_424800.000000.vtu buildup_test_ts_107_t_424800.000000.vtu pressure pressure 1.6e-5 0.0
    standard_solution_buildup_test_ts_211_t_720000.000000.vtu buildup_test_ts_211_t_720000.000000.vtu pressure pressure 5e-5 0.0
)

AddTest(
    NAME LiquidFlow_Anisotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    RUNTIME 115
    EXECUTABLE ogs
    EXECUTABLE_ARGS anisotropic_gravity_driven3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hex.vtu anisotropic_gravity_driven3D_ts_1_t_1.000000.vtu analytic_pressure pressure 1e-6 1e-6
)

AddTest(
    NAME LiquidFlow_Isotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    RUNTIME 130
    EXECUTABLE ogs
    EXECUTABLE_ARGS isotropic_gravity_driven3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hex.vtu isotropic_gravity_driven3D_ts_1_t_1.000000.vtu analytic_pressure pressure 1e-6 1e-6
)

AddTest(
    NAME LiquidFlowDirichletBCWithinTimeInterval
    PATH Parabolic/LiquidFlow/TimeIntervalDirichletBC
    EXECUTABLE ogs
    EXECUTABLE_ARGS TimeIntervalDirichletBC.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    mesh2D.vtu dirichlet_bc_wihin_interval_ts_2_t_10.000000.vtu analytical_solution_t_lt_10 pressure 1e-6 1e-12
    mesh2D.vtu dirichlet_bc_wihin_interval_ts_4_t_20.000000.vtu analytical_solution_t_gt_10 pressure 1e-6 1e-12
)

AddTest(
    NAME LiquidFlow_h1_1Dsource
    PATH Parabolic/LiquidFlow/Verification/h1_1Dsource
    EXECUTABLE ogs
    EXECUTABLE_ARGS h1_1Dsource.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h1_1Dsource_ts_1_t_1.000000.vtu h1_1Dsource_ts_1_t_1.000000.vtu pressure pressure 5e-8 0.0
)

AddTest(
    NAME LiquidFlow_h1_1Dsteady
    PATH Parabolic/LiquidFlow/Verification/h1_1Dsteady
    EXECUTABLE ogs
    EXECUTABLE_ARGS h1_1Dsteady.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h1_1Dsteady_ts_1_t_1.000000.vtu h1_1Dsteady_ts_1_t_1.000000.vtu pressure pressure 1e-9 0.0
)

AddTest(
    NAME LiquidFlow_h1_3Dhydstat
    PATH Parabolic/LiquidFlow/Verification/h1_3Dhydstat
    EXECUTABLE ogs
    EXECUTABLE_ARGS h1_3Dhydstat.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h1_3Dhydstat_ts_1_t_1.000000.vtu h1_3Dhydstat_ts_1_t_1.000000.vtu pressure pressure 1e-5 0.0
)

AddTest(
    NAME LiquidFlow_h2_1D1bt
    PATH Parabolic/LiquidFlow/Verification/h2_1D1bt
    EXECUTABLE ogs
    EXECUTABLE_ARGS h2_1D1bt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h2_1D1bt_ts_500_t_21600.000000.vtu h2_1D1bt_ts_500_t_21600.000000.vtu pressure pressure 1e-5 0.0
    h2_1D1bt_ts_1000_t_43200.000000.vtu h2_1D1bt_ts_1000_t_43200.000000.vtu pressure pressure 1e-5 0.0
)

AddTest(
    NAME LiquidFlow_h2_1D2bt
    PATH Parabolic/LiquidFlow/Verification/h2_1D2bt
    EXECUTABLE ogs
    EXECUTABLE_ARGS h2_1D2bt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h2_1D2bt_ts_1500_t_3888.000000.vtu h2_1D2bt_ts_1500_t_3888.000000.vtu pressure pressure 1e-5 0.0
    h2_1D2bt_ts_3000_t_7776.000000.vtu h2_1D2bt_ts_3000_t_7776.000000.vtu pressure pressure 1e-5 0.0
)

#===============================================================================
# PETSc/MPI
AddTest(
    NAME LiquidFlow_LineDirichletNeumannBC
    PATH Parabolic/LiquidFlow/LineDirichletNeumannBC
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_dirichlet_neumannBC.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    sat1D.vtu sat_1D_ts_1_t_1_000000_0.vtu AnalyticPressure pressure 1e-8 1e-8
#    sat1D.vtu sat_1D_ts_1_t_1_000000_0.vtu AnalyticVec v 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_GravityDriven
    PATH Parabolic/LiquidFlow/GravityDriven
    EXECUTABLE ogs
    EXECUTABLE_ARGS gravity_driven.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    mesh2D.vtu gravity_driven_ts_1_t_1_000000_0.vtu AnalyticPressure pressure 1e-8 1e-8
#    mesh2D.vtu gravity_driven_ts_1_t_1_000000_0.vtu v_ref v 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_PressureBCatCornerOfAnisotropicSquare
    PATH Parabolic/LiquidFlow/PressureBCatCornerOfAnisotropicSquare
    EXECUTABLE ogs
    EXECUTABLE_ARGS pressureBC_at_corner_of_anisotropic_square.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    mesh2D.vtu sat_2D_lflow_ts_1_t_1_000000_0.vtu OGS5_Results pressure 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_AxisymTheis
    PATH Parabolic/LiquidFlow/AxiSymTheis
    EXECUTABLE ogs
    EXECUTABLE_ARGS axisym_theis.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    axisym_theis.vtu liquid_pcs_ts_30_t_1728_000000_0.vtu OGS5_pressure pressure 1e-8 1e-8
)
AddTest(
    NAME LiquidFlow_Anisotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    RUNTIME 115
    EXECUTABLE ogs
    EXECUTABLE_ARGS anisotropic_gravity_driven3D.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    hex.vtu anisotropic_gravity_driven3D_ts_1_t_1_000000_0.vtu analytic_pressure pressure 1e-6 1e-6
)

AddTest(
    NAME LiquidFlow_Isotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    RUNTIME 130
    EXECUTABLE ogs
    EXECUTABLE_ARGS isotropic_gravity_driven3D.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    hex.vtu isotropic_gravity_driven3D_ts_1_t_1_000000_0.vtu analytic_pressure pressure 1e-6 1e-6
)

AddTest(
    NAME LiquidFlowDirichletBCWithinTimeInterval
    PATH Parabolic/LiquidFlow/TimeIntervalDirichletBC
    EXECUTABLE ogs
    EXECUTABLE_ARGS TimeIntervalDirichletBC.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    mesh2D.vtu dirichlet_bc_wihin_interval_ts_2_t_10_000000_0.vtu analytical_solution_t_lt_10 pressure 1e-6 1e-12
    mesh2D.vtu dirichlet_bc_wihin_interval_ts_4_t_20_000000_0.vtu analytical_solution_t_gt_10 pressure 1e-6 1e-12
)

# Dupuit
AddTest(
    NAME LiquidFlow_Dupuit_BC_BC
    PATH Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC
    EXECUTABLE ogs
    EXECUTABLE_ARGS TestSet_01.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    Dupuit_TestSet_01_ts_1_t_8640.000000.vtu TestSet_01_ts_1_t_8640.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_2_t_17280.000000.vtu TestSet_01_ts_2_t_17280.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_3_t_25920.000000.vtu TestSet_01_ts_3_t_25920.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_4_t_34560.000000.vtu TestSet_01_ts_4_t_34560.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_5_t_43200.000000.vtu TestSet_01_ts_5_t_43200.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_6_t_51840.000000.vtu TestSet_01_ts_6_t_51840.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_7_t_60480.000000.vtu TestSet_01_ts_7_t_60480.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_8_t_69120.000000.vtu TestSet_01_ts_8_t_69120.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_9_t_77760.000000.vtu TestSet_01_ts_9_t_77760.000000.vtu pressure pressure 5e-10 1e-11
    Dupuit_TestSet_01_ts_10_t_86400.000000.vtu TestSet_01_ts_10_t_86400.000000.vtu pressure pressure 5e-10 1e-11
)

AddTest(
    NAME LiquidFlow_Dupuit_BC_BC_Recharge
    PATH Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC_RECHARGE
    EXECUTABLE ogs
    EXECUTABLE_ARGS TestSet_01.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    Dupuit_TestSet_01_ts_1_t_8640.000000.vtu TestSet_01_ts_1_t_8640.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_2_t_17280.000000.vtu TestSet_01_ts_2_t_17280.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_3_t_25920.000000.vtu TestSet_01_ts_3_t_25920.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_4_t_34560.000000.vtu TestSet_01_ts_4_t_34560.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_5_t_43200.000000.vtu TestSet_01_ts_5_t_43200.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_6_t_51840.000000.vtu TestSet_01_ts_6_t_51840.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_7_t_60480.000000.vtu TestSet_01_ts_7_t_60480.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_8_t_69120.000000.vtu TestSet_01_ts_8_t_69120.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_9_t_77760.000000.vtu TestSet_01_ts_9_t_77760.000000.vtu pressure pressure 2e-10 1e-11
    Dupuit_TestSet_01_ts_10_t_86400.000000.vtu TestSet_01_ts_10_t_86400.000000.vtu pressure pressure 2e-10 1e-11
)

AddTest(
    NAME LiquidFlow_Dupuit_BC_BC_Recharge2
    PATH Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC_RECHARGE2
    EXECUTABLE ogs
    EXECUTABLE_ARGS TestSet_01.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    Dupuit_TestSet_01_ts_1_t_8640.000000.vtu TestSet_01_ts_1_t_8640.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_2_t_17280.000000.vtu TestSet_01_ts_2_t_17280.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_3_t_25920.000000.vtu TestSet_01_ts_3_t_25920.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_4_t_34560.000000.vtu TestSet_01_ts_4_t_34560.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_5_t_43200.000000.vtu TestSet_01_ts_5_t_43200.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_6_t_51840.000000.vtu TestSet_01_ts_6_t_51840.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_7_t_60480.000000.vtu TestSet_01_ts_7_t_60480.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_8_t_69120.000000.vtu TestSet_01_ts_8_t_69120.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_9_t_77760.000000.vtu TestSet_01_ts_9_t_77760.000000.vtu pressure pressure 5e-10 1.8e-11
    Dupuit_TestSet_01_ts_10_t_86400.000000.vtu TestSet_01_ts_10_t_86400.000000.vtu pressure pressure 5e-10 1.8e-11
)

AddTest(
    NAME LiquidFlow_Dupuit_BC_BC_Storage
    PATH Parabolic/LiquidFlow/Unconfined_Aquifer/BC_BC_STORAGE
    EXECUTABLE ogs
    EXECUTABLE_ARGS TestSet_01.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    Dupuit_TestSet_01_ts_1_t_86400.000000.vtu TestSet_01_ts_1_t_86400.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_2_t_172800.000000.vtu TestSet_01_ts_2_t_172800.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_3_t_259200.000000.vtu TestSet_01_ts_3_t_259200.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_4_t_345600.000000.vtu TestSet_01_ts_4_t_345600.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_5_t_432000.000000.vtu TestSet_01_ts_5_t_432000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_6_t_518400.000000.vtu TestSet_01_ts_6_t_518400.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_7_t_604800.000000.vtu TestSet_01_ts_7_t_604800.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_8_t_691200.000000.vtu TestSet_01_ts_8_t_691200.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_9_t_777600.000000.vtu TestSet_01_ts_9_t_777600.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_10_t_864000.000000.vtu TestSet_01_ts_10_t_864000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_20_t_1728000.000000.vtu TestSet_01_ts_20_t_1728000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_30_t_2592000.000000.vtu TestSet_01_ts_30_t_2592000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_40_t_3456000.000000.vtu TestSet_01_ts_40_t_3456000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_50_t_4320000.000000.vtu TestSet_01_ts_50_t_4320000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_60_t_5184000.000000.vtu TestSet_01_ts_60_t_5184000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_70_t_6048000.000000.vtu TestSet_01_ts_70_t_6048000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_80_t_6912000.000000.vtu TestSet_01_ts_80_t_6912000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_90_t_7776000.000000.vtu TestSet_01_ts_90_t_7776000.000000.vtu pressure pressure 3e-12 2e-13
    Dupuit_TestSet_01_ts_100_t_8640000.000000.vtu TestSet_01_ts_100_t_8640000.000000.vtu pressure pressure 3e-12 2e-13
)

AddTest(
    NAME LiquidFlow_TimeDependentHeterogeneousDirichletBCs
    PATH Parabolic/LiquidFlow/TimeDependentHeterogeneousBoundaryConditions
    EXECUTABLE ogs
    EXECUTABLE_ARGS TimeDependentHeterogeneousBoundaryConditions.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    time_dependent_heterogeneous_bcs_ts_1_t_10.000000.vtu time_dependent_heterogeneous_bcs_ts_1_t_10.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_bcs_ts_118_t_1180.000000.vtu time_dependent_heterogeneous_bcs_ts_118_t_1180.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_bcs_ts_119_t_1190.000000.vtu time_dependent_heterogeneous_bcs_ts_119_t_1190.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_bcs_ts_120_t_1200.000000.vtu time_dependent_heterogeneous_bcs_ts_120_t_1200.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_bcs_ts_200_t_2000.000000.vtu time_dependent_heterogeneous_bcs_ts_200_t_2000.000000.vtu pressure pressure 1e-7 1e-13
)

AddTest(
    NAME LiquidFlow_TimeDependentHeterogeneousSourceTerm
    PATH Parabolic/LiquidFlow/TimeDependentHeterogeneousSourceTerm
    EXECUTABLE ogs
    EXECUTABLE_ARGS TimeDependentHeterogeneousSourceTerm.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    time_dependent_heterogeneous_source_term_ts_1_t_10.000000.vtu time_dependent_heterogeneous_source_term_ts_1_t_10.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_source_term_ts_118_t_1180.000000.vtu time_dependent_heterogeneous_source_term_ts_118_t_1180.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_source_term_ts_119_t_1190.000000.vtu time_dependent_heterogeneous_source_term_ts_119_t_1190.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_source_term_ts_120_t_1200.000000.vtu time_dependent_heterogeneous_source_term_ts_120_t_1200.000000.vtu pressure pressure 1e-7 1e-13
    time_dependent_heterogeneous_source_term_ts_200_t_2000.000000.vtu time_dependent_heterogeneous_source_term_ts_200_t_2000.000000.vtu pressure pressure 1e-7 1e-13
)

AddTest(
    NAME LiquidFlow_Flux_3D_Hex
    PATH Parabolic/LiquidFlow/Flux
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3_calculatesurfaceflux.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_0_t_0.000000.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_0_t_0.000000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_1_t_0.432000.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_1_t_0.432000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_2_t_0.864000.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface_ts_2_t_0.864000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_cube_1e3_calculatesurfaceflux_ts_1_t_0.432000_expected.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_ts_1_t_0.432000.vtu pressure pressure 1e-7 1e-13
    LF_cube_1e3_calculatesurfaceflux_ts_2_t_0.864000_expected.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_ts_2_t_0.864000.vtu pressure pressure 1e-7 1e-13
    LF_cube_1e3_calculatesurfaceflux_ts_1_t_0.432000_expected.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_ts_1_t_0.432000.vtu HydraulicFlow HydraulicFlow 1e-13 0
    LF_cube_1e3_calculatesurfaceflux_ts_2_t_0.864000_expected.vtu LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_ts_2_t_0.864000.vtu HydraulicFlow HydraulicFlow 1e-13 0
)

AddTest(
    NAME LiquidFlow_Flux_3D_Pyramid
    PATH Parabolic/LiquidFlow/Flux/3D/Pyramid
    EXECUTABLE ogs
    EXECUTABLE_ARGS cuboid_1x1x1_pyramid_6000_calculatesurfaceflux.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_entire_boundary_ts_0_t_0.000000.vtu LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_entire_boundary_ts_0_t_0.000000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_entire_boundary_ts_1_t_0.432000.vtu LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_entire_boundary_ts_1_t_0.432000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_entire_boundary_ts_2_t_0.864000.vtu LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_entire_boundary_ts_2_t_0.864000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_ts_1_t_0.432000_expected.vtu LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_ts_1_t_0.432000.vtu pressure pressure 1e-7 1e-13
    LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_ts_2_t_0.864000_expected.vtu LF_cuboid_1x1x1_pyramid_6000_calculatesurfaceflux_cuboid_1x1x1_pyramid_6000_ts_2_t_0.864000.vtu pressure pressure 1e-7 1e-13
)

AddTest(
    NAME LiquidFlow_Flux_2D_Quads
    PATH Parabolic/LiquidFlow/Flux/2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e1_calculatesurfaceflux.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    LF_square_1e1_surfaceflux_square_1x1_quad_1e1_complete_surface_ts_0_t_0.000000.vtu LF_square_1e1_surfaceflux_square_1x1_quad_1e1_complete_surface_ts_0_t_0.000000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_square_1e1_surfaceflux_square_1x1_quad_1e1_complete_surface_ts_1_t_0.432000.vtu LF_square_1e1_surfaceflux_square_1x1_quad_1e1_complete_surface_ts_1_t_0.432000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_square_1e1_surfaceflux_square_1x1_quad_1e1_complete_surface_ts_2_t_0.864000.vtu LF_square_1e1_surfaceflux_square_1x1_quad_1e1_complete_surface_ts_2_t_0.864000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_square_1e1_surfaceflux_ts_1_t_0.432000_expected.vtu LF_square_1e1_surfaceflux_square_1x1_quad_1e1_ts_1_t_0.432000.vtu pressure pressure 1e-7 1e-13
    LF_square_1e1_surfaceflux_ts_2_t_0.864000_expected.vtu LF_square_1e1_surfaceflux_square_1x1_quad_1e1_ts_2_t_0.864000.vtu pressure pressure 1e-7 1e-13
)

AddTest(
    NAME LiquidFlow_Flux_2D_Tris
    PATH Parabolic/LiquidFlow/Flux/2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1.8e1_calculatesurfaceflux.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_complete_boundary_ts_0_t_0.000000.vtu LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_complete_boundary_ts_0_t_0.000000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_complete_boundary_ts_1_t_0.432000.vtu LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_complete_boundary_ts_1_t_0.432000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_complete_boundary_ts_2_t_0.864000.vtu LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_complete_boundary_ts_2_t_0.864000.vtu specific_flux specific_flux 1e-7 1e-13
    LF_square_1x1_tri_1.8e1_surfaceflux_ts_1_t_0.432000_expected.vtu LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_ts_1_t_0.432000.vtu pressure pressure 1e-7 1e-13
    LF_square_1x1_tri_1.8e1_surfaceflux_ts_2_t_0.864000_expected.vtu LF_square_1x1_tri_1.8e1_surfaceflux_square_1x1_tri_1.8e1_ts_2_t_0.864000.vtu pressure pressure 1e-7 1e-13
)

AddTest(
    NAME LiquidFlow_Flux_3D_HEX_Parallel_2
    PATH Parabolic/LiquidFlow/Flux/3D/Hex/Parallel
    EXECUTABLE ogs
    EXECUTABLE_ARGS cuboid_1x1x1_hex_27_Dirichlet_Dirichlet.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_ts_2_t_86400_000000_0.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_ts_2_t_86400_000000_0.vtu pressure pressure 1e-10 1e-15
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_ts_2_t_86400_000000_1.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_ts_2_t_86400_000000_1.vtu pressure pressure 1e-10 1e-15
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_bottom_boundary_ts_2_t_86400_000000_0.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_bottom_boundary_ts_2_t_86400_000000_0.vtu pressure pressure  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_bottom_boundary_ts_2_t_86400_000000_1.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_bottom_boundary_ts_2_t_86400_000000_1.vtu pressure pressure  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_top_boundary_ts_2_t_86400_000000_0.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_top_boundary_ts_2_t_86400_000000_0.vtu pressure pressure  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_top_boundary_ts_2_t_86400_000000_1.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_top_boundary_ts_2_t_86400_000000_1.vtu pressure pressure  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_0_t_0_000000_0.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_0_t_0_000000_0.vtu specific_flux specific_flux  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_0_t_0_000000_1.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_0_t_0_000000_1.vtu specific_flux specific_flux  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_1_t_43200_000000_0.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_1_t_43200_000000_0.vtu specific_flux specific_flux  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_1_t_43200_000000_1.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_1_t_43200_000000_1.vtu specific_flux specific_flux  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_2_t_86400_000000_0.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_2_t_86400_000000_0.vtu specific_flux specific_flux  1e-7 1e-13
    top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_2_t_86400_000000_1.vtu top_boundary_to_bottom_boundary_cuboid_1x1x1_hex_27_entire_boundary_ts_2_t_86400_000000_1.vtu specific_flux specific_flux  1e-7 1e-13
)

AddTest(
    NAME SimpleSynthetics_XDMF
    PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF
    EXECUTABLE ogs
    EXECUTABLE_ARGS FunctionParameterTest_XDMF.prj
    WRAPPER time
    TESTER xdmfdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf pressure pressure 1e-7 1e-13
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf HydraulicFlow HydraulicFlow 1e-7 1e-13
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf MaterialIDs MaterialIDs 1e-7 1e-13
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf v v 1e-7 1e-13
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf pressure pressure 1e-7 1e-13
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf pressure pressure 1e-7 1e-13
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13
)

AddTest(
    NAME SimpleSynthetics_XDMF_compression_off
    PATH Parabolic/LiquidFlow/SimpleSynthetics/XDMF_compression_off
    EXECUTABLE ogs
    EXECUTABLE_ARGS FunctionParameterTest_XDMF.prj
    WRAPPER time
    TESTER xdmfdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf pressure pressure 1e-7 1e-13
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf HydraulicFlow HydraulicFlow 1e-7 1e-13
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf MaterialIDs MaterialIDs 1e-7 1e-13
    square_5x5_tris_32.xdmf square_5x5_tris_32_square_5x5_tris_32.xdmf v v 1e-7 1e-13
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf pressure pressure 1e-7 1e-13
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13
    square_5x5_tris_32_right_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_right_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf pressure pressure 1e-7 1e-13
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_element_ids bulk_element_ids 1e-7 1e-13
    square_5x5_tris_32_left_boundary.xdmf square_5x5_tris_32_square_5x5_tris_32_left_boundary.xdmf bulk_node_ids bulk_node_ids 1e-7 1e-13
)



#AddTest(
#    NAME LiquidFlow_SimpleSynthetics_constraint_dirichlet_bc
#    PATH Parabolic/LiquidFlow/SimpleSynthetics
#    EXECUTABLE ogs
#    EXECUTABLE_ARGS constraint_bc_1e3.prj
#    TESTER vtkdiff
#    REQUIREMENTS NOT OGS_USE_MPI
#    DIFF_DATA
#    GLOB LF_constraint_bc_1e3_ts_*.vtu p p 1e-15 1e-14
#)

if (NOT (OGS_USE_MPI))
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/PrimaryVariableConstraintDirichletBC/cuboid_1x1x1_hex_1000_Dirichlet_Dirichlet_1.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/PrimaryVariableConstraintDirichletBC/cuboid_1x1x1_hex_1000_Dirichlet_Dirichlet_2.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/PrimaryVariableConstraintDirichletBC/cuboid_1x1x1_hex_1000_Dirichlet_Dirichlet_3.prj)
    OgsTest(PROJECTFILE Parabolic/LiquidFlow/SimpleSynthetics/FunctionParameterTest.prj)
endif()

# inclined mesh
AddTest(
    NAME LiquidFlow_HydrostaticFlowInInclined_2D_Plane
    PATH Parabolic/LiquidFlow/InclinedMeshElements/Inclined2DMesh
    EXECUTABLE ogs
    EXECUTABLE_ARGS hydrostatic_flow_in_inclined_2D_plane.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hydrostatic_flow_in_inclined_2D_plane_ts_t_1.000000.vtu hydrostatic_flow_in_inclined_2D_plane_ts_t_1.000000.vtu pressure pressure 1e-07 1e-13
    hydrostatic_flow_in_inclined_2D_plane_ts_t_1.000000.vtu hydrostatic_flow_in_inclined_2D_plane_ts_t_1.000000.vtu v v 1e-14 1e-14
)

AddTest(
    NAME LiquidFlow_TransientFlowInInclined_2D_Plane
    PATH Parabolic/LiquidFlow/InclinedMeshElements/Inclined2DMesh
    EXECUTABLE ogs
    EXECUTABLE_ARGS transient_flow_in_inclined_2D_plane.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    transient_flow_in_inclined_2D_plane_ts_t_864000.000000.vtu transient_flow_in_inclined_2D_plane_ts_t_864000.000000.vtu pressure pressure 1e-14 1e-11
    transient_flow_in_inclined_2D_plane_ts_t_864000.000000.vtu transient_flow_in_inclined_2D_plane_ts_t_864000.000000.vtu v v 1e-14 1e-14
)

AddTest(
    NAME LiquidFlow_fractures_in_3D
    PATH Parabolic/LiquidFlow/InclinedMeshElements/FractureIn3D
    EXECUTABLE ogs
    EXECUTABLE_ARGS fractures_in_3D.prj
    RUNTIME 10
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    fractures_in_3D_ts_100_t_8640000.000000.vtu fractures_in_3D_ts_100_t_8640000.000000.vtu pressure pressure 1e-14 1e-14
    fractures_in_3D_ts_100_t_8640000.000000.vtu fractures_in_3D_ts_100_t_8640000.000000.vtu v v 1e-14 1e-14
)

AddTest(
    NAME LiquidFlow_line_fractures_in_3D
    PATH Parabolic/LiquidFlow/InclinedMeshElements/1Din3D
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_fractures_in_3D.prj
    RUNTIME 4
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    line_fractures_in_3D_ts_100_t_8640000.000000.vtu line_fractures_in_3D_ts_100_t_8640000.000000.vtu pressure pressure 1e-14 1e-14
    line_fractures_in_3D_ts_100_t_8640000.000000.vtu line_fractures_in_3D_ts_100_t_8640000.000000.vtu v v 1e-14 1e-14
)

AddTest(
    NAME GMSH2OGS_quadratic_mesh_assembly_test
    PATH Utils/GMSH2OGS
    EXECUTABLE ogs
    EXECUTABLE_ARGS quadratic_mesh_assembly_test.prj
    RUNTIME 1
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    quadratic_mesh_assembly_test_ts_1_t_1.000000.vtu quadratic_mesh_assembly_test_ts_1_t_1.000000.vtu pressure pressure 1e-8 1e-12
)
