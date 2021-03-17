# HydroMechanics; Small deformations, linear poroelastic (HML)

### With monolithic scheme
if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE HydroMechanics/Linear/Confined_Compression/square_1e2.prj RUNTIME 9)
endif()
if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE HydroMechanics/Linear/Confined_Compression/square_1e2_linear.prj)
endif()
if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE HydroMechanics/Linear/DrainageEexcavation/HMdrainage.prj RUNTIME 330)
endif()
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad9_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_quad9.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_quad9_ts_1_t_5.000000.vtu square_1e2_quad9_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_quad9_ts_1_t_5.000000.vtu square_1e2_quad9_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
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
    expected_square_1e2_tri_ts_1_t_5.000000.vtu square_1e2_tri_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_tri_ts_1_t_5.000000.vtu square_1e2_tri_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
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
    expected_cube_1e3_ts_1_t_5.000000.vtu cube_1e3_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_cube_1e3_ts_1_t_5.000000.vtu cube_1e3_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
)

AddTest(
    NAME HydroMechanics_hm2_1D1bt
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1D1bt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
    hm2_1D1bt_ts_50_t_5.000000.vtu hm2_1D1bt_ts_50_t_5.000000.vtu pressure pressure 1.0e-11 0.0
    hm2_1D1bt_ts_50_t_5.000000.vtu hm2_1D1bt_ts_50_t_5.000000.vtu displacement displacement 1.0e-11 0.0
    hm2_1D1bt_ts_50_t_5.000000.vtu hm2_1D1bt_ts_50_t_5.000000.vtu sigma sigma 1.0e-10 0.0
    hm2_1D1bt_ts_50_t_5.000000.vtu hm2_1D1bt_ts_50_t_5.000000.vtu epsilon epsilon 1.0e-11 0.0
    hm2_1D1bt_ts_100_t_10.000000.vtu hm2_1D1bt_ts_100_t_10.000000.vtu pressure pressure 2.0e-11 0.0
    hm2_1D1bt_ts_100_t_10.000000.vtu hm2_1D1bt_ts_100_t_10.000000.vtu displacement displacement 1.0e-11 0.0
    hm2_1D1bt_ts_100_t_10.000000.vtu hm2_1D1bt_ts_100_t_10.000000.vtu sigma sigma 5.0e-10 0.0
    hm2_1D1bt_ts_100_t_10.000000.vtu hm2_1D1bt_ts_100_t_10.000000.vtu epsilon epsilon 1.0e-11 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1D2bt
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1D2bt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 38
    DIFF_DATA
    hm2_1D2bt_ts_125_t_5.000000.vtu hm2_1D2bt_ts_125_t_5.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1D2bt_ts_125_t_5.000000.vtu hm2_1D2bt_ts_125_t_5.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1D2bt_ts_125_t_5.000000.vtu hm2_1D2bt_ts_125_t_5.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1D2bt_ts_125_t_5.000000.vtu hm2_1D2bt_ts_125_t_5.000000.vtu epsilon epsilon 1.0e-10 0.0
    hm2_1D2bt_ts_250_t_10.000000.vtu hm2_1D2bt_ts_250_t_10.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1D2bt_ts_250_t_10.000000.vtu hm2_1D2bt_ts_250_t_10.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1D2bt_ts_250_t_10.000000.vtu hm2_1D2bt_ts_250_t_10.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1D2bt_ts_250_t_10.000000.vtu hm2_1D2bt_ts_250_t_10.000000.vtu epsilon epsilon 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1Dbiot
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1Dbiot.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 9
    DIFF_DATA
    hm2_1Dbiot_ts_20_t_5.000000.vtu hm2_1Dbiot_ts_20_t_5.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dbiot_ts_20_t_5.000000.vtu hm2_1Dbiot_ts_20_t_5.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dbiot_ts_20_t_5.000000.vtu hm2_1Dbiot_ts_20_t_5.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1Dbiot_ts_20_t_5.000000.vtu hm2_1Dbiot_ts_20_t_5.000000.vtu epsilon epsilon 1.0e-10 0.0
    hm2_1Dbiot_ts_40_t_10.000000.vtu hm2_1Dbiot_ts_40_t_10.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dbiot_ts_40_t_10.000000.vtu hm2_1Dbiot_ts_40_t_10.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dbiot_ts_40_t_10.000000.vtu hm2_1Dbiot_ts_40_t_10.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1Dbiot_ts_40_t_10.000000.vtu hm2_1Dbiot_ts_40_t_10.000000.vtu epsilon epsilon 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1Dcolumn1
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1Dcolumn1.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 28
    DIFF_DATA
    hm2_1Dcolumn1_ts_120_t_6.000000.vtu hm2_1Dcolumn1_ts_120_t_6.000000.vtu pressure pressure 5.0e-10 0.0
    hm2_1Dcolumn1_ts_120_t_6.000000.vtu hm2_1Dcolumn1_ts_120_t_6.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn1_ts_120_t_6.000000.vtu hm2_1Dcolumn1_ts_120_t_6.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1Dcolumn1_ts_120_t_6.000000.vtu hm2_1Dcolumn1_ts_120_t_6.000000.vtu epsilon epsilon 1.0e-10 0.0
    hm2_1Dcolumn1_ts_240_t_12.000000.vtu hm2_1Dcolumn1_ts_240_t_12.000000.vtu pressure pressure 5.0e-10 0.0
    hm2_1Dcolumn1_ts_240_t_12.000000.vtu hm2_1Dcolumn1_ts_240_t_12.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn1_ts_240_t_12.000000.vtu hm2_1Dcolumn1_ts_240_t_12.000000.vtu sigma sigma 5.0e-9 0.0
    hm2_1Dcolumn1_ts_240_t_12.000000.vtu hm2_1Dcolumn1_ts_240_t_12.000000.vtu epsilon epsilon 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1Dcolumn2
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1Dcolumn2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 9
    DIFF_DATA
    hm2_1Dcolumn2_ts_16_t_4.000000.vtu hm2_1Dcolumn2_ts_16_t_4.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dcolumn2_ts_16_t_4.000000.vtu hm2_1Dcolumn2_ts_16_t_4.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn2_ts_16_t_4.000000.vtu hm2_1Dcolumn2_ts_16_t_4.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1Dcolumn2_ts_16_t_4.000000.vtu hm2_1Dcolumn2_ts_16_t_4.000000.vtu epsilon epsilon 1.0e-10 0.0
    hm2_1Dcolumn2_ts_40_t_10.000000.vtu hm2_1Dcolumn2_ts_40_t_10.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dcolumn2_ts_40_t_10.000000.vtu hm2_1Dcolumn2_ts_40_t_10.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn2_ts_40_t_10.000000.vtu hm2_1Dcolumn2_ts_40_t_10.000000.vtu sigma sigma 1.0e-9 0.0
    hm2_1Dcolumn2_ts_40_t_10.000000.vtu hm2_1Dcolumn2_ts_40_t_10.000000.vtu epsilon epsilon 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_2Dmandel
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_2Dmandel.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 166
    DIFF_DATA
    hm2_2Dmandel_ts_40_t_2.000000.vtu hm2_2Dmandel_ts_40_t_2.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_2Dmandel_ts_40_t_2.000000.vtu hm2_2Dmandel_ts_40_t_2.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_2Dmandel_ts_40_t_2.000000.vtu hm2_2Dmandel_ts_40_t_2.000000.vtu sigma sigma 3.0e-9 0.0
    hm2_2Dmandel_ts_40_t_2.000000.vtu hm2_2Dmandel_ts_40_t_2.000000.vtu epsilon epsilon 1.0e-10 0.0
    hm2_2Dmandel_ts_160_t_8.000000.vtu hm2_2Dmandel_ts_160_t_8.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_2Dmandel_ts_160_t_8.000000.vtu hm2_2Dmandel_ts_160_t_8.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_2Dmandel_ts_160_t_8.000000.vtu hm2_2Dmandel_ts_160_t_8.000000.vtu sigma sigma 3.0e-9 0.0
    hm2_2Dmandel_ts_160_t_8.000000.vtu hm2_2Dmandel_ts_160_t_8.000000.vtu epsilon epsilon 1.0e-10 0.0
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
    expected_square_1e2_UC_early_ts_10_t_1.000000.vtu square_1e2_UC_early_ts_10_t_1.000000.vtu displacement displacement 1e-11 1e-16
    expected_square_1e2_UC_early_ts_10_t_1.000000.vtu square_1e2_UC_early_ts_10_t_1.000000.vtu pressure pressure 1e-10 1e-16
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
    expected_square_1e2_UC_late_ts_10_t_1000.000000.vtu square_1e2_UC_late_ts_10_t_1000.000000.vtu displacement displacement 1e-13 1e-16
    expected_square_1e2_UC_late_ts_10_t_1000.000000.vtu square_1e2_UC_late_ts_10_t_1000.000000.vtu pressure pressure 1e-13 1e-16
)

AddTest(
    NAME HydroMechanics_HML_flow_gravity
    PATH HydroMechanics/Linear/Gravity
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_gravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flow_gravity_ts_16_t_40000000.000000.vtu flow_gravity_ts_16_t_40000000.000000.vtu displacement displacement 1e-14 0
    flow_gravity_ts_16_t_40000000.000000.vtu flow_gravity_ts_16_t_40000000.000000.vtu pressure pressure 1e-10 0
    flow_gravity_ts_16_t_40000000.000000.vtu flow_gravity_ts_16_t_40000000.000000.vtu velocity velocity 1e-10 0
)

# Tests for Principal Stress Output
AddTest(
    NAME HydroMechanics_hollow_sphere
    PATH HydroMechanics/Principal_Stress/Hollow_Sphere
    EXECUTABLE ogs
    EXECUTABLE_ARGS sphere.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    output_ts_1_t_1.000000.vtu output_ts_1_t_1.000000.vtu displacement displacement 0 1e-12
    output_ts_1_t_1.000000.vtu output_ts_1_t_1.000000.vtu principal_stress_values principal_stress_values 0 1e-10
)

AddTest(
    NAME HydroMechanics_tube
    PATH HydroMechanics/Principal_Stress/Tube
    EXECUTABLE ogs
    EXECUTABLE_ARGS tube.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    output_ts_1_t_1.000000.vtu output_ts_1_t_1.000000.vtu displacement displacement 1e-12 0
    output_ts_1_t_1.000000.vtu output_ts_1_t_1.000000.vtu principal_stress_values principal_stress_values 0 1e-10
)

## Tests for Ideal gas
# flow_no_strain
AddTest(
    NAME HydroMechanics_IdealGas_flow_no_strain
    PATH HydroMechanics/IdealGas/flow_no_strain
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_no_strain.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 40
    DIFF_DATA
    # tolerance could be lowered with smaller time steps, gut would increase runtime
    # TODO (FZill) another solver might provide better results, same goes for "flow_free_expansion"
    flow_no_strain_ts_1000_t_100.000000.vtu flow_no_strain_ts_1000_t_100.000000.vtu pressure pressure 0 1e-6
    flow_no_strain_ts_1000_t_100.000000.vtu flow_no_strain_ts_1000_t_100.000000.vtu displacement displacement 1e-12 0
)

# flow_free_expansion
AddTest(
    NAME HydroMechanics_IdealGas_flow_free_expansion
    PATH HydroMechanics/IdealGas/flow_free_expansion
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_free_expansion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flow_free_expansion_ts_1000_t_10000.000000.vtu flow_free_expansion_ts_1000_t_10000.000000.vtu pressure pressure 0 5e-7
    flow_free_expansion_ts_1000_t_10000.000000.vtu flow_free_expansion_ts_1000_t_10000.000000.vtu displacement displacement 1e-11 0
)

# flow_pressure_boundary
AddTest(
    NAME HydroMechanics_IdealGas_flow_pressure_boundary
    PATH HydroMechanics/IdealGas/flow_pressure_boundary
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_pressure_boundary.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flow_pressure_boundary_ts_100_t_4000.000000.vtu flow_pressure_boundary_ts_100_t_4000.000000.vtu pressure pressure 0 1e-13
    flow_pressure_boundary_ts_100_t_4000.000000.vtu flow_pressure_boundary_ts_100_t_4000.000000.vtu displacement displacement 1e-12 0
)

# Permeability models
AddTest(
    NAME HydroMechanics_Permeability_EmbeddedFracture
    PATH HydroMechanics/EmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_ts_1_t_1.000000.vtu cube_ts_1_t_1.000000.vtu pressure pressure 0 1e-14
    cube_ts_1_t_1.000000.vtu cube_ts_1_t_1.000000.vtu velocity velocity 1e-15 0
    cube_ts_1_t_1.000000.vtu cube_ts_1_t_1.000000.vtu displacement displacement 1e-15 0
)

AddTest(
    NAME FailureIndexDependentPermeability
    PATH HydroMechanics/FailureIndexDependentPermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS quad_with_half_hole.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    quad_with_half_hole_ts_1_t_1.000000.vtu quad_with_half_hole_ts_1_t_1.000000.vtu pressure pressure 1.0e-10 1.e-10
    quad_with_half_hole_ts_1_t_1.000000.vtu quad_with_half_hole_ts_1_t_1.000000.vtu displacement displacement 1.0e-10 1.e-10
    quad_with_half_hole_ts_1_t_1.000000.vtu quad_with_half_hole_ts_1_t_1.000000.vtu sigma sigma 1.0e-6 1.e-9
    quad_with_half_hole_ts_1_t_1.000000.vtu quad_with_half_hole_ts_1_t_1.000000.vtu permeability permeability 1.0e-6 1.e-9
)

AddTest(
    NAME HydroMechanicsStrainDependentPermeability
    PATH HydroMechanics/StrainDependentPermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS gas_loading.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 21
    DIFF_DATA
    output_curve_ts_2_t_400.000000.vtu  output_curve_ts_2_t_400.000000.vtu pressure pressure 1.0e-10 1.e-10
    output_curve_ts_2_t_400.000000.vtu  output_curve_ts_2_t_400.000000.vtu displacement displacement 1.0e-10 1.e-10
    output_curve_ts_2_t_400.000000.vtu  output_curve_ts_2_t_400.000000.vtu sigma sigma 1.0e-6 1.e-9
    output_curve_ts_2_t_400.000000.vtu  output_curve_ts_2_t_400.000000.vtu permeability permeability 1.0e-6 1.e-9
)

AddTest(
    NAME MonolithicInjectionProduction1D
    PATH HydroMechanics/StaggeredScheme/InjectionProduction1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS InjectionProduction1DMono.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu displacement displacement 1e-14 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu pressure pressure 1e-7 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu velocity velocity 1e-19 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu epsilon epsilon 1e-14 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu sigma sigma 1e-7 0
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
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu displacement displacement 1e-13 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu pressure pressure 1e-6 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu velocity velocity 1e-18 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu epsilon epsilon 1e-13 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu sigma sigma 1e-6 0
)

AddTest(
    NAME HydroMechanics_hm1_1Dbeam
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_1Dbeam.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_1Dbeam_ts_1_t_1.000000.vtu hm1_1Dbeam_ts_1_t_1.000000.vtu pressure pressure 1.0e-6 0.0
    hm1_1Dbeam_ts_1_t_1.000000.vtu hm1_1Dbeam_ts_1_t_1.000000.vtu displacement displacement 1.0e-16 0.0
    hm1_1Dbeam_ts_1_t_1.000000.vtu hm1_1Dbeam_ts_1_t_1.000000.vtu sigma sigma 1.0e-6 0.0
    hm1_1Dbeam_ts_1_t_1.000000.vtu hm1_1Dbeam_ts_1_t_1.000000.vtu epsilon epsilon 1.0e-16 0.0
)

AddTest(
    NAME HydroMechanics_hm1_2Dsquare
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_2Dsquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_2Dsquare_ts_1_t_1.000000.vtu hm1_2Dsquare_ts_1_t_1.000000.vtu pressure pressure 1.0e-6 0.0
    hm1_2Dsquare_ts_1_t_1.000000.vtu hm1_2Dsquare_ts_1_t_1.000000.vtu displacement displacement 1.0e-16 0.0
    hm1_2Dsquare_ts_1_t_1.000000.vtu hm1_2Dsquare_ts_1_t_1.000000.vtu sigma sigma 1.0e-6 0.0
    hm1_2Dsquare_ts_1_t_1.000000.vtu hm1_2Dsquare_ts_1_t_1.000000.vtu epsilon epsilon 1.0e-16 0.0
)

AddTest(
    NAME HydroMechanics_hm1_3Dcube
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_3Dcube.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_3Dcube_ts_1_t_1.000000.vtu hm1_3Dcube_ts_1_t_1.000000.vtu pressure pressure 1.0e-5 0.0
    hm1_3Dcube_ts_1_t_1.000000.vtu hm1_3Dcube_ts_1_t_1.000000.vtu displacement displacement 1.0e-15 0.0
    hm1_3Dcube_ts_1_t_1.000000.vtu hm1_3Dcube_ts_1_t_1.000000.vtu sigma sigma 1.0e-5 0.0
    hm1_3Dcube_ts_1_t_1.000000.vtu hm1_3Dcube_ts_1_t_1.000000.vtu epsilon epsilon 1.0e-15 0.0
)

AddTest(
    NAME HydroMechanics_hm1_3Dgravity
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_3Dgravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_3Dgravity_ts_1_t_1.000000.vtu hm1_3Dgravity_ts_1_t_1.000000.vtu pressure pressure 1.0e-6 0.0
    hm1_3Dgravity_ts_1_t_1.000000.vtu hm1_3Dgravity_ts_1_t_1.000000.vtu displacement displacement 1.0e-15 0.0
    hm1_3Dgravity_ts_1_t_1.000000.vtu hm1_3Dgravity_ts_1_t_1.000000.vtu sigma sigma 1.0e-6 0.0
    hm1_3Dgravity_ts_1_t_1.000000.vtu hm1_3Dgravity_ts_1_t_1.000000.vtu epsilon epsilon 1.0e-15 0.0
)
