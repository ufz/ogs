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
    OgsTest(PROJECTFILE HydroMechanics/A2/A2.prj RUNTIME 10)
    OgsTest(PROJECTFILE HydroMechanics/ExcavationNiches/excavation_niches.prj RUNTIME 60)
endif()

# Ground equilibrium
AddTest(
    NAME GroundEquilibrium
    PATH HydroMechanics/GroundEquilibrium
    EXECUTABLE ogs
    EXECUTABLE_ARGS simHM_ground.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_ts_1_t_1000000.000000.vtu pressure_interpolated pressure_interpolated 1e-10 0
    simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_ts_1_t_1000000.000000.vtu displacement displacement 1e-14 0
    simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_ts_1_t_1000000.000000.vtu sigma sigma 1e-7 0
    simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_ts_1_t_1000000.000000.vtu DarcyVelocity DarcyVelocity 1e-7 0
)
if (OGS_USE_PYTHON)
    # same, but using Python BCs
    AddTest(
        NAME GroundEquilibriumPython
        PATH HydroMechanics/GroundEquilibrium
        EXECUTABLE ogs
        EXECUTABLE_ARGS simHM_ground_python.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
		simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_python_ts_1_t_1000000.000000.vtu pressure_interpolated pressure_interpolated 1e-10 0
		simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_python_ts_1_t_1000000.000000.vtu displacement displacement 1e-14 0
		simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_python_ts_1_t_1000000.000000.vtu sigma sigma 1e-7 0
		simHM_ground_ts_1_t_1000000.000000.vtu simHM_ground_python_ts_1_t_1000000.000000.vtu DarcyVelocity DarcyVelocity 1e-7 0
    )
endif()
AddTest(
    NAME GroundEquilibriumQuadBCu
    PATH HydroMechanics/GroundEquilibrium
    EXECUTABLE ogs
    EXECUTABLE_ARGS simHM_ground_quadBCu.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu pressure_interpolated pressure_interpolated 1e-10 0
    simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu displacement displacement 1e-14 0
    simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu sigma sigma 1e-7 0
    simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu DarcyVelocity DarcyVelocity 1e-7 0
)
if (OGS_USE_PYTHON)
    # same, but using Python BCs
    AddTest(
        NAME GroundEquilibriumQuadBCuPython
        PATH HydroMechanics/GroundEquilibrium
        EXECUTABLE ogs
        EXECUTABLE_ARGS simHM_ground_quadBCu_python.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
		simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu pressure_interpolated pressure_interpolated 1e-10 0
		simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu displacement displacement 1e-14 0
		simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu sigma sigma 1e-7 0
		simHM_ground_quadBCu_ts_10_t_1000000.000000.vtu simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu DarcyVelocity DarcyVelocity 1e-7 0
    )
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
if (OGS_USE_PYTHON)
    # same, but using Python BCs
    AddTest(
        NAME HydroMechanics_HML_square_1e2_unconfined_compression_early_python
        PATH HydroMechanics/Linear/Unconfined_Compression_early
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e2_UC_early_python.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        expected_square_1e2_UC_early_ts_10_t_1.000000.vtu square_1e2_UC_early_python_ts_10_t_1.000000.vtu displacement displacement 1e-11 1e-16
        expected_square_1e2_UC_early_ts_10_t_1.000000.vtu square_1e2_UC_early_python_ts_10_t_1.000000.vtu pressure pressure 1e-10 1e-16
    )
endif()

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
    RUNTIME 7
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
if (OGS_USE_PYTHON)
    # same, but using Python BCs
    AddTest(
        NAME HydroMechanics_IdealGas_flow_pressure_boundary_python
        PATH HydroMechanics/IdealGas/flow_pressure_boundary
        EXECUTABLE ogs
        EXECUTABLE_ARGS flow_pressure_boundary_python.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        flow_pressure_boundary_ts_100_t_4000.000000.vtu flow_pressure_boundary_python_ts_100_t_4000.000000.vtu pressure pressure 0 1e-13
        flow_pressure_boundary_ts_100_t_4000.000000.vtu flow_pressure_boundary_python_ts_100_t_4000.000000.vtu displacement displacement 1e-12 0
    )
endif()

# Permeability models
AddTest(
    NAME HydroMechanics_Permeability_EmbeddedFracture_square
    PATH HydroMechanics/EmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS square.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_ts_1_t_1.000000.vtu square_ts_1_t_1.000000.vtu pressure pressure 0 1e-14
    square_ts_1_t_1.000000.vtu square_ts_1_t_1.000000.vtu velocity velocity 1e-15 0
    square_ts_1_t_1.000000.vtu square_ts_1_t_1.000000.vtu displacement displacement 1e-15 0
)

AddTest(
    NAME HydroMechanics_Permeability_EmbeddedFracture_cube
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
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_x_strain_y_flow
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS x_strain_y_flow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    x_strain_y_flow_ts_10_t_1000000000.000000.vtu x_strain_y_flow_ts_10_t_1000000000.000000.vtu pressure pressure 0 1e-14
    x_strain_y_flow_ts_10_t_1000000000.000000.vtu x_strain_y_flow_ts_10_t_1000000000.000000.vtu permeability permeability 0 1e-14
    x_strain_y_flow_ts_10_t_1000000000.000000.vtu x_strain_y_flow_ts_10_t_1000000000.000000.vtu displacement displacement 1e-15 0
    x_strain_y_flow_ts_10_t_1000000000.000000.vtu x_strain_y_flow_ts_10_t_1000000000.000000.vtu epsilon epsilon 1e-16 0
)

AddTest(
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_y_strain_z_flow
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS y_strain_z_flow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    y_strain_z_flow_ts_10_t_1000000000.000000.vtu y_strain_z_flow_ts_10_t_1000000000.000000.vtu pressure pressure 0 1e-14
    y_strain_z_flow_ts_10_t_1000000000.000000.vtu y_strain_z_flow_ts_10_t_1000000000.000000.vtu permeability permeability 0 1e-14
    y_strain_z_flow_ts_10_t_1000000000.000000.vtu y_strain_z_flow_ts_10_t_1000000000.000000.vtu displacement displacement 1e-15 0
    y_strain_z_flow_ts_10_t_1000000000.000000.vtu y_strain_z_flow_ts_10_t_1000000000.000000.vtu epsilon epsilon 1e-16 0
)

AddTest(
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_z_strain_x_flow
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS z_strain_x_flow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    z_strain_x_flow_ts_10_t_1000000000.000000.vtu z_strain_x_flow_ts_10_t_1000000000.000000.vtu pressure pressure 0 1e-14
    z_strain_x_flow_ts_10_t_1000000000.000000.vtu z_strain_x_flow_ts_10_t_1000000000.000000.vtu permeability permeability 0 1e-14
    z_strain_x_flow_ts_10_t_1000000000.000000.vtu z_strain_x_flow_ts_10_t_1000000000.000000.vtu displacement displacement 1e-15 0
    z_strain_x_flow_ts_10_t_1000000000.000000.vtu z_strain_x_flow_ts_10_t_1000000000.000000.vtu epsilon epsilon 1e-16 0
)

AddTest(
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_disc_anisotropic
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole_anisotropic.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu pressure pressure 0 1e-13
    disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu permeability permeability 0 5e-12
    disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu velocity velocity 0 1e-10
    disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_ts_1_t_1000000000.000000.vtu displacement displacement 1e-15 0
)

AddTest(
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_disc_anisotropic_rotated
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole_anisotropic_rotated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu pressure pressure 0 1e-13
    disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu permeability permeability 0 5e-12
    disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu velocity velocity 0 1e-9
    disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu disc_with_hole_anisotropic_rotated_ts_1_t_1000000000.000000.vtu displacement displacement 1e-15 0
)

AddTest(
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_disc_quasiisotropic
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole_quasiisotropic.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu pressure pressure 0 1e-13
    disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu permeability permeability 0 1e-11
    disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu velocity velocity 0 1e-10
    disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu disc_with_hole_quasiisotropic_ts_1_t_1000000000.000000.vtu displacement displacement 1e-15 0
)

AddTest(
    NAME HydroMechanics_Permeability_OrthotropicEmbeddedFracture_unconfined_biot
    PATH HydroMechanics/OrthotropicEmbeddedFracturePermeability
    EXECUTABLE ogs
    EXECUTABLE_ARGS unconfined_biot.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    unconfined_biot_ts_10_t_10000000.000000.vtu unconfined_biot_ts_10_t_10000000.000000.vtu pressure pressure 0 1e-14
    unconfined_biot_ts_10_t_10000000.000000.vtu unconfined_biot_ts_10_t_10000000.000000.vtu permeability permeability 0 1e-14
    unconfined_biot_ts_10_t_10000000.000000.vtu unconfined_biot_ts_10_t_10000000.000000.vtu displacement displacement 1e-15 0
    unconfined_biot_ts_10_t_10000000.000000.vtu unconfined_biot_ts_10_t_10000000.000000.vtu epsilon epsilon 1e-16 0
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
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu displacement displacement 1e-14 1.e-9
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu pressure pressure 1e-7 1.e-9
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu velocity velocity 1e-19 1.e-9
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu epsilon epsilon 1e-14 1.e-9
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_ts_1_t_86400.000000.vtu sigma sigma 1e-7 1.e-9
)
if (OGS_USE_PYTHON)
    # same but using Python BCs
    # TODO are output file name the same as above? Will they be overwritten during the CI run?
    AddTest(
        NAME MonolithicInjectionProduction1DPython
        PATH HydroMechanics/StaggeredScheme/InjectionProduction1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS InjectionProduction1DMono_python.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_Python_ts_1_t_86400.000000.vtu displacement displacement 1e-14 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_Python_ts_1_t_86400.000000.vtu pressure pressure 1e-7 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_Python_ts_1_t_86400.000000.vtu velocity velocity 1e-19 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_Python_ts_1_t_86400.000000.vtu epsilon epsilon 1e-14 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Mono_Python_ts_1_t_86400.000000.vtu sigma sigma 1e-7 0
    )
endif()


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
if (OGS_USE_PYTHON)
    # same but using Python BCs
    # TODO are output file name the same as above? Will they be overwritten during the CI run?
    AddTest(
        NAME StaggeredInjectionProduction1DPython
        PATH HydroMechanics/StaggeredScheme/InjectionProduction1D
        EXECUTABLE ogs
        EXECUTABLE_ARGS InjectionProduction1D_python.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        DIFF_DATA
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Python_ts_1_t_86400.000000.vtu displacement displacement 1e-13 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Python_ts_1_t_86400.000000.vtu pressure pressure 1e-6 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Python_ts_1_t_86400.000000.vtu velocity velocity 1e-18 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Python_ts_1_t_86400.000000.vtu epsilon epsilon 1e-13 0
        InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_Python_ts_1_t_86400.000000.vtu sigma sigma 1e-6 0
    )
endif()

AddTest(
    NAME StaggeredInjectionProduction1D_MFront
    PATH HydroMechanics/StaggeredScheme/InjectionProduction1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS InjectionProduction1D_MFront.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI AND OGS_USE_MFRONT
    DIFF_DATA
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu displacement displacement 1e-13 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu pressure pressure 1e-6 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu velocity velocity 1e-18 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu epsilon epsilon 1e-13 0
    InjectionProduction1D_Reference_ts_1_t_86400.000000.vtu InjectionProduction1D_ts_1_t_86400.000000.vtu sigma sigma 1e-6 0
)

AddTest(
    NAME Staggered_MandelCryer
    PATH HydroMechanics/StaggeredScheme/MandelCryer
    EXECUTABLE ogs
    EXECUTABLE_ARGS MandelCryerStaggered.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 170
    DIFF_DATA
    reference_results_MandelCryerStaggered_ts_1_t_0.010000.vtu results_MandelCryerStaggered_ts_1_t_0.010000.vtu displacement displacement 1e-4 0
    reference_results_MandelCryerStaggered_ts_1_t_0.010000.vtu results_MandelCryerStaggered_ts_1_t_0.010000.vtu pressure pressure 1e-1 0
    reference_results_MandelCryerStaggered_ts_1_t_0.010000.vtu results_MandelCryerStaggered_ts_1_t_0.010000.vtu velocity velocity 1e-1 0
    reference_results_MandelCryerStaggered_ts_1_t_0.010000.vtu results_MandelCryerStaggered_ts_1_t_0.010000.vtu epsilon epsilon 1e-4 0
    reference_results_MandelCryerStaggered_ts_1_t_0.010000.vtu results_MandelCryerStaggered_ts_1_t_0.010000.vtu sigma sigma 1e-1 0
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

AddTest(
    NAME HydroMechanics_nodal_source_test
    PATH HydroMechanics/NodalSourceTerm
    EXECUTABLE ogs
    EXECUTABLE_ARGS nodal_source_test.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    HM_NodalSourceTem_ts_100_t_86400.000000.vtu HM_NodalSourceTem_ts_100_t_86400.000000.vtu pressure pressure 1.0e-8 0.0
    HM_NodalSourceTem_ts_100_t_86400.000000.vtu HM_NodalSourceTem_ts_100_t_86400.000000.vtu displacement displacement 1.0e-15 0.0
    HM_NodalSourceTem_ts_100_t_86400.000000.vtu HM_NodalSourceTem_ts_100_t_86400.000000.vtu sigma sigma 1.5e-8 0.0
    HM_NodalSourceTem_ts_100_t_86400.000000.vtu HM_NodalSourceTem_ts_100_t_86400.000000.vtu epsilon epsilon 1.0e-15 0.0
)

# Rotated 2D mesh
AddTest(
    NAME HydroMechanics_HML_flow_gravity_rotated_2D_mesh
    PATH HydroMechanics/Linear/Gravity/RotatedAroundVerticalAxis
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_gravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flow_gravity_ts_16_t_40000000.000000.vtu flow_gravity_ts_16_t_40000000.000000.vtu displacement displacement 1e-9 0
    flow_gravity_ts_16_t_40000000.000000.vtu flow_gravity_ts_16_t_40000000.000000.vtu pressure pressure 1e-9 0
    flow_gravity_ts_16_t_40000000.000000.vtu flow_gravity_ts_16_t_40000000.000000.vtu velocity velocity 1e-9 0
)

#Parallel computing with PETSc enabled.
AddTest(
    NAME ParallelFEM_SimpleHM_SingleMesh
    PATH HydroMechanics/ParallelComputing/SimpleHM/SingleMesh
    EXECUTABLE ogs
    EXECUTABLE_ARGS drainage.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu pressure pressure 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu pressure pressure 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu q q 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu q q 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu displacement displacement 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu displacement displacement 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu sigma sigma 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu sigma sigma 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu epsilon epsilon 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu epsilon epsilon 1e-10 1e-9
)
AddTest(
    NAME ParallelFEM_SimpleHM_MultiMesh
    PATH HydroMechanics/ParallelComputing/SimpleHM/MultiMesh
    EXECUTABLE ogs
    EXECUTABLE_ARGS drainage.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu pressure pressure 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu pressure pressure 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu q q 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu q q 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu displacement displacement 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu displacement displacement 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu sigma sigma 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu sigma sigma 1e-10 1e-9
    drainage_ts_10_t_10_000000_0.vtu drainage_ts_10_t_10_000000_0.vtu epsilon epsilon 1e-10 1e-9
    drainage_ts_10_t_10_000000_1.vtu drainage_ts_10_t_10_000000_1.vtu epsilon epsilon 1e-10 1e-9
)

AddTest(
    NAME ParallelFEM_SimpleHM_SingleMesh_StaggeredScheme
    PATH HydroMechanics/ParallelComputing/SimpleHM/SingleMesh
    EXECUTABLE ogs
    EXECUTABLE_ARGS drainage_staggered.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu pressure pressure 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu pressure pressure 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu q q 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu q q 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu displacement displacement 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu displacement displacement 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu sigma sigma 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu sigma sigma 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu epsilon epsilon 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu epsilon epsilon 1e-10 1e-9
)

AddTest(
    NAME ParallelFEM_SimpleHM_MultiMesh_StaggeredScheme
    PATH HydroMechanics/ParallelComputing/SimpleHM/MultiMesh
    EXECUTABLE ogs
    EXECUTABLE_ARGS drainage_staggered.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu pressure pressure 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu pressure pressure 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu q q 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu q q 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu displacement displacement 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu displacement displacement 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu sigma sigma 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu sigma sigma 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_0.vtu drainage_staggered_ts_10_t_10_000000_0.vtu epsilon epsilon 1e-10 1e-9
    drainage_staggered_ts_10_t_10_000000_1.vtu drainage_staggered_ts_10_t_10_000000_1.vtu epsilon epsilon 1e-10 1e-9
)
