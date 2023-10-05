# ThermoHydroMechanics; Small deformation, linear poroelastic, homogeneous
if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/verification/thm2_1Dfixd/thm2_1Dfixd.prj RUNTIME 60)
    OgsTest(PROJECTFILE ThermoHydroMechanics/A2/A2.prj RUNTIME 23)
    OgsTest(PROJECTFILE ThermoHydroMechanics/A2/A2_heating.prj RUNTIME 23)
    OgsTest(PROJECTFILE ThermoHydroMechanics/1D_freezing_column_Stefan/Stefan_problem.prj RUNTIME 15)
    OgsTest(PROJECTFILE ThermoHydroMechanics/ColumnDeformationFreezing/TM.prj RUNTIME 13)
    OgsTest(PROJECTFILE ThermoHydroMechanics/HeatingHomogeneousDomain/hex_THM.prj RUNTIME 30)
    OgsTest(PROJECTFILE ThermoHydroMechanics/9percentWaterFreezingExpansion/UnitSquare.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow/WithFreezingPhase.prj RUNTIME 1)
endif()

AddTest(
    NAME ThermoHydroMechanics_square_1e0
    PATH ThermoHydroMechanics/Linear/Square_sealed_homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e0_ts_10_t_1000.000000.vtu square_1e0_ts_10_t_1000.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e0_ts_10_t_1000.000000.vtu square_1e0_ts_10_t_1000.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e0_ts_10_t_1000.000000.vtu square_1e0_ts_10_t_1000.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e0_ts_10_t_1000.000000.vtu square_1e0_ts_10_t_1000.000000.vtu epsilon epsilon 1e-8 1e-8
    expected_square_1e0_ts_10_t_1000.000000.vtu square_1e0_ts_10_t_1000.000000.vtu sigma sigma 1e-8 1e-8
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
    expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_ts_10_t_100.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_ts_10_t_100.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_ts_10_t_100.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_ts_10_t_100.000000.vtu epsilon epsilon 1e-8 1e-8
    expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_ts_10_t_100.000000.vtu sigma sigma 1e-8 1e-8
)

# Same as above, but with function instead of group based parameter for Young's modulus
AddTest(
NAME ThermoHydroMechanics_square_1e2_sealed_bimaterial_function
PATH ThermoHydroMechanics/Linear/Beam_sealed_bimaterial
RUNTIME 5
EXECUTABLE ogs
EXECUTABLE_ARGS square_1e2_function.xml
WRAPPER time
TESTER vtkdiff
REQUIREMENTS NOT OGS_USE_MPI
DIFF_DATA
expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_function_ts_10_t_100.000000.vtu displacement displacement 1e-8 1e-8
expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_function_ts_10_t_100.000000.vtu pressure pressure 1e-8 1e-8
expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_function_ts_10_t_100.000000.vtu temperature temperature 1e-8 1e-8
expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_function_ts_10_t_100.000000.vtu epsilon epsilon 1e-8 1e-8
expected_square_1e2_ts_10_t_100.000000.vtu square_1e2_function_ts_10_t_100.000000.vtu sigma sigma 1e-8 1e-8
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
    expected_square_1e2_ts_10_t_1000.000000.vtu square_1e2_ts_10_t_1000.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e2_ts_10_t_1000.000000.vtu square_1e2_ts_10_t_1000.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e2_ts_10_t_1000.000000.vtu square_1e2_ts_10_t_1000.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e2_ts_10_t_1000.000000.vtu square_1e2_ts_10_t_1000.000000.vtu epsilon epsilon 1e-8 1e-8
    expected_square_1e2_ts_10_t_1000.000000.vtu square_1e2_ts_10_t_1000.000000.vtu sigma sigma 1e-8 1e-8
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation
AddTest(
    NAME ThermoHydroMechanics_point_heat_injection
    PATH ThermoHydroMechanics/Linear/Point_injection
    RUNTIME 45
    EXECUTABLE ogs
    EXECUTABLE_ARGS pointheatsource_quadratic-mesh.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu displacement displacement 1e-5 1e-5
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu pressure pressure 1e-5 1e-5
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu temperature temperature 1e-5 1e-5
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu epsilon epsilon 1e-5 1e-5
    expected_pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu pointheatsource_quadratic-mesh_ts_10_t_50000.000000.vtu sigma sigma 1e-5 1e-5
)
# ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation, linear elements for displacement
AddTest(
    NAME ThermoHydroMechanics_point_heat_injection_lin
    PATH ThermoHydroMechanics/Linear/Point_injection
    RUNTIME 15
    EXECUTABLE ogs
    EXECUTABLE_ARGS pointheatsource_linear-mesh.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu displacement displacement 1e-5 1e-5
    expected_pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu pressure pressure 1e-5 1e-5
    expected_pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu temperature temperature 1e-5 1e-5
    expected_pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu epsilon epsilon 1e-5 1e-5
    expected_pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu pointheatsource_linear-mesh_ts_10_t_50000.000000.vtu sigma sigma 1e-5 1e-5
)
# ThermoHydroMechanics; Small deformation, linear elastic, porosity=0, anisotropic thermal expansion
AddTest(
    NAME ThermoHydroMechanics_cube_ortho-thermal-expansion-phi0
    PATH ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_ortho_phi0.0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
    expected_cube_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
    expected_cube_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
    expected_cube_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
    expected_cube_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
# ThermoHydroMechanics; Small deformation, linear elastic, porosity!=0, anisotropic thermal expansion
AddTest(
    NAME ThermoHydroMechanics_cube_ortho-thermal-expansion
    PATH ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_ortho_phi0.183.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
AddTest(
    NAME ThermoHydroMechanics_cube_ortho-thermal-expansion_petsc
    PATH ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_ortho_phi0.183_petsc.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
    expected_cube_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
AddTest(
     NAME ThermoHydroMechanics_square_ortho-thermal-expansion-phi0
     PATH ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity
     RUNTIME 5
     EXECUTABLE ogs
     EXECUTABLE_ARGS square_ortho_phi0.0.prj
     WRAPPER time
     TESTER vtkdiff
     REQUIREMENTS NOT OGS_USE_MPI
     DIFF_DATA
     expected_square_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
     expected_square_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
     expected_square_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
     expected_square_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
     expected_square_ortho-thermal-expansion_phi0ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
AddTest(
     NAME ThermoHydroMechanics_square_ortho-thermal-expansion
     PATH ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity
     RUNTIME 5
     EXECUTABLE ogs
     EXECUTABLE_ARGS square_ortho_phi0.183.prj
     WRAPPER time
     TESTER vtkdiff
     REQUIREMENTS NOT OGS_USE_MPI
     DIFF_DATA
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
AddTest(
     NAME ThermoHydroMechanics_cube_storage_incompressible_fluid
     PATH ThermoHydroMechanics/Linear/Storage
     RUNTIME 5
     EXECUTABLE ogs
     EXECUTABLE_ARGS cube_incompressible_fluid.prj
     WRAPPER time
     TESTER vtkdiff
     REQUIREMENTS NOT OGS_USE_MPI
     DIFF_DATA
     expected_THM_incompressible_fluidts_10_t_1.000000.vtu THM_incompressible_fluid_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
     expected_THM_incompressible_fluidts_10_t_1.000000.vtu THM_incompressible_fluid_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
     expected_THM_incompressible_fluidts_10_t_1.000000.vtu THM_incompressible_fluid_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
     expected_THM_incompressible_fluidts_10_t_1.000000.vtu THM_incompressible_fluid_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
     expected_THM_incompressible_fluidts_10_t_1.000000.vtu THM_incompressible_fluid_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
AddTest(
     NAME ThermoHydroMechanics_cube_storage_isochoric_heat-up
     PATH ThermoHydroMechanics/Linear/Storage
     RUNTIME 5
     EXECUTABLE ogs
     EXECUTABLE_ARGS cube_isochoric_heat-up.prj
     WRAPPER time
     TESTER vtkdiff
     REQUIREMENTS NOT OGS_USE_MPI
     DIFF_DATA
     expected_THM_isochoric_heat-upts_10_t_1.000000.vtu THM_isochoric_heat-up_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
     expected_THM_isochoric_heat-upts_10_t_1.000000.vtu THM_isochoric_heat-up_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
     expected_THM_isochoric_heat-upts_10_t_1.000000.vtu THM_isochoric_heat-up_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
     expected_THM_isochoric_heat-upts_10_t_1.000000.vtu THM_isochoric_heat-up_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
     expected_THM_isochoric_heat-upts_10_t_1.000000.vtu THM_isochoric_heat-up_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
AddTest(
    NAME ThermoHydroMechanics_HeatTransportInStationaryFlow
    PATH ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS HeatTransportInStationaryFlow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature_interpolated 5e-3 1e-10
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure_interpolated 1e-10 1e-10
)
AddTest(
    NAME ThermoHydroMechanics_BGRaCreepAndInitialStressAtIP_AREHS
    PATH ThermoHydroMechanics/BGRaCreepAndInitialStressAtIP_AREHS
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS arehs-salt-THM01_0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    arehs-salt-THM0_ts_5_t_157680000000.000000.vtu arehs-salt-THM0_ts_5_t_157680000000.000000.vtu displacement displacement 1e-8 1e-8
    arehs-salt-THM0_ts_5_t_157680000000.000000.vtu arehs-salt-THM0_ts_5_t_157680000000.000000.vtu pressure pressure 1e-8 1e-8
    arehs-salt-THM0_ts_5_t_157680000000.000000.vtu arehs-salt-THM0_ts_5_t_157680000000.000000.vtu temperature temperature 1e-8 1e-8
    arehs-salt-THM0_ts_5_t_157680000000.000000.vtu arehs-salt-THM0_ts_5_t_157680000000.000000.vtu epsilon epsilon 1e-8 1e-8
    arehs-salt-THM0_ts_5_t_157680000000.000000.vtu arehs-salt-THM0_ts_5_t_157680000000.000000.vtu sigma sigma 1.0e-5 1e-6
)
# ThermoHydroMechanics; thermo_osmosis and thermo_filtration effects, linear poroelastic, column consolidation
AddTest(
    NAME ThermoHydroMechanics_thermo_osmosis_filtration_effects_Column
    PATH ThermoHydroMechanics/Linear/ThermoOsmosis
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

AddTest(
    NAME THM_HT_ClassicalTransportExample
    PATH ThermoHydroMechanics/Linear/TH_ClassicalTransportExample
    EXECUTABLE ogs
    EXECUTABLE_ARGS classical_transport_example.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    classical_transport_example_t_4800.00.vtu classical_transport_example_t_4800.00.vtu C C 1.e-9 1.0e-12
    classical_transport_example_t_4800.00.vtu classical_transport_example_t_4800.00.vtu pressure pressure 1.e-9 1.0e-12
    classical_transport_example_t_4800.00.vtu classical_transport_example_t_4800.00.vtu velocity velocity 1.e-12 1.0e-12
    classical_transport_example_t_7200.00.vtu classical_transport_example_t_7200.00.vtu C C 1.e-9 1.0e-12
    classical_transport_example_t_7200.00.vtu classical_transport_example_t_7200.00.vtu pressure pressure 1.e-9 1.0e-12
    classical_transport_example_t_7200.00.vtu classical_transport_example_t_7200.00.vtu velocity velocity 1.e-12 1.0e-12
)

AddTest(
    NAME THM_HT_ClassicalTransportExample_full_upwind
    PATH ThermoHydroMechanics/Linear/TH_ClassicalTransportExample
    EXECUTABLE ogs
    EXECUTABLE_ARGS classical_transport_example_full_upwind.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    classical_transport_example_full_upwind_t_4800.00.vtu classical_transport_example_full_upwind_t_4800.00.vtu C C 1.e-9 1.0e-12
    classical_transport_example_full_upwind_t_4800.00.vtu classical_transport_example_full_upwind_t_4800.00.vtu pressure pressure 1.e-9 1.0e-12
    classical_transport_example_full_upwind_t_4800.00.vtu classical_transport_example_full_upwind_t_4800.00.vtu velocity velocity 1.e-12 1.0e-12
    classical_transport_example_full_upwind_t_7200.00.vtu classical_transport_example_full_upwind_t_7200.00.vtu C C 1.e-9 1.0e-12
    classical_transport_example_full_upwind_t_7200.00.vtu classical_transport_example_full_upwind_t_7200.00.vtu pressure pressure 1.e-9 1.0e-12
    classical_transport_example_full_upwind_t_7200.00.vtu classical_transport_example_full_upwind_t_7200.00.vtu velocity velocity 1.e-12 1.0e-12
)

AddTest(
    NAME ParallelFEM_HeatingHomogeneousDomain
    PATH ThermoHydroMechanics/HeatingHomogeneousDomain
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS hex_THM.prj -m partitioned_mesh
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
#
    hex_ts_5_t_864000_000000_0.vtu hex_ts_5_t_864000_000000_0.vtu displacement displacement 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_0.vtu hex_ts_5_t_864000_000000_0.vtu pressure pressure 1e-10 1.0e-7
    hex_ts_5_t_864000_000000_0.vtu hex_ts_5_t_864000_000000_0.vtu temperature temperature 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_0.vtu hex_ts_5_t_864000_000000_0.vtu epsilon epsilon 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_0.vtu hex_ts_5_t_864000_000000_0.vtu sigma sigma 1e-8 1.0e-7
#
    hex_ts_5_t_864000_000000_1.vtu hex_ts_5_t_864000_000000_1.vtu displacement displacement 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_1.vtu hex_ts_5_t_864000_000000_1.vtu pressure pressure 1e-10 1.0e-7
    hex_ts_5_t_864000_000000_1.vtu hex_ts_5_t_864000_000000_1.vtu temperature temperature 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_1.vtu hex_ts_5_t_864000_000000_1.vtu epsilon epsilon 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_1.vtu hex_ts_5_t_864000_000000_1.vtu sigma sigma 1e-8 1.0e-7
#
    hex_ts_5_t_864000_000000_2.vtu hex_ts_5_t_864000_000000_2.vtu displacement displacement 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_2.vtu hex_ts_5_t_864000_000000_2.vtu pressure pressure 1e-10 1.0e-7
    hex_ts_5_t_864000_000000_2.vtu hex_ts_5_t_864000_000000_2.vtu temperature temperature 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_2.vtu hex_ts_5_t_864000_000000_2.vtu epsilon epsilon 1e-10 1.0e-9
    hex_ts_5_t_864000_000000_2.vtu hex_ts_5_t_864000_000000_2.vtu sigma sigma 1e-8 1.0e-7
)
