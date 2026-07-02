# ThermoHydroMechanics; Small deformation, linear poroelastic, homogeneous
if (NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/verification/thm2_1Dfixd/thm2_1Dfixd.prj RUNTIME 17)
    OgsTest(PROJECTFILE ThermoHydroMechanics/A2/A2.prj RUNTIME 7)
    OgsTest(PROJECTFILE ThermoHydroMechanics/A2/A2_heating.prj RUNTIME 6)
    OgsTest(PROJECTFILE ThermoHydroMechanics/1D_freezing_column_Stefan/Stefan_problem.prj RUNTIME 5)
    OgsTest(PROJECTFILE ThermoHydroMechanics/1D_freezing_column_Stefan/Stefan_problem_Function.xml RUNTIME 5)
    OgsTest(PROJECTFILE ThermoHydroMechanics/ColumnDeformationFreezing/TM.prj RUNTIME 13)
    OgsTest(PROJECTFILE ThermoHydroMechanics/9percentWaterFreezingExpansion/UnitSquare.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow/WithFreezingPhase.prj RUNTIME 1)

    if (OGS_USE_MFRONT)
        OgsTest(PROJECTFILE ThermoHydroMechanics/MultiMaterial/DP_Ehlers/M/square_1e1_2_matIDs.prj RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/MultiMaterial/DP_Ehlers/M/square_1e1_2_matIDs_restart.xml RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/MultiMaterial/DP_Ehlers/TM/square_1e1_2_matIDs.prj RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/MultiMaterial/DP_Ehlers/TM/square_1e1_2_matIDs_restart.xml RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/MultiMaterial/DP_MCC/TM/square_1e1_2_matIDs.prj RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/MultiMaterial/DP_MCC/TM/square_1e1_2_matIDs_restart.xml RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/RestartMCC/mfront_restart_part1.prj RUNTIME 1)
        OgsTest(PROJECTFILE ThermoHydroMechanics/RestartMCC/mfront_restart_part2.xml RUNTIME 1)
    endif()
endif()
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoHydroMechanics/HeatingHomogeneousDomain/hex_THM.prj RUNTIME 15)
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE ThermoHydroMechanics/HeatingHomogeneousDomain/heating_homogenous_vol.py RUNTIME 20)
endif()
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Square_sealed_homogeneous/square_1e0.prj)
    # ThermoHydroMechanics; Small deformation, linear poroelastic, sealed, bimaterial
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Beam_sealed_bimaterial/square_1e2.prj RUNTIME 5)
    # Same as above, but with function instead of group based parameter for Young's modulus
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Beam_sealed_bimaterial/square_1e2_function.xml RUNTIME 5)
    # ThermoHydroMechanics; Small deformation, linear poroelastic, unsealed, bimaterial
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Beam_unsealed_bimaterial/square_1e2.prj RUNTIME 2)
    # ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Point_injection/pointheatsource_quadratic-mesh.prj RUNTIME 6)
    # ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation, linear elements for displacement
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Point_injection/pointheatsource_linear-mesh.prj RUNTIME 3)
    # ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation, variation with a volumetric source term
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/Linear/Point_injection/with-volumetric-source-term/pointheatsource_quadratic-mesh.xml
        RUNTIME 6
    )
    # ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation, variation with a Python source term
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/Linear/Point_injection/with-python-source-term/pointheatsource_quadratic-mesh.xml
        RUNTIME 6
    )
endif()
AddTest(
    NAME ThermoHydroMechanics_cube_ortho-thermal-expansion-phi0
    PATH ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_ortho_phi0.0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
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
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_ortho_phi0.183.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
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
     RUNTIME 1
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
     RUNTIME 1
     EXECUTABLE ogs
     EXECUTABLE_ARGS square_ortho_phi0.183.prj
     WRAPPER time
     TESTER vtkdiff
     REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
     DIFF_DATA
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu displacement displacement 1e-8 1e-8
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu pressure pressure 1e-5 1e-5
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu temperature temperature 1e-8 1e-8
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu sigma sigma 1e-5 1e-5
     expected_square_ortho-thermal-expansion_phi0.183ts_10_t_1.000000.vtu THM_square_ortho-thermal-expansion-phi0.183_ts_10_t_1.000000.vtu epsilon epsilon 1e-8 1e-8
)
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Storage/cube_incompressible_fluid.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Storage/cube_isochoric_heat-up.prj RUNTIME 1)
    OgsTest(
        PROJECTFILE ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj
        RUNTIME 1
    )
    OgsTest(PROJECTFILE ThermoHydroMechanics/BGRaCreepAndInitialStressAtIP_AREHS/arehs-salt-THM01_0.prj RUNTIME 29)
    # ThermoHydroMechanics; thermo_osmosis and thermo_filtration effects, linear poroelastic, column consolidation
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/ThermoOsmosis/Column.prj RUNTIME 9)
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/Linear/TH_ClassicalTransportExample/classical_transport_example.prj
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/Linear/TH_ClassicalTransportExample/classical_transport_example_full_upwind.prj
        RUNTIME 1
    )
endif()

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

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoHydroMechanics/TotalInitialStress/total_initial_stress_HM.prj RUNTIME 1)
endif()
