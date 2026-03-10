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
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/StaggeredSchemeClonedFromHydroMechanics/ConsolidationBenchmark/consolidation_benchmark_THM_staggered.xml
        EXECUTABLE_ARGS -m ${Data_SOURCE_DIR}/HydroMechanics/StaggeredScheme/ConsolidationBenchmark
        RUNTIME 2
    )
    OgsTest(PROJECTFILE ThermoHydroMechanics/StaggeredSchemeClonedFromHydroMechanics/InjectionProduction1D/InjectionProduction1D_THM_staggered.xml RUNTIME 1)
    OgsTest(PROJECTFILE ThermoHydroMechanics/StaggeredSchemeClonedFromHydroMechanics/MandelCryer/MandelCryerStaggered_THM_staggered.xml RUNTIME 18)
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/StaggeredSchemeClonedFromHT/ADecovalexTHMCBasedHTExample/th_decovalex_THM.xml
        EXECUTABLE_ARGS -m ${Data_SOURCE_DIR}/Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample
        RUNTIME 5
    )
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/StaggeredSchemeClonedFromHT/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow_THM.xml
        EXECUTABLE_ARGS -m ${Data_SOURCE_DIR}/Parabolic/HT/StaggeredCoupling/HeatTransportInStationaryFlow
        RUNTIME 1
    )

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
    NotebookTest(NOTEBOOKFILE ThermoHydroMechanics/MontTerriDeepBorehole/MontTerriDeepBorehole.py RUNTIME 20)
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
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity/cube_ortho_phi0_0.xml RUNTIME 1)
endif()
# ThermoHydroMechanics; Small deformation, linear elastic, porosity!=0, anisotropic thermal expansion
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity/cube_ortho_phi0_183.xml RUNTIME 1)
endif()
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity/cube_ortho_phi0_183_petsc.xml RUNTIME 5)
endif()
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity/square_ortho_phi0_0.xml RUNTIME 1)
endif()
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/anisotropic_thermal_expansivity/square_ortho_phi0_183.xml RUNTIME 1)
endif()
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Storage/cube_incompressible_fluid.prj RUNTIME 1)
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/Storage/cube_isochoric_heat-up.prj RUNTIME 1)
    OgsTest(
        PROJECTFILE ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj
        RUNTIME 1
    )
    # Rejected coupling scheme options must abort with the fatal message.
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/RejectedCouplingSchemeOptions/fixed_stress_over_time_step.xml
        EXECUTABLE_ARGS
            -m
            ${Data_SOURCE_DIR}/ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow
        NO_TEST_DEFINITION
        PROPERTIES
            PASS_REGULAR_EXPRESSION
            "fixed_stress_over_time_step is not supported by the ThermoHydroMechanics process"
    )
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/RejectedCouplingSchemeOptions/missing_coupling_scheme_type.xml
        EXECUTABLE_ARGS
            -m
            ${Data_SOURCE_DIR}/ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow
        NO_TEST_DEFINITION
        PROPERTIES
            PASS_REGULAR_EXPRESSION
            "coupling scheme configuration is missing its <type> element"
    )
    # Out-of-range stabilization parameters must warn about the range.
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/RejectedCouplingSchemeOptions/out_of_range_fixed_stress_stabilization_parameter.xml
        EXECUTABLE_ARGS
            -m
            ${Data_SOURCE_DIR}/ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow
        NO_TEST_DEFINITION
        PROPERTIES
            PASS_REGULAR_EXPRESSION
            "is out of reasonable range"
    )
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/RejectedCouplingSchemeOptions/unknown_coupling_scheme_type.xml
        EXECUTABLE_ARGS
            -m
            ${Data_SOURCE_DIR}/ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow
        NO_TEST_DEFINITION
        PROPERTIES PASS_REGULAR_EXPRESSION "Unknown coupling scheme type"
    )
    OgsTest(PROJECTFILE ThermoHydroMechanics/BGRaCreepAndInitialStressAtIP_AREHS/arehs-salt-THM01_0.prj RUNTIME 29)
    OgsTest(
        PROJECTFILE
            ThermoHydroMechanics/BGRaCreepAndInitialStressAtIP_AREHS/arehs-salt-THM01_0_staggered.xml
        RUNTIME 95
    )
    # ThermoHydroMechanics; thermo_osmosis and thermo_filtration effects, linear poroelastic, column consolidation
    OgsTest(PROJECTFILE ThermoHydroMechanics/Linear/ThermoOsmosis/Column.prj RUNTIME 9)
    OgsTest(
        PROJECTFILE ThermoHydroMechanics/Linear/ThermoOsmosis/Column.prj
        PATCH_FILES Column_with_eT.xml
        NAME_SUFFIX with_eT
        RUNTIME 9
    )
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

if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE ThermoHydroMechanics/HeatingHomogeneousDomain/partitioned_mesh/hex_THM.prj
        WRAPPER mpirun -np 3
        RUNTIME 5
        NAME_SUFFIX ParallelFEM_HeatingHomogeneousDomain
    )
endif()

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoHydroMechanics/TotalInitialStress/total_initial_stress_HM.prj RUNTIME 1)
endif()
