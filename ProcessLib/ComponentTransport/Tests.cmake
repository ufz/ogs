if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/InclinedElements/Inclined2DMesh/inclined_2D_mesh_HC.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ClassicalTransportExample/classical_transport_example_full_upwind.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/Theis_Axisymmetric/axisym_theis_CT.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/HTCWithFracture/2D_single_fracture_HTC.prj RUNTIME 29)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/HTCWithFracture/2D_single_fracture_HTC_Monolithic.prj RUNTIME 1
        PROPERTIES WILL_FAIL true)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/ConcentrationDiffusionOnly.prj)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/ConcentrationDiffusionOnly.prj)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/ConcentrationDiffusionOnly_3Components.prj)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/ConcentrationDiffusionAndStorage.prj)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/ConcentrationDiffusionAndStorage.prj RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/DiffusionAndStorageAndAdvection.prj RUNTIME 10)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/MassConservation/mass_conservation.prj RUNTIME 18)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/OpenBoundaryWithTets/box_flow.prj RUNTIME 4)
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/DiffusionAndStorageAndGravityAndDispersionHalf.prj
        RUNTIME 17
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/DiffusionAndStorageAndAdvectionAndDispersion.prj
        RUNTIME 10
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/open_boundary_component-transport_cube_1e3.prj
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/open_boundary_component-transport_cube_1e3_advective_form.prj
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/DiffusionAndStorageAndAdvectionAndDecay.prj
        RUNTIME 9
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/DiffusionAndStorageAndAdvectionAndDispersionHalf.prj
        RUNTIME 9
    )
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/SimpleSynthetics/surfaceflux_component-transport_cube_1e3.prj)
endif()

AddTest(
    NAME 3D_ComponentTransport_surfaceflux_pvd
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS surfaceflux_component-transport_cube_1e3.prj
    WRAPPER time
    TESTER diff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    cube_1x1x1_hex_1e3_right.pvd
    cube_1x1x1_hex_1e3_left.pvd
    cube_1x1x1_hex_1e3.pvd
    cube_1x1x1_hex_1e3_complete_surface.pvd
)

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/surfaceflux_component-transport_cube_1e3.prj)
endif()

AddTest(
    NAME 2D_StaggeredScheme_ComponentTransport_TracerSimulation
    PATH Parabolic/ComponentTransport/TracerSimulation
    EXECUTABLE ogs
    RUNTIME 3
    EXECUTABLE_ARGS TracerSimulation.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    TracerSimulation_ts_20_t_20000.000000_expected_ogs5.vtu TracerSimulation_ts_20_t_20000.000000.vtu Cs Cs 4e-7 1e-10
    TracerSimulation_ts_40_t_40000.000000_expected_ogs5.vtu TracerSimulation_ts_40_t_40000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_ts_60_t_60000.000000_expected_ogs5.vtu TracerSimulation_ts_60_t_60000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_ts_80_t_80000.000000_expected_ogs5.vtu TracerSimulation_ts_80_t_80000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_ts_100_t_100000.000000_expected_ogs5.vtu TracerSimulation_ts_100_t_100000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_ts_20_t_20000.000000_expected_ogs5.vtu TracerSimulation_ts_20_t_20000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_ts_40_t_40000.000000_expected_ogs5.vtu TracerSimulation_ts_40_t_40000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_ts_60_t_60000.000000_expected_ogs5.vtu TracerSimulation_ts_60_t_60000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_ts_80_t_80000.000000_expected_ogs5.vtu TracerSimulation_ts_80_t_80000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_ts_100_t_100000.000000_expected_ogs5.vtu TracerSimulation_ts_100_t_100000.000000.vtu pressure pressure 1.3e-7 1e-10
)

AddTest(
    NAME Parallel_2D_StaggeredScheme_ComponentTransport_TracerSimulation
    PATH Parabolic/ComponentTransport/TracerSimulation
    EXECUTABLE ogs
    EXECUTABLE_ARGS TracerSimulation.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    TracerSimulation_ts_100_t_100000_000000_0_expected.vtu TracerSimulation_ts_100_t_100000_000000_0.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_ts_100_t_100000_000000_0_expected.vtu TracerSimulation_ts_100_t_100000_000000_0.vtu pressure pressure 4.0e-5 2.7e-9
    TracerSimulation_ts_100_t_100000_000000_1_expected.vtu TracerSimulation_ts_100_t_100000_000000_1.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_ts_100_t_100000_000000_1_expected.vtu TracerSimulation_ts_100_t_100000_000000_1.vtu pressure pressure 4.0e-5 3.7e-9
    TracerSimulation_ts_100_t_100000_000000_2_expected.vtu TracerSimulation_ts_100_t_100000_000000_2.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_ts_100_t_100000_000000_2_expected.vtu TracerSimulation_ts_100_t_100000_000000_2.vtu pressure pressure 3.5e-5 2.5e-9
    TracerSimulation_ts_100_t_100000_000000_3_expected.vtu TracerSimulation_ts_100_t_100000_000000_3.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_ts_100_t_100000_000000_3_expected.vtu TracerSimulation_ts_100_t_100000_000000_3.vtu pressure pressure 3.6e-5 4.7e-9
)

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/goswami/goswami_input.prj RUNTIME 291)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/elder/elder.prj RUNTIME 1068)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/elder/elder-python.prj RUNTIME 1094)
    OgsTest(PROJECTFILE Elliptic/square_100x100_ComponentTransport/square_1e4_heterogeneity.prj)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/heterogeneous/ogs5_H_2D/ogs5_H_2d.prj)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/heterogeneous/ogs5_H_3D/ogs5_H_3d.prj RUNTIME 20)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/VariableNeumannBoundary/vdbc_input.prj RUNTIME 7)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/Theis/theis.prj RUNTIME 45)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ConTracer/ConTracer_2d.prj RUNTIME 12)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ConTracer/ConTracer_1d.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/DiffusionSorptionDecay/1D_DiffusionSorptionDecay.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/AdvectionDiffusionSorptionDecay/1D_AdvectionDiffusionSorptionDecay.prj RUNTIME 6)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/MultiLayerDiffusion/1D_MultiLayerDiffusion.prj RUNTIME 21)

    # variation that will fail because material ids are missing in the input
    # mesh
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/MultiLayerDiffusion/1D_MultiLayerDiffusion_fail_no_mat_ids.xml
        PROPERTIES
        PASS_REGULAR_EXPRESSION "More than one porous medium definition.*but no MaterialIDs are present in the bulk mesh"
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/failing_tests/non_existent_python_script.xml
        PROPERTIES
        PASS_REGULAR_EXPRESSION "File \".*[/\\]Tests[/\\]Data[/\\]Parabolic[/\\]ComponentTransport[/\\]ReactiveTransport[/\\]DecayChain[/\\]GlobalImplicitApproach[/\\]failing_tests[/\\]nonexistent.py\" could not be opened!"
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/failing_tests/python_import_error.xml
        PROPERTIES
        PASS_REGULAR_EXPRESSION "Error evaluating python script .*[/\\]Tests[/\\]Data[/\\]Parabolic[/\\]ComponentTransport[/\\]ReactiveTransport[/\\]DecayChain[/\\]GlobalImplicitApproach[/\\]failing_tests[/\\]import_error.py: ImportError: this script raises an error on purpose"
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/failing_tests/python_error_Dirichlet.xml
        PROPERTIES
        PASS_REGULAR_EXPRESSION "error: Error evaluating Dirichlet BC in Python: RuntimeError: this exception is thrown on purpose"
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/failing_tests/python_error_Neumann.xml
        PROPERTIES
        PASS_REGULAR_EXPRESSION "error: Error evaluating natural BC in Python: RuntimeError: this exception is thrown on purpose"
        RUNTIME 1
    )
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/failing_tests/python_error_source_term.xml
        PROPERTIES
        PASS_REGULAR_EXPRESSION "error: Error evaluating source term in Python: RuntimeError: this exception is thrown on purpose"
        RUNTIME 1
    )

    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/DiffusionAndStorageAndAdvection.prj RUNTIME 11)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/DiffusionAndStorageAndAdvectionAndDecay.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/DiffusionAndStorageAndAdvectionAndDispersionHalf.prj RUNTIME 20)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/DiffusionAndStorageAndAdvectionAndDispersion.prj RUNTIME 15)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/DiffusionAndStorageAndAdvectionAndDispersion_3Components.prj RUNTIME 29)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/StaggeredScheme/DiffusionAndStorageAndGravityAndDispersionHalf.prj RUNTIME 37)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/EquilibriumPhase/calcite.prj RUNTIME 6)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant/1d_isofrac.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant/1d_isofrac_small_domain.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant/1d_isofrac_flag_formula.prj RUNTIME 8)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant_AllAsComponents/KineticReactant2.prj RUNTIME 10)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/RadionuclidesMigration/Cs.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/RadionuclidesMigration/U.prj RUNTIME 17)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/SurfaceComplexation/RadionuclideSorption.prj RUNTIME 11)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/SurfaceComplexation/RadionuclideSorption_fixed_pe.prj RUNTIME 10)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/SurfaceComplexation/RadionuclideSorption_SitesMoles.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/SurfaceComplexation/LookupTable/RadionuclideSorption.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/EquilibriumPhase/calciteDissolvePrecipitateOnly.prj RUNTIME 5)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/CationExchange/exchange.prj RUNTIME 7)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/CationExchange/exchangeAndSurface.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/1d_decay_chain_OS.prj RUNTIME 808)

    # several variations of 1d_decay_chain_GIA
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/1d_decay_chain_GIA.prj RUNTIME 12)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/linear/1d_decay_chain_GIA.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/linear_compute_only_on_dt_change/1d_decay_chain_GIA.xml RUNTIME 4)

    # further variations of 1d_decay_chain_GIA with Eigen's SparseLU solver
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/SparseLU/1d_decay_chain_GIA.xml RUNTIME 12)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/SparseLU_linear/1d_decay_chain_GIA.xml RUNTIME 4)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/SparseLU_linear_compute_only_on_dt_change/1d_decay_chain_GIA.xml RUNTIME 4)

    # variation with changing timestep size
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/varying_dt/1d_decay_chain_GIA.xml RUNTIME 12)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/varying_dt_linear/1d_decay_chain_GIA.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/varying_dt_linear_compute_only_on_dt_change/1d_decay_chain_GIA.xml RUNTIME 4)

    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/SolidPhasePositionTest/1d_vertical_test.prj RUNTIME 1)

    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ThermalDiffusion/TemperatureField_transport.prj RUNTIME 15)
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/FCT_test/1d_step_func.prj RUNTIME 2)
endif()
if(OGS_USE_PETSC)
    OgsTest(WRAPPER mpirun -np 1 PROJECTFILE Parabolic/ComponentTransport/FCT_test/1d_step_func.prj RUNTIME 16)
    OgsTest(WRAPPER mpirun -np 1 PROJECTFILE Parabolic/ComponentTransport/HTCWithFracture/2D_single_fracture_HTC.prj RUNTIME 1400)
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS) AND NOT WIN32)
    NotebookTest(
        NOTEBOOKFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/GlobalImplicitApproach/performance_measurements.py
        RUNTIME 60
        SKIP_WEB
    )
endif()

# TODO: does not work on m1-mac (mac-bilke-623)
# see https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/5307
# and https://github.com/pmgbergen/porepy/issues/1473
if(NOT (OGS_USE_PETSC OR OGS_USE_LIS) AND
   NOT WIN32 AND
   NOT "${HOSTNAME}" MATCHES "mac-bilke-623")
    NotebookTest(
        NOTEBOOKFILE Parabolic/ComponentTransport/DFN_PorePy/DFNbyPorePy.py
        RUNTIME 22
        PYTHON_PACKAGES porepy@git+https://github.com/pmgbergen/porepy.git@v1.12
    )
    NotebookTest(
        NOTEBOOKFILE Parabolic/ComponentTransport/DFN_PorePy/DFNbyPorePy_to_OGS.py
        RUNTIME 66
        PYTHON_PACKAGES porepy@git+https://github.com/pmgbergen/porepy.git@v1.12
    )
endif()
if(NOT (OGS_USE_MPI OR OGS_USE_LIS) AND ("${HOSTNAME}" MATCHES "envinf1" OR APPLE OR MSVC))
    OgsTest(
        PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant_AllAsComponents/KineticReactant2_2d.prj
        RUNTIME 3300
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/Wetland/Wetland_1d.prj RUNTIME 15)
endif()

if(OGS_USE_MPI)
    OgsTest(WRAPPER mpirun -np 1 PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/EquilibriumPhase/calcitePorosityChange.prj RUNTIME 25)
    OgsTest(WRAPPER mpirun -np 2 PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/EquilibriumPhase/calcite_mpi.xml RUNTIME 60)
    OgsTest(WRAPPER mpirun -np 2 PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/EquilibriumPhase/calcite_mpi_chem_threads.xml RUNTIME 60)
    OgsTest(WRAPPER mpirun -np 2 PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/SurfaceComplexation/ParallelTest/RadionuclideSorption.prj RUNTIME 60)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE
            Parabolic/ComponentTransport/ClassicalTransportExample/classical_transport_example.prj
        RUNTIME 1
    )
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE Parabolic/LiquidFlow/AxiSymTheis/axisym_theis.py RUNTIME 7)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/ReactiveTransport/RadionuclidesMigration/RadionuclidesMigration.py RUNTIME 31)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/ReactiveTransport/CO2Injection/CO2Injection.py RUNTIME 5)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/MultiLayerDiffusion/MultiLayerDiffusion.py RUNTIME 25)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/DiffusionSorptionDecay/DiffusionSorptionDecay.py RUNTIME 12)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/elder_jupyter/elder_jupyter.py RUNTIME 11)
    NotebookTest(NOTEBOOKFILE Parabolic/ThermalTwoPhaseFlowPP/HeatPipe/heatpipe.py RUNTIME 10)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/ThermalDiffusion/ThermalDiffusion.py RUNTIME 33)

endif()

if(OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/ReactiveTransport/DecayChain/DecayChain.py RUNTIME 160)
    NotebookTest(NOTEBOOKFILE Parabolic/ComponentTransport/ReactiveTransport/Acidification/PorosityIncrease.py RUNTIME 250)
endif()

if(NOT OGS_USE_LIS)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/MassFlux/only_grad_c.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/MassFlux/grad_c_and_grad_p.prj RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/MassFlux/only_grad_p.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/MassFlux/grad_c_and_grad_p_and_r.prj RUNTIME 3)
endif()
