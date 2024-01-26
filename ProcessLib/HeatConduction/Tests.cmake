if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/T/1D_freezing_column_Stefan/Stefan_problem.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/T/1D_freezing_column_Stefan/Stefan_problem_homogen.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/2D_freezing_disk/circle_disk.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/T/2D_Robin/square_1e4_robin.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/2D_Robin/square_1e4_robin_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/2D_Ice_melting-forming_manuf_solution/ManSol3_IceWaterMix_Scaled.prj RUNTIME 5)
    OgsTest(PROJECTFILE Parabolic/T/1D_Two-phase_Stefan_problem_for_ice_melting/Two-phase_Stefan_problem.prj RUNTIME 5)
    OgsTest(PROJECTFILE Parabolic/T/2D_Soil_freezing_round_BHE/m16m15projectB.prj RUNTIME 35)
endif()

AddTest(
        NAME 1D_HeatConduction_dirichlet
        PATH Parabolic/T/1D_dirichlet
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.prj
        TESTER vtkdiff
        DIFF_DATA
        # analytical solution
        temperature_analytical.vtu line_60_heat_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature 1e-5 1e-5
        temperature_analytical.vtu line_60_heat_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature 1e-5 1e-5
        # numerical solution – ts 65
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        temperature temperature 2.1e-12 0
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        heat_flux heat_flux 2.8e-12 0
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        HeatFlowRate HeatFlowRate 6.5e-12 0
        # numerical solution – ts 405
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        temperature temperature 1.2e-11 0
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        heat_flux heat_flux 6.2e-12 0
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        HeatFlowRate HeatFlowRate 5.9e-12 0
        # numerical solution – ts 500
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        temperature temperature 1.5e-11 0
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        heat_flux heat_flux 7.5e-12 0
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        HeatFlowRate HeatFlowRate 7.9e-12 0
    REQUIREMENTS NOT OGS_USE_MPI
)

# variation of the test above with assembly optimization
AddTest(
        NAME 1D_HeatConduction_dirichlet_linear
        PATH Parabolic/T/1D_dirichlet/linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.xml
        TESTER vtkdiff
        DIFF_DATA
        # numerical solution – ts 65
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        temperature temperature 2.2e-12 0
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        heat_flux heat_flux 2.1e-12 0
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        HeatFlowRate HeatFlowRate 4.4e-12 0
        # numerical solution – ts 405
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        temperature temperature 1.3e-11 0
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        heat_flux heat_flux 5.7e-12 0
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        HeatFlowRate HeatFlowRate 5.9e-12 0
        # numerical solution – ts 500
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        temperature temperature 1.5e-11 0
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        heat_flux heat_flux 7.5e-12 0
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        HeatFlowRate HeatFlowRate 9.8e-12 0
    REQUIREMENTS NOT OGS_USE_MPI
)

# another variation with more optimizations
AddTest(
        NAME 1D_HeatConduction_dirichlet_linear_dt_change
        PATH Parabolic/T/1D_dirichlet/linear_compute_only_on_dt_change
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.xml
        TESTER vtkdiff
        DIFF_DATA
        # numerical solution – ts 65
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        temperature temperature 2.2e-12 0
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        heat_flux heat_flux 2.1e-12 0
        line_60_heat_ts_65_t_5078125.000000.vtu
        line_60_heat_ts_65_t_5078125.000000.vtu
        HeatFlowRate HeatFlowRate 4.4e-12 0
        # numerical solution – ts 405
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        temperature temperature 1.3e-11 0
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        heat_flux heat_flux 5.7e-12 0
        line_60_heat_ts_405_t_31640625.000000.vtu
        line_60_heat_ts_405_t_31640625.000000.vtu
        HeatFlowRate HeatFlowRate 5.9e-12 0
        # numerical solution – ts 500
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        temperature temperature 1.5e-11 0
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        heat_flux heat_flux 7.5e-12 0
        line_60_heat_ts_500_t_39062500.000000.vtu
        line_60_heat_ts_500_t_39062500.000000.vtu
        HeatFlowRate HeatFlowRate 7.9e-12 0
    REQUIREMENTS NOT OGS_USE_MPI
)

# a variation with a varying time step size
AddTest(
        NAME 1D_HeatConduction_dirichlet_varying_dt
        PATH Parabolic/T/1D_dirichlet/varying_dt
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.xml
        TESTER vtkdiff
        DIFF_DATA
        # numerical solution – ts 130
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        temperature temperature 6.3e-13 0
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        heat_flux heat_flux 7.8e-13 0
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        HeatFlowRate HeatFlowRate 7.5e-12 0
        # numerical solution – ts 505
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        temperature temperature 9.9e-12 0
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        heat_flux heat_flux 5.7e-12 0
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        HeatFlowRate HeatFlowRate 5.7e-12 0
        # numerical solution – ts 600
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        temperature temperature 1.3e-11 0
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        heat_flux heat_flux 6.8e-12 0
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        HeatFlowRate HeatFlowRate 6.8e-12 0
    REQUIREMENTS NOT OGS_USE_MPI
)

# a variation with a varying time step size and <linear> optimization
AddTest(
        NAME 1D_HeatConduction_dirichlet_varying_dt_linear
        PATH Parabolic/T/1D_dirichlet/varying_dt_linear
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.xml
        TESTER vtkdiff
        DIFF_DATA
        # numerical solution – ts 130
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        temperature temperature 8.6e-13 0
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        heat_flux heat_flux 1.1e-12 0
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        HeatFlowRate HeatFlowRate 8.5e-12 0
        # numerical solution – ts 505
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        temperature temperature 1.1e-11 0
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        heat_flux heat_flux 5.2e-12 0
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        HeatFlowRate HeatFlowRate 5.2e-12 0
        # numerical solution – ts 600
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        temperature temperature 1.3e-11 0
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        heat_flux heat_flux 5.9e-12 0
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        HeatFlowRate HeatFlowRate 6.2e-12 0
    REQUIREMENTS NOT OGS_USE_MPI
)

# a variation with a varying time step size and more optimizations
AddTest(
        NAME 1D_HeatConduction_dirichlet_varying_dt_linear_compute_only_on_dt_change
        PATH Parabolic/T/1D_dirichlet/varying_dt_linear_compute_only_on_dt_change
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.xml
        TESTER vtkdiff
        DIFF_DATA
        # numerical solution – ts 130
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        temperature temperature 8.6e-13 0
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        heat_flux heat_flux 1.1e-12 0
        line_60_heat_ts_130_t_5078125.000000.vtu
        line_60_heat_ts_130_t_5078125.000000.vtu
        HeatFlowRate HeatFlowRate 8.5e-12 0
        # numerical solution – ts 505
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        temperature temperature 1.1e-11 0
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        heat_flux heat_flux 5.2e-12 0
        line_60_heat_ts_505_t_31640625.000000.vtu
        line_60_heat_ts_505_t_31640625.000000.vtu
        HeatFlowRate HeatFlowRate 5.2e-12 0
        # numerical solution – ts 600
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        temperature temperature 1.3e-11 0
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        heat_flux heat_flux 5.9e-12 0
        line_60_heat_ts_600_t_39062500.000000.vtu
        line_60_heat_ts_600_t_39062500.000000.vtu
        HeatFlowRate HeatFlowRate 6.2e-12 0
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME 1D_HeatConduction_neumann_picard
        PATH Parabolic/T/1D_neumann
        EXECUTABLE ogs
        EXECUTABLE_ARGS picard.prj
        TESTER vtkdiff
        DIFF_DATA
        picard_ts_1_t_78125.000000.vtu picard_ts_1_t_78125.000000.vtu temperature temperature 5e-12 1e-14
        picard_ts_3_t_234375.000000.vtu picard_ts_3_t_234375.000000.vtu temperature temperature 5e-12 1e-14
        picard_ts_65_t_5078125.000000.vtu picard_ts_65_t_5078125.000000.vtu temperature temperature 5e-12 1e-14
        picard_ts_405_t_31640625.000000.vtu picard_ts_405_t_31640625.000000.vtu temperature temperature 2e-11 1e-14
        picard_ts_500_t_39062500.000000.vtu picard_ts_500_t_39062500.000000.vtu temperature temperature 2e-11 1e-14
        picard_ts_1_t_78125.000000.vtu picard_ts_1_t_78125.000000.vtu HeatFlowRate HeatFlowRate 5e-12 0
        picard_ts_3_t_234375.000000.vtu picard_ts_3_t_234375.000000.vtu HeatFlowRate HeatFlowRate 5e-12 0
        picard_ts_65_t_5078125.000000.vtu picard_ts_65_t_5078125.000000.vtu HeatFlowRate HeatFlowRate 8e-12 0
        picard_ts_405_t_31640625.000000.vtu picard_ts_405_t_31640625.000000.vtu HeatFlowRate HeatFlowRate 5e-12 1e-16
        picard_ts_500_t_39062500.000000.vtu picard_ts_500_t_39062500.000000.vtu HeatFlowRate HeatFlowRate 6e-12 1e-16
        temperature_analytical.vtu picard_ts_1_t_78125.000000.vtu temperature_78125s temperature 8e-2 1e-4
        temperature_analytical.vtu picard_ts_3_t_234375.000000.vtu temperature_234375s temperature 6e-2 1e-4
        temperature_analytical.vtu picard_ts_65_t_5078125.000000.vtu temperature_5078125s temperature 1e-4 1e-4
        temperature_analytical.vtu picard_ts_405_t_31640625.000000.vtu temperature_31640625s temperature 1e-4 1e-4
        temperature_analytical.vtu picard_ts_500_t_39062500.000000.vtu temperature_39062500s temperature 1e-4 1e-4
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 1D_HeatConduction_neumann_newton
    PATH Parabolic/T/1D_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS newton.prj
    TESTER vtkdiff
    DIFF_DATA
    newton_ts_1_t_78125.000000.vtu newton_ts_1_t_78125.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_3_t_234375.000000.vtu newton_ts_3_t_234375.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_65_t_5078125.000000.vtu newton_ts_65_t_5078125.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_405_t_31640625.000000.vtu newton_ts_405_t_31640625.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_500_t_39062500.000000.vtu newton_ts_500_t_39062500.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_1_t_78125.000000.vtu newton_ts_1_t_78125.000000.vtu HeatFlowRate HeatFlowRate 1e-12 1e-16
    newton_ts_3_t_234375.000000.vtu newton_ts_3_t_234375.000000.vtu HeatFlowRate HeatFlowRate 3e-12 0
    newton_ts_65_t_5078125.000000.vtu newton_ts_65_t_5078125.000000.vtu HeatFlowRate HeatFlowRate 2e-12 0
    newton_ts_405_t_31640625.000000.vtu newton_ts_405_t_31640625.000000.vtu HeatFlowRate HeatFlowRate 3e-12 0
    newton_ts_500_t_39062500.000000.vtu newton_ts_500_t_39062500.000000.vtu HeatFlowRate HeatFlowRate 2e-12 0
    temperature_analytical.vtu newton_ts_1_t_78125.000000.vtu temperature_78125s temperature 8e-2 1e-4
    temperature_analytical.vtu newton_ts_3_t_234375.000000.vtu temperature_234375s temperature 6e-2 1e-4
    temperature_analytical.vtu newton_ts_65_t_5078125.000000.vtu temperature_5078125s temperature 1e-4 1e-4
    temperature_analytical.vtu newton_ts_405_t_31640625.000000.vtu temperature_31640625s temperature 1e-4 1e-4
    temperature_analytical.vtu newton_ts_500_t_39062500.000000.vtu temperature_39062500s temperature 1e-4 1e-4
    # TODO: Fix on Apple M1
    REQUIREMENTS NOT OGS_USE_MPI AND NOT ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64"
)

AddTest(
    NAME 1D_HeatConduction_neumann_picard_masslumping
    PATH Parabolic/T/1D_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS picard_masslumping.prj
    TESTER vtkdiff
    DIFF_DATA
    picard_masslumping_ts_1_t_78125.000000.vtu picard_masslumping_ts_1_t_78125.000000.vtu temperature temperature 1e-12 1e-16
    picard_masslumping_ts_3_t_234375.000000.vtu picard_masslumping_ts_3_t_234375.000000.vtu temperature temperature 1e-12 1e-16
    picard_masslumping_ts_65_t_5078125.000000.vtu picard_masslumping_ts_65_t_5078125.000000.vtu temperature temperature 1e-12 1e-16
    picard_masslumping_ts_405_t_31640625.000000.vtu picard_masslumping_ts_405_t_31640625.000000.vtu temperature temperature 1e-12 1e-16
    picard_masslumping_ts_500_t_39062500.000000.vtu picard_masslumping_ts_500_t_39062500.000000.vtu temperature temperature 1e-12 1e-16
    temperature_analytical.vtu picard_masslumping_ts_1_t_78125.000000.vtu temperature_78125s temperature 2e-1 1e-4
    temperature_analytical.vtu picard_masslumping_ts_3_t_234375.000000.vtu temperature_234375s temperature 2e-1 1e-4
    temperature_analytical.vtu picard_masslumping_ts_65_t_5078125.000000.vtu temperature_5078125s temperature 1e-4 1e-4
    temperature_analytical.vtu picard_masslumping_ts_405_t_31640625.000000.vtu temperature_31640625s temperature 1e-4 1e-4
    temperature_analytical.vtu picard_masslumping_ts_500_t_39062500.000000.vtu temperature_39062500s temperature 1e-4 1e-4
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 1D_HeatConduction_neumann_newton_masslumping
    PATH Parabolic/T/1D_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS newton_masslumping.prj
    TESTER vtkdiff
    DIFF_DATA
    newton_masslumping_ts_1_t_78125.000000.vtu newton_masslumping_ts_1_t_78125.000000.vtu temperature temperature 1e-12 1e-16
    newton_masslumping_ts_3_t_234375.000000.vtu newton_masslumping_ts_3_t_234375.000000.vtu temperature temperature 1e-12 1e-16
    newton_masslumping_ts_65_t_5078125.000000.vtu newton_masslumping_ts_65_t_5078125.000000.vtu temperature temperature 1e-12 1e-16
    newton_masslumping_ts_405_t_31640625.000000.vtu newton_masslumping_ts_405_t_31640625.000000.vtu temperature temperature 1e-12 1e-16
    newton_masslumping_ts_500_t_39062500.000000.vtu newton_masslumping_ts_500_t_39062500.000000.vtu temperature temperature 1e-12 1e-16
    temperature_analytical.vtu newton_masslumping_ts_1_t_78125.000000.vtu temperature_78125s temperature 2e-1 1e-4
    temperature_analytical.vtu newton_masslumping_ts_3_t_234375.000000.vtu temperature_234375s temperature 2e-1 1e-4
    temperature_analytical.vtu newton_masslumping_ts_65_t_5078125.000000.vtu temperature_5078125s temperature 1e-4 1e-4
    temperature_analytical.vtu newton_masslumping_ts_405_t_31640625.000000.vtu temperature_31640625s temperature 1e-4 1e-4
    temperature_analytical.vtu newton_masslumping_ts_500_t_39062500.000000.vtu temperature_39062500s temperature 1e-4 1e-4
    # TODO: Fix on Apple M1
    REQUIREMENTS NOT OGS_USE_MPI AND NOT ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64"
)

AddTest(
    NAME 1D_HeatConduction_neumann_petsc_newtonls
    PATH Parabolic/T/1D_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS petsc_newtonls.prj
    TESTER vtkdiff
    DIFF_DATA
    newton_ts_1_t_78125.000000.vtu petsc_newtonls_ts_1_t_78125.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_3_t_234375.000000.vtu petsc_newtonls_ts_3_t_234375.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_65_t_5078125.000000.vtu petsc_newtonls_ts_65_t_5078125.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_405_t_31640625.000000.vtu petsc_newtonls_ts_405_t_31640625.000000.vtu temperature temperature 1e-12 1e-16
    newton_ts_500_t_39062500.000000.vtu petsc_newtonls_ts_500_t_39062500.000000.vtu temperature temperature 1e-12 1e-16
    temperature_analytical.vtu petsc_newtonls_ts_1_t_78125.000000.vtu temperature_78125s temperature 8e-2 1e-4
    temperature_analytical.vtu petsc_newtonls_ts_3_t_234375.000000.vtu temperature_234375s temperature 6e-2 1e-4
    temperature_analytical.vtu petsc_newtonls_ts_65_t_5078125.000000.vtu temperature_5078125s temperature 1e-4 1e-4
    temperature_analytical.vtu petsc_newtonls_ts_405_t_31640625.000000.vtu temperature_31640625s temperature 1e-4 1e-4
    temperature_analytical.vtu petsc_newtonls_ts_500_t_39062500.000000.vtu temperature_39062500s temperature 1e-4 1e-4
    REQUIREMENTS OGS_USE_MPI
)
# SQUARE 1x1 HEAT CONDUCTION TEST -- AXIALLY SYMMETRIC
# The results were compared to an analytical solution (method of manufactured
# solutions). The vtkdiff comparison is against the numerical solution.
AddTest(
    NAME 2D_HeatConduction_axi
    PATH Parabolic/T/2D_axially_symmetric
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_axi.prj
    TESTER vtkdiff
    DIFF_DATA
    square_1e2_axi_ts_10_t_1.000000.vtu square_1e2_axi_ts_10_t_1.000000.vtu temperature temperature 2e-15 0
    square_1e2_axi_ts_10_t_1.000000.vtu square_1e2_axi_ts_10_t_1.000000.vtu heat_flux heat_flux 1e-14 0
    REQUIREMENTS NOT OGS_USE_MPI
)
# WEDGE 1x1 HEAT CONDUCTION TEST -- same setup as above test but in cartesian coordinates
AddTest(
    NAME 2D_HeatConduction_wedge
    PATH Parabolic/T/2D_axially_symmetric
    EXECUTABLE ogs
    EXECUTABLE_ARGS wedge_1e2_axi_ang_0.02.prj
    TESTER vtkdiff
    DIFF_DATA
    wedge_ang_0.02_ts_10_t_1.000000.vtu wedge_ang_0.02_ts_10_t_1.000000.vtu temperature temperature 2e-12 0
    wedge_ang_0.02_ts_10_t_1.000000.vtu wedge_ang_0.02_ts_10_t_1.000000.vtu heat_flux heat_flux 1e-11 0
    REQUIREMENTS NOT OGS_USE_MPI
)

# The 25 BHE array benchmark
# test results are compared to 2D simulation result
AddTest(
    NAME BHE_Array_2D
    PATH Parabolic/T/2D_BHE_array
    RUNTIME 90
    EXECUTABLE ogs
    EXECUTABLE_ARGS bhe2d.prj
    TESTER vtkdiff
    DIFF_DATA
    standard_solution_bhe2d_ts_840_t_72576000.000000.vtu bhe2d_ts_840_t_72576000.000000.vtu temperature temperature 1e-12 0.0
    REQUIREMENTS NOT OGS_USE_MPI
)

# test the source term on a subdomain
AddTest(
    NAME 1D_HeatConduction_dirichlet_SourceTerm
    PATH Parabolic/T/1D_dirichlet_source-term
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_1_line_1e2_source_term.prj
    TESTER vtkdiff
    DIFF_DATA
    line_1_line_1e2_ts_500_t_39062500.000000_reference.vtu line_1_line_1e2_ts_500_t_39062500.000000.vtu temperature temperature 1.4e-11 0.0
    REQUIREMENTS NOT OGS_USE_MPI
)

# test the source term on a subdomain - parallel version
AddTest(
    NAME 1D_HeatConduction_dirichlet_SourceTerm_Parallel2
    PATH Parabolic/T/1D_dirichlet_source-term/2
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_1_line_1e2_source_term.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    REQUIREMENTS OGS_USE_MPI
    TESTER vtkdiff
    DIFF_DATA
    line_1_line_1e2_ts_0_t_0_000000_0.vtu line_1_line_1e2_ts_0_t_0_000000_0.vtu temperature temperature 1e-15 0.0
    line_1_line_1e2_ts_0_t_0_000000_1.vtu line_1_line_1e2_ts_0_t_0_000000_1.vtu temperature temperature 1e-15 0.0
    line_1_line_1e2_ts_500_t_39062500_000000_0.vtu line_1_line_1e2_ts_500_t_39062500_000000_0.vtu temperature temperature 1.4e-9 0.0
    line_1_line_1e2_ts_500_t_39062500_000000_1.vtu line_1_line_1e2_ts_500_t_39062500_000000_1.vtu temperature temperature 1.4e-9 0.0
)

# test the source term on a subdomain with the PETSc embedded executable file
AddTest(
    NAME 1D_HeatConduction_dirichlet_SourceTerm_PETSc
    PATH Parabolic/T/1D_dirichlet_source-term
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_1_line_1e2_source_term.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    line_1_line_1e2_ts_500_t_39062500.000000_reference.vtu line_1_line_1e2_ts_500_t_39062500.000000.vtu temperature temperature 1e-10 0.0
)

AddTest(
        NAME HeatConduction_t1_1Dsource
        PATH Parabolic/T/t1_1Dsource
        EXECUTABLE ogs
        EXECUTABLE_ARGS t1_1Dsource.prj
        TESTER vtkdiff
        DIFF_DATA
        t1_1Dsource_ts_1_t_1.000000.vtu t1_1Dsource_ts_1_t_1.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_t1_1Dsteady
        PATH Parabolic/T/t1_1Dsteady
        EXECUTABLE ogs
        EXECUTABLE_ARGS t1_1Dsteady.prj
        TESTER vtkdiff
        DIFF_DATA
        t1_1Dsteady_ts_1_t_1.000000.vtu t1_1Dsteady_ts_1_t_1.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_t2_1D1bt
        PATH Parabolic/T/t2_1D1bt
        EXECUTABLE ogs
        EXECUTABLE_ARGS t2_1D1bt.prj
        TESTER vtkdiff
        DIFF_DATA
        t2_1D1bt_ts_500_t_21600.000000.vtu t2_1D1bt_ts_500_t_21600.000000.vtu temperature temperature 10e-12 0.0
        t2_1D1bt_ts_1000_t_43200.000000.vtu t2_1D1bt_ts_1000_t_43200.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_t2_1D2bt
        PATH Parabolic/T/t2_1D2bt
        EXECUTABLE ogs
        EXECUTABLE_ARGS t2_1D2bt.prj
        TESTER vtkdiff
        DIFF_DATA
        t2_1D2bt_ts_3000_t_7776.000000.vtu t2_1D2bt_ts_3000_t_7776.000000.vtu temperature temperature 10e-12 0.0
        t2_1D2bt_ts_1500_t_3888.000000.vtu t2_1D2bt_ts_1500_t_3888.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_1D_LineSourceTerm
        PATH Parabolic/T/1D_line_source_term_tests
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term.prj
        TESTER vtkdiff
        DIFF_DATA
        mesh_1_line_100.vtu source_term_at_entire_line_ts_1_t_1.000000.vtu analytical_temperature temperature 1e-13 1e-13
        REQUIREMENTS NOT OGS_USE_MPI
)

if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/T/1D_line_source_term_tests/moving_source_term.prj)
endif()

# tests for line source term implementation
AddTest(
        NAME HeatConduction_2D_LineSourceTermLeft
        PATH Parabolic/T/2D_source_term_tests/line_source_term_left
        EXECUTABLE ogs
        EXECUTABLE_ARGS source_term_left.prj
        TESTER vtkdiff
        DIFF_DATA
        source_term_left_ts_1_t_1.000000.vtu source_term_left_ts_1_t_1.000000.vtu temperature temperature 1e-15 0.0
        source_term_left_ts_1_t_1.000000.vtu source_term_left_ts_1_t_1.000000.vtu heat_flux heat_flux 1e-15 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

# tests for line source term implementation with inclined elements
AddTest(
        NAME HeatConduction_2D_LineSourceTermLeft_inclined_elements
        PATH Parabolic/T/2D_source_term_tests/line_source_term_left
        EXECUTABLE ogs
        EXECUTABLE_ARGS source_term_left_r.prj
        TESTER vtkdiff
        DIFF_DATA
        source_term_left_r_ts_1_t_1.000000.vtu source_term_left_r_ts_1_t_1.000000.vtu temperature temperature 1e-15 0.0
        source_term_left_r_ts_1_t_1.000000.vtu source_term_left_r_ts_1_t_1.000000.vtu heat_flux heat_flux 1e-15 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

# For the special setup with a 'dirac' line source term at x=0.5 the
# analytical solution in 2 dimensions is valid:
# u(x,y) = -ln(sqrt((x-0.5)^2))/(2 * Pi)
AddTest(
        NAME HeatConduction_2D_LineSourceTermMiddle
        PATH Parabolic/T/2D_source_term_tests/line_source_term_x_0.5
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term_x_0.5.prj
        TESTER vtkdiff
        DIFF_DATA
        source_term_middle_ts_1_t_1.000000.vtu source_term_middle_ts_1_t_1.000000.vtu temperature temperature 7e-15 2e-14
        source_term_middle_ts_1_t_1.000000.vtu source_term_middle_ts_1_t_1.000000.vtu heat_flux heat_flux 7e-14 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_2D_LineSourceTermMiddle_Restricted
        PATH
        Parabolic/T/2D_source_term_tests/line_source_term_x_0.5_restricted_to_middle
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term_x_0.5.prj
        TESTER vtkdiff
        DIFF_DATA
        source_term_middle_restricted_ts_1_t_1.000000.vtu source_term_middle_restricted_ts_1_t_1.000000.vtu temperature temperature 1e-15 0.0
        source_term_middle_restricted_ts_1_t_1.000000.vtu source_term_middle_restricted_ts_1_t_1.000000.vtu heat_flux heat_flux 3e-15 4e-7
        REQUIREMENTS NOT OGS_USE_MPI
)

# tests for line source term implementation on a cubic domain
AddTest(
        NAME HeatConduction_3D_LineSourceTermMiddle
        PATH
        Parabolic/T/3D_line_source_term_tests/3D_line_source_term_middle
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term_x_0.5_y_0.5.prj
        TESTER vtkdiff
        DIFF_DATA
        3D_line_source_term_ts_1_t_1.000000.vtu 3D_line_source_term_ts_1_t_1.000000.vtu temperature temperature 2e-15 0.0
        3D_line_source_term_ts_1_t_1.000000.vtu 3D_line_source_term_ts_1_t_1.000000.vtu heat_flux heat_flux 7e-15 7e-13
        REQUIREMENTS NOT OGS_USE_MPI
)

# tests for line source term implementation on a cubic domain
AddTest(
        NAME HeatConduction_3D_LineSourceTermMiddle_Restricted
        PATH
        Parabolic/T/3D_line_source_term_tests/3D_line_source_term_middle_restricted
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term_x_0.5_y_0.5_restricted.prj
        TESTER vtkdiff
        DIFF_DATA
        3D_line_source_term_restricted_ts_1_t_1.000000.vtu 3D_line_source_term_restricted_ts_1_t_1.000000.vtu temperature temperature 1e-15 0.0
        3D_line_source_term_restricted_ts_1_t_1.000000.vtu 3D_line_source_term_restricted_ts_1_t_1.000000.vtu heat_flux heat_flux 1.1e-14 5e-12
        REQUIREMENTS NOT OGS_USE_MPI
)

# tests for line source term implementation on a cylindrical domain
# For the special setup with a line source term at position (xi, eta) the
# analytical solution in 2 dimensions is valid:
# u(x,y) = -ln(sqrt((x-xi)^2+(y-eta)^2))/(2 * Pi)
AddTest(
        NAME HeatConduction_3D_LineSourceTermInMiddleOfCylinder_49k_prisms
        PATH
        Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/49k_prisms
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term_in_cylinder.prj
        TESTER vtkdiff
        DIFF_DATA
        Cylinder_r_1_h_1_prism_49k.vtu 3D_line_source_term_in_cylinder_49k_ts_1_t_1.000000.vtu analytical_solution_temperature temperature 0.2 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_3D_LineSourceTermInMiddleOfCylinder_286k_prisms
        PATH
        Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/286k_prisms
        RUNTIME 32
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_source_term_in_cylinder.prj
        TESTER vtkdiff
        DIFF_DATA
        Cylinder_r_1_h_1_prism_286k.vtu 3D_line_source_term_in_cylinder_286k_ts_1_t_1.000000.vtu analytical_solution_temperature temperature 4e-3 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

if(NOT OGS_USE_PETSC)
    # Both tests above are executed in this notebook (without diff check). Maybe remove regular tests later?
    NotebookTest(NOTEBOOKFILE Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/heatconduction-line_source_term.md RUNTIME 10)
endif()
