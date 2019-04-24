AddTest(
        NAME 1D_HeatConduction_dirichlet
        PATH Parabolic/T/1D_dirichlet
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.prj
        TESTER vtkdiff
        DIFF_DATA
        temperature_analytical.vtu line_60_heat_pcs_0_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature 1e-5 1e-5
        temperature_analytical.vtu line_60_heat_pcs_0_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature 1e-5 1e-5
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME 1D_HeatConduction_neumann
        PATH Parabolic/T/1D_neumann
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.prj
        TESTER vtkdiff
        DIFF_DATA
        temperature_analytical.vtu line_60_heat_pcs_0_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature 1e-4 1e-4
        temperature_analytical.vtu line_60_heat_pcs_0_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature 1e-4 1e-4
    REQUIREMENTS NOT OGS_USE_MPI
)

# SQUARE 1x1 HEAT CONDUCTION TEST -- AXIALLY SYMMETRIC
# test results are compared to 3D simulation on a wedge-shaped domain
AddTest(
        NAME 2D_HeatConduction_axi
        PATH Parabolic/T/2D_axially_symmetric
        EXECUTABLE ogs
        EXECUTABLE_ARGS square_1e2_axi.prj
        TESTER vtkdiff
        DIFF_DATA
        wedge_1e2_axi_ang_0.02_t_2s_extracted_surface.vtu square_1e2_axi_pcs_0_ts_2_t_2.000000.vtu temperature temperature 1.7e-5 1e-5
        wedge_1e2_axi_ang_0.02_t_2s_extracted_surface.vtu square_1e2_axi_pcs_0_ts_2_t_2.000000.vtu heat_flux_x heat_flux_x 1.7e-5 1e-5
    REQUIREMENTS NOT OGS_USE_MPI
)
# # WEDGE 1x1 HEATCONDUCTION TEST -- computes reference results for the above test
# AddTest(
#      NAME 2D_HeatConduction_wedge
#      PATH Parabolic/T/2D_axially_symmetric
#      EXECUTABLE ogs
#      EXECUTABLE_ARGS wedge_1e2_axi_ang_0.02.prj
# )

# The 25 BHE array benchmark
# test results are compared to 2D simulation result
AddTest(
    NAME LARGE_BHE_Array_2D
    PATH Parabolic/T/2D_BHE_array
    RUNTIME 90
    EXECUTABLE ogs
    EXECUTABLE_ARGS bhe2d.prj
    TESTER vtkdiff
    DIFF_DATA
    standard_solution_bhe2d_pcs_0_ts_840_t_72576000.000000.vtu bhe2d_pcs_0_ts_840_t_72576000.000000.vtu temperature temperature 1e-12 0.0
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
    line_1_line_1e2_pcs_0_ts_500_t_39062500.000000_reference.vtu line_1_line_1e2_pcs_0_ts_500_t_39062500.000000.vtu temperature temperature 1e-11 0.0
    REQUIREMENTS NOT OGS_USE_MPI
)

# test the source term on a subdomain with the PETSc embedded executable file
AddTest(
    NAME 1D_HeatConduction_dirichlet_SourceTerm_PETSc
    PATH Parabolic/T/1D_dirichlet_source-term
    EXECUTABLE_ARGS line_1_line_1e2_source_term.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    line_1_line_1e2_pcs_0_ts_500_t_39062500.000000_reference.vtu line_1_line_1e2_pcs_0_ts_500_t_39062500_000000_0.vtu temperature temperature 1e-10 0.0
)

AddTest(
        NAME HeatConduction_t1_1Dsource
        PATH Parabolic/T/t1_1Dsource
        EXECUTABLE ogs
        EXECUTABLE_ARGS t1_1Dsource.prj
        TESTER vtkdiff
        DIFF_DATA
        t1_1Dsource_pcs_0_ts_1_t_1.000000.vtu t1_1Dsource_pcs_0_ts_1_t_1.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_t1_1Dsteady
        PATH Parabolic/T/t1_1Dsteady
        EXECUTABLE ogs
        EXECUTABLE_ARGS t1_1Dsteady.prj
        TESTER vtkdiff
        DIFF_DATA
        t1_1Dsteady_pcs_0_ts_1_t_1.000000.vtu t1_1Dsteady_pcs_0_ts_1_t_1.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_t2_1D1bt
        PATH Parabolic/T/t2_1D1bt
        EXECUTABLE ogs
        EXECUTABLE_ARGS t2_1D1bt.prj
        TESTER vtkdiff
        DIFF_DATA
        t2_1D1bt_pcs_0_ts_500_t_21600.000000.vtu t2_1D1bt_pcs_0_ts_500_t_21600.000000.vtu temperature temperature 10e-12 0.0
        t2_1D1bt_pcs_0_ts_1000_t_43200.000000.vtu t2_1D1bt_pcs_0_ts_1000_t_43200.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME HeatConduction_t2_1D2bt
        PATH Parabolic/T/t2_1D2bt
        EXECUTABLE ogs
        EXECUTABLE_ARGS t2_1D2bt.prj
        TESTER vtkdiff
        DIFF_DATA
        t2_1D2bt_pcs_0_ts_3000_t_7776.000000.vtu t2_1D2bt_pcs_0_ts_3000_t_7776.000000.vtu temperature temperature 10e-12 0.0
        t2_1D2bt_pcs_0_ts_1500_t_3888.000000.vtu t2_1D2bt_pcs_0_ts_1500_t_3888.000000.vtu temperature temperature 10e-12 0.0
        REQUIREMENTS NOT OGS_USE_MPI
)
