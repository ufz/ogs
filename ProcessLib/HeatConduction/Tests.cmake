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
    NAME BHE_Array_2D
    PATH Parabolic/T/2D_BHE_array
    EXECUTABLE ogs
    EXECUTABLE_ARGS bhe2d.prj
    TESTER vtkdiff
    DIFF_DATA
    standard_solution_bhe2d_pcs_0_ts_840_t_72576000.000000.vtu bhe2d_pcs_0_ts_840_t_72576000.000000.vtu temperature temperature 1e-12 0.0
    REQUIREMENTS NOT OGS_USE_MPI
)
