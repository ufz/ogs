AddTest(
        NAME 1D_HeatConduction_dirichlet
        PATH Parabolic/T/1D_dirichlet
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.prj
        TESTER vtkdiff
        ABSTOL 1e-5 RELTOL 1e-5
        DIFF_DATA
        temperature_analytical.vtu line_60_heat_pcs_0_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature
        temperature_analytical.vtu line_60_heat_pcs_0_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME 1D_HeatConduction_neumann
        PATH Parabolic/T/1D_neumann
        EXECUTABLE ogs
        EXECUTABLE_ARGS line_60_heat.prj
        TESTER vtkdiff
        ABSTOL 1e-4 RELTOL 1e-4
        DIFF_DATA
        temperature_analytical.vtu line_60_heat_pcs_0_ts_65_t_5078125.000000.vtu Temperature_Analytical_2months temperature
        temperature_analytical.vtu line_60_heat_pcs_0_ts_405_t_31640625.000000.vtu Temperature_Analytical_1year temperature
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
        ABSTOL 1.7e-5 RELTOL 1e-5
        DIFF_DATA
        wedge_1e2_axi_ang_0.02_t_2s_extracted_surface.vtu square_1e2_axi_pcs_0_ts_2_t_2.000000.vtu temperature temperature
        wedge_1e2_axi_ang_0.02_t_2s_extracted_surface.vtu square_1e2_axi_pcs_0_ts_2_t_2.000000.vtu heat_flux_x heat_flux_x
    REQUIREMENTS NOT OGS_USE_MPI
)
# # WEDGE 1x1 HEATCONDUCTION TEST -- computes reference results for the above test
# AddTest(
#      NAME 2D_HeatConduction_wedge
#      PATH Parabolic/T/2D_axially_symmetric
#      EXECUTABLE ogs
#      EXECUTABLE_ARGS wedge_1e2_axi_ang_0.02.prj
# )
