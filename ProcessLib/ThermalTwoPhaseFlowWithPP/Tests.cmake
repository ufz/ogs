AddTest(
    NAME 2D_Thermal_TwoPhase_heatpipe
    PATH Parabolic/ThermalTwoPhaseFlowPP/HeatPipe
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_HeatPipe_quad_curve_small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_10000.000000.vtu thermaltwophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu capillary_pressure capillary_pressure 1e-8 1e-10
    ref_t_10000.000000.vtu thermaltwophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu gas_pressure gas_pressure 1e-8 1e-10
    ref_t_10000.000000.vtu thermaltwophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu saturation saturation 1e-8 1e-10
    ref_t_10000.000000.vtu thermaltwophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu temperature temperature 1e-8 1e-10
)
AddTest(
    NAME LARGE_2D_Thermal_TwoPhase_heatpipe
    PATH Parabolic/ThermalTwoPhaseFlowPP/HeatPipe
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_HeatPipe_quad_curve_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_1100000.000000.vtu thermaltwophaseflow_large_pcs_0_ts_2000_t_1100000.000000.vtu capillary_pressure capillary_pressure 1e-8 1e-10
    ref_t_1100000.000000.vtu thermaltwophaseflow_large_pcs_0_ts_2000_t_1100000.000000.vtu gas_pressure gas_pressure 1e-8 1e-10
    ref_t_1100000.000000.vtu thermaltwophaseflow_large_pcs_0_ts_2000_t_1100000.000000.vtu saturation saturation 1e-8 1e-10
    ref_t_1100000.000000.vtu thermaltwophaseflow_large_pcs_0_ts_2000_t_1100000.000000.vtu temperature temperature 1e-8 1e-10
)
