AddTest(
    NAME 2D_Thermal_TwoPhase_heatpipe_small
    PATH Parabolic/ThermalTwoPhaseFlowPP/HeatPipe
    RUNTIME 10
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_HeatPipe_quad_curve_small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_10000.000000.vtu thermaltwophaseflow_small_ts_100_t_10000.000000.vtu capillary_pressure capillary_pressure 1e-8 1e-10
    ref_t_10000.000000.vtu thermaltwophaseflow_small_ts_100_t_10000.000000.vtu gas_pressure gas_pressure 1e-8 1e-10
    ref_t_10000.000000.vtu thermaltwophaseflow_small_ts_100_t_10000.000000.vtu saturation saturation 1e-8 1e-10
    ref_t_10000.000000.vtu thermaltwophaseflow_small_ts_100_t_10000.000000.vtu temperature temperature 1e-8 1e-10
)
AddTest(
    NAME 2D_Thermal_TwoPhase_heatpipe_large
    PATH Parabolic/ThermalTwoPhaseFlowPP/HeatPipe
    RUNTIME 170
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_HeatPipe_quad_curve_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_1400000.000000.vtu thermaltwophaseflow_large_ts_1400_t_1400000.000000.vtu capillary_pressure capillary_pressure 9e-5 5e-10
    ref_t_1400000.000000.vtu thermaltwophaseflow_large_ts_1400_t_1400000.000000.vtu gas_pressure gas_pressure 5e-7 1e-10
    ref_t_1400000.000000.vtu thermaltwophaseflow_large_ts_1400_t_1400000.000000.vtu saturation saturation 1e-8 1e-10
    ref_t_1400000.000000.vtu thermaltwophaseflow_large_ts_1400_t_1400000.000000.vtu temperature temperature 1e-8 1e-10
)
AddTest(
    NAME 1D_TwoPhase_TCE_diffusion_small
    PATH Parabolic/ThermalTwoPhaseFlowPP/TCEDiffusion
    RUNTIME 12
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_TCE_diffusion_1D_small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_1000000.000000.vtu Twophase_TCE_diffusion_1D_ts_100_t_1000000.000000.vtu capillary_pressure capillary_pressure 1e-8 1e-8
    ref_t_1000000.000000.vtu Twophase_TCE_diffusion_1D_ts_100_t_1000000.000000.vtu gas_pressure gas_pressure 1e-8 1e-8
    ref_t_1000000.000000.vtu Twophase_TCE_diffusion_1D_ts_100_t_1000000.000000.vtu Xc Xc 1e-8 1e-8
    ref_t_1000000.000000.vtu Twophase_TCE_diffusion_1D_ts_100_t_1000000.000000.vtu temperature temperature 1e-8 1e-8
)
