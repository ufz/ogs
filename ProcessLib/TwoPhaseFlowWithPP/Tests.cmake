AddTest(
    NAME 2D_TwoPhase_PP_Lia_quad_short
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_Lia_quad_short.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 11
    DIFF_DATA
    h2_Liako_20.vtu twophaseflow_ts_290_t_20.000000.vtu saturation saturation 1e-2 1e-4
    h2_Liako_20.vtu twophaseflow_ts_290_t_20.000000.vtu capillary_pressure capillary_pressure 20 1e-3
    h2_Liako_20.vtu twophaseflow_ts_290_t_20.000000.vtu gas_pressure gas_pressure 20 1e-3
)
AddTest(
    NAME 2D_TwoPhase_PP_Lia_quad_large
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_Lia_quad_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 16
    DIFF_DATA
    h2_Liako_1198.vtu twophaseflow_ts_1198_t_1000.000000.vtu saturation saturation 1e-2 1e-3
    h2_Liako_1198.vtu twophaseflow_ts_1198_t_1000.000000.vtu capillary_pressure capillary_pressure 20 1e-2
    h2_Liako_1198.vtu twophaseflow_ts_1198_t_1000.000000.vtu gas_pressure gas_pressure 20 1e-2
)
AddTest(
    NAME 1D_TwoPhase_PP_mcwt
    PATH Parabolic/TwoPhaseFlowPP/McWhorter
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_mcwt_line.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 12
    DIFF_DATA
    twophaseflow_ts_627_t_1000.000000.vtu twophaseflow_ts_627_t_1000.000000.vtu saturation saturation 1e-3 1e-4
    twophaseflow_ts_627_t_1000.000000.vtu twophaseflow_ts_627_t_1000.000000.vtu capillary_pressure capillary_pressure 1e-3 1e-4
    twophaseflow_ts_627_t_1000.000000.vtu twophaseflow_ts_627_t_1000.000000.vtu gas_pressure gas_pressure 1e-3 1e-4
)
