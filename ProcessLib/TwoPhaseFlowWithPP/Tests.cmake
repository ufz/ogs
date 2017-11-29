AddTest(
    NAME 2D_TwoPhase_PP_Lia_quad_1
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_Lia_quad_short.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h2_Liako_20.vtu twophaseflow_pcs_0_ts_218_t_20.000000.vtu saturation saturation 1e-2 1e-4
)
AddTest(
    NAME 2D_TwoPhase_PP_Lia_quad_2
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_Lia_quad_short.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h2_Liako_20.vtu twophaseflow_pcs_0_ts_218_t_20.000000.vtu capillary_pressure capillary_pressure 20 1e-3
    h2_Liako_20.vtu twophaseflow_pcs_0_ts_218_t_20.000000.vtu gas_pressure gas_pressure 20 1e-3
)
AddTest(
    NAME LARGE_2D_TwoPhase_PP_Lia_quad_1
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_Lia_quad_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h2_Liako_1198.vtu twophaseflow_pcs_0_ts_1198_t_1000.000000.vtu SATURATION1 saturation 1e-2 1e-3
)
AddTest(
    NAME LARGE_2D_TwoPhase_PP_Lia_quad_2
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_Lia_quad_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    h2_Liako_1198.vtu twophaseflow_pcs_0_ts_1198_t_1000.000000.vtu PRESSURE1 capillary_pressure 20 1e-2
    h2_Liako_1198.vtu twophaseflow_pcs_0_ts_1198_t_1000.000000.vtu PRESSURE2 gas_pressure 20 1e-2
)
AddTest(
    NAME 1D_TwoPhase_PP_mcwt_1
    PATH Parabolic/TwoPhaseFlowPP/McWorts
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_mcwt_line.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    mcwt_1000.vtu twophaseflow_pcs_0_ts_519_t_1000.000000.vtu SATURATION1 saturation 1e-2 1e-4
)
AddTest(
    NAME 1D_TwoPhase_PP_mcwt_2
    PATH Parabolic/TwoPhaseFlowPP/McWorts
    EXECUTABLE ogs
    EXECUTABLE_ARGS TwoPhase_mcwt_line.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    mcwt_1000.vtu twophaseflow_pcs_0_ts_519_t_1000.000000.vtu PRESSURE1 capillary_pressure 10 1e-3
    mcwt_1000.vtu twophaseflow_pcs_0_ts_519_t_1000.000000.vtu PRESSURE2 gas_pressure 10 1e-3
)
