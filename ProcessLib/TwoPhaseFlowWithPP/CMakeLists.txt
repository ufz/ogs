AddTest(
    NAME 2D_TwoPhase_PP_Lia_quad
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_Lia_quad2_short.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-12
    DIFF_DATA
    Lia_20.vtu twophaseflow_pcs_0_ts_119_t_20.000000.vtu capillary_pressure capillary_pressure
    Lia_20.vtu twophaseflow_pcs_0_ts_119_t_20.000000.vtu gas_pressure gas_pressure
    Lia_20.vtu twophaseflow_pcs_0_ts_119_t_20.000000.vtu saturation saturation
)
AddTest(
    NAME LARGE_2D_TwoPhase_PP_Lia_quad
    PATH Parabolic/TwoPhaseFlowPP/Liakopoulos
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_Lia_quad2_large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-12
    DIFF_DATA
    Lia_1000.vtu twophaseflow_pcs_0_ts_1180_t_1000.000000.vtu capillary_pressure capillary_pressure
    Lia_1000.vtu twophaseflow_pcs_0_ts_1180_t_1000.000000.vtu gas_pressure gas_pressure
    Lia_1000.vtu twophaseflow_pcs_0_ts_1180_t_1000.000000.vtu saturation saturation
)
