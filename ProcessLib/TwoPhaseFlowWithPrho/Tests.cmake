AddTest(
    NAME LARGE_2D_TwoPhase_Prho_MoMaS
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_quad.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_ts_10_t_10000.000000.vtu twophaseflow_pcs_0_ts_10_t_10000.000000.vtu liquid_pressure liquid_pressure 1e-8 1e-12
    ref_ts_10_t_10000.000000.vtu twophaseflow_pcs_0_ts_10_t_10000.000000.vtu overall_mass_density overall_mass_density 1e-8 1e-12
    ref_ts_10_t_10000.000000.vtu twophaseflow_pcs_0_ts_10_t_10000.000000.vtu saturation saturation 1e-8 1e-12
)

AddTest(
    NAME LARGE_2D_TwoPhase_Prho_MoMaS_Adaptive_dt
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_quad_adaptive_dt.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_twophaseflow_pcs_0_ts_19_t_100000.000000.vtu twophaseflow_adaptive_dt_pcs_0_ts_19_t_100000.000000.vtu liquid_pressure liquid_pressure 1e-7 1e-6
    ref_twophaseflow_pcs_0_ts_19_t_100000.000000.vtu twophaseflow_adaptive_dt_pcs_0_ts_19_t_100000.000000.vtu overall_mass_density overall_mass_density 1e-7 1e-6
    ref_twophaseflow_pcs_0_ts_19_t_100000.000000.vtu twophaseflow_adaptive_dt_pcs_0_ts_19_t_100000.000000.vtu saturation saturation 1e-7 1e-6
)

AddTest(
    NAME 2D_TwoPhase_Prho_MoMaS
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_10000.000000.vtu twophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu liquid_pressure liquid_pressure 1e-6 1e-10
    ref_t_10000.000000.vtu twophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu overall_mass_density overall_mass_density 1e-6 1e-10
    ref_t_10000.000000.vtu twophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu saturation saturation 1e-6 1e-10
)
