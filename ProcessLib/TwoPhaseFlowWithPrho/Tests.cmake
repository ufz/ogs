AddTest(
    NAME LARGE_2D_TwoPhase_Prho_MoMaS
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_quad.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-12
    DIFF_DATA
    ref_ts_10_t_10000.000000.vtu twophaseflow_pcs_0_ts_10_t_10000.000000.vtu liquid_pressure liquid_pressure
    ref_ts_10_t_10000.000000.vtu twophaseflow_pcs_0_ts_10_t_10000.000000.vtu overall_mass_density overall_mass_density
    ref_ts_10_t_10000.000000.vtu twophaseflow_pcs_0_ts_10_t_10000.000000.vtu saturation saturation
)

AddTest(
    NAME LARGE_2D_TwoPhase_Prho_MoMaS_Adaptive_dt
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_quad_adaptive_dt.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-7 RELTOL 1e-6
    DIFF_DATA
    ref_twophaseflow_pcs_0_ts_19_t_100000.000000.vtu twophaseflow_adaptive_dt_pcs_0_ts_19_t_100000.000000.vtu liquid_pressure liquid_pressure
    ref_twophaseflow_pcs_0_ts_19_t_100000.000000.vtu twophaseflow_adaptive_dt_pcs_0_ts_19_t_100000.000000.vtu overall_mass_density overall_mass_density
    ref_twophaseflow_pcs_0_ts_19_t_100000.000000.vtu twophaseflow_adaptive_dt_pcs_0_ts_19_t_100000.000000.vtu saturation saturation
)

AddTest(
    NAME 2D_TwoPhase_Prho_MoMaS
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-6 RELTOL 1e-10
    DIFF_DATA
    ref_t_10000.000000.vtu twophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu liquid_pressure liquid_pressure
    ref_t_10000.000000.vtu twophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu overall_mass_density overall_mass_density
    ref_t_10000.000000.vtu twophaseflow_small_pcs_0_ts_10_t_10000.000000.vtu saturation saturation
)
