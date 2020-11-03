AddTest(
    NAME 2D_TwoPhase_Prho_MoMaS
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    RUNTIME 30
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_quad.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    reference_t_10000.000000.vtu twophaseflow_ts_10_t_10000.000000.vtu liquid_pressure liquid_pressure 1e-6 1e-12
    reference_t_10000.000000.vtu twophaseflow_ts_10_t_10000.000000.vtu overall_mass_density overall_mass_density 0 1e-12
    reference_t_10000.000000.vtu twophaseflow_ts_10_t_10000.000000.vtu saturation saturation 1e-8 0
)

AddTest(
    NAME 2D_TwoPhase_Prho_MoMaS_Adaptive_dt
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    RUNTIME 550
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_quad_adaptive_dt.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    reference_t_10000.000000.vtu
    twophaseflow_adaptive_dt_ts_10_t_10000.000000.vtu liquid_pressure liquid_pressure 1e-6 1e-12
    reference_t_10000.000000.vtu twophaseflow_adaptive_dt_ts_10_t_10000.000000.vtu overall_mass_density overall_mass_density 1e-10 1e-16
    reference_t_10000.000000.vtu twophaseflow_adaptive_dt_ts_10_t_10000.000000.vtu saturation saturation 1e-10 0
    reference_t_100000.000000.vtu twophaseflow_adaptive_dt_ts_108_t_100000.000000.vtu liquid_pressure liquid_pressure 1e-6 1e-12
    reference_t_100000.000000.vtu twophaseflow_adaptive_dt_ts_108_t_100000.000000.vtu overall_mass_density overall_mass_density 1e-10 1e-16
    reference_t_100000.000000.vtu twophaseflow_adaptive_dt_ts_108_t_100000.000000.vtu saturation saturation 1e-10 0.0
)

AddTest(
    NAME 2D_TwoPhase_Prho_MoMaS_small
    PATH Parabolic/TwoPhaseFlowPrho/MoMaS
    RUNTIME 40
    EXECUTABLE ogs
    EXECUTABLE_ARGS Twophase_MoMaS_small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_t_10000.000000.vtu twophaseflow_small_ts_50_t_10000.000000.vtu liquid_pressure liquid_pressure 1 1e-6
    ref_t_10000.000000.vtu twophaseflow_small_ts_50_t_10000.000000.vtu overall_mass_density overall_mass_density 1e-4 1e-10
    ref_t_10000.000000.vtu twophaseflow_small_ts_50_t_10000.000000.vtu saturation saturation 1e-6 1e-10
)
