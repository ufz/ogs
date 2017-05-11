AddTest(
    NAME LARGE_3D_Crack_Beam_Tension
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-1 RELTOL 1e-2
    DIFF_DATA
    analytical_results.vtu beam3d_pcs_0_ts_10000_t_20.000000.vtu displacement displacement
    analytical_results.vtu beam3d_pcs_0_ts_10000_t_20.000000.vtu phasefield phasefield