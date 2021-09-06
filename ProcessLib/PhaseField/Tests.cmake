AddTest(
    NAME PhaseField_3D_beam_tens_AT2_iso
    PATH PhaseField/beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT2_iso_tensile.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_AT2_iso_tension_ts_10_t_1_000000_0.vtu AT2_iso_tension_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_AT2_iso_tension_ts_10_t_1_000000_0.vtu AT2_iso_tension_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
)
