AddTest(
    NAME PhaseField_3D_beam_tens_AT1_iso
    PATH PhaseField/beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT1_iso_tensile.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_AT1_iso_tension_ts_10_t_1_000000_0.vtu AT1_iso_tension_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_AT1_iso_tension_ts_10_t_1_000000_0.vtu AT1_iso_tension_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
)

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

AddTest(
    NAME PhaseField_3D_beam_tens_AT1_vd
    PATH PhaseField/beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT1_vd_tensile.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_AT1_vd_tension_ts_10_t_1_000000_0.vtu AT1_vd_tension_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_AT1_vd_tension_ts_10_t_1_000000_0.vtu AT1_vd_tension_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_3D_beam_tens_AT1_vd_2core
    PATH PhaseField/beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT1_vd_tensile_2core.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu displacement displacement 1e-5 0
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_3D_beam_tens_COHESIVE_linear_es
    PATH PhaseField/beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar_COHESIVE_linear.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 300
    DIFF_DATA
        expected_bar_COHESIVE_linear_ts_10_t_1_000000_0.vtu bar_COHESIVE_linear_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_bar_COHESIVE_linear_ts_10_t_1_000000_0.vtu bar_COHESIVE_linear_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_3D_beam_tens_COHESIVE_exponential_es
    PATH PhaseField/beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar_COHESIVE_exponential.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 300
    DIFF_DATA
        expected_bar_COHESIVE_exponential_ts_10_t_1_000000_0.vtu bar_COHESIVE_exponential_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_bar_COHESIVE_exponential_ts_10_t_1_000000_0.vtu bar_COHESIVE_exponential_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
)
