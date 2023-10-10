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
        expected_AT1_iso_tension_ts_10_t_1_000000_0.vtu AT1_iso_tension_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
        expected_AT1_iso_tension_ts_10_t_1_000000_0.vtu AT1_iso_tension_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
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
    RUNTIME 260
    DIFF_DATA
        expected_AT2_iso_tension_ts_10_t_1_000000_0.vtu AT2_iso_tension_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
        expected_AT2_iso_tension_ts_10_t_1_000000_0.vtu AT2_iso_tension_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
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
        expected_AT1_vd_tension_ts_10_t_1_000000_0.vtu AT1_vd_tension_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
        expected_AT1_vd_tension_ts_10_t_1_000000_0.vtu AT1_vd_tension_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
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
    RUNTIME 140
    DIFF_DATA
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu displacement displacement 1e-5 0
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu displacement displacement 1e-5 0
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu AT1_vd_tension_2core_ts_10_t_1_000000_0.vtu phasefield phasefield 1e-6 0
        expected_AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu AT1_vd_tension_2core_ts_10_t_1_000000_1.vtu phasefield phasefield 1e-6 0
)

if(OGS_USE_MPI)
    OgsTest(PROJECTFILE PhaseField/beam/bar_COHESIVE_linear.prj
        RUNTIME 300
        WRAPPER mpirun -np 1
    )
endif()

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
        expected_bar_COHESIVE_exponential_ts_11_t_1.100000.vtu bar_COHESIVE_exponential_ts_11_t_1.100000.vtu displacement displacement 1e-5 0
        expected_bar_COHESIVE_exponential_ts_11_t_1.100000.vtu bar_COHESIVE_exponential_ts_11_t_1.100000.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_2D_surfing_AT1_vd
    PATH PhaseField/surfing
    EXECUTABLE ogs
    EXECUTABLE_ARGS surfing.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_surfing_ts_20_t_1_000000_0.vtu surfing_ts_20_t_1.000000.vtu displacement displacement 1e-5 0
        expected_surfing_ts_20_t_1_000000_0.vtu surfing_ts_20_t_1.000000.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_2D_K_regime_HF_2cores
    PATH PhaseField/k_regime_HF
    EXECUTABLE ogs
    EXECUTABLE_ARGS 2D_bm_0p01.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 2
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_0.vtu 2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_0.vtu displacement displacement 1e-5 0
        expected_2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_1.vtu 2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_1.vtu displacement displacement 1e-5 0
        expected_2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_0.vtu 2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_0.vtu phasefield phasefield 1e-6 0
        expected_2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_1.vtu 2D_PropagatingCrack_AT1_h0p01_ts_2_t_0_020000_1.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_3D_beam_tens_AT2_vd_ortho
    PATH PhaseField/beam/voldev-ortho
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT2_vd_tensile_VZ.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 120
    DIFF_DATA
    expected_AT2_vd_tensile_ts_10_t_1.000000.vtu AT2_vd_tensile_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
    expected_AT2_vd_tensile_ts_10_t_1.000000.vtu AT2_vd_tensile_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME PhaseField_2D_square_tens_AT2_masonry_ortho
    PATH PhaseField/single_edge_notched/masonry-ortho
    EXECUTABLE ogs
    EXECUTABLE_ARGS shear.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 600
    DIFF_DATA
    expected_AT2_OrthoMasonry_ts_100_t_1.000000.vtu AT2_OrthoMasonry_ts_100_t_1.000000.vtu displacement displacement 1e-5 0
    expected_AT2_OrthoMasonry_ts_100_t_1.000000.vtu AT2_OrthoMasonry_ts_100_t_1.000000.vtu phasefield phasefield 1e-6 0
)

if(OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE PhaseField/surfing_jupyter_notebook/surfing_pyvista.ipynb RUNTIME 25)
    NotebookTest(NOTEBOOKFILE PhaseField/beam_jupyter_notebook/beam.ipynb RUNTIME 500 PROPERTIES PROCESSORS 3)
    NotebookTest(NOTEBOOKFILE PhaseField/tpb_jupyter_notebook/TPB.ipynb RUNTIME 110 PROPERTIES PROCESSORS 4)
    NotebookTest(NOTEBOOKFILE PhaseField/kregime_jupyter_notebook/Kregime_Static_jupyter.ipynb RUNTIME 40)
    NotebookTest(NOTEBOOKFILE PhaseField/PForthotropy_jupyter_notebook/sen_shear.ipynb RUNTIME 500 PROPERTIES PROCESSORS 4)
    NotebookTest(NOTEBOOKFILE PhaseField/Kregime_Propagating_jupyter_notebook/Kregime_Propagating_jupyter.ipynb RUNTIME 550)
endif()
