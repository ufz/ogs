# HMPhaseField; Phase-field fracture, Small deformation, Hydromechanics,
# Hydraulic fracture

# Beam fracture for M and PhaseField processes
AddTest(
    NAME HMPhaseField_3D_beam_tens_AT1_iso
    PATH HMPhaseField/Beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_AT1_iso_tension_ts_10_t_1.000000.vtu AT1_iso_tension_ts_10_t_1.000000.vtu displacement displacement 1e-10 1e-8
        expected_AT1_iso_tension_ts_10_t_1.000000.vtu AT1_iso_tension_ts_10_t_1.000000.vtu phasefield phasefield 1e-10 1e-8
)

AddTest(
    NAME HMPhaseField_3D_beam_tens_AT2_iso
    PATH HMPhaseField/Beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT2_iso_tensile.xml
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 120
    DIFF_DATA
        expected_AT2_iso_tension_ts_10_t_1.000000.vtu AT2_iso_tension_ts_10_t_1.000000.vtu displacement displacement 1e-9 1e-6
        expected_AT2_iso_tension_ts_10_t_1.000000.vtu AT2_iso_tension_ts_10_t_1.000000.vtu phasefield phasefield 1e-9 1e-6
)

AddTest(
    NAME HMPhaseField_3D_beam_tens_AT1_vd
    PATH HMPhaseField/Beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS AT1_vd_tensile.xml
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 18
    DIFF_DATA
        expected_AT1_vd_tension_ts_10_t_1.000000.vtu AT1_vd_tension_ts_10_t_1.000000.vtu displacement displacement 1e-10 1e-8
        expected_AT1_vd_tension_ts_10_t_1.000000.vtu AT1_vd_tension_ts_10_t_1.000000.vtu phasefield phasefield 1e-10 1e-8
)

if(OGS_USE_MPI)
    OgsTest(PROJECTFILE HMPhaseField/Beam/bar_COHESIVE_linear.xml RUNTIME 100
            WRAPPER mpirun -np 1
    )
endif()

AddTest(
    NAME HMPhaseField_3D_beam_tens_COHESIVE_exponential_es
    PATH HMPhaseField/Beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar_COHESIVE_exponential.xml
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 80
    DIFF_DATA
        expected_bar_COHESIVE_exponential_ts_11_t_1.100000.vtu bar_COHESIVE_exponential_ts_11_t_1.100000.vtu displacement displacement 1e-9 1e-6
        expected_bar_COHESIVE_exponential_ts_11_t_1.100000.vtu bar_COHESIVE_exponential_ts_11_t_1.100000.vtu phasefield phasefield 1e-9 1e-6
)

# Small deformation for M process
OgsTest(PROJECTFILE HMPhaseField/Mechanics/bar.prj)
AddTest(
        NAME HMPhaseField_LinearElasatics_3D
        PATH HMPhaseField/Mechanics
        EXECUTABLE ogs
        EXECUTABLE_ARGS bar_3D.xml
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        ../../Mechanics/Linear/MaterialForces/bar_3D_out_ts_2_t_1.000000.vtu bar_3D_out_ts_2_t_1.000000.vtu displacement displacement 1e-14 1e-15
)
AddTest(
        NAME HMPhaseField_SDL_disc_with_hole_mfront
        PATH HMPhaseField/Mechanics
        EXECUTABLE ogs
        EXECUTABLE_ARGS disc_with_hole.xml
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_MPI
        DIFF_DATA
        disc_with_hole_expected_ts_4_t_1.000000.vtu disc_with_hole_ts_4_t_1.000000.vtu displacement displacement 2e-16 1e-16
)

# Fluid flow in a fracture for H process
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE HMPhaseField/LiquidFlowInFracture/fractureflow.prj
            RUNTIME 10 WRAPPER mpirun -np 1
    )
endif()

# Hydromechanics (1D Consolidation) for H and M processes
AddTest(
    NAME HMPhaseField_consolidation
    PATH HMPhaseField/Consolidation
    EXECUTABLE ogs
    EXECUTABLE_ARGS consolidation.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 20
    DIFF_DATA
        expected_consolidation_ts_10_t_100.000000.vtu consolidation_ts_10_t_100.000000.vtu displacement displacement 1e-12 1e-10
        expected_consolidation_ts_10_t_100.000000.vtu consolidation_ts_10_t_100.000000.vtu phasefield phasefield 1e-14 1e-14
        expected_consolidation_ts_10_t_100.000000.vtu consolidation_ts_10_t_100.000000.vtu pressure pressure 1e-12 1e-10
)

# KGD (toughness-dominated hydraulic fracture) for H, M and PhaseField processes
AddTest(
    NAME HMPhaseField_KGD_4core
    PATH HMPhaseField/KGD
    EXECUTABLE ogs
    EXECUTABLE_ARGS KGD.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 300
    DIFF_DATA
        expected_KGD_4core_ts_5_t_200_000000_0.vtu KGD_4core_ts_5_t_200_000000_0.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_1.vtu KGD_4core_ts_5_t_200_000000_1.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_2.vtu KGD_4core_ts_5_t_200_000000_2.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_3.vtu KGD_4core_ts_5_t_200_000000_3.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_0.vtu KGD_4core_ts_5_t_200_000000_0.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_1.vtu KGD_4core_ts_5_t_200_000000_1.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_2.vtu KGD_4core_ts_5_t_200_000000_2.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_3.vtu KGD_4core_ts_5_t_200_000000_3.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_0.vtu KGD_4core_ts_5_t_200_000000_0.vtu pressure pressure 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_1.vtu KGD_4core_ts_5_t_200_000000_1.vtu pressure pressure 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_2.vtu KGD_4core_ts_5_t_200_000000_2.vtu pressure pressure 1e-8 1e-5
        expected_KGD_4core_ts_5_t_200_000000_3.vtu KGD_4core_ts_5_t_200_000000_3.vtu pressure pressure 1e-8 1e-5
)

AddTest(
    NAME HMPhaseField_KGD_4core_cohesive_linear_softening
    PATH HMPhaseField/KGD
    EXECUTABLE ogs
    EXECUTABLE_ARGS KGD_cohesive_linear_softening.xml
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    RUNTIME 300
    DIFF_DATA
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_0.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_0.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_1.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_1.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_2.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_2.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_3.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_3.vtu displacement displacement 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_0.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_0.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_1.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_1.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_2.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_2.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_3.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_3.vtu phasefield phasefield 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_0.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_0.vtu pressure pressure 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_1.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_1.vtu pressure pressure 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_2.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_2.vtu pressure pressure 1e-8 1e-5
        expected_KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_3.vtu KGD_4core_cohesive_linear_softening_ts_5_t_200_000000_3.vtu pressure pressure 1e-8 1e-5
)

if( OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE HMPhaseField/GreatCell/GreatCellHM_VPF.py RUNTIME 1300)
endif()
