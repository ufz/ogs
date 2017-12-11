AddTest(
    NAME 2D_ComponentTransport_ConcentrationDiffusionOnly
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionOnly.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu linear_1_to_0 concentration 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero pressure 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero_vector_2d darcy_velocity 1e-7 1e-10
    VIS DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_ConcentrationDiffusionAndStorage
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionAndStorage.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero pressure 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero_vector_2d darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvection
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvection.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndGravityAndDispersionHalf
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndGravityAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu concentration concentration 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    VIS DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDispersion
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDecay
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDecay.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDispersionHalf
    PATH Parabolic/ComponentTransport/SimpleSynthetics/
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu concentration concentration 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu
)

AddTest(
    NAME LARGE_2D_ComponentTransport_Goswami
    PATH Parabolic/ComponentTransport/goswami/
    EXECUTABLE ogs
    EXECUTABLE_ARGS goswami_input.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu concentration concentration 1e-1 1e-5
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu pressure pressure 1e-1 1e-5
    ComponentTransport_pcs_0_ts_100_t_150.000000_expected.vtu ComponentTransport_pcs_0_ts_100_t_150.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_150_t_350.000000_expected.vtu ComponentTransport_pcs_0_ts_150_t_350.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_200_t_750.000000_expected.vtu ComponentTransport_pcs_0_ts_200_t_750.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_250_t_1150.000000_expected.vtu ComponentTransport_pcs_0_ts_250_t_1150.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_300_t_1950.000000_expected.vtu ComponentTransport_pcs_0_ts_300_t_1950.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_350_t_2750.000000_expected.vtu ComponentTransport_pcs_0_ts_350_t_2750.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_366_t_3000.000000_expected.vtu ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    ComponentTransport_pcs_0_ts_50_t_50.000000_expected.vtu ComponentTransport_pcs_0_ts_50_t_50.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    VIS ComponentTransport_pcs_0_ts_366_t_3000.000000.vtu
)

AddTest(
    NAME LARGE_2D_ComponentTransport_Elder
    PATH Parabolic/ComponentTransport/elder/
    EXECUTABLE ogs
    EXECUTABLE_ARGS elder.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    elder_pcs_0_ts_0_t_0.000000_reference.vtu elder_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_100_t_26298000.000000_reference.vtu elder_pcs_0_ts_100_t_26298000.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_120_t_31557600.000000_reference.vtu elder_pcs_0_ts_120_t_31557600.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_140_t_36817200.000000_reference.vtu elder_pcs_0_ts_140_t_36817200.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_160_t_42076800.000000_reference.vtu elder_pcs_0_ts_160_t_42076800.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_180_t_47336400.000000_reference.vtu elder_pcs_0_ts_180_t_47336400.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_200_t_52596000.000000_reference.vtu elder_pcs_0_ts_200_t_52596000.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_20_t_5259600.000000_reference.vtu elder_pcs_0_ts_20_t_5259600.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_220_t_57855600.000000_reference.vtu elder_pcs_0_ts_220_t_57855600.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_240_t_63115200.000000_reference.vtu elder_pcs_0_ts_240_t_63115200.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_260_t_68374800.000000_reference.vtu elder_pcs_0_ts_260_t_68374800.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_280_t_73634400.000000_reference.vtu elder_pcs_0_ts_280_t_73634400.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_300_t_78894000.000000_reference.vtu elder_pcs_0_ts_300_t_78894000.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_320_t_84153600.000000_reference.vtu elder_pcs_0_ts_320_t_84153600.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_340_t_89413200.000000_reference.vtu elder_pcs_0_ts_340_t_89413200.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_360_t_94672800.000000_reference.vtu elder_pcs_0_ts_360_t_94672800.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_380_t_99932400.000000_reference.vtu elder_pcs_0_ts_380_t_99932400.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_400_t_105192000.000000_reference.vtu elder_pcs_0_ts_400_t_105192000.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_40_t_10519200.000000_reference.vtu elder_pcs_0_ts_40_t_10519200.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_420_t_110451600.000000_reference.vtu elder_pcs_0_ts_420_t_110451600.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_440_t_115711200.000000_reference.vtu elder_pcs_0_ts_440_t_115711200.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_460_t_120970800.000000_reference.vtu elder_pcs_0_ts_460_t_120970800.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_480_t_126230400.000000_reference.vtu elder_pcs_0_ts_480_t_126230400.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_500_t_131490000.000000_reference.vtu elder_pcs_0_ts_500_t_131490000.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_60_t_15778800.000000_reference.vtu elder_pcs_0_ts_60_t_15778800.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_80_t_21038400.000000_reference.vtu elder_pcs_0_ts_80_t_21038400.000000.vtu pressure pressure 1e-1 1e-5
    elder_pcs_0_ts_0_t_0.000000_reference.vtu elder_pcs_0_ts_0_t_0.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_100_t_26298000.000000_reference.vtu elder_pcs_0_ts_100_t_26298000.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_120_t_31557600.000000_reference.vtu elder_pcs_0_ts_120_t_31557600.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_140_t_36817200.000000_reference.vtu elder_pcs_0_ts_140_t_36817200.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_160_t_42076800.000000_reference.vtu elder_pcs_0_ts_160_t_42076800.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_180_t_47336400.000000_reference.vtu elder_pcs_0_ts_180_t_47336400.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_200_t_52596000.000000_reference.vtu elder_pcs_0_ts_200_t_52596000.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_20_t_5259600.000000_reference.vtu elder_pcs_0_ts_20_t_5259600.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_220_t_57855600.000000_reference.vtu elder_pcs_0_ts_220_t_57855600.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_240_t_63115200.000000_reference.vtu elder_pcs_0_ts_240_t_63115200.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_260_t_68374800.000000_reference.vtu elder_pcs_0_ts_260_t_68374800.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_280_t_73634400.000000_reference.vtu elder_pcs_0_ts_280_t_73634400.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_300_t_78894000.000000_reference.vtu elder_pcs_0_ts_300_t_78894000.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_320_t_84153600.000000_reference.vtu elder_pcs_0_ts_320_t_84153600.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_340_t_89413200.000000_reference.vtu elder_pcs_0_ts_340_t_89413200.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_360_t_94672800.000000_reference.vtu elder_pcs_0_ts_360_t_94672800.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_380_t_99932400.000000_reference.vtu elder_pcs_0_ts_380_t_99932400.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_400_t_105192000.000000_reference.vtu elder_pcs_0_ts_400_t_105192000.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_40_t_10519200.000000_reference.vtu elder_pcs_0_ts_40_t_10519200.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_420_t_110451600.000000_reference.vtu elder_pcs_0_ts_420_t_110451600.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_440_t_115711200.000000_reference.vtu elder_pcs_0_ts_440_t_115711200.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_460_t_120970800.000000_reference.vtu elder_pcs_0_ts_460_t_120970800.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_480_t_126230400.000000_reference.vtu elder_pcs_0_ts_480_t_126230400.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_500_t_131490000.000000_reference.vtu elder_pcs_0_ts_500_t_131490000.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_60_t_15778800.000000_reference.vtu elder_pcs_0_ts_60_t_15778800.000000.vtu conc conc 1e-1 1e-5
    elder_pcs_0_ts_80_t_21038400.000000_reference.vtu elder_pcs_0_ts_80_t_21038400.000000.vtu conc conc 1e-1 1e-5
)

AddTest(
    NAME 2D_ComponentTransport_HeterogeneousPermeability
    PATH Elliptic/square_100x100_ComponentTransport/
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e4_heterogeneity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu pressure_expected pressure 1e-7 1e-10
    square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu darcy_velocity_expected darcy_velocity 1e-7 1e-10
    VIS square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_HeterogeneousPermeability_Comparison_OGS5
    PATH Parabolic/ComponentTransport/heterogeneous/ogs5_H_2D/
    EXECUTABLE ogs
    EXECUTABLE_ARGS ogs5_H_2d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    2D1P-GWFlow_1_reference.vtu out_ogs5_H_pcs_0_ts_1_t_10000000.000000.vtu pressure_OGS5 pressure 1e-1 1e-3
    2D1P-GWFlow_1_reference.vtu out_ogs5_H_pcs_0_ts_1_t_10000000.000000.vtu NODAL_VELOCITY1 darcy_velocity 2e-11 0
    VIS out_ogs5_H_pcs_0_ts_10_t_100000000.000000.vtu
)

AddTest(
    NAME 3D_ComponentTransport_HeterogeneousPermeability_Comparison_OGS5
    PATH Parabolic/ComponentTransport/heterogeneous/ogs5_H_3D/
    EXECUTABLE ogs
    EXECUTABLE_ARGS ogs5_H_3d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D1P-GWFlow_1_reference.vtu out_ogs5_H_pcs_0_ts_1_t_10000000.000000.vtu pressure_ogs5 pressure 2.4e1 1.4e-2
    3D1P-GWFlow_1_reference.vtu out_ogs5_H_pcs_0_ts_1_t_10000000.000000.vtu NODAL_VELOCITY1 darcy_velocity 1e-10 1.4e-2
    VIS out_ogs5_H_pcs_0_ts_1_t_10000000.000000.vtu
)
