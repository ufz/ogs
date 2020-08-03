AddTest(
    NAME 2D_ComponentTransport_ConcentrationDiffusionOnly
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionOnly.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu linear_1_to_0 Si 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero pressure 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu zero_vector_2d darcy_velocity 1e-7 1e-10
    VIS DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_StaggeredScheme_ConcentrationDiffusionOnly
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionOnly.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionOnly_pcs_1_ts_1_t_1.000000_expected.vtu DiffusionOnly_ts_1_t_1.000000.vtu Si Si 1e-7 1e-10
    DiffusionOnly_pcs_1_ts_1_t_1.000000_expected.vtu DiffusionOnly_ts_1_t_1.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionOnly_pcs_1_ts_1_t_1.000000_expected.vtu DiffusionOnly_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionOnly_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_MultiComponentTransport_ConcentrationDiffusionOnly
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionOnly_3Components.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_3Components_pcs_0_ts_1_t_1.000000.vtu linear_1_to_0 Si 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_3Components_pcs_0_ts_1_t_1.000000.vtu linear_1_to_0 Al 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_3Components_pcs_0_ts_1_t_1.000000.vtu linear_1_to_0 Cl 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_3Components_pcs_0_ts_1_t_1.000000.vtu zero pressure 1e-7 1e-10
    DiffusionOnly_pcs_0_ts_1_t_1.000000_expected.vtu DiffusionOnly_3Components_pcs_0_ts_1_t_1.000000.vtu zero_vector_2d darcy_velocity 1e-7 1e-10
    VIS DiffusionOnly_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_ConcentrationDiffusionAndStorage
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionAndStorage.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu concentration Si 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_pcs_0_ts_100_t_0.150000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu concentration Si 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero pressure 1e-7 1e-10
    DiffusionAndStorage_pcs_0_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu zero_vector_2d darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorage_pcs_0_ts_134_t_1.500000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_StaggeredScheme_ConcentrationDiffusionAndStorage
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConcentrationDiffusionAndStorage.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    DiffusionAndStorage_pcs_1_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_ts_100_t_0.150000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorage_pcs_1_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_ts_100_t_0.150000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorage_pcs_1_ts_100_t_0.150000_expected.vtu DiffusionAndStorage_ts_100_t_0.150000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorage_pcs_1_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_ts_134_t_1.500000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorage_pcs_1_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_ts_134_t_1.500000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorage_pcs_1_ts_134_t_1.500000_expected.vtu DiffusionAndStorage_ts_134_t_1.500000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorage_ts_134_t_1.500000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvection
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvection.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 23
    DIFF_DATA
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_StaggeredScheme_DiffusionAndStorageAndAdvection
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvection.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 23
    DIFF_DATA
    DiffusionAndStorageAndAdvection_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvection_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvection_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvection_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_ImpermeableBoundaries
    PATH Parabolic/ComponentTransport/MassConservation
    EXECUTABLE ogs
    EXECUTABLE_ARGS mass_conservation.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 23
    DIFF_DATA
    mass_conservation_ogsOutput_pcs_0_ts_0_t_0.000000_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_0_t_0.000000.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_300_t_34895.986246_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_300_t_34895.986246.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_600_t_81993.310506_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_600_t_81993.310506.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_900_t_145558.519328_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_900_t_145558.519328.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1200_t_231349.715241_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1200_t_231349.715241.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1500_t_347138.358629_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1500_t_347138.358629.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1800_t_503413.251350_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1800_t_503413.251350.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_2100_t_714330.672785_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_2100_t_714330.672785.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000.vtu concentration concentration 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_0_t_0.000000_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_300_t_34895.986246_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_300_t_34895.986246.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_600_t_81993.310506_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_600_t_81993.310506.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_900_t_145558.519328_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_900_t_145558.519328.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1200_t_231349.715241_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1200_t_231349.715241.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1500_t_347138.358629_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1500_t_347138.358629.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1800_t_503413.251350_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1800_t_503413.251350.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_2100_t_714330.672785_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_2100_t_714330.672785.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000.vtu pressure pressure 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_0_t_0.000000_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_0_t_0.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_300_t_34895.986246_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_300_t_34895.986246.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_600_t_81993.310506_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_600_t_81993.310506.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_900_t_145558.519328_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_900_t_145558.519328.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1200_t_231349.715241_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1200_t_231349.715241.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1500_t_347138.358629_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1500_t_347138.358629.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_1800_t_503413.251350_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_1800_t_503413.251350.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_2100_t_714330.672785_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_2100_t_714330.672785.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000_expected.vtu mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS mass_conservation_ogsOutput_pcs_0_ts_2323_t_1000000.000000.vtu
)

AddTest(
    NAME 3D_ComponentTransport_NonAdvective_OpenBoundary
    PATH Parabolic/ComponentTransport/OpenBoundaryWithTets
    EXECUTABLE ogs
    EXECUTABLE_ARGS box_flow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 4
    DIFF_DATA
    box_ogsOutput_pcs_0_ts_0_t_0.000000_expected.vtu box_ogsOutput_pcs_0_ts_0_t_0.000000.vtu concentration concentration 5e-7 5e-10
    box_ogsOutput_pcs_0_ts_20_t_100000.000000_expected.vtu box_ogsOutput_pcs_0_ts_20_t_100000.000000.vtu concentration concentration 5e-7 5e-10
    box_ogsOutput_pcs_0_ts_0_t_0.000000_expected.vtu box_ogsOutput_pcs_0_ts_0_t_0.000000.vtu pressure pressure 5e-7 5e-10
    box_ogsOutput_pcs_0_ts_20_t_100000.000000_expected.vtu box_ogsOutput_pcs_0_ts_20_t_100000.000000.vtu pressure pressure 5e-7 5e-10
    box_ogsOutput_pcs_0_ts_0_t_0.000000_expected.vtu box_ogsOutput_pcs_0_ts_0_t_0.000000.vtu darcy_velocity darcy_velocity 5e-7 5e-10
    box_ogsOutput_pcs_0_ts_20_t_100000.000000_expected.vtu box_ogsOutput_pcs_0_ts_20_t_100000.000000.vtu darcy_velocity darcy_velocity 5e-7 5e-10
    VIS box_ogsOutput_pcs_0_ts_20_t_100000.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndGravityAndDispersionHalf
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndGravityAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 46
    DIFF_DATA
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1000_t_2500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1000_t_2500.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1100_t_5000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1100_t_5000.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1200_t_7500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1200_t_7500.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1300_t_10000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1300_t_10000.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1400_t_12500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1400_t_12500.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1000_t_2500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1000_t_2500.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1100_t_5000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1100_t_5000.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1200_t_7500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1200_t_7500.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1300_t_10000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1300_t_10000.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1400_t_12500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1400_t_12500.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1000_t_2500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1000_t_2500.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1100_t_5000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1100_t_5000.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1200_t_7500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1200_t_7500.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1300_t_10000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1300_t_10000.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1400_t_12500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1400_t_12500.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    VIS DiffusionAndStorageAndGravityAndDispersionHalf_pcs_0_ts_1500_t_15000.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_StaggeredScheme_DiffusionAndStorageAndGravityAndDispersionHalf
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndGravityAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 46
    DIFF_DATA
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1000_t_2500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1000_t_2500.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1100_t_5000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1100_t_5000.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1200_t_7500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1200_t_7500.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1300_t_10000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1300_t_10000.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1400_t_12500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1400_t_12500.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1500_t_15000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1500_t_15000.000000.vtu Si Si 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1000_t_2500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1000_t_2500.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1100_t_5000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1100_t_5000.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1200_t_7500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1200_t_7500.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1300_t_10000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1300_t_10000.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1400_t_12500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1400_t_12500.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1500_t_15000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1500_t_15000.000000.vtu pressure pressure 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1000_t_2500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1000_t_2500.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1100_t_5000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1100_t_5000.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1200_t_7500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1200_t_7500.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1300_t_10000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1300_t_10000.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1400_t_12500.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1400_t_12500.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    DiffusionAndStorageAndGravityAndDispersionHalf_pcs_1_ts_1500_t_15000.000000_expected.vtu DiffusionAndStorageAndGravityAndDispersionHalf_ts_1500_t_15000.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-10
    VIS DiffusionAndStorageAndGravityAndDispersionHalf_ts_1500_t_15000.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDispersion
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 26
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_pcs_0_ts_672_t_900.000000.vtu pressure pressure 2e-7 1e-10
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
    NAME 2D_ComponentTransport_OpenBoundaryBC
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS open_boundary_component-transport_cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 1
    DIFF_DATA
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_0_t_0.000000_expected.vtu OpenBoundaryBC_pcs_0_ts_0_t_0.000000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_1_t_0.020000_expected.vtu OpenBoundaryBC_pcs_0_ts_1_t_0.020000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_2_t_0.040000_expected.vtu OpenBoundaryBC_pcs_0_ts_2_t_0.040000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_3_t_0.060000_expected.vtu OpenBoundaryBC_pcs_0_ts_3_t_0.060000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_4_t_0.080000_expected.vtu OpenBoundaryBC_pcs_0_ts_4_t_0.080000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_5_t_0.100000_expected.vtu OpenBoundaryBC_pcs_0_ts_5_t_0.100000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_0_t_0.000000_expected.vtu OpenBoundaryBC_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_1_t_0.020000_expected.vtu OpenBoundaryBC_pcs_0_ts_1_t_0.020000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_2_t_0.040000_expected.vtu OpenBoundaryBC_pcs_0_ts_2_t_0.040000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_3_t_0.060000_expected.vtu OpenBoundaryBC_pcs_0_ts_3_t_0.060000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_4_t_0.080000_expected.vtu OpenBoundaryBC_pcs_0_ts_4_t_0.080000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_5_t_0.100000_expected.vtu OpenBoundaryBC_pcs_0_ts_5_t_0.100000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_0_t_0.000000_expected.vtu OpenBoundaryBC_pcs_0_ts_0_t_0.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_1_t_0.020000_expected.vtu OpenBoundaryBC_pcs_0_ts_1_t_0.020000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_2_t_0.040000_expected.vtu OpenBoundaryBC_pcs_0_ts_2_t_0.040000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_3_t_0.060000_expected.vtu OpenBoundaryBC_pcs_0_ts_3_t_0.060000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_4_t_0.080000_expected.vtu OpenBoundaryBC_pcs_0_ts_4_t_0.080000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_5_t_0.100000_expected.vtu OpenBoundaryBC_pcs_0_ts_5_t_0.100000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS OpenBoundaryBC_pcs_0_ts_5_t_0.100000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_Advective_and_NonAdvective_comparison
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS open_boundary_component-transport_cube_1e3_advective_form.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 1
    DIFF_DATA
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_0_t_0.000000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_0_t_0.000000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_1_t_0.020000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_1_t_0.020000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_2_t_0.040000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_2_t_0.040000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_3_t_0.060000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_3_t_0.060000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_4_t_0.080000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_4_t_0.080000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_5_t_0.100000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_5_t_0.100000.vtu Si Si 5e-6 5e-6
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_0_t_0.000000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_1_t_0.020000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_1_t_0.020000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_2_t_0.040000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_2_t_0.040000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_3_t_0.060000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_3_t_0.060000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_4_t_0.080000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_4_t_0.080000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_5_t_0.100000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_5_t_0.100000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_0_t_0.000000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_0_t_0.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_1_t_0.020000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_1_t_0.020000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_2_t_0.040000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_2_t_0.040000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_3_t_0.060000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_3_t_0.060000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_4_t_0.080000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_4_t_0.080000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndAdvection_surfaceflux_pcs_0_ts_5_t_0.100000_expected.vtu AdvectiveNonAdvectiveComparison_pcs_0_ts_5_t_0.100000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS AdvectiveNonAdvectiveComparison_pcs_0_ts_5_t_0.100000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_StaggeredScheme_DiffusionAndStorageAndAdvectionAndDispersion
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 26
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDispersion_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_MultiComponentTransport_StaggeredScheme_DiffusionAndStorageAndAdvectionAndDispersion
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersion_3Components.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 26
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_100_t_5.700000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_200_t_35.700000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_300_t_155.700000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_400_t_315.700000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_500_t_495.700000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_600_t_720.700000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_672_t_900.000000.vtu Al Al 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_100_t_5.700000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_200_t_35.700000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_300_t_155.700000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_400_t_315.700000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_500_t_495.700000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_600_t_720.700000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_672_t_900.000000.vtu Cl Cl 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersion_3Components_pcs_3_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDispersion_3Components_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDecay
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDecay.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 23
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_pcs_0_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
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
    NAME 2D_ComponentTransport_StaggeredScheme_DiffusionAndStorageAndAdvectionAndDecay
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDecay.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 23
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDecay_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDecay_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDecay_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_DiffusionAndStorageAndAdvectionAndDispersionHalf
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 25
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_0_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
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
    NAME 2D_ComponentTransport_StaggeredScheme_DiffusionAndStorageAndAdvectionAndDispersionHalf
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS DiffusionAndStorageAndAdvectionAndDispersionHalf.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 25
    DIFF_DATA
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_100_t_5.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_200_t_35.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_300_t_155.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_400_t_315.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_500_t_495.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_600_t_720.700000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_672_t_900.000000.vtu Si Si 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_100_t_5.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_200_t_35.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_300_t_155.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_400_t_315.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_500_t_495.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_600_t_720.700000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_672_t_900.000000.vtu pressure pressure 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_100_t_5.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_100_t_5.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_200_t_35.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_200_t_35.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_300_t_155.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_300_t_155.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_400_t_315.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_400_t_315.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_500_t_495.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_500_t_495.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_600_t_720.700000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_600_t_720.700000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    DiffusionAndStorageAndAdvectionAndDispersionHalf_pcs_1_ts_672_t_900.000000_expected.vtu DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_672_t_900.000000.vtu darcy_velocity darcy_velocity 1e-7 1e-10
    VIS DiffusionAndStorageAndAdvectionAndDispersionHalf_ts_672_t_900.000000.vtu
)

AddTest(
    NAME 3D_ComponentTransport_surfaceflux
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS surfaceflux_component-transport_cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_complete_surface_expected_specificflux.vtu flux_1e3_t_1.000000.vtu specific_flux specific_flux 1e-10 1e-16
    cube_1x1x1_hex_1e3_left_pcs_0_ts_0_t_0.000000.vtu cube_1x1x1_hex_1e3_left_pcs_0_ts_0_t_0.000000.vtu concentration Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_left_pcs_0_ts_1_t_1.000000.vtu cube_1x1x1_hex_1e3_left_pcs_0_ts_1_t_1.000000.vtu concentration Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_left_pcs_0_ts_0_t_0.000000.vtu cube_1x1x1_hex_1e3_left_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-10 1e-16
    cube_1x1x1_hex_1e3_left_pcs_0_ts_1_t_1.000000.vtu cube_1x1x1_hex_1e3_left_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-10 1e-16
    cube_1x1x1_hex_1e3_right_pcs_0_ts_0_t_0.000000.vtu cube_1x1x1_hex_1e3_right_pcs_0_ts_0_t_0.000000.vtu concentration Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_right_pcs_0_ts_1_t_1.000000.vtu cube_1x1x1_hex_1e3_right_pcs_0_ts_1_t_1.000000.vtu concentration Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_right_pcs_0_ts_0_t_0.000000.vtu cube_1x1x1_hex_1e3_right_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-10 1e-16
    cube_1x1x1_hex_1e3_right_pcs_0_ts_1_t_1.000000.vtu cube_1x1x1_hex_1e3_right_pcs_0_ts_1_t_1.000000.vtu pressure pressure 1e-10 1e-16
)

AddTest(
    NAME 3D_ComponentTransport_surfaceflux_pvd
    PATH Parabolic/ComponentTransport/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS surfaceflux_component-transport_cube_1e3.prj
    WRAPPER time
    TESTER diff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1x1x1_hex_1e3_right_pcs_0.pvd
    cube_1x1x1_hex_1e3_left_pcs_0.pvd
    cube_1x1x1_hex_1e3_pcs_0.pvd
)

AddTest(
    NAME 3D_StaggeredScheme_ComponentTransport_surfaceflux
    PATH Parabolic/ComponentTransport/StaggeredScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS surfaceflux_component-transport_cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flux_1e3_t_1.000000.vtu flux_1e3_t_1.000000.vtu specific_flux specific_flux 1e-10 1e-16
    cube_1x1x1_hex_1e3_left_pcs_1_ts_0_t_0.000000.vtu cube_1x1x1_hex_1e3_left_ts_0_t_0.000000.vtu Si Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_left_pcs_1_ts_1_t_1.000000.vtu cube_1x1x1_hex_1e3_left_ts_1_t_1.000000.vtu Si Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_right_pcs_1_ts_0_t_0.000000.vtu cube_1x1x1_hex_1e3_right_ts_0_t_0.000000.vtu Si Si 1e-10 1e-16
    cube_1x1x1_hex_1e3_right_pcs_1_ts_1_t_1.000000.vtu cube_1x1x1_hex_1e3_right_ts_1_t_1.000000.vtu Si Si 1e-10 1e-16
)

AddTest(
    NAME 2D_StaggeredScheme_ComponentTransport_TracerSimulation
    PATH Parabolic/ComponentTransport/TracerSimulation
    EXECUTABLE ogs
    EXECUTABLE_ARGS TracerSimulation.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    TracerSimulation_pcs_1_ts_20_t_20000.000000_expected_ogs5.vtu TracerSimulation_ts_20_t_20000.000000.vtu Cs Cs 4e-7 1e-10
    TracerSimulation_pcs_1_ts_40_t_40000.000000_expected_ogs5.vtu TracerSimulation_ts_40_t_40000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_60_t_60000.000000_expected_ogs5.vtu TracerSimulation_ts_60_t_60000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_80_t_80000.000000_expected_ogs5.vtu TracerSimulation_ts_80_t_80000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_100_t_100000.000000_expected_ogs5.vtu TracerSimulation_ts_100_t_100000.000000.vtu Cs Cs 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_20_t_20000.000000_expected_ogs5.vtu TracerSimulation_ts_20_t_20000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_40_t_40000.000000_expected_ogs5.vtu TracerSimulation_ts_40_t_40000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_60_t_60000.000000_expected_ogs5.vtu TracerSimulation_ts_60_t_60000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_80_t_80000.000000_expected_ogs5.vtu TracerSimulation_ts_80_t_80000.000000.vtu pressure pressure 1.3e-7 1e-10
    TracerSimulation_pcs_1_ts_100_t_100000.000000_expected_ogs5.vtu TracerSimulation_ts_100_t_100000.000000.vtu pressure pressure 1.3e-7 1e-10
)

AddTest(
    NAME Parallel_2D_StaggeredScheme_ComponentTransport_TracerSimulation
    PATH Parabolic/ComponentTransport/TracerSimulation
    EXECUTABLE_ARGS TracerSimulation.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    TracerSimulation_pcs_1_ts_100_t_100000_000000_0_expected.vtu TracerSimulation_ts_100_t_100000_000000_0.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_pcs_1_ts_100_t_100000_000000_0_expected.vtu TracerSimulation_ts_100_t_100000_000000_0.vtu pressure pressure 2.5e-5 1.7e-9
    TracerSimulation_pcs_1_ts_100_t_100000_000000_1_expected.vtu TracerSimulation_ts_100_t_100000_000000_1.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_pcs_1_ts_100_t_100000_000000_1_expected.vtu TracerSimulation_ts_100_t_100000_000000_1.vtu pressure pressure 2.5e-5 1.7e-9
    TracerSimulation_pcs_1_ts_100_t_100000_000000_2_expected.vtu TracerSimulation_ts_100_t_100000_000000_2.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_pcs_1_ts_100_t_100000_000000_2_expected.vtu TracerSimulation_ts_100_t_100000_000000_2.vtu pressure pressure 2.5e-5 1.7e-9
    TracerSimulation_pcs_1_ts_100_t_100000_000000_3_expected.vtu TracerSimulation_ts_100_t_100000_000000_3.vtu Cs Cs 1e-10 1e-16
    TracerSimulation_pcs_1_ts_100_t_100000_000000_3_expected.vtu TracerSimulation_ts_100_t_100000_000000_3.vtu pressure pressure 2.5e-5 1.7e-9
)

AddTest(
    NAME LARGE_2D_ComponentTransport_Goswami
    PATH Parabolic/ComponentTransport/goswami
    RUNTIME 900
    EXECUTABLE ogs
    EXECUTABLE_ARGS goswami_input.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    Goswami_Component_Transport_pcs_0_ts_1185_t_600.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_1185_t_600.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_1785_t_1200.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_1785_t_1200.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_2385_t_1800.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_2385_t_1800.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_2985_t_2400.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_2985_t_2400.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_3585_t_3000.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_3585_t_3000.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_4185_t_3600.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_4185_t_3600.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_4785_t_4200.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_4785_t_4200.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_5385_t_4800.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_5385_t_4800.000000.vtu Si Si 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_1185_t_600.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_1185_t_600.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_1785_t_1200.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_1785_t_1200.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_2385_t_1800.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_2385_t_1800.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_2985_t_2400.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_2985_t_2400.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_3585_t_3000.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_3585_t_3000.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_4185_t_3600.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_4185_t_3600.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_4785_t_4200.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_4785_t_4200.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_5385_t_4800.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_5385_t_4800.000000.vtu pressure pressure 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_1185_t_600.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_1185_t_600.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_1785_t_1200.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_1785_t_1200.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_2385_t_1800.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_2385_t_1800.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_2985_t_2400.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_2985_t_2400.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_3585_t_3000.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_3585_t_3000.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_4185_t_3600.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_4185_t_3600.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_4785_t_4200.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_4785_t_4200.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
    Goswami_Component_Transport_pcs_0_ts_5385_t_4800.000000_expected.vtu Goswami_Component_Transport_pcs_0_ts_5385_t_4800.000000.vtu darcy_velocity darcy_velocity 1e-1 1e-5
)

AddTest(
    NAME LARGE_2D_ComponentTransport_Elder
    PATH Parabolic/ComponentTransport/elder
    RUNTIME 2700
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
    elder_pcs_0_ts_0_t_0.000000_reference.vtu elder_pcs_0_ts_0_t_0.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_100_t_26298000.000000_reference.vtu elder_pcs_0_ts_100_t_26298000.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_120_t_31557600.000000_reference.vtu elder_pcs_0_ts_120_t_31557600.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_140_t_36817200.000000_reference.vtu elder_pcs_0_ts_140_t_36817200.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_160_t_42076800.000000_reference.vtu elder_pcs_0_ts_160_t_42076800.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_180_t_47336400.000000_reference.vtu elder_pcs_0_ts_180_t_47336400.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_200_t_52596000.000000_reference.vtu elder_pcs_0_ts_200_t_52596000.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_20_t_5259600.000000_reference.vtu elder_pcs_0_ts_20_t_5259600.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_220_t_57855600.000000_reference.vtu elder_pcs_0_ts_220_t_57855600.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_240_t_63115200.000000_reference.vtu elder_pcs_0_ts_240_t_63115200.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_260_t_68374800.000000_reference.vtu elder_pcs_0_ts_260_t_68374800.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_280_t_73634400.000000_reference.vtu elder_pcs_0_ts_280_t_73634400.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_300_t_78894000.000000_reference.vtu elder_pcs_0_ts_300_t_78894000.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_320_t_84153600.000000_reference.vtu elder_pcs_0_ts_320_t_84153600.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_340_t_89413200.000000_reference.vtu elder_pcs_0_ts_340_t_89413200.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_360_t_94672800.000000_reference.vtu elder_pcs_0_ts_360_t_94672800.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_380_t_99932400.000000_reference.vtu elder_pcs_0_ts_380_t_99932400.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_400_t_105192000.000000_reference.vtu elder_pcs_0_ts_400_t_105192000.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_40_t_10519200.000000_reference.vtu elder_pcs_0_ts_40_t_10519200.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_420_t_110451600.000000_reference.vtu elder_pcs_0_ts_420_t_110451600.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_440_t_115711200.000000_reference.vtu elder_pcs_0_ts_440_t_115711200.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_460_t_120970800.000000_reference.vtu elder_pcs_0_ts_460_t_120970800.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_480_t_126230400.000000_reference.vtu elder_pcs_0_ts_480_t_126230400.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_500_t_131490000.000000_reference.vtu elder_pcs_0_ts_500_t_131490000.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_60_t_15778800.000000_reference.vtu elder_pcs_0_ts_60_t_15778800.000000.vtu Si Si 1e-1 1e-5
    elder_pcs_0_ts_80_t_21038400.000000_reference.vtu elder_pcs_0_ts_80_t_21038400.000000.vtu Si Si 1e-1 1e-5
)

AddTest(
    NAME LARGE_2D_ComponentTransport_ElderPython
    PATH Parabolic/ComponentTransport/elder
    EXECUTABLE ogs
    EXECUTABLE_ARGS elder-python.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_PYTHON AND NOT OGS_USE_MPI
    DIFF_DATA
    elder_pcs_0_ts_0_t_0.000000_reference.vtu            elder_python_pcs_0_ts_0_t_0.000000.vtu            pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_100_t_26298000.000000_reference.vtu   elder_python_pcs_0_ts_100_t_26298000.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_120_t_31557600.000000_reference.vtu   elder_python_pcs_0_ts_120_t_31557600.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_140_t_36817200.000000_reference.vtu   elder_python_pcs_0_ts_140_t_36817200.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_160_t_42076800.000000_reference.vtu   elder_python_pcs_0_ts_160_t_42076800.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_180_t_47336400.000000_reference.vtu   elder_python_pcs_0_ts_180_t_47336400.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_200_t_52596000.000000_reference.vtu   elder_python_pcs_0_ts_200_t_52596000.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_20_t_5259600.000000_reference.vtu     elder_python_pcs_0_ts_20_t_5259600.000000.vtu     pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_220_t_57855600.000000_reference.vtu   elder_python_pcs_0_ts_220_t_57855600.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_240_t_63115200.000000_reference.vtu   elder_python_pcs_0_ts_240_t_63115200.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_260_t_68374800.000000_reference.vtu   elder_python_pcs_0_ts_260_t_68374800.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_280_t_73634400.000000_reference.vtu   elder_python_pcs_0_ts_280_t_73634400.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_300_t_78894000.000000_reference.vtu   elder_python_pcs_0_ts_300_t_78894000.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_320_t_84153600.000000_reference.vtu   elder_python_pcs_0_ts_320_t_84153600.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_340_t_89413200.000000_reference.vtu   elder_python_pcs_0_ts_340_t_89413200.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_360_t_94672800.000000_reference.vtu   elder_python_pcs_0_ts_360_t_94672800.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_380_t_99932400.000000_reference.vtu   elder_python_pcs_0_ts_380_t_99932400.000000.vtu   pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_400_t_105192000.000000_reference.vtu  elder_python_pcs_0_ts_400_t_105192000.000000.vtu  pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_40_t_10519200.000000_reference.vtu    elder_python_pcs_0_ts_40_t_10519200.000000.vtu    pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_420_t_110451600.000000_reference.vtu  elder_python_pcs_0_ts_420_t_110451600.000000.vtu  pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_440_t_115711200.000000_reference.vtu  elder_python_pcs_0_ts_440_t_115711200.000000.vtu  pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_460_t_120970800.000000_reference.vtu  elder_python_pcs_0_ts_460_t_120970800.000000.vtu  pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_480_t_126230400.000000_reference.vtu  elder_python_pcs_0_ts_480_t_126230400.000000.vtu  pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_500_t_131490000.000000_reference.vtu  elder_python_pcs_0_ts_500_t_131490000.000000.vtu  pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_60_t_15778800.000000_reference.vtu    elder_python_pcs_0_ts_60_t_15778800.000000.vtu    pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_80_t_21038400.000000_reference.vtu    elder_python_pcs_0_ts_80_t_21038400.000000.vtu    pressure  pressure  1e-1  1e-5
    elder_pcs_0_ts_0_t_0.000000_reference.vtu            elder_python_pcs_0_ts_0_t_0.000000.vtu            conc      Si      1e-1  1e-5
    elder_pcs_0_ts_100_t_26298000.000000_reference.vtu   elder_python_pcs_0_ts_100_t_26298000.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_120_t_31557600.000000_reference.vtu   elder_python_pcs_0_ts_120_t_31557600.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_140_t_36817200.000000_reference.vtu   elder_python_pcs_0_ts_140_t_36817200.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_160_t_42076800.000000_reference.vtu   elder_python_pcs_0_ts_160_t_42076800.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_180_t_47336400.000000_reference.vtu   elder_python_pcs_0_ts_180_t_47336400.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_200_t_52596000.000000_reference.vtu   elder_python_pcs_0_ts_200_t_52596000.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_20_t_5259600.000000_reference.vtu     elder_python_pcs_0_ts_20_t_5259600.000000.vtu     conc      Si      1e-1  1e-5
    elder_pcs_0_ts_220_t_57855600.000000_reference.vtu   elder_python_pcs_0_ts_220_t_57855600.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_240_t_63115200.000000_reference.vtu   elder_python_pcs_0_ts_240_t_63115200.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_260_t_68374800.000000_reference.vtu   elder_python_pcs_0_ts_260_t_68374800.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_280_t_73634400.000000_reference.vtu   elder_python_pcs_0_ts_280_t_73634400.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_300_t_78894000.000000_reference.vtu   elder_python_pcs_0_ts_300_t_78894000.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_320_t_84153600.000000_reference.vtu   elder_python_pcs_0_ts_320_t_84153600.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_340_t_89413200.000000_reference.vtu   elder_python_pcs_0_ts_340_t_89413200.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_360_t_94672800.000000_reference.vtu   elder_python_pcs_0_ts_360_t_94672800.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_380_t_99932400.000000_reference.vtu   elder_python_pcs_0_ts_380_t_99932400.000000.vtu   conc      Si      1e-1  1e-5
    elder_pcs_0_ts_400_t_105192000.000000_reference.vtu  elder_python_pcs_0_ts_400_t_105192000.000000.vtu  conc      Si      1e-1  1e-5
    elder_pcs_0_ts_40_t_10519200.000000_reference.vtu    elder_python_pcs_0_ts_40_t_10519200.000000.vtu    conc      Si      1e-1  1e-5
    elder_pcs_0_ts_420_t_110451600.000000_reference.vtu  elder_python_pcs_0_ts_420_t_110451600.000000.vtu  conc      Si      1e-1  1e-5
    elder_pcs_0_ts_440_t_115711200.000000_reference.vtu  elder_python_pcs_0_ts_440_t_115711200.000000.vtu  conc      Si      1e-1  1e-5
    elder_pcs_0_ts_460_t_120970800.000000_reference.vtu  elder_python_pcs_0_ts_460_t_120970800.000000.vtu  conc      Si      1e-1  1e-5
    elder_pcs_0_ts_480_t_126230400.000000_reference.vtu  elder_python_pcs_0_ts_480_t_126230400.000000.vtu  conc      Si      1e-1  1e-5
    elder_pcs_0_ts_500_t_131490000.000000_reference.vtu  elder_python_pcs_0_ts_500_t_131490000.000000.vtu  conc      Si      1e-1  1e-5
    elder_pcs_0_ts_60_t_15778800.000000_reference.vtu    elder_python_pcs_0_ts_60_t_15778800.000000.vtu    conc      Si      1e-1  1e-5
    elder_pcs_0_ts_80_t_21038400.000000_reference.vtu    elder_python_pcs_0_ts_80_t_21038400.000000.vtu    conc      Si      1e-1  1e-5
)

AddTest(
    NAME 2D_ComponentTransport_HeterogeneousPermeability
    PATH Elliptic/square_100x100_ComponentTransport
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e4_heterogeneity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu pressure_expected pressure 2e-2 1e-10
    square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu darcy_velocity_expected darcy_velocity 1e-7 1e-10
    VIS square_100x100_quad_1e4_ComponentTransport_pcs_0_ts_1_t_1.000000.vtu
)

AddTest(
    NAME 2D_ComponentTransport_HeterogeneousPermeability_Comparison_OGS5
    PATH Parabolic/ComponentTransport/heterogeneous/ogs5_H_2D
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
    PATH Parabolic/ComponentTransport/heterogeneous/ogs5_H_3D
    EXECUTABLE ogs
    EXECUTABLE_ARGS ogs5_H_3d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 15
    DIFF_DATA
    3D1P-GWFlow_1_reference.vtu out_ogs5_H_pcs_0_ts_10_t_10000000.000000.vtu pressure_ogs5 pressure 2.4e1 1.4e-2
    3D1P-GWFlow_1_reference.vtu out_ogs5_H_pcs_0_ts_10_t_10000000.000000.vtu NODAL_VELOCITY1 darcy_velocity 1e-10 1.4e-2
    VIS out_ogs5_H_pcs_0_ts_1_t_10000000.000000.vtu
)


AddTest(
    NAME 1D_ComponentTransport_VariableDependentBoundary
    PATH Parabolic/ComponentTransport/VariableNeumannBoundary
    EXECUTABLE ogs
    EXECUTABLE_ARGS vdbc_input.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    vdbc_pcs_0_ts_0_t_0.000000_expected.vtu vdbc_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_1590_t_6000.000000_expected.vtu vdbc_pcs_0_ts_1590_t_6000.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_3990_t_30000.000000_expected.vtu vdbc_pcs_0_ts_3990_t_30000.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_9990_t_90000.000000_expected.vtu vdbc_pcs_0_ts_9990_t_90000.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_15990_t_150000.000000_expected.vtu vdbc_pcs_0_ts_15990_t_150000.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_21990_t_210000.000000_expected.vtu vdbc_pcs_0_ts_21990_t_210000.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_25990_t_250000.000000_expected.vtu vdbc_pcs_0_ts_25990_t_250000.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_0_t_0.000000_expected.vtu vdbc_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-5 1e-4
    vdbc_pcs_0_ts_1590_t_6000.000000_expected.vtu vdbc_pcs_0_ts_1590_t_6000.000000.vtu concentration Si 1e-5 1e-4
    vdbc_pcs_0_ts_3990_t_30000.000000_expected.vtu vdbc_pcs_0_ts_3990_t_30000.000000.vtu concentration Si 1e-5 1e-4
    vdbc_pcs_0_ts_9990_t_90000.000000_expected.vtu vdbc_pcs_0_ts_9990_t_90000.000000.vtu concentration Si 1e-5 1e-4
    vdbc_pcs_0_ts_15990_t_150000.000000_expected.vtu vdbc_pcs_0_ts_15990_t_150000.000000.vtu concentration Si 1e-5 1e-4
    vdbc_pcs_0_ts_21990_t_210000.000000_expected.vtu vdbc_pcs_0_ts_21990_t_210000.000000.vtu concentration Si 1e-5 1e-4
    vdbc_pcs_0_ts_25990_t_250000.000000_expected.vtu vdbc_pcs_0_ts_25990_t_250000.000000.vtu concentration Si 1e-5 1e-4
)

AddTest(
    NAME ComponentTransport_Theis
    PATH Parabolic/ComponentTransport/Theis
    EXECUTABLE ogs
    EXECUTABLE_ARGS theis.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 90
    DIFF_DATA
    theis_pcs_0_ts_0_t_0.000000.vtu
    theis_pcs_0_ts_0_t_0.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_0_t_0.000000.vtu
    theis_pcs_0_ts_0_t_0.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_10_t_100.000000.vtu
    theis_pcs_0_ts_10_t_100.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_10_t_100.000000.vtu
    theis_pcs_0_ts_10_t_100.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_20_t_400.000000.vtu
    theis_pcs_0_ts_20_t_400.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_20_t_400.000000.vtu
    theis_pcs_0_ts_20_t_400.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_30_t_1000.000000.vtu
    theis_pcs_0_ts_30_t_1000.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_30_t_1000.000000.vtu
    theis_pcs_0_ts_30_t_1000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_40_t_2000.000000.vtu
    theis_pcs_0_ts_40_t_2000.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_40_t_2000.000000.vtu
    theis_pcs_0_ts_40_t_2000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_60_t_20000.000000.vtu
    theis_pcs_0_ts_60_t_20000.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_60_t_20000.000000.vtu
    theis_pcs_0_ts_60_t_20000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_70_t_70000.000000.vtu
    theis_pcs_0_ts_70_t_70000.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_70_t_70000.000000.vtu
    theis_pcs_0_ts_70_t_70000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
    theis_pcs_0_ts_73_t_100000.000000.vtu
    theis_pcs_0_ts_73_t_100000.000000.vtu pressure pressure 1e-3 1e-6
    theis_pcs_0_ts_73_t_100000.000000.vtu
    theis_pcs_0_ts_73_t_100000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-6
)

AddTest(
    NAME ComponentTransport_ConTracer_1d
    PATH Parabolic/ComponentTransport/ConTracer
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConTracer_1d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 10
    DIFF_DATA
    ConTracer_1d_pcs_1_ts_33_t_118800.000000_expected.vtu ConTracer_1d_ts_33_t_118800.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_66_t_237600.000000_expected.vtu ConTracer_1d_ts_66_t_237600.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_99_t_356400.000000_expected.vtu ConTracer_1d_ts_99_t_356400.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_132_t_475200.000000_expected.vtu ConTracer_1d_ts_132_t_475200.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_165_t_594000.000000_expected.vtu ConTracer_1d_ts_165_t_594000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_198_t_712800.000000_expected.vtu ConTracer_1d_ts_198_t_712800.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_231_t_831600.000000_expected.vtu ConTracer_1d_ts_231_t_831600.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_264_t_950400.000000_expected.vtu ConTracer_1d_ts_264_t_950400.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_297_t_1069200.000000_expected.vtu ConTracer_1d_ts_297_t_1069200.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_329_t_1184400.000000_expected.vtu ConTracer_1d_ts_329_t_1184400.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_1d_pcs_1_ts_33_t_118800.000000_expected.vtu ConTracer_1d_ts_33_t_118800.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_66_t_237600.000000_expected.vtu ConTracer_1d_ts_66_t_237600.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_99_t_356400.000000_expected.vtu ConTracer_1d_ts_99_t_356400.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_132_t_475200.000000_expected.vtu ConTracer_1d_ts_132_t_475200.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_165_t_594000.000000_expected.vtu ConTracer_1d_ts_165_t_594000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_198_t_712800.000000_expected.vtu ConTracer_1d_ts_198_t_712800.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_231_t_831600.000000_expected.vtu ConTracer_1d_ts_231_t_831600.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_264_t_950400.000000_expected.vtu ConTracer_1d_ts_264_t_950400.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_297_t_1069200.000000_expected.vtu ConTracer_1d_ts_297_t_1069200.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_1d_pcs_1_ts_329_t_1184400.000000_expected.vtu ConTracer_1d_ts_329_t_1184400.000000.vtu Cs Cs 1e-10 1e-16
)

AddTest(
    NAME ComponentTransport_ConTracer_2d
    PATH Parabolic/ComponentTransport/ConTracer
    EXECUTABLE ogs
    EXECUTABLE_ARGS ConTracer_2d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 40
    DIFF_DATA
    ConTracer_2d_pcs_1_ts_30_t_108000.000000_expected.vtu ConTracer_2d_ts_30_t_108000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_60_t_216000.000000_expected.vtu ConTracer_2d_ts_60_t_216000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_90_t_324000.000000_expected.vtu ConTracer_2d_ts_90_t_324000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_120_t_432000.000000_expected.vtu ConTracer_2d_ts_120_t_432000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_150_t_540000.000000_expected.vtu ConTracer_2d_ts_150_t_540000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_180_t_648000.000000_expected.vtu ConTracer_2d_ts_180_t_648000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_210_t_756000.000000_expected.vtu ConTracer_2d_ts_210_t_756000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_240_t_864000.000000_expected.vtu ConTracer_2d_ts_240_t_864000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_270_t_972000.000000_expected.vtu ConTracer_2d_ts_270_t_972000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_300_t_1080000.000000_expected.vtu ConTracer_2d_ts_300_t_1080000.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_329_t_1184400.000000_expected.vtu ConTracer_2d_ts_329_t_1184400.000000.vtu pressure pressure 1e-6 1e-10
    ConTracer_2d_pcs_1_ts_30_t_108000.000000_expected.vtu ConTracer_2d_ts_30_t_108000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_60_t_216000.000000_expected.vtu ConTracer_2d_ts_60_t_216000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_90_t_324000.000000_expected.vtu ConTracer_2d_ts_90_t_324000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_120_t_432000.000000_expected.vtu ConTracer_2d_ts_120_t_432000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_150_t_540000.000000_expected.vtu ConTracer_2d_ts_150_t_540000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_180_t_648000.000000_expected.vtu ConTracer_2d_ts_180_t_648000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_210_t_756000.000000_expected.vtu ConTracer_2d_ts_210_t_756000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_240_t_864000.000000_expected.vtu ConTracer_2d_ts_240_t_864000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_270_t_972000.000000_expected.vtu ConTracer_2d_ts_270_t_972000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_300_t_1080000.000000_expected.vtu ConTracer_2d_ts_300_t_1080000.000000.vtu Cs Cs 1e-10 1e-16
    ConTracer_2d_pcs_1_ts_329_t_1184400.000000_expected.vtu ConTracer_2d_ts_329_t_1184400.000000.vtu Cs Cs 1e-10 1e-16
)

if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/EquilibriumPhase/calcite.prj RUNTIME 25)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant/1d_isofrac.prj RUNTIME 20)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant/1d_isofrac_flag_formula.prj RUNTIME 20)
    OgsTest(PROJECTFILE Parabolic/ComponentTransport/ReactiveTransport/KineticReactant_AllAsComponents/KineticReactant2.prj RUNTIME 60)
endif()

AddTest(
    NAME LARGE_2D_ReactiveMassTransport_Phreeqc_KineticReactantBlockTest_AllAsComponents
    PATH Parabolic/ComponentTransport/ReactiveTransport/KineticReactant_AllAsComponents
    EXECUTABLE ogs
    EXECUTABLE_ARGS KineticReactant2_2d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    KineticReactant2_2d_ts_4_t_400.000000_expected.vtu KineticReactant2_2d_ts_4_t_400.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_8_t_800.000000_expected.vtu KineticReactant2_2d_ts_8_t_800.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_12_t_1200.000000_expected.vtu KineticReactant2_2d_ts_12_t_1200.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_16_t_1600.000000_expected.vtu KineticReactant2_2d_ts_16_t_1600.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_20_t_2000.000000_expected.vtu KineticReactant2_2d_ts_20_t_2000.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_4_t_400.000000_expected.vtu KineticReactant2_2d_ts_4_t_400.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_8_t_800.000000_expected.vtu KineticReactant2_2d_ts_8_t_800.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_12_t_1200.000000_expected.vtu KineticReactant2_2d_ts_12_t_1200.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_16_t_1600.000000_expected.vtu KineticReactant2_2d_ts_16_t_1600.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_20_t_2000.000000_expected.vtu KineticReactant2_2d_ts_20_t_2000.000000.vtu pressure pressure 1e-6 1e-10
    KineticReactant2_2d_ts_4_t_400.000000_expected.vtu KineticReactant2_2d_ts_4_t_400.000000.vtu Synthetica Synthetica 1e-10 1e-16
    KineticReactant2_2d_ts_8_t_800.000000_expected.vtu KineticReactant2_2d_ts_8_t_800.000000.vtu Synthetica Synthetica 1e-10 1e-16
    KineticReactant2_2d_ts_12_t_1200.000000_expected.vtu KineticReactant2_2d_ts_12_t_1200.000000.vtu Synthetica Synthetica 1e-10 1e-16
    KineticReactant2_2d_ts_16_t_1600.000000_expected.vtu KineticReactant2_2d_ts_16_t_1600.000000.vtu Synthetica Synthetica 1e-10 1e-16
    KineticReactant2_2d_ts_20_t_2000.000000_expected.vtu KineticReactant2_2d_ts_20_t_2000.000000.vtu Synthetica Synthetica 1e-10 1e-16
    KineticReactant2_2d_ts_4_t_400.000000_expected.vtu KineticReactant2_2d_ts_4_t_400.000000.vtu Syntheticb Syntheticb 1e-10 1e-16
    KineticReactant2_2d_ts_8_t_800.000000_expected.vtu KineticReactant2_2d_ts_8_t_800.000000.vtu Syntheticb Syntheticb 1e-10 1e-16
    KineticReactant2_2d_ts_12_t_1200.000000_expected.vtu KineticReactant2_2d_ts_12_t_1200.000000.vtu Syntheticb Syntheticb 1e-10 1e-16
    KineticReactant2_2d_ts_16_t_1600.000000_expected.vtu KineticReactant2_2d_ts_16_t_1600.000000.vtu Syntheticb Syntheticb 1e-10 1e-16
    KineticReactant2_2d_ts_20_t_2000.000000_expected.vtu KineticReactant2_2d_ts_20_t_2000.000000.vtu Syntheticb Syntheticb 1e-10 1e-16
    KineticReactant2_2d_ts_4_t_400.000000_expected.vtu KineticReactant2_2d_ts_4_t_400.000000.vtu Productd Productd 1e-10 1e-16
    KineticReactant2_2d_ts_8_t_800.000000_expected.vtu KineticReactant2_2d_ts_8_t_800.000000.vtu Productd Productd 1e-10 1e-16
    KineticReactant2_2d_ts_12_t_1200.000000_expected.vtu KineticReactant2_2d_ts_12_t_1200.000000.vtu Productd Productd 1e-10 1e-16
    KineticReactant2_2d_ts_16_t_1600.000000_expected.vtu KineticReactant2_2d_ts_16_t_1600.000000.vtu Productd Productd 1e-10 1e-16
    KineticReactant2_2d_ts_20_t_2000.000000_expected.vtu KineticReactant2_2d_ts_20_t_2000.000000.vtu Productd Productd 1e-10 1e-16
    KineticReactant2_2d_ts_4_t_400.000000_expected.vtu KineticReactant2_2d_ts_4_t_400.000000.vtu H H 1e-10 1e-16
    KineticReactant2_2d_ts_8_t_800.000000_expected.vtu KineticReactant2_2d_ts_8_t_800.000000.vtu H H 1e-10 1e-16
    KineticReactant2_2d_ts_12_t_1200.000000_expected.vtu KineticReactant2_2d_ts_12_t_1200.000000.vtu H H 1e-10 1e-16
    KineticReactant2_2d_ts_16_t_1600.000000_expected.vtu KineticReactant2_2d_ts_16_t_1600.000000.vtu H H 1e-10 1e-16
    KineticReactant2_2d_ts_20_t_2000.000000_expected.vtu KineticReactant2_2d_ts_20_t_2000.000000.vtu H H 1e-10 1e-16
    RUNTIME 1600
)

AddTest(
    NAME 1D_ReactiveMassTransport_Wetland
    PATH Parabolic/ComponentTransport/ReactiveTransport/Wetland
    EXECUTABLE ogs
    EXECUTABLE_ARGS Wetland_1d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu pressure pressure 1e-6 1e-10
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu H H 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Do Do 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Sf Sf 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Sa Sa 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Sin Sin 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Xs_m Xs_m 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Xi_m Xi_m 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Snh Snh 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Sno Sno 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Sso Sso 1e-10 1e-16
    Wetland_1d_ts_4_t_28800.000000_expected.vtu Wetland_1d_ts_4_t_28800.000000.vtu Sulphide Sulphide 1e-10 1e-16
    RUNTIME 24
)
