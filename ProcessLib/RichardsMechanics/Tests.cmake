if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/gravity.prj)
    OgsTest(PROJECTFILE RichardsMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/confined_compression_fully_saturated.prj RUNTIME 7)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_anisotropic.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_coordinate_system.prj)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small.prj RUNTIME 9)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small_masslumping.prj RUNTIME 10)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_quasinewton.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/double_porosity_swelling.prj RUNTIME 20)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity.prj RUNTIME 8)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity_swelling.prj RUNTIME 11)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_power_law_permeability_xyz.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xyz.prj)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xy.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_power_law.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_saturation_cutoff.prj)
    OgsTest(PROJECTFILE RichardsMechanics/alternative_mass_balance_anzInterval_10.prj)
    if(NOT ENABLE_ASAN)
        OgsTest(PROJECTFILE RichardsMechanics/rotated_consolidation.prj RUNTIME 2)
    endif()
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos.prj RUNTIME 17)
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos_restart.xml RUNTIME 17)
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos_QN.prj RUNTIME 50)
    OgsTest(PROJECTFILE RichardsMechanics/A2.prj RUNTIME 20)
endif()

if (NOT OGS_USE_MPI AND OGS_USE_MFRONT)
    OgsTest(PROJECTFILE RichardsMechanics/mfront_restart_part1.prj RUNTIME 1)
    OgsTest(PROJECTFILE RichardsMechanics/mfront_restart_part2.prj RUNTIME 1)
endif()

AddTest(
    NAME RichardsMechanics_square_1e2_confined_compression_restart
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 8
    EXECUTABLE_ARGS confined_compression_fully_saturated_restart.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    # Does not exist?
    # PROPERTIES DEPENDS ogs-RichardsMechanics_square_1e2_confined_compression-time-vtkdiff
    DIFF_DATA
    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu displacement displacement 1e-16 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu displacement displacement 1e-16 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu displacement displacement 1e-16 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu pressure pressure 1e-16 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu pressure pressure 1e-16 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu pressure pressure 1e-16 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu sigma sigma 5e-14 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu sigma sigma 5e-14 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu sigma sigma 5e-14 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu epsilon epsilon 5e-14 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu epsilon epsilon 5e-14 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu epsilon epsilon 5e-14 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu saturation saturation 2.5e-15 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu saturation saturation 2.5e-15 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu saturation saturation 2.5e-15 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu velocity velocity 1e-16 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu velocity velocity 1e-16 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu velocity velocity 1e-16 0
)
