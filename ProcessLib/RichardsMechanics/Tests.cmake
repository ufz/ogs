if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/gravity.prj)
    OgsTest(PROJECTFILE RichardsMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/confined_compression_fully_saturated.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_anisotropic.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_coordinate_system.prj)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small.prj)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small_masslumping.prj)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_quasinewton.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity.prj)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity_swelling.prj)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_power_law_permeability_xyz.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xyz.prj)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xy.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_power_law.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_saturation_cutoff.prj)
    OgsTest(PROJECTFILE RichardsMechanics/alternative_mass_balance_anzInterval_10.prj)
endif()



AddTest(
    NAME RichardsMechanics_square_1e2_confined_compression_restart
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS confined_compression_fully_saturated_restart.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    confined_compression_fully_saturated_pcs_0_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_0_t_100.000000.vtu displacement displacement 1e-16 0
    confined_compression_fully_saturated_pcs_0_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_100_t_1000.000000.vtu displacement displacement 1e-16 0
    confined_compression_fully_saturated_pcs_0_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_400_t_4000.000000.vtu displacement displacement 1e-16 0

    confined_compression_fully_saturated_pcs_0_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_0_t_100.000000.vtu pressure pressure 1e-16 0
    confined_compression_fully_saturated_pcs_0_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_100_t_1000.000000.vtu pressure pressure 1e-16 0
    confined_compression_fully_saturated_pcs_0_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_400_t_4000.000000.vtu pressure pressure 1e-16 0

    confined_compression_fully_saturated_pcs_0_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_0_t_100.000000.vtu sigma sigma 5e-14 0
    confined_compression_fully_saturated_pcs_0_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_100_t_1000.000000.vtu sigma sigma 5e-14 0
    confined_compression_fully_saturated_pcs_0_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_400_t_4000.000000.vtu sigma sigma 5e-14 0

    confined_compression_fully_saturated_pcs_0_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_0_t_100.000000.vtu epsilon epsilon 5e-14 0
    confined_compression_fully_saturated_pcs_0_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_100_t_1000.000000.vtu epsilon epsilon 5e-14 0
    confined_compression_fully_saturated_pcs_0_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_400_t_4000.000000.vtu epsilon epsilon 5e-14 0

    confined_compression_fully_saturated_pcs_0_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_0_t_100.000000.vtu saturation saturation 2e-15 0
    confined_compression_fully_saturated_pcs_0_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_100_t_1000.000000.vtu saturation saturation 2e-15 0
    confined_compression_fully_saturated_pcs_0_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_400_t_4000.000000.vtu saturation saturation 2e-15 0

    confined_compression_fully_saturated_pcs_0_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_0_t_100.000000.vtu velocity velocity 1e-16 0
    confined_compression_fully_saturated_pcs_0_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_100_t_1000.000000.vtu velocity velocity 1e-16 0
    confined_compression_fully_saturated_pcs_0_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_pcs_0_ts_400_t_4000.000000.vtu velocity velocity 1e-16 0
)
if(TEST ogs-RichardsMechanics_square_1e2_confined_compression_restart-time)
    set_tests_properties(ogs-RichardsMechanics_square_1e2_confined_compression_restart-time PROPERTIES
        DEPENDS ogs-RichardsMechanics_square_1e2_confined_compression-time-vtkdiff)
endif()

AddTest(
    NAME RichardsMechanics_liakopoulosHM
    PATH RichardsMechanics/LiakopoulosHM
    EXECUTABLE ogs
    EXECUTABLE_ARGS liakopoulos.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    liakopoulos_t_300.vtu liakopoulos_t_300.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_300.vtu liakopoulos_t_300.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_300.vtu liakopoulos_t_300.vtu saturation saturation 1e-10 1e-12
    liakopoulos_t_600.vtu liakopoulos_t_600.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_600.vtu liakopoulos_t_600.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_600.vtu liakopoulos_t_600.vtu saturation saturation 1e-10 1e-12
    liakopoulos_t_7200.vtu liakopoulos_t_7200.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_7200.vtu liakopoulos_t_7200.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_7200.vtu liakopoulos_t_7200.vtu saturation saturation 1e-10 1e-12
)