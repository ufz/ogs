# Monolithic scheme
AddTest(
    NAME 2D_ThermalConvection_constviscosityMonolithic
    PATH Parabolic/HT/ConstViscosity
    RUNTIME 66
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu T T 1e-7 1e-9
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu p p 1e-9 1e-3
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_IsothermalFluidFlow
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS IsothermalFluidFlow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlow_ts_1_t_1.000000.vtu T T 1e-10 1e-16
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlow_ts_1_t_1.000000.vtu p p 1e-10 1e-16
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlow_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS IsothermalFluidFlow_ts_1_t_1.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_PressureDiffusionTemperatureDiffusion
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS PressureDiffusionTemperatureDiffusion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu linear_top2_to_bottom1 T 1e-10 1e-16
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu Linear_1_to_minus1 p 1e-10 1e-16
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_IsothermalFluidFlowWithGravity
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS IsothermalFluidFlowWithGravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravity_ts_1_t_1.000000.vtu T T 1e-10 1e-16
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravity_ts_1_t_1.000000.vtu p p 1e-10 1e-16
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravity_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS IsothermalFluidFlow_ts_1_t_1.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_PressureParabolicTemperatureParabolic
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS PressureParabolicTemperatureParabolic.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolic_ts_1_t_0.100000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolic_ts_1_t_0.100000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolic_ts_1_t_0.100000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolic_ts_2_t_0.200000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolic_ts_2_t_0.200000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolic_ts_2_t_0.200000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolic_ts_3_t_0.300000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolic_ts_3_t_0.300000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolic_ts_3_t_0.300000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolic_ts_4_t_0.400000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolic_ts_4_t_0.400000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolic_ts_4_t_0.400000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolic_ts_5_t_0.500000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolic_ts_5_t_0.500000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolic_ts_5_t_0.500000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolic_ts_6_t_0.600000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolic_ts_6_t_0.600000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolic_ts_6_t_0.600000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolic_ts_7_t_0.700000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolic_ts_7_t_0.700000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolic_ts_7_t_0.700000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolic_ts_8_t_0.800000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolic_ts_8_t_0.800000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolic_ts_8_t_0.800000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolic_ts_9_t_0.900000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolic_ts_9_t_0.900000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolic_ts_9_t_0.900000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolic_ts_10_t_1.000000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolic_ts_10_t_1.000000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolic_ts_10_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
)

AddTest(
    NAME HT_SimpleSynthetics_CoupledPressureParabolicTemperatureParabolic
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS CoupledPressureParabolicTemperatureParabolic.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_1_t_0.100000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_1_t_0.100000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_1_t_0.100000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_2_t_0.200000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_2_t_0.200000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_2_t_0.200000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_3_t_0.300000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_3_t_0.300000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_3_t_0.300000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_4_t_0.400000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_4_t_0.400000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_4_t_0.400000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_5_t_0.500000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_5_t_0.500000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_5_t_0.500000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_6_t_0.600000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_6_t_0.600000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_6_t_0.600000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_7_t_0.700000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_7_t_0.700000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_7_t_0.700000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_8_t_0.800000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_8_t_0.800000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_8_t_0.800000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_9_t_0.900000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_9_t_0.900000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_9_t_0.900000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_10_t_1.000000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_10_t_1.000000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolic_ts_10_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
)

AddTest(
    NAME HT_calculatesurfaceflux
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS calculatesurfaceflux_ht_cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flux_1e3_t_0.000000.vtu ThermalConvection_cube_1x1x1_hex_1e3_complete_surface_ts_1_t_0.000000.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_0.000010.vtu ThermalConvection_cube_1x1x1_hex_1e3_complete_surface_ts_2_t_0.000010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_0.001010.vtu ThermalConvection_cube_1x1x1_hex_1e3_complete_surface_ts_3_t_0.001010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_0.101010.vtu ThermalConvection_cube_1x1x1_hex_1e3_complete_surface_ts_4_t_0.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_1.101010.vtu ThermalConvection_cube_1x1x1_hex_1e3_complete_surface_ts_5_t_1.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_10.000000.vtu ThermalConvection_cube_1x1x1_hex_1e3_complete_surface_ts_6_t_10.000000.vtu specific_flux specific_flux 1e-10 1e-16
)

AddTest(
    NAME HT_calculatesurfaceflux
    PATH Parabolic/HT/SimpleSynthetics
    RUNTIME 190
    EXECUTABLE ogs
    EXECUTABLE_ARGS calculatesurfaceflux_ht_cube_1e4.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flux_1e4_t_0.000000.vtu ThermalConvection_cube_1x1x1_hex_1e4_complete_surface_ts_1_t_0.000000.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_0.000010.vtu ThermalConvection_cube_1x1x1_hex_1e4_complete_surface_ts_2_t_0.000010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_0.001010.vtu ThermalConvection_cube_1x1x1_hex_1e4_complete_surface_ts_3_t_0.001010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_0.101010.vtu ThermalConvection_cube_1x1x1_hex_1e4_complete_surface_ts_4_t_0.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_1.101010.vtu ThermalConvection_cube_1x1x1_hex_1e4_complete_surface_ts_5_t_1.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_10.000000.vtu ThermalConvection_cube_1x1x1_hex_1e4_complete_surface_ts_6_t_10.000000.vtu specific_flux specific_flux 1e-10 1e-16
)

# Staggered scheme
AddTest(
    NAME HT_a_DECOVALEX_THMC_based_Example
    PATH Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample
    EXECUTABLE ogs
    EXECUTABLE_ARGS th_decovalex.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    th_decovalex.vtu th_decovalex_ts_40_t_18.000000.vtu T_ref T 6e-12 1.e-14
    th_decovalex.vtu th_decovalex_ts_40_t_18.000000.vtu p_ref p 1e-7 1.e-14
    VIS th_decovalex_ts_78_t_1000.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_IsothermalFluidFlowStaggered
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS IsothermalFluidFlowStaggered.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlowStaggered_ts_1_t_1.000000.vtu T T 1e-10 1e-10
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlowStaggered_ts_1_t_1.000000.vtu p p 1e-10 1e-10
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlowStaggered_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    VIS IsothermalFluidFlowStaggered_ts_1_t_1.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_PressureDiffusionTemperatureDiffusionStaggered
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS PressureDiffusionTemperatureDiffusionStaggered.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusionStaggered_ts_1_t_1.000000.vtu linear_top2_to_bottom1 T 1e-10 1e-10
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusionStaggered_ts_1_t_1.000000.vtu Linear_1_to_minus1 p 1e-10 1e-10
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusionStaggered_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    VIS PressureDiffusionTemperatureDiffusionStaggered_ts_1_t_1.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_IsothermalFluidFlowWithGravityStaggered
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS IsothermalFluidFlowWithGravityStaggered.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravityStaggered_ts_1_t_1.000000.vtu T T 1e-10 1e-10
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravityStaggered_ts_1_t_1.000000.vtu p p 0.04 1e-10
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravityStaggered_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    VIS IsothermalFluidFlowWithGravityStaggered_ts_1_t_1.000000.vtu
)

AddTest(
    NAME HT_SimpleSynthetics_PressureParabolicTemperatureParabolicStaggered
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS PressureParabolicTemperatureParabolicStaggered.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_1_t_0.100000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_1_t_0.100000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_1_t_0.100000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_2_t_0.200000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_2_t_0.200000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_2_t_0.200000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_3_t_0.300000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_3_t_0.300000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_3_t_0.300000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_4_t_0.400000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_4_t_0.400000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_4_t_0.400000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_5_t_0.500000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_5_t_0.500000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_5_t_0.500000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_6_t_0.600000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_6_t_0.600000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_6_t_0.600000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_7_t_0.700000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_7_t_0.700000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_7_t_0.700000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_8_t_0.800000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_8_t_0.800000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_8_t_0.800000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_9_t_0.900000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_9_t_0.900000.vtu p p 1e-7 1e-7
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_9_t_0.900000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_10_t_1.000000.vtu T T 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_10_t_1.000000.vtu p p 1e-10 1e-10
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolicStaggered_ts_10_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
)

AddTest(
    NAME HT_SimpleSynthetics_CoupledPressureParabolicTemperatureParabolicStaggered
    PATH Parabolic/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS CoupledPressureParabolicTemperatureParabolicStaggered.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_1_t_0.100000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_1_t_0.100000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_1_t_0.100000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_2_t_0.200000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_2_t_0.200000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_2_t_0.200000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_3_t_0.300000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_3_t_0.300000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_3_t_0.300000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_4_t_0.400000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_4_t_0.400000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_4_t_0.400000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_5_t_0.500000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_5_t_0.500000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_5_t_0.500000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_6_t_0.600000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_6_t_0.600000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_6_t_0.600000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_7_t_0.700000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_7_t_0.700000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_7_t_0.700000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_8_t_0.800000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_8_t_0.800000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_8_t_0.800000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_9_t_0.900000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_9_t_0.900000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_9_t_0.900000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_10_t_1.000000.vtu T T 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_10_t_1.000000.vtu p p 1e-10 1e-10
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolicStaggered_ts_10_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
)

# 2019-05-09 TF disable the test until the MPL can deal with parameters as properties
# AddTest(
#     NAME HT_a_DECOVALEX_THMC_based_Example
#     PATH Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample
#     EXECUTABLE ogs
#     EXECUTABLE_ARGS th_decovalex.prj
#     WRAPPER mpirun
#     WRAPPER_ARGS -np 1
#     TESTER vtkdiff
#     REQUIREMENTS OGS_USE_MPI
#     RUNTIME 186
#     DIFF_DATA
#     th_decovalex.vtu th_decovalex_pcs_1_ts_40_t_18_000000_0.vtu T_ref T 1e-10  1.e-10
#     th_decovalex.vtu th_decovalex_pcs_1_ts_40_t_18_000000_0.vtu p_ref p 1e-10  1.e-10
#     VIS th_decovalex_pcs_1_ts_78_t_1000_000000_0.vtu
# )

AddTest(
    NAME HT_FaultedCube_rev0
    PATH Parabolic/HT/FaultedCube
    EXECUTABLE_ARGS Ra_795_fault_bcgs_jacobi.prj
    EXECUTABLE ogs
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ThermalConvection_ts_1_t_0.000000_expected.vtu ThermalConvection_ts_1_t_0.000000.vtu T T 1e-10 1e-16
    ThermalConvection_ts_1_t_0.000000_expected.vtu ThermalConvection_ts_1_t_0.000000.vtu p p 7e-7 1e-12
    ThermalConvection_ts_1_t_0.000000_expected.vtu ThermalConvection_ts_1_t_0.000000.vtu darcy_velocity darcy_velocity 1e-8 1e-13
    VIS ThermalConvection_ts_1_t_0.000000.vtu
)

if(NOT OGS_USE_MPI AND OGS_BUILD_TESTING AND Python3_FOUND)
    add_custom_target(generate_invalid_project_files_ht
        ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/ThirdParty/ogs6py/generateInvalidMediaForHT.py
                                                generateInvalidMediaForHT.py
        WORKING_DIRECTORY ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/)
    file(GLOB HT_INVALID_PRJ_FILES ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/*.prj)
    foreach(ht_invalid_prj_file ${HT_INVALID_PRJ_FILES})
        string(REPLACE ${Data_SOURCE_DIR}/Parabolic/HT/InvalidProjectFiles/HT "invalid" ht_invalid_prj_file_short ${ht_invalid_prj_file})
        AddTest(
            NAME HT_${ht_invalid_prj_file_short}
            PATH Parabolic/HT/InvalidProjectFiles
            EXECUTABLE ogs
            EXECUTABLE_ARGS ${ht_invalid_prj_file}
            RUNTIME 1
            DEPENDS generate_invalid_project_files_ht
        )
        set_tests_properties(ogs-HT_${ht_invalid_prj_file_short} PROPERTIES WILL_FAIL TRUE)
    endforeach()
endif()

if (NOT (OGS_USE_MPI))
    OgsTest(PROJECTFILE Parabolic/HT/SimpleSynthetics/deactivated_subdomain/HT_DeactivatedSubdomain.prj)
endif()

AddTest(
    NAME HT_HeatTransportInStationaryFlow
    PATH Parabolic/HT/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS HeatTransportInStationaryFlow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature 1.e-6 1e-10
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure 1e-10 1e-10
)

AddTest(
    NAME HT_ComponentTransport_ThermalDiffusion_TemperatureField
    PATH Parabolic/ComponentTransport/ThermalDiffusion
    EXECUTABLE ogs
    EXECUTABLE_ARGS TemperatureField.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 29
    DIFF_DATA
    TemperatureFieldts_0_t_0.000000_expected.vtu TemperatureField_ts_0_t_0.000000.vtu T T 1e-10 1e-10
    TemperatureFieldts_73_t_6307200.000000_expected.vtu TemperatureField_ts_73_t_6307200.000000.vtu T T 1e-10 1e-10
    TemperatureFieldts_146_t_12614400.000000_expected.vtu TemperatureField_ts_146_t_12614400.000000.vtu T T 1e-10 1e-10
    TemperatureFieldts_219_t_18921600.000000_expected.vtu TemperatureField_ts_219_t_18921600.000000.vtu T T 1e-10 1e-10
    TemperatureFieldts_292_t_25228800.000000_expected.vtu TemperatureField_ts_292_t_25228800.000000.vtu T T 1e-10 1e-10
    TemperatureFieldts_365_t_31536000.000000_expected.vtu TemperatureField_ts_365_t_31536000.000000.vtu T T 1e-10 1e-10
    TemperatureFieldts_0_t_0.000000_expected.vtu TemperatureField_ts_0_t_0.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    TemperatureFieldts_73_t_6307200.000000_expected.vtu TemperatureField_ts_73_t_6307200.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    TemperatureFieldts_146_t_12614400.000000_expected.vtu TemperatureField_ts_146_t_12614400.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    TemperatureFieldts_219_t_18921600.000000_expected.vtu TemperatureField_ts_219_t_18921600.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    TemperatureFieldts_292_t_25228800.000000_expected.vtu TemperatureField_ts_292_t_25228800.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    TemperatureFieldts_365_t_31536000.000000_expected.vtu TemperatureField_ts_365_t_31536000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-10
    TemperatureFieldts_0_t_0.000000_expected.vtu TemperatureField_ts_0_t_0.000000.vtu p p 1e-10 1e-10
    TemperatureFieldts_73_t_6307200.000000_expected.vtu TemperatureField_ts_73_t_6307200.000000.vtu p p 1e-10 1e-10
    TemperatureFieldts_146_t_12614400.000000_expected.vtu TemperatureField_ts_146_t_12614400.000000.vtu p p 1e-10 1e-10
    TemperatureFieldts_219_t_18921600.000000_expected.vtu TemperatureField_ts_219_t_18921600.000000.vtu p p 1e-10 1e-10
    TemperatureFieldts_292_t_25228800.000000_expected.vtu TemperatureField_ts_292_t_25228800.000000.vtu p p 1e-10 1e-10
    TemperatureFieldts_365_t_31536000.000000_expected.vtu TemperatureField_ts_365_t_31536000.000000.vtu p p 1e-10 1e-10
)

AddTest(
    NAME HT_HeatTransportInStationaryFlow_Staggered_Scheme
    PATH Parabolic/HT/StaggeredCoupling/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS HeatTransportInStationaryFlow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 4
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000_mono.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature 5.e-5 1.0e-4
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000_mono.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure 2e-5 1e-5
)

#MPI/PETSc
AddTest(
    NAME HT_ParallelComputing_HeatTransportInStationaryFlow_Staggered_Scheme
    PATH Parabolic/HT/StaggeredCoupling/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    RUNTIME 2
    EXECUTABLE_ARGS HeatTransportInStationaryFlow.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 3
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    HeatTransportInStationaryFlow_ts_50_t_50000_000000_0.vtu HeatTransportInStationaryFlow_ts_50_t_50000_000000_0.vtu temperature temperature 1.e-9 1.0e-8
    HeatTransportInStationaryFlow_ts_50_t_50000_000000_0.vtu HeatTransportInStationaryFlow_ts_50_t_50000_000000_0.vtu pressure pressure 1.e-9 1.0e-8
    HeatTransportInStationaryFlow_ts_50_t_50000_000000_1.vtu HeatTransportInStationaryFlow_ts_50_t_50000_000000_1.vtu temperature temperature 1.e-9 1.0e-8
    HeatTransportInStationaryFlow_ts_50_t_50000_000000_1.vtu HeatTransportInStationaryFlow_ts_50_t_50000_000000_1.vtu pressure pressure 1.e-9 1.0e-8
    HeatTransportInStationaryFlow_ts_50_t_50000_000000_2.vtu HeatTransportInStationaryFlow_ts_50_t_50000_000000_2.vtu temperature temperature 1.e-9 1.0e-8
    HeatTransportInStationaryFlow_ts_50_t_50000_000000_2.vtu HeatTransportInStationaryFlow_ts_50_t_50000_000000_2.vtu pressure pressure 1.e-9 1.0e-8
)
