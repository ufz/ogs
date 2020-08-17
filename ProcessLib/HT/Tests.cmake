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
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu T T 1e-7 1e-9
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu p p 1e-9 1e-3
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu
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
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlow_pcs_0_ts_1_t_1.000000.vtu T T 1e-10 1e-16
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlow_pcs_0_ts_1_t_1.000000.vtu p p 1e-10 1e-16
    IsothermalFluidFlow_expected.vtu IsothermalFluidFlow_pcs_0_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS IsothermalFluidFlow_pcs_0_ts_1_t_1.000000.vtu
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
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_pcs_0_ts_1_t_1.000000.vtu linear_top2_to_bottom1 T 1e-10 1e-16
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_pcs_0_ts_1_t_1.000000.vtu Linear_1_to_minus1 p 1e-10 1e-16
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_pcs_0_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS PressureDiffusionTemperatureDiffusion_pcs_0_ts_1_t_1.000000.vtu
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
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravity_pcs_0_ts_1_t_1.000000.vtu T T 1e-10 1e-16
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravity_pcs_0_ts_1_t_1.000000.vtu p p 1e-10 1e-16
    IsothermalFluidFlowWithGravity_expected.vtu IsothermalFluidFlowWithGravity_pcs_0_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    VIS IsothermalFluidFlow_pcs_0_ts_1_t_1.000000.vtu
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
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_1_t_0.100000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_1_t_0.100000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_1_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_1_t_0.100000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_2_t_0.200000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_2_t_0.200000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_2_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_2_t_0.200000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_3_t_0.300000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_3_t_0.300000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_3_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_3_t_0.300000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_4_t_0.400000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_4_t_0.400000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_4_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_4_t_0.400000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_5_t_0.500000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_5_t_0.500000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_5_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_5_t_0.500000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_6_t_0.600000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_6_t_0.600000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_6_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_6_t_0.600000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_7_t_0.700000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_7_t_0.700000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_7_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_7_t_0.700000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_8_t_0.800000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_8_t_0.800000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_8_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_8_t_0.800000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_9_t_0.900000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_9_t_0.900000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_9_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_9_t_0.900000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_10_t_1.000000.vtu T T 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_10_t_1.000000.vtu p p 1e-10 1e-16
    PressureParabolicTemperatureParabolic_ts_10_expected.vtu PressureParabolicTemperatureParabolic_pcs_0_ts_10_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
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
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_1_t_0.100000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_1_t_0.100000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_1_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_1_t_0.100000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_2_t_0.200000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_2_t_0.200000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_2_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_2_t_0.200000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_3_t_0.300000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_3_t_0.300000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_3_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_3_t_0.300000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_4_t_0.400000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_4_t_0.400000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_4_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_4_t_0.400000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_5_t_0.500000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_5_t_0.500000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_5_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_5_t_0.500000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_6_t_0.600000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_6_t_0.600000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_6_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_6_t_0.600000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_7_t_0.700000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_7_t_0.700000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_7_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_7_t_0.700000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_8_t_0.800000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_8_t_0.800000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_8_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_8_t_0.800000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_9_t_0.900000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_9_t_0.900000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_9_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_9_t_0.900000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_10_t_1.000000.vtu T T 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_10_t_1.000000.vtu p p 1e-10 1e-16
    CoupledPressureParabolicTemperatureParabolic_ts_10_expected.vtu CoupledPressureParabolicTemperatureParabolic_pcs_0_ts_10_t_1.000000.vtu darcy_velocity darcy_velocity 1e-10 1e-16
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
    flux_1e3_t_0.000000.vtu flux_1e3_t_0.000000.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_0.000010.vtu flux_1e3_t_0.000010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_0.001010.vtu flux_1e3_t_0.001010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_0.101010.vtu flux_1e3_t_0.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_1.101010.vtu flux_1e3_t_1.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e3_t_10.000000.vtu flux_1e3_t_10.000000.vtu specific_flux specific_flux 1e-10 1e-16
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
    flux_1e4_t_0.000000.vtu flux_1e4_t_0.000000.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_0.000010.vtu flux_1e4_t_0.000010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_0.001010.vtu flux_1e4_t_0.001010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_0.101010.vtu flux_1e4_t_0.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_1.101010.vtu flux_1e4_t_1.101010.vtu specific_flux specific_flux 1e-10 1e-16
    flux_1e4_t_10.000000.vtu flux_1e4_t_10.000000.vtu specific_flux specific_flux 1e-10 1e-16
)

# Staggered scheme
AddTest(
    NAME 2D_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu T_ref T 1e-1  1.e-3
    square_5500x5500.vtu ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu p_ref p  2e+4  1.e-2
    square_5500x5500.vtu ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu darcy_velocity_ref darcy_velocity  1e-1  1.e-3
    VIS ConstViscosityThermalConvection_ts_149_t_50000000000.000000.vtu
)

AddTest(
    NAME 2D_Adaptive_dt_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme_adaptive_dt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_ts_141_t_50000000000.000000.vtu T_ref T 1e-3  1.e-3
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_ts_141_t_50000000000.000000.vtu p_ref p  1e-3  2.e-3
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_ts_141_t_50000000000.000000.vtu darcy_velocity_ref darcy_velocity  1e-3  1.e-3
    VIS ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_141_t_50000000000.000000.vtu
)
# Workaround sporadic timeouts on macOS
if(APPLE AND TEST ogs-2D_Adaptive_dt_ThermalConvection_constviscosityStaggeredScheme-time)
    set_tests_properties(
        ogs-2D_Adaptive_dt_ThermalConvection_constviscosityStaggeredScheme-time
        PROPERTIES TIMEOUT 1800)
endif()

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

# MPI/PETSc tests
AddTest(
    DISABLED
    NAME Parallel_2D_ThermalConvection_constviscosity
    PATH Parabolic/HT/ConstViscosity
    RUNTIME 61 # Actual RUNTIME?
    EXECUTABLE_ARGS square_5500x5500.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu p p 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu p p 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu p p 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu p p 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu T T 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu T T 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu T T 1e-15 1e-14
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu T T 1e-15 1e-14
)

AddTest(
    DISABLED
    NAME Parallel_2D_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    RUNTIME 61 # Actual RUNTIME?
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000_000000_0.vtu T_ref T  1e-14  1.e-14
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000_000000_0.vtu p_ref p  1e-14  1.e-14
    VIS ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000_0.vtu
)

AddTest(
    DISABLED
    NAME 2D_Adaptive_dt_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    RUNTIME 61 # Actual RUNTIME?
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme_adaptive_dt.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_141_t_50000000000_000000_0.vtu T_ref T 1e-3  1.e-3
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_141_t_50000000000_000000_0.vtu p_ref p  1e-3  1.e-3
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_141_t_50000000000_000000_0.vtu darcy_velocity_ref darcy_velocity  1e-3  1.e-3
    VIS ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_141_t_50000000000_000000_0.vtu
)

# 2019-05-09 TF disable the test until the MPL can deal with parameters as properties
# AddTest(
#     NAME HT_a_DECOVALEX_THMC_based_Example
#     PATH Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample
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
    ThermalConvection_pcs_0_ts_1_t_0.000000_expected.vtu ThermalConvection_pcs_0_ts_1_t_0.000000.vtu T T 1e-10 1e-16
    ThermalConvection_pcs_0_ts_1_t_0.000000_expected.vtu ThermalConvection_pcs_0_ts_1_t_0.000000.vtu p p 7e-7 1e-12
    ThermalConvection_pcs_0_ts_1_t_0.000000_expected.vtu ThermalConvection_pcs_0_ts_1_t_0.000000.vtu darcy_velocity darcy_velocity 1e-8 1e-13
    VIS ThermalConvection_pcs_0_ts_1_t_0.000000.vtu
)

if(NOT OGS_USE_MPI AND BUILD_TESTING AND Python3_FOUND)
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
