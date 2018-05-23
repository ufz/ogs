AddTest(
    NAME LARGE_2D_ThermalConvection_constviscosity
    PATH Parabolic/HT/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu temperature T 1e-10 1e-16
    VIS ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu
)

AddTest(
    NAME LARGE_2D_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000.vtu T_ref T 1e-16  1.e-16
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000.vtu p_ref p  1e-16  1.e-16
    VIS ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000.vtu
)

AddTest(
    NAME LARGE_2D_Adaptive_dt_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme_adaptive_dt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_135_t_50000000000.000000.vtu T_ref T 1e-3  1.e-3
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_135_t_50000000000.000000.vtu p_ref p  1e-3  1.e-3
    VIS ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_135_t_50000000000.000000.vtu
)

AddTest(
    NAME HT_a_DECOVALEX_THMC_based_Example
    PATH Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample
    EXECUTABLE ogs
    EXECUTABLE_ARGS th_decovalex.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    th_decovalex.vtu th_decovalex_pcs_1_ts_40_t_18.000000.vtu T_ref T 1e-14  1.e-14
    th_decovalex.vtu th_decovalex_pcs_1_ts_40_t_18.000000.vtu p_ref p 1e-14  1.e-14
    VIS th_decovalex_pcs_1_ts_78_t_1000.000000.vtu
)

# MPI/PETSc tests
AddTest(
    NAME Parallel_LARGE_2D_ThermalConvection_constviscosity
    PATH Parabolic/HT/ConstViscosity
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
    NAME LARGE_2D_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
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
    NAME LARGE_2D_Adaptive_dt_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme_adaptive_dt.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_135_t_50000000000_000000_0.vtu T_ref T 1e-3  1.e-3
    square_5500x5500.vtu ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_135_t_50000000000_000000_0.vtu p_ref p  1e-3  1.e-3
    VIS ConstViscosityThermalConvectionStaggeredAdaptive_dt_pcs_1_ts_135_t_50000000000_000000_0.vtu
)

AddTest(
    NAME HT_a_DECOVALEX_THMC_based_Example
    PATH Parabolic/HT/StaggeredCoupling/ADecovalexTHMCBasedHTExample
    EXECUTABLE_ARGS th_decovalex.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    th_decovalex.vtu th_decovalex_pcs_1_ts_40_t_18_000000_0.vtu T_ref T 1e-10  1.e-10
    th_decovalex.vtu th_decovalex_pcs_1_ts_40_t_18_000000_0.vtu p_ref p 1e-10  1.e-10
    VIS th_decovalex_pcs_1_ts_78_t_1000_000000_0.vtu
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
