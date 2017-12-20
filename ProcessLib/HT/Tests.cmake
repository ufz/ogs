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
