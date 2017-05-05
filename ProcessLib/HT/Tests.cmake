AddTest(
    NAME LARGE_2D_ThermalConvection_constviscosity
    PATH Parabolic/HT/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-10 RELTOL 1e-16
    DIFF_DATA
    ThermalConvection_const_viscosity_expected.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu temperature T
    VIS ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu
)
