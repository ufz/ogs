# Liquid flow
AddTest(
    NAME LiquidFlow_LineDirichletNeumannBC
    PATH Parabolic/LiquidFlow/LineDirichletNeumannBC
    EXECUTABLE ogs
    EXECUTABLE_ARGS line_dirichlet_neumannBC.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    sat1D.vtu sat_1D_pcs_0_ts_1_t_1.000000.vtu AnalyticPressure pressure
    sat1D.vtu sat_1D_pcs_0_ts_1_t_1.000000.vtu AnalyticVec v_x
)
AddTest(
    NAME LiquidFlow_PressureBCatCornerOfAnisotropicSquare
    PATH Parabolic/LiquidFlow/PressureBCatCornerOfAnisotropicSquare
    EXECUTABLE ogs
    EXECUTABLE_ARGS pressureBC_at_corner_of_anisotropic_square.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    mesh2D.vtu sat_2D_lflow_pcs_0_ts_1_t_1.000000.vtu OGS5_Results pressure
)
AddTest(
    NAME LiquidFlow_GravityDriven
    PATH Parabolic/LiquidFlow/GravityDriven
    EXECUTABLE ogs
    EXECUTABLE_ARGS gravity_driven.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    mesh2D.vtu gravity_driven_pcs_0_ts_1_t_1.000000.vtu AnalyticPressure pressure
    mesh2D.vtu gravity_driven_pcs_0_ts_1_t_1.000000.vtu AnalyticVx v_x
    mesh2D.vtu gravity_driven_pcs_0_ts_1_t_1.000000.vtu AnalyticVy v_y
)
AddTest(
    NAME LiquidFlow_AxisymTheis
    PATH Parabolic/LiquidFlow/AxiSymTheis
    EXECUTABLE ogs
    EXECUTABLE_ARGS axisym_theis.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    axisym_theis.vtu liquid_pcs_pcs_0_ts_30_t_1728.000000.vtu OGS5_pressure pressure
)

AddTest(
    NAME LARGE_LiquidFlow_Anisotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    EXECUTABLE ogs
    EXECUTABLE_ARGS anisotropic_gravity_driven3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-6 RELTOL 1e-6
    DIFF_DATA
    hex.vtu anisotropic_gravity_driven3D_pcs_0_ts_1_t_1.000000.vtu analytic_pressure pressure
)

AddTest(
    NAME LARGE_LiquidFlow_Isotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    EXECUTABLE ogs
    EXECUTABLE_ARGS isotropic_gravity_driven3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-6 RELTOL 1e-6
    DIFF_DATA
    hex.vtu isotropic_gravity_driven3D_pcs_0_ts_1_t_1.000000.vtu analytic_pressure pressure
)

# Coupling
AddTest(
    NAME StaggeredTH_ThermalDensityDrivenFlow2D
    PATH StaggeredCoupledProcesses/TH/ThermalDensityDrivenFlow2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS thermal_gravity_driven_flow2d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1.1 RELTOL 1e-7
    DIFF_DATA
    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300.000000.vtu OGS5_PRESSURE1 pressure
    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300.000000.vtu OGS5_TEMPERATURE1 temperature
    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300.000000.vtu v_x_ref v_x
    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300.000000.vtu v_y_ref v_y
)

AddTest(
    NAME Adaptive_dt_StaggeredTH_ThermalDensityDrivenFlow2D
    PATH StaggeredCoupledProcesses/TH/ThermalDensityDrivenFlow2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS thermal_gravity_driven_flow2d_adaptive_dt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1.1 RELTOL 1e-6
    DIFF_DATA
    mesh2D.vtu gravity_driven_adaptive_dt_pcs_1_ts_11_t_300.000000.vtu OGS5_PRESSURE1 pressure
    mesh2D.vtu gravity_driven_adaptive_dt_pcs_1_ts_11_t_300.000000.vtu OGS5_TEMPERATURE1 temperature
    mesh2D.vtu gravity_driven_adaptive_dt_pcs_1_ts_11_t_300.000000.vtu v_x_ref v_x
    mesh2D.vtu gravity_driven_adaptive_dt_pcs_1_ts_11_t_300.000000.vtu v_y_ref v_y
)

AddTest(
    NAME Adaptive_dt_ThermalConvectionFlow2D
    PATH StaggeredCoupledProcesses/TH/ThermalConvection2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS quad_5500x5500_adaptive_dt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1.e-16 RELTOL 1e-16
    DIFF_DATA
    ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000_non_const_mu.vtu  ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000.vtu  pressure pressure
    ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000_non_const_mu.vtu  ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000.vtu  temperature temperature
)

AddTest(
    NAME Adaptive_dt_ThermalConvectionFlow2D_Constant_Viscosity
    PATH StaggeredCoupledProcesses/TH/ThermalConvection2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS quad_5500x5500_adaptive_dt_constant_viscosity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1.e-16 RELTOL 1e-16
    DIFF_DATA
    ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000_const_mu.vtu  ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000.vtu  pressure pressure
    ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000_const_mu-vtu  ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000.vtu  temperature temperature
)

#===============================================================================
# PETSc/MPI
AddTest(
    NAME LiquidFlow_LineDirichletNeumannBC
    PATH Parabolic/LiquidFlow/LineDirichletNeumannBC
    EXECUTABLE_ARGS line_dirichlet_neumannBC.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    sat1D.vtu sat_1D_pcs_0_ts_1_t_1_000000_0.vtu AnalyticPressure pressure
#    sat1D.vtu sat_1D_pcs_0_ts_1_t_1_000000_0.vtu AnalyticVec v_x
)
AddTest(
    NAME LiquidFlow_GravityDriven
    PATH Parabolic/LiquidFlow/GravityDriven
    EXECUTABLE_ARGS gravity_driven.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    mesh2D.vtu gravity_driven_pcs_0_ts_1_t_1_000000_0.vtu AnalyticPressure pressure
#    mesh2D.vtu gravity_driven_pcs_0_ts_1_t_1_000000_0.vtu AnalyticVx v_x
#    mesh2D.vtu gravity_driven_pcs_0_ts_1_t_1_000000_0.vtu AnalyticVy v_y
)
AddTest(
    NAME LiquidFlow_PressureBCatCornerOfAnisotropicSquare
    PATH Parabolic/LiquidFlow/PressureBCatCornerOfAnisotropicSquare
    EXECUTABLE_ARGS pressureBC_at_corner_of_anisotropic_square.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    mesh2D.vtu sat_2D_lflow_pcs_0_ts_1_t_1_000000_0.vtu OGS5_Results pressure
)
AddTest(
    NAME LiquidFlow_AxisymTheis
    PATH Parabolic/LiquidFlow/AxiSymTheis
    EXECUTABLE_ARGS axisym_theis.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-8 RELTOL 1e-8
    DIFF_DATA
    axisym_theis.vtu liquid_pcs_pcs_0_ts_30_t_1728_000000_0.vtu OGS5_pressure pressure
)
AddTest(
    NAME LARGE_LiquidFlow_Anisotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    EXECUTABLE_ARGS anisotropic_gravity_driven3D.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-6 RELTOL 1e-6
    DIFF_DATA
    hex.vtu anisotropic_gravity_driven3D_pcs_0_ts_1_t_1_000000_0.vtu analytic_pressure pressure
)

AddTest(
    NAME LARGE_LiquidFlow_Isotropic_GravityDriven3D
    PATH Parabolic/LiquidFlow/GravityDriven3D
    EXECUTABLE_ARGS isotropic_gravity_driven3D.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-6 RELTOL 1e-6
    DIFF_DATA
    hex.vtu isotropic_gravity_driven3D_pcs_0_ts_1_t_1_000000_0.vtu analytic_pressure pressure
)

# Coupling
AddTest(
    NAME StaggeredTH_ThermalDensityDrivenFlow2D
    PATH StaggeredCoupledProcesses/TH/ThermalDensityDrivenFlow2D
    EXECUTABLE_ARGS thermal_gravity_driven_flow2d.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1.1 RELTOL 1e-7
    DIFF_DATA
    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300_000000_0.vtu OGS5_PRESSURE1  pressure
    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300_000000_0.vtu OGS5_TEMPERATURE1 temperature
# To be activated when the output of velocity is correct under PETSc version.
#    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300_000000_0.vtu v_x_ref v_x
#    mesh2D.vtu gravity_driven_pcs_1_ts_10_t_300_000000_0.vtu v_y_ref v_y

)

AddTest(
    NAME Adaptive_dt_StaggeredTH_ThermalDensityDrivenFlow2D
    PATH StaggeredCoupledProcesses/TH/ThermalDensityDrivenFlow2D
    EXECUTABLE_ARGS thermal_gravity_driven_flow2d_adaptive_dt.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1.1 RELTOL 1e-7
    DIFF_DATA
    mesh2D.vtu gravity_driven_adaptive_dt_pcs_1_ts_11_t_300_000000_0.vtu OGS5_PRESSURE1 pressure
    mesh2D.vtu gravity_driven_adaptive_dt_pcs_1_ts_11_t_300_000000_0.vtu OGS5_TEMPERATURE1 temperature
# To be activated when the output of velocity is correct under PETSc version.
#    mesh2D.vtu gravity_driven_pcs_1_ts_5_t_300_000000_0.vtu v_x_ref v_x
#    mesh2D.vtu gravity_driven_pcs_1_ts_5_t_300_000000_0.vtu v_y_ref v_y

)
AddTest(
    NAME Adaptive_dt_ThermalConvectionFlow2D
    PATH StaggeredCoupledProcesses/TH/ThermalConvection2D
    EXECUTABLE_ARGS quad_5500x5500_adaptive_dt.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1.e-16 RELTOL 1e-16
    DIFF_DATA
    ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000_non_const_mu  ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000_0.vtu  pressure_mu pressure
    ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000_non_const_mu  ConstViscosityThermalConvection_pcs_1_ts_231_t_50000000000.000000_0.vtu  temperature temperature
)

AddTest(
    NAME Adaptive_dt_ThermalConvectionFlow2D_Constant_Viscosity
    PATH StaggeredCoupledProcesses/TH/ThermalConvection2D
    EXECUTABLE_ARGS quad_5500x5500_adaptive_dt_constant_viscosity.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1.e-16 RELTOL 1e-16
    DIFF_DATA
    ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000_const_mu.vtu  ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000_0.vtu  pressure pressure
    ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000_const_mu.vtu  ConstViscosityThermalConvection_pcs_1_ts_137_t_50000000000.000000_0.vtu  temperature temperature
)
