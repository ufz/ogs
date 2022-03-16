# LIE; HydroMechanics
if (NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_constK.prj RUNTIME 17)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_CZ_kf_const.prj RUNTIME 250)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_CZ_kf_cubic.prj RUNTIME 120)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_LE_kf_const.prj RUNTIME 90)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_LE_kf_cubic.prj RUNTIME 125)
endif()

AddTest(
    NAME LIE_HM_single_fracture
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 20
    DIFF_DATA
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu pressure pressure 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu displacement displacement 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu nodal_w nodal_w 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu nodal_aperture nodal_aperture 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu strain_xx strain_xx 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu strain_yy strain_yy 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu strain_xy strain_xy 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu stress_xx stress_xx 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu stress_yy stress_yy 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu velocity velocity 1e-12 1e-12
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu MassFlowRate MassFlowRate 1e-12 0
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu NodalForces NodalForces 5e-6 0
    expected_single_fracture_ts_10_t_100.000000.vtu single_fracture_ts_10_t_100.000000.vtu NodalForcesJump NodalForcesJump 1e-6 0
)

AddTest(
    NAME LIE_HM_single_fracture_3D
    PATH LIE/HydroMechanics
    RUNTIME 320
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu pressure pressure 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu displacement displacement 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu nodal_w nodal_w 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu nodal_aperture nodal_aperture 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu strain_xx strain_xx 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu strain_yy strain_yy 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu strain_zz strain_zz 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu strain_xy strain_xy 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu strain_xz strain_xz 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu strain_yz strain_yz 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu stress_xx stress_xx 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu stress_yy stress_yy 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu stress_zz stress_zz 1e-12 1e-12
    expected_single_fracture_3D_ts_10_t_100.000000.vtu single_fracture_3D_ts_10_t_100.000000.vtu velocity velocity 1e-12 1e-12
)

AddTest(
    NAME LIE_HM_TaskB
    PATH LIE/HydroMechanics
    RUNTIME 80
    EXECUTABLE ogs
    EXECUTABLE_ARGS TaskB.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu pressure pressure 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu displacement displacement 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu nodal_w nodal_w 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu nodal_aperture nodal_aperture 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu strain_xx strain_xx 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu strain_yy strain_yy 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu strain_xy strain_xy 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu stress_xx stress_xx 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu stress_yy stress_yy 1e-12 1e-12
    expected_TaskB_ts_4_t_18.000000.vtu TaskB_ts_4_t_18.000000.vtu velocity velocity 1e-12 1e-12
)

AddTest(
    NAME LIE_HM_single_fracture_3compartments_flow
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3compartments_flow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu pressure pressure 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu displacement displacement 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu displacement_jump1 displacement_jump1 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu nodal_w nodal_w 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu nodal_aperture nodal_aperture 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu strain_xx strain_xx 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu strain_yy strain_yy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu strain_xy strain_xy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu stress_xx stress_xx 1e-12 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu stress_yy stress_yy 1e-12 1e-15
    expected_single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_ts_1108_t_10000.000000.vtu velocity velocity 1e-15 1e-15
)

AddTest(
    NAME LIE_HM_single_fracture_3compartments_flow_linear_aperture0
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3compartments_flow_linear_aperture0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu pressure pressure 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu displacement displacement 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu displacement_jump1 displacement_jump1 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu nodal_w nodal_w 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu nodal_aperture nodal_aperture 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu strain_xx strain_xx 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu strain_yy strain_yy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu strain_xy strain_xy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu stress_xx stress_xx 1e-12 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu stress_yy stress_yy 1e-12 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_ts_1108_t_10000.000000.vtu velocity velocity 1e-15 1e-15
)

# Same as the LIE_HM_single_fracture_3compartments_flow_linear_aperture0 but with
# aperture0 defined on the elements (and discontinuous on the nodes).
AddTest(
    NAME LIE_HM_single_fracture_3compartments_flow_linear_aperture0_e
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3compartments_flow_linear_aperture0_e.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu pressure pressure 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu displacement displacement 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu displacement_jump1 displacement_jump1 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu nodal_w nodal_w 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu nodal_aperture nodal_aperture 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu strain_xx strain_xx 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu strain_yy strain_yy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu strain_xy strain_xy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu stress_xx stress_xx 1e-12 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu stress_yy stress_yy 1e-12 1e-15
    expected_single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu single_fracture_3compartments_flow_linear_aperture0_e_ts_1108_t_10000.000000.vtu velocity velocity 1e-15 1e-15
)

AddTest(
    NAME LIE_HM_single_fracture_3compartments_flow_CHZ
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3compartments_flow_CHZ.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu pressure pressure 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu displacement displacement 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu displacement_jump1 displacement_jump1 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu nodal_w nodal_w 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu nodal_aperture nodal_aperture 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu strain_xx strain_xx 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu strain_yy strain_yy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu strain_xy strain_xy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu stress_xx stress_xx 1e-10 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu stress_yy stress_yy 2e-10 1e-15
    expected_single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_ts_356_t_500.000000.vtu velocity velocity 1e-15 1e-15
)

AddTest(
    NAME LIE_HM_single_fracture_3compartments_flow_CHZ_sigma0
    PATH LIE/HydroMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_fracture_3compartments_flow_CHZ_sigma0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu pressure pressure 1e-10 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu pressure_interpolated pressure_interpolated 1e-10 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu displacement displacement 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu displacement_jump1 displacement_jump1 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu nodal_w nodal_w 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu nodal_aperture nodal_aperture 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu strain_xx strain_xx 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu strain_yy strain_yy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu strain_xy strain_xy 1e-15 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu stress_xx stress_xx 3e-10 1e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu stress_yy stress_yy 1e-12 1.5e-15
    expected_single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu single_fracture_3compartments_flow_CHZ_sigma0_ts_356_t_500.000000.vtu velocity velocity 1e-15 1e-15
)
