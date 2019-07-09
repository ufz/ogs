if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE HydroMechanics/InitialStates/displacement/into_initial_state.prj)
    OgsTest(PROJECTFILE HydroMechanics/InitialStates/displacement/non_equilibrium_initial_state.prj)
    OgsTest(PROJECTFILE HydroMechanics/InitialStates/pressure/into_initial_state.prj)
    OgsTest(PROJECTFILE HydroMechanics/InitialStates/pressure/non_equilibrium_initial_state.prj)
endif()

# HydroMechanics; Small deformations, linear poroelastic (HML)

### With monolithic scheme
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad8_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB square_1e2_pcs_0_ts_*.vtu displacement displacement 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu pressure pressure 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-15 1e-15
    GLOB square_1e2_pcs_0_ts_*.vtu NodalForces NodalForces 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_HML_square_1e2_quad9_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_quad9.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu square_1e2_quad9_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_HML_square_1e2_tri6_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_tri.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu square_1e2_tri_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_HML_cube_1e3_hex20_confined_compression
    PATH HydroMechanics/Linear/Confined_Compression
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e3_pcs_0_ts_1_t_5.000000.vtu cube_1e3_pcs_0_ts_1_t_5.000000.vtu displacement displacement 1e-15 1e-15
    expected_cube_1e3_pcs_0_ts_1_t_5.000000.vtu cube_1e3_pcs_0_ts_1_t_5.000000.vtu pressure pressure 1e-15 1e-15
)
AddTest(
    NAME HydroMechanics_hm1_1Dbeam
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_1Dbeam.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu pressure pressure 2.0e-9 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu displacement displacement 1.0e-9 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu sigma_xx sigma_xx 5.0e-8 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu sigma_yy sigma_yy 5.0e-8 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu sigma_zz sigma_zz 5.0e-8 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu sigma_xy sigma_xy 5.0e-8 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu sigma_xz sigma_xz 5.0e-8 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu sigma_yz sigma_yz 5.0e-8 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu epsilon_xx epsilon_xx 1.0e-9 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu epsilon_yy epsilon_yy 1.0e-9 0.0
    hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu hm1_1Dbeam_pcs_1_ts_1_t_1.000000.vtu epsilon_zz epsilon_zz 1.0e-9 0.0
)

AddTest(
    NAME HydroMechanics_hm1_2Dsquare
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_2Dsquare.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu pressure pressure 1.0e-9 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu displacement displacement 1.0e-9 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu sigma_xx sigma_xx 5.0e-8 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu sigma_yy sigma_yy 5.0e-8 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu sigma_zz sigma_zz 5.0e-8 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu sigma_xy sigma_xy 5.0e-8 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu sigma_xz sigma_xz 5.0e-8 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu sigma_yz sigma_yz 5.0e-8 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu epsilon_xx epsilon_xx 1.0e-9 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu epsilon_yy epsilon_yy 1.0e-9 0.0
    hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu hm1_2Dsquare_pcs_1_ts_1_t_1.000000.vtu epsilon_zz epsilon_zz 1.0e-9 0.0
)

AddTest(
    NAME HydroMechanics_hm1_3Dcube
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_3Dcube.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu pressure pressure 2.0e-8 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu displacement displacement 1.0e-8 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu sigma_xx sigma_xx 5.0e-7 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu sigma_yy sigma_yy 5.0e-7 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu sigma_zz sigma_zz 5.0e-7 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu sigma_xy sigma_xy 5.0e-7 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu sigma_xz sigma_xz 5.0e-7 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu sigma_yz sigma_yz 5.0e-7 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu epsilon_xx epsilon_xx 1.0e-8 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu epsilon_yy epsilon_yy 1.0e-8 0.0
    hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu hm1_3Dcube_pcs_1_ts_1_t_1.000000.vtu epsilon_zz epsilon_zz 1.0e-8 0.0
)

AddTest(
    NAME HydroMechanics_hm1_3Dgravity
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm1_3Dgravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu pressure pressure 1.0e-9 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu displacement displacement 1.0e-10 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu sigma_xx sigma_xx 1.0e-8 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu sigma_yy sigma_yy 1.0e-8 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu sigma_zz sigma_zz 2.0e-8 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu sigma_xy sigma_xy 1.0e-8 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu sigma_xz sigma_xz 1.0e-8 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu sigma_yz sigma_yz 1.0e-8 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu epsilon_xx epsilon_xx 1.0e-9 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu epsilon_yy epsilon_yy 1.0e-9 0.0
    hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu hm1_3Dgravity_pcs_1_ts_1_t_1.000000.vtu epsilon_zz epsilon_zz 1.0e-9 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1D1bt
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1D1bt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu pressure pressure 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu displacement displacement 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu sigma_xx sigma_xx 1.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu sigma_yy sigma_yy 1.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu sigma_zz sigma_zz 1.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu sigma_xy sigma_xy 1.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu sigma_xz sigma_xz 1.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu sigma_yz sigma_yz 1.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu epsilon_xx epsilon_xx 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu epsilon_yy epsilon_yy 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu hm2_1D1bt_pcs_0_ts_50_t_5.000000.vtu epsilon_zz epsilon_zz 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu pressure pressure 2.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu displacement displacement 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu sigma_xx sigma_xx 5.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu sigma_yy sigma_yy 5.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu sigma_zz sigma_zz 5.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu sigma_xy sigma_xy 5.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu sigma_xz sigma_xz 5.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu sigma_yz sigma_yz 5.0e-10 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu epsilon_xx epsilon_xx 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu epsilon_yy epsilon_yy 1.0e-11 0.0
    hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu hm2_1D1bt_pcs_0_ts_100_t_10.000000.vtu epsilon_zz epsilon_zz 1.0e-11 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1D2bt
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1D2bt.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu hm2_1D2bt_pcs_0_ts_125_t_5.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu hm2_1D2bt_pcs_0_ts_250_t_10.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1Dbiot
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1Dbiot.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu hm2_1Dbiot_pcs_0_ts_20_t_5.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu hm2_1Dbiot_pcs_0_ts_40_t_10.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1Dcolumn1
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1Dcolumn1.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu pressure pressure 5.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu hm2_1Dcolumn1_pcs_0_ts_120_t_6.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu pressure pressure 5.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu sigma_xx sigma_xx 5.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu sigma_yy sigma_yy 5.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu sigma_zz sigma_zz 5.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu sigma_xy sigma_xy 5.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu sigma_xz sigma_xz 5.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu sigma_yz sigma_yz 5.0e-9 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu hm2_1Dcolumn1_pcs_0_ts_240_t_12.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_1Dcolumn2
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_1Dcolumn2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu hm2_1Dcolumn2_pcs_0_ts_16_t_4.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu sigma_xx sigma_xx 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu sigma_yy sigma_yy 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu sigma_zz sigma_zz 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu sigma_xy sigma_xy 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu sigma_xz sigma_xz 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu sigma_yz sigma_yz 1.0e-9 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu hm2_1Dcolumn2_pcs_0_ts_40_t_10.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
)

AddTest(
    NAME HydroMechanics_hm2_2Dmandel
    PATH HydroMechanics/Verification
    EXECUTABLE ogs
    EXECUTABLE_ARGS hm2_2Dmandel.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 166
    DIFF_DATA
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu sigma_xx sigma_xx 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu sigma_yy sigma_yy 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu sigma_zz sigma_zz 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu sigma_xy sigma_xy 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu sigma_xz sigma_xz 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu sigma_yz sigma_yz 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu hm2_2Dmandel_pcs_0_ts_40_t_2.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu pressure pressure 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu displacement displacement 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu sigma_xx sigma_xx 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu sigma_yy sigma_yy 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu sigma_zz sigma_zz 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu sigma_xy sigma_xy 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu sigma_xz sigma_xz 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu sigma_yz sigma_yz 2.0e-9 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu epsilon_xx epsilon_xx 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu epsilon_yy epsilon_yy 1.0e-10 0.0
    hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu hm2_2Dmandel_pcs_0_ts_160_t_8.000000.vtu epsilon_zz epsilon_zz 1.0e-10 0.0
)

# HydroMechanics; Small deformation, linear poroelastic (unconfined compression early) The drainage process is ongoing and the displacement behaviour is related to water pressure and solid properties.
AddTest(
    NAME HydroMechanics_HML_square_1e2_unconfined_compression_early
    PATH HydroMechanics/Linear/Unconfined_Compression_early
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_UC_early.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu displacement displacement 1e-11 1e-16
    expected_square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu square_1e2_UC_early_pcs_0_ts_10_t_1.000000.vtu pressure pressure 1e-10 1e-16
)
# HydroMechanics; Small deformation, linear poroelastic (unconfined compression late) the drainage process is finished and the displacement of the porous media is only a result of solid properties.
AddTest(
    NAME HydroMechanics_HML_square_1e2_unconfined_compression_late
    PATH HydroMechanics/Linear/Unconfined_Compression_late
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_UC_late.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu displacement displacement 1e-13 1e-16
    expected_square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu square_1e2_UC_late_pcs_0_ts_10_t_1000.000000.vtu pressure pressure 1e-13 1e-16
)

AddTest(
    NAME HydroMechanics_HML_flow_gravity
    PATH HydroMechanics/Linear/Gravity
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_gravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    flow_gravity_pcs_0_ts_16_t_40000000.000000.vtu flow_gravity_pcs_0_ts_16_t_40000000.000000.vtu displacement displacement 1e-14 0
    flow_gravity_pcs_0_ts_16_t_40000000.000000.vtu flow_gravity_pcs_0_ts_16_t_40000000.000000.vtu pressure pressure 1e-10 0
    flow_gravity_pcs_0_ts_16_t_40000000.000000.vtu flow_gravity_pcs_0_ts_16_t_40000000.000000.vtu velocity velocity 1e-10 0
)

## Test as the reference of InjectionProduction1D
AddTest(
    NAME MonolithicInjectionProduction1D
    PATH HydroMechanics/StaggeredScheme/InjectionProduction1D/RerefenceSolutionByMonolithicScheme
    EXECUTABLE ogs
    EXECUTABLE_ARGS InjectionProduction1DMono.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu displacement displacement 1e-11 1e-11
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu pressure pressure 1e-11 1e-11
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu velocity velocity 1e-11 1e-11
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu epsilon_yy epsilon_yy 1e-11 1e-11
    InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu sigma_yy sigma_yy 1e-11 1e-11
)

### With staggered scheme
AddTest(
    NAME StaggeredInjectionProduction1D
    PATH HydroMechanics/StaggeredScheme/InjectionProduction1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS InjectionProduction1D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu displacement displacement 1e-11 1e-11
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu pressure pressure 1e-11 1e-11
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu velocity velocity 1e-11 1e-11
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu epsilon_yy epsilon_yy 1e-11 1e-11
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu sigma_yy sigma_yy 1e-11 1e-11
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu HydraulicFlow HydraulicFlow 1e-11 0
    InjectionProduction1D_Mono_pcs_0_ts_100_t_8640000.000000.vtu InjectionProduction1D_pcs_1_ts_100_t_8640000.000000.vtu NodalForces NodalForces 3e-7 0
)
