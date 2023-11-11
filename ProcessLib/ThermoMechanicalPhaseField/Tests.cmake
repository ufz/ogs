AddTest(
    NAME ThermoMechanicalPhaseField_3D_isothermal
    PATH ThermoMechanicalPhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu cube_1e0_tm_ts_1_t_100.000000.vtu displacement displacement 1e-6 1e-6
    expected_cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu cube_1e0_tm_ts_1_t_100.000000.vtu temperature temperature 1e-6 1e-6
    analytical_cube_1e0_tm_pcs_2_ts_1_t_100.000000.vtu cube_1e0_tm_ts_1_t_100.000000.vtu phasefield phasefield 1e-5 1e-5
    mechanical_cube_1e0_pcs_1_ts_1_t_100.000000.vtu cube_1e0_tm_ts_1_t_100.000000.vtu phasefield phasefield 1e-5 1e-5
)

AddTest(
    NAME ThermoMechanicalPhaseField_3D_beam
    PATH ThermoMechanicalPhaseField
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_ts_1_t_10.000000.vtu displacement displacement 1e-5 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_ts_1_t_10.000000.vtu phasefield phasefield 1e-6 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_ts_1_t_10.000000.vtu temperature temperature 1e-6 0
)

AddTest(
    NAME ThermoMechanicalPhaseField_TES_IGLU
    PATH ThermoMechanicalPhaseField
    RUNTIME 1000
    EXECUTABLE ogs
    EXECUTABLE_ARGS tes_hx3_iglu.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_tes_hx3_pcs_2_ts_1_t_20.000000.vtu tes_hx3_ts_1_t_20.000000.vtu phasefield phasefield 1e-6 0
    expected_tes_hx3_pcs_2_ts_1_t_20.000000.vtu tes_hx3_ts_1_t_20.000000.vtu heat_flux heat_flux 1e-4 0
)

AddTest(
    NAME ThermoMechanicalPhaseField_THERMAL_SHOCK
    PATH ThermoMechanicalPhaseField
    RUNTIME 2200
    EXECUTABLE ogs
    EXECUTABLE_ARGS slab_5.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_slab_5_pcs_2_ts_1_t_1.000000.vtu slab_5_ts_1_t_1.000000.vtu phasefield phasefield 1e-6 0
)

AddTest(
    NAME ThermoMechanicalPhaseField_3D_beam_local_coupling_temperature_displacement
    PATH ThermoMechanicalPhaseField
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d_local_coupling_temperature_displacement.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_temperature_displacement_ts_1_t_10.000000.vtu displacement displacement 1e-5 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_temperature_displacement_ts_1_t_10.000000.vtu phasefield phasefield 1e-6 1e-4
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_temperature_displacement_ts_1_t_10.000000.vtu temperature temperature 1e-6 0
)

AddTest(
    NAME ThermoMechanicalPhaseField_3D_beam_local_coupling_temperature_phasefield
    PATH ThermoMechanicalPhaseField
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d_local_coupling_temperature_phasefield.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_temperature_phasefield_ts_1_t_10.000000.vtu displacement displacement 1e-5 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_temperature_phasefield_ts_1_t_10.000000.vtu phasefield phasefield 1e-6 1e-4
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_temperature_phasefield_ts_1_t_10.000000.vtu temperature temperature 1e-6 0
)

AddTest(
    NAME ThermoMechanicalPhaseField_3D_beam_local_coupling_displacement_phasefield
    PATH ThermoMechanicalPhaseField
    RUNTIME 5
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d_local_coupling_displacement_phasefield.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_displacement_phasefield_ts_1_t_10.000000.vtu displacement displacement 1e-5 0
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_displacement_phasefield_ts_1_t_10.000000.vtu phasefield phasefield 1e-6 1e-4
    expected_beam3d_pcs_2_ts_1_t_10.000000.vtu beam3d_local_coupling_displacement_phasefield_ts_1_t_10.000000.vtu temperature temperature 1e-6 0
)
