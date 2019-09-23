AddTest(
    NAME HeatTransportBHE_1U_3D_beier_sandbox
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 20
    DIFF_DATA
    beier_sandbox_pcs_0_ts_10_t_600.000000.vtu beier_sandbox_pcs_0_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    beier_sandbox_pcs_0_ts_10_t_600.000000.vtu beier_sandbox_pcs_0_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_beier_sandbox_fixed_power_constant_flow
    PATH Parabolic/T/3D_Beier_sandbox
    RUNTIME 220
    EXECUTABLE ogs
    EXECUTABLE_ARGS fixed_power_constant_flow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    fixed_power_constant_flow_pcs_0_ts_10_t_600.000000.vtu fixed_power_constant_flow_pcs_0_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    fixed_power_constant_flow_pcs_0_ts_10_t_600.000000.vtu fixed_power_constant_flow_pcs_0_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_coaxial_pipe_3D_deep_BHE_CXA
    PATH Parabolic/T/3D_deep_BHE
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3D_deep_BHE_CXA.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D_deep_BHE_CXA_pcs_0_ts_10_t_600.000000.vtu 3D_deep_BHE_CXA_pcs_0_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    3D_deep_BHE_CXA_pcs_0_ts_10_t_600.000000.vtu 3D_deep_BHE_CXA_pcs_0_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_coaxial_pipe_3D_deep_BHE_CXC
    PATH Parabolic/T/3D_deep_BHE
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3D_deep_BHE_CXC.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D_deep_BHE_CXC_pcs_0_ts_10_t_600.000000.vtu 3D_deep_BHE_CXC_pcs_0_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    3D_deep_BHE_CXC_pcs_0_ts_10_t_600.000000.vtu 3D_deep_BHE_CXC_pcs_0_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_3D_2U_BHE
    PATH Parabolic/T/3D_2U_BHE
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3D_2U_BHE.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D_2U_BHE_pcs_0_ts_10_t_600.000000.vtu 3D_2U_BHE_pcs_0_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-12 1e-14
    3D_2U_BHE_pcs_0_ts_10_t_600.000000.vtu 3D_2U_BHE_pcs_0_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
)
