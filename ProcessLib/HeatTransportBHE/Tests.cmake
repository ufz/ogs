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
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
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
    fixed_power_constant_flow_ts_10_t_600.000000.vtu fixed_power_constant_flow_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    fixed_power_constant_flow_ts_10_t_600.000000.vtu fixed_power_constant_flow_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_coaxial_pipe_3D_deep_BHE_CXA
    PATH Parabolic/T/3D_deep_BHE
    RUNTIME 32
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3D_deep_BHE_CXA.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D_deep_BHE_CXA_ts_10_t_600.000000.vtu 3D_deep_BHE_CXA_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    3D_deep_BHE_CXA_ts_10_t_600.000000.vtu 3D_deep_BHE_CXA_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_coaxial_pipe_3D_deep_BHE_CXC
    PATH Parabolic/T/3D_deep_BHE
    RUNTIME 32
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3D_deep_BHE_CXC.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D_deep_BHE_CXC_ts_10_t_600.000000.vtu 3D_deep_BHE_CXC_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    3D_deep_BHE_CXC_ts_10_t_600.000000.vtu 3D_deep_BHE_CXC_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_3D_2U_BHE
    PATH Parabolic/T/3D_2U_BHE
    RUNTIME 14
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3D_2U_BHE.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3D_2U_BHE_ts_10_t_600.000000.vtu 3D_2U_BHE_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-12 1e-14
    3D_2U_BHE_ts_10_t_600.000000.vtu 3D_2U_BHE_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
)

AddTest(
    NAME HeatTransportBHE_3D_BHE_groundwater_advection
    PATH Parabolic/T/3D_BHE_GW_advection
    RUNTIME 8
    EXECUTABLE ogs
    EXECUTABLE_ARGS BHE_GW_advection.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    BHE_GW_advection_ts_10_t_500.000000.vtu BHE_GW_advection_ts_10_t_500.000000.vtu temperature_BHE1 temperature_BHE1 1e-12 1e-14
    BHE_GW_advection_ts_10_t_500.000000.vtu BHE_GW_advection_ts_10_t_500.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
)

if("${Python3_VERSION}" VERSION_LESS 3.9)
    AddTest(
        NAME HeatTransportBHE_3D_3BHEs_array
        PATH Parabolic/T/3D_3BHEs_array
        RUNTIME 50
        EXECUTABLE ogs
        EXECUTABLE_ARGS 3bhes_1U.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_PYTHON AND NOT OGS_USE_MPI
        PYTHON_PACKAGES "TESPy=0.3.2"
        DIFF_DATA
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-10 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 1e-10 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 1e-10 1e-13
    )

    AddTest(
        NAME HeatTransportBHE_1U_3D_beier_sandbox_server_communication
        PATH Parabolic/T/3D_Beier_sandbox_SimX
        EXECUTABLE ogs
        EXECUTABLE_ARGS beier_sandbox.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_PYTHON AND NOT OGS_USE_MPI
        RUNTIME 50
        PYTHON_PACKAGES "pandas"
        DIFF_DATA
        beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
        beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
    )

    AddTest(
        NAME HeatTransportBHE_3D_3BHEs_array_server_communication
        PATH Parabolic/T/3D_3BHEs_array_SimX
        RUNTIME 50
        EXECUTABLE ogs
        EXECUTABLE_ARGS 3bhes_1U.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS OGS_USE_PYTHON AND NOT OGS_USE_MPI
        PYTHON_PACKAGES "TESPy=0.3.2"
        DIFF_DATA
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-9 1e-12
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 1e-9 1e-12
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 1e-9 1e-12
    )
endif()

AddTest(
    NAME HeatTransportBHE_single_pipe_flow_EUBHE
    PATH Parabolic/T/BHE_1P
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS BHE_1P.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    BHE_1P_ts_10_t_600.000000.vtu BHE_1P_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-12 1e-14
    BHE_1P_ts_10_t_600.000000.vtu BHE_1P_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
)

AddTest(
    NAME HeatTransportBHE_single_pipe_flow_EUBHE_newton
    PATH Parabolic/T/BHE_1P
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS BHE_1P_newton.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    BHE_1P_newton_ts_10_t_600.000000.vtu BHE_1P_newton_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-12 1e-14
    BHE_1P_newton_ts_10_t_600.000000.vtu BHE_1P_newton_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
)
