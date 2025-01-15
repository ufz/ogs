AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich
    PATH Parabolic/T/3D_BHE_Sandwich
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 10
    DIFF_DATA
    sandwich_ts_10_t_600.000000.vtu sandwich_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    sandwich_ts_10_t_600.000000.vtu sandwich_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_linear
    PATH Parabolic/T/3D_BHE_Sandwich
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_linear.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 10
    DIFF_DATA
    sandwich_ts_10_t_600.000000.vtu sandwich_linear_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    sandwich_ts_10_t_600.000000.vtu sandwich_linear_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_Newton
    PATH Parabolic/T/3D_BHE_Sandwich
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 10
    DIFF_DATA
    sandwich_ts_10_t_600.000000.vtu sandwich_newton_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-14
    sandwich_ts_10_t_600.000000.vtu sandwich_newton_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_algebraicBC
    PATH Parabolic/T/3D_BHE_Sandwich
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_algebraicBC.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 3
    DIFF_DATA
    sandwich_ts_10_t_600.000000.vtu sandwich_algebraic_bc_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-6
    sandwich_ts_10_t_600.000000.vtu sandwich_algebraic_bc_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-9
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_algebraicBC_LSCG
    PATH Parabolic/T/3D_BHE_Sandwich
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_algebraicBC_LSCG.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 3
    DIFF_DATA
    sandwich_ts_10_t_600.000000.vtu sandwich_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-6
    sandwich_ts_10_t_600.000000.vtu sandwich_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-9
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_fixed_power
    PATH Parabolic/T/3D_BHE_Sandwich
    RUNTIME 100
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_fixed_power.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    sandwich_fixed_power_ts_10_t_600.000000.vtu sandwich_fixed_power_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    sandwich_fixed_power_ts_10_t_600.000000.vtu sandwich_fixed_power_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_fixed_power_algebraicBC
    PATH Parabolic/T/3D_BHE_Sandwich
    RUNTIME 3
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_fixed_power_algebraicBC.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    sandwich_fixed_power_ts_10_t_600.000000.vtu sandwich_fixed_power_algebraic_bc_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-3
    sandwich_fixed_power_ts_10_t_600.000000.vtu sandwich_fixed_power_algebraic_bc_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-6
)

AddTest(
    NAME HeatTransportBHE_1U_3D_bhe_sandwich_fixed_power_algebraicBC_LSCG
    PATH Parabolic/T/3D_BHE_Sandwich
    RUNTIME 3
    EXECUTABLE ogs
    EXECUTABLE_ARGS sandwich_fixed_power_algebraicBC_LSCG.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    sandwich_fixed_power_ts_10_t_600.000000.vtu sandwich_fixed_power_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-3
    sandwich_fixed_power_ts_10_t_600.000000.vtu sandwich_fixed_power_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-6
)

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
    NAME HeatTransportBHE_1U_3D_beier_sandbox_linear
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox_linear.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 20
    DIFF_DATA
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_linear_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_linear_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_3D_beier_sandbox_Newton
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 20
    DIFF_DATA
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_newton_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_newton_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
)

AddTest(
    NAME HeatTransportBHE_1U_3D_MassLumping
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox_MassLumping.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 20
    DIFF_DATA
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_mass_lumping_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-6
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_mass_lumping_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-4
)

AddTest(
    NAME HeatTransportBHE_1U_3D_beier_sandbox_binary_curve
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox_binary_curve.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 20
    DIFF_DATA
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_binary_curve_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-14
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_binary_curve_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-12
)

AddTest(
    NAME HeatTransportBHE_1U_3D_beier_sandbox_algebraicBC
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox_algebraicBC.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 3
    DIFF_DATA
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_algebraic_bc_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-7
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_algebraic_bc_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-10
)

AddTest(
    NAME HeatTransportBHE_1U_3D_beier_sandbox_algebraicBC_LSCG
    PATH Parabolic/T/3D_Beier_sandbox
    EXECUTABLE ogs
    EXECUTABLE_ARGS beier_sandbox_algebraicBC_LSCG.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 3
    DIFF_DATA
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-7
    beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-10
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
    NAME HeatTransportBHE_1U_beier_sandbox_fixed_power_constant_flow_algebraicBC
    PATH Parabolic/T/3D_Beier_sandbox
    RUNTIME 3
    EXECUTABLE ogs
    EXECUTABLE_ARGS fixed_power_constant_flow_algebraicBC.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    fixed_power_constant_flow_ts_10_t_600.000000.vtu fixed_power_constant_flow_algebraic_bc_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-4
    fixed_power_constant_flow_ts_10_t_600.000000.vtu fixed_power_constant_flow_algebraic_bc_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-9
)

AddTest(
    NAME HeatTransportBHE_1U_beier_sandbox_fixed_power_constant_flow_algebraicBC_LSCG
    PATH Parabolic/T/3D_Beier_sandbox
    RUNTIME 3
    EXECUTABLE ogs
    EXECUTABLE_ARGS fixed_power_constant_flow_algebraicBC_LSCG.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    fixed_power_constant_flow_ts_10_t_600.000000.vtu fixed_power_constant_flow_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-4
    fixed_power_constant_flow_ts_10_t_600.000000.vtu fixed_power_constant_flow_algebraic_bc_LSCG_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 5e-9
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
    NAME HeatTransportBHE_3D_3BHEs
    PATH Parabolic/T/3D_3BHEs
    RUNTIME 14
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3bhes.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3bhes_ts_10_t_600.000000.vtu 3bhes_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-12
    3bhes_ts_10_t_600.000000.vtu 3bhes_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 0 1e-12
    3bhes_ts_10_t_600.000000.vtu 3bhes_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 0 1e-12
    3bhes_ts_10_t_600.000000.vtu 3bhes_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-12
    3bhes_1_ts_10_t_600.000000.vtu 3bhes_1_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 2e-13 1e-12
    3bhes_1_ts_10_t_600.000000.vtu 3bhes_1U_BHE_1_mesh_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 2e-13 1e-12
)

AddTest(
    NAME HeatTransportBHE_3D_3BHEs_BHE_id
    PATH Parabolic/T/3D_3BHEs
    RUNTIME 14
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3bhes_id_1U_2U_1U.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3bhes_ts_10_t_600.000000.vtu 3bhes_id_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-12
    3bhes_ts_10_t_600.000000.vtu 3bhes_id_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 0 1e-12
    3bhes_ts_10_t_600.000000.vtu 3bhes_id_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 0 1e-12
    3bhes_ts_10_t_600.000000.vtu 3bhes_id_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-12
)

AddTest(
    NAME HeatTransportBHE_3D_3BHEs_BHE_id_failcase_id_out_of_range
    PATH Parabolic/T/3D_3BHEs
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS id_out_of_range.xml
    REQUIREMENTS NOT OGS_USE_MPI
    PROPERTIES
    PASS_REGULAR_EXPRESSION "The maximum given BHE id '100' did not match the number of given BHE definitions '3'. The BHE ids needs to be defined starting from 0, so the maximum BHE id needs to be number of BHE definitions minus 1. After all definitions there are no gaps allowed between the given ids."
)

AddTest(
    NAME HeatTransportBHE_3D_3BHEs_BHE_id_failcase_duplicate_id
    PATH Parabolic/T/3D_3BHEs
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS duplicate_id.xml
    REQUIREMENTS NOT OGS_USE_MPI
    PROPERTIES
    PASS_REGULAR_EXPRESSION "BHE with id '1' is already present in the list! Check for duplicate definitions of BHE ids."
)

AddTest(
    NAME HeatTransportBHE_3D_3BHEs_BHE_id_1U
    PATH Parabolic/T/3D_3BHEs
    RUNTIME 14
    EXECUTABLE ogs
    EXECUTABLE_ARGS 3bhes_id_1U.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_id_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 1e-12
    3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_id_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 0 1e-12
    3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_id_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 0 1e-12
    3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_id_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-12
)

AddTest(
    NAME HeatTransportBHE_3D_BHE_groundwater_advection
    PATH Parabolic/T/3D_BHE_GW_advection
    RUNTIME 4
    EXECUTABLE ogs
    EXECUTABLE_ARGS BHE_GW_advection.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    BHE_GW_advection_ts_10_t_500.000000.vtu BHE_GW_advection_ts_10_t_500.000000.vtu temperature_BHE1 temperature_BHE1 1e-12 1e-14
    BHE_GW_advection_ts_10_t_500.000000.vtu BHE_GW_advection_ts_10_t_500.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
)

if("${Python_VERSION}" VERSION_LESS 3.9)
    AddTest(
        NAME HeatTransportBHE_3D_3BHEs_array
        PATH Parabolic/T/3D_3BHEs_array
        RUNTIME 50
        EXECUTABLE ogs
        EXECUTABLE_ARGS 3bhes_1U.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        PYTHON_PACKAGES "TESPy==0.3.2"
        DIFF_DATA
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-10 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 1e-10 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 1e-10 1e-13
    )

    AddTest(
        NAME HeatTransportBHE_1U_3D_beier_sandbox_python_interface
        PATH Parabolic/T/3D_Beier_sandbox_python_interface
        EXECUTABLE ogs
        EXECUTABLE_ARGS beier_sandbox.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        RUNTIME 50
        PYTHON_PACKAGES "pandas==1.4.2"
        DIFF_DATA
        beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 0 5e-15
        beier_sandbox_ts_10_t_600.000000.vtu beier_sandbox_ts_10_t_600.000000.vtu temperature_soil temperature_soil 0 1e-13
    )

    AddTest(
        NAME HeatTransportBHE_3D_3BHEs_array_python_interface
        PATH Parabolic/T/3D_3BHEs_array_python_interface
        RUNTIME 50
        EXECUTABLE ogs
        EXECUTABLE_ARGS 3bhes_1U.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT OGS_USE_MPI
        PYTHON_PACKAGES "TESPy==0.3.2"
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

if(NOT OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE Parabolic/T/BHE_1P/pipe_flow_ebhe.py RUNTIME 200)
endif()
