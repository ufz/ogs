if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_linear.xml RUNTIME 1)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_algebraicBC.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_algebraicBC_LSCG.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_fixed_power.prj
            RUNTIME 1
    )
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_fixed_power_algebraicBC.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_Sandwich/sandwich_fixed_power_algebraicBC_LSCG.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox.prj
            RUNTIME 10
    )
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox_linear.xml RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox_newton.xml RUNTIME 11)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox_MassLumping.xml RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox_binary_curve.xml RUNTIME 10)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox_algebraicBC.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/beier_sandbox_algebraicBC_LSCG.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/fixed_power_constant_flow.prj RUNTIME 125)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/fixed_power_constant_flow_algebraicBC.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_Beier_sandbox/fixed_power_constant_flow_algebraicBC_LSCG.xml RUNTIME 3)
    OgsTest(PROJECTFILE Parabolic/T/3D_deep_BHE/3D_deep_BHE_CXA.prj RUNTIME 19)
    OgsTest(PROJECTFILE Parabolic/T/3D_deep_BHE/3D_deep_BHE_CXC.prj RUNTIME 19)
    OgsTest(PROJECTFILE Parabolic/T/3D_2U_BHE/3D_2U_BHE.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/T/3D_2U_BHE/3D_2U_BHE_sections.prj RUNTIME 14)
    OgsTest(PROJECTFILE Parabolic/T/3D_3BHEs/3bhes.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/T/3D_3BHEs/3bhes_id_1U_2U_1U.prj RUNTIME 9)
    OgsTest(
        PROJECTFILE Parabolic/T/3D_3BHEs/id_out_of_range.xml
        RUNTIME 1
        PROPERTIES PASS_REGULAR_EXPRESSION "BHE id 100 is out of range. The mesh contains 3 BHE"
    )
    OgsTest(
        PROJECTFILE Parabolic/T/3D_3BHEs/duplicate_id.xml
        RUNTIME 1
        PROPERTIES PASS_REGULAR_EXPRESSION "BHE with id '1' is already present in the list! Check for duplicate definitions of BHE ids."
    )
    OgsTest(PROJECTFILE Parabolic/T/3D_3BHEs/3bhes_id_1U.prj RUNTIME 9)
    OgsTest(PROJECTFILE Parabolic/T/3D_3BHEs/3bhes_id_1U_sections.prj
            RUNTIME 14
    )
    OgsTest(PROJECTFILE Parabolic/T/3D_BHE_GW_advection/BHE_GW_advection.prj
            RUNTIME 4
    )
endif()

# TODO: update to newer Python or remove!
if("${Python_VERSION}" VERSION_LESS 3.9)
    AddTest(
        NAME HeatTransportBHE_3D_3BHEs_array
        PATH Parabolic/T/3D_3BHEs_array
        RUNTIME 50
        EXECUTABLE ogs
        EXECUTABLE_ARGS 3bhes_1U.prj
        WRAPPER time
        TESTER vtkdiff
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
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
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
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
        REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
        PYTHON_PACKAGES "TESPy==0.3.2"
        DIFF_DATA
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_soil temperature_soil 1e-12 1e-13
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE1 temperature_BHE1 1e-9 1e-12
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE2 temperature_BHE2 1e-9 1e-12
        3bhes_1U_ts_10_t_600.000000.vtu 3bhes_1U_ts_10_t_600.000000.vtu temperature_BHE3 temperature_BHE3 1e-9 1e-12
    )
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/T/BHE_1P/BHE_1P.prj RUNTIME 27)
    OgsTest(PROJECTFILE Parabolic/T/BHE_1P/BHE_1P_newton.prj RUNTIME 34)
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE Parabolic/T/BHE_1P/pipe_flow_ebhe.py RUNTIME 107)
endif()
