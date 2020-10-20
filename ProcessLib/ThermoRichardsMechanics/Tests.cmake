if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/LinearMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/FullySaturatedFlowMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE ThermoRichardsMechanics/RichardsFlow2D/RichardsFlow_2d_small.prj)
endif()

AddTest(
    NAME ThermoRichardsMechanics_liakopoulosHM
    PATH ThermoRichardsMechanics/LiakopoulosHM
    EXECUTABLE ogs
    EXECUTABLE_ARGS liakopoulos.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    liakopoulos_t_300.vtu liakopoulos_t_300.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_300.vtu liakopoulos_t_300.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_300.vtu liakopoulos_t_300.vtu saturation saturation 1e-10 1e-12
    liakopoulos_t_600.vtu liakopoulos_t_600.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_600.vtu liakopoulos_t_600.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_600.vtu liakopoulos_t_600.vtu saturation saturation 1e-10 1e-12
    liakopoulos_t_7200.vtu liakopoulos_t_7200.vtu sigma sigma 1e-9 1e-12
    liakopoulos_t_7200.vtu liakopoulos_t_7200.vtu displacement displacement 1e-10 1e-12
    liakopoulos_t_7200.vtu liakopoulos_t_7200.vtu saturation saturation 1e-10 1e-12
)

AddTest(
    NAME ThermoRichardsMechanics_3D_ThermoElastic_Stress_Analysis
    PATH ThermoRichardsMechanics/Simple3DThermoMechanicsFromTM
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    stress_analytical.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu sigma sigma 1e-5 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu displacement displacement 1e-10 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu temperature temperature 1e-10 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu sigma sigma 1e-6 1e-12
    expected_cube_1e3_tm_ts_17_t_72000.000000.vtu cube_1e3_tm_ts_17_t_72000.000000.vtu epsilon epsilon 1e-16 0
)
