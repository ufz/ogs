if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated_gas.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated_gas_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/Confined_Compression/HM_confined_compression_gas.prj RUNTIME 50)
    OgsTest(PROJECTFILE TH2M/HM/Confined_Compression/HM_confined_compression_liquid.prj RUNTIME 50)
    OgsTest(PROJECTFILE TH2M/THM/Confined_Compression/THM_confined_compression_gas.prj RUNTIME 55)
    OgsTest(PROJECTFILE TH2M/THM/Confined_Compression/THM_confined_compression_liquid.prj RUNTIME 55)
    OgsTest(PROJECTFILE TH2M/TH/idealGasLaw/compression_gas.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/H2M/Liakopoulos/liakopoulos_TH2M.prj RUNTIME 15)
    OgsTest(PROJECTFILE TH2M/H2M/Liakopoulos/liakopoulos_newton.xml RUNTIME 5)
    OgsTest(PROJECTFILE TH2M/H2/mcWorther/mcWorther_h2.prj RUNTIME 55)
endif()

# TH2M 1d heat diffusion w/ Dirichlet-BC
AddTest(
    NAME TH2M_T_1d_dirichlet
    PATH TH2M/T/T_1d_dirichlet
    EXECUTABLE ogs
    EXECUTABLE_ARGS T_1d_dirichlet.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    displacement displacement 1e-8 1e-8

    # secondary variables
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    velocity_gas velocity_gas 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    velocity_liquid velocity_liquid 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    sigma sigma 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    epsilon epsilon 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    liquid_density liquid_density 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    gas_density gas_density 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    porosity porosity 1e-8 1e-8

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    saturation saturation 1e-8 1e-8
)
AddTest(
    NAME TH2M_T_1d_dirichlet_newton
    PATH TH2M/T/T_1d_dirichlet
    EXECUTABLE ogs
    EXECUTABLE_ARGS T_1d_dirichlet_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 3e-6 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu displacement displacement 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 3e-6 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu sigma sigma 4e-6 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu epsilon epsilon 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu gas_density gas_density 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu porosity porosity 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M 2d linear elastic mechanics w/ neumann BC
AddTest(
    NAME TH2M_M_2d_neumann
    PATH TH2M/M/M_2d_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS M_2d_neumann.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu displacement displacement 1e-8 1e-8

    # secondary variables
    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu sigma sigma 6e-7 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu epsilon epsilon 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu liquid_density liquid_density 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu gas_density gas_density 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu porosity porosity 1e-8 1e-8

    result_TH2M_M_ts_2_t_2.000000.vtu result_TH2M_M_ts_2_t_2.000000.vtu saturation saturation 1e-8 1e-8
)
AddTest(
    NAME TH2M_M_2d_neumann_newton
    PATH TH2M/M/M_2d_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS M_2d_neumann_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu displacement displacement 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu sigma sigma 6e-7 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu epsilon epsilon 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu gas_density gas_density 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu porosity porosity 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M THM point_heatsource benchmark
AddTest(
    NAME TH2M_THM_point_heatsource
    PATH TH2M/THM/sphere
    RUNTIME 40
    EXECUTABLE ogs
    EXECUTABLE_ARGS point_heatsource.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 6e-7 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu displacement displacement 1e-8 1e-8

    # secondary variables
    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 6e-7 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu sigma sigma 3e-7 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu porosity porosity 1e-8 1e-8

    result_point_heatsource_ts_5_t_100000.000000.vtu result_point_heatsource_ts_5_t_100000.000000.vtu saturation saturation 1e-8 1e-8
)
AddTest(
    NAME TH2M_THM_point_heatsource_newton
    PATH TH2M/THM/sphere
    RUNTIME 40
    EXECUTABLE ogs
    EXECUTABLE_ARGS point_heatsource_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 7e-7 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu displacement displacement 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 7e-7 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu sigma sigma 4e-7 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu porosity porosity 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M Thermohydromechanics in a slab
AddTest(
    NAME TH2M_THM_THM_1d_dirichlet
    PATH TH2M/THM/slab
    RUNTIME 15
    EXECUTABLE ogs
    EXECUTABLE_ARGS THM_1d_dirichlet.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 2e-6 2e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu displacement displacement 1e-8 1e-8

    # secondary variables
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 2e-6 2e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu sigma sigma 3e-6 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu porosity porosity 1e-8 1e-8

    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu saturation saturation 1e-8 1e-8
)
AddTest(
    NAME TH2M_THM_THM_1d_dirichlet_newton
    PATH TH2M/THM/slab
    RUNTIME 15
    EXECUTABLE ogs
    EXECUTABLE_ARGS THM_1d_dirichlet_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 3e-6 2e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu displacement displacement 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 3e-6 2e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu sigma sigma 4e-6 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu porosity porosity 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M Heatpipe w/ static gas phase in radial domain
AddTest(
    NAME TH2M_TH_unsaturated_heatpipe_radial
    PATH TH2M/TH_unsaturated/heatpipe_radial_static_gas
    RUNTIME 20
    EXECUTABLE ogs
    EXECUTABLE_ARGS heatpipe_radial_static_gas.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu results_heatpipe_radial_static_gas_ts_0_t_0.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu results_heatpipe_radial_static_gas_ts_4_t_400000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu results_heatpipe_radial_static_gas_ts_8_t_800000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    # secondary variables
    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu sigma sigma 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu epsilon epsilon 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu liquid_density liquid_density 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu gas_density gas_density 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu porosity porosity 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu saturation saturation 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu solid_density solid_density 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu vapour_pressure vapour_pressure 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu xnCG xnCG 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu xmCG xmCG 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu xmWL xmWL 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu k_rel_G k_rel_G 1e-8 1e-8

    results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu results_heatpipe_radial_static_gas_ts_10_t_1000000.000000.vtu k_rel_L k_rel_L 1e-8 1e-8
)

# TH2M Heatpipe w/ static gas phase in a slab
AddTest(
    NAME TH2M_TH_unsaturated_heatpipe_slab
    PATH TH2M/TH_unsaturated/heatpipe_slab_static_gas
    EXECUTABLE ogs
    EXECUTABLE_ARGS heatpipe_slab_static_gas.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu results_heatpipe_slab_static_gas_ts_0_t_0.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu results_heatpipe_slab_static_gas_ts_10_t_10000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    # secondary variables
    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu sigma sigma 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu porosity porosity 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu saturation saturation 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu solid_density solid_density 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu vapour_pressure vapour_pressure 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu xnCG xnCG 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu xmCG xmCG 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu xmWL xmWL 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu k_rel_G k_rel_G 1e-8 1e-8

    results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu results_heatpipe_slab_static_gas_ts_19_t_100000.000000.vtu k_rel_L k_rel_L 1e-8 1e-8
)

# TH2M Heatpipe multiphase in radial domain
AddTest(
    NAME TH2M_TH2_heatpipe_radial
    PATH TH2M/TH2/heatpipe_radial
    RUNTIME 25
    EXECUTABLE ogs
    EXECUTABLE_ARGS heatpipe_radial.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    results_heatpipe_radial_ts_0_t_0.000000.vtu results_heatpipe_radial_ts_0_t_0.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_ts_0_t_0.000000.vtu results_heatpipe_radial_ts_0_t_0.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_0_t_0.000000.vtu results_heatpipe_radial_ts_0_t_0.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_0_t_0.000000.vtu results_heatpipe_radial_ts_0_t_0.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_radial_ts_4_t_400000.000000.vtu results_heatpipe_radial_ts_4_t_400000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_ts_4_t_400000.000000.vtu results_heatpipe_radial_ts_4_t_400000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_4_t_400000.000000.vtu results_heatpipe_radial_ts_4_t_400000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_4_t_400000.000000.vtu results_heatpipe_radial_ts_4_t_400000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_radial_ts_8_t_800000.000000.vtu results_heatpipe_radial_ts_8_t_800000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_ts_8_t_800000.000000.vtu results_heatpipe_radial_ts_8_t_800000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_8_t_800000.000000.vtu results_heatpipe_radial_ts_8_t_800000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_8_t_800000.000000.vtu results_heatpipe_radial_ts_8_t_800000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    # secondary variables
    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu sigma sigma 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu epsilon epsilon 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu liquid_density liquid_density 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu gas_density gas_density 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu porosity porosity 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu saturation saturation 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu solid_density solid_density 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu vapour_pressure vapour_pressure 1e-8 2e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu xnCG xnCG 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu xmCG xmCG 2e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu xmWL xmWL 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu k_rel_G k_rel_G 1e-8 1e-8

    results_heatpipe_radial_ts_10_t_1000000.000000.vtu results_heatpipe_radial_ts_10_t_1000000.000000.vtu k_rel_L k_rel_L 1e-8 1e-8
)

# TH2M Heatpipe multiphase in a slab
AddTest(
    NAME TH2M_TH2_heatpipe_slab
    PATH TH2M/TH2/heatpipe_slab
    RUNTIME 45
    EXECUTABLE ogs
    EXECUTABLE_ARGS heatpipe_slab.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA

    # primary variables
    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8

    # secondary variables
    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu sigma sigma 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu porosity porosity 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu saturation saturation 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu solid_density solid_density 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu vapour_pressure vapour_pressure 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu xnCG xnCG 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu xmCG xmCG 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu xmWL xmWL 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu k_rel_G k_rel_G 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu k_rel_L k_rel_L 1e-8 1e-8
)
