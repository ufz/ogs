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

    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu
    temperature_interpolated temperature_interpolated 1e-8 1e-8
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

    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_0_t_0.000000.vtu results_heatpipe_slab_ts_0_t_0.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8


    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_10_t_10000.000000.vtu results_heatpipe_slab_ts_10_t_10000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8


    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu displacement displacement 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8

    results_heatpipe_slab_ts_19_t_100000.000000.vtu results_heatpipe_slab_ts_19_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
)
