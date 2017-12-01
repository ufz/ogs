# TES tests
AddTest(
    NAME TES_zeolite_discharge_small
    PATH Parabolic/TES/1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS tes-1D-zeolite-discharge-small.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu pressure pressure 1e-7 5e-9
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu temperature temperature 1e-7 5e-9
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu vapour_partial_pressure vapour_partial_pressure 1e-7 5e-9
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0.100000.vtu solid_density solid_density 1e-7 5e-9
)

AddTest(
    NAME LARGE_TES_zeolite_discharge
    PATH Parabolic/TES/1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS tes-1D-zeolite-discharge-large.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu pressure pressure 1e-3 1e-4
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu temperature temperature 1e-3 1e-4
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu vapour_partial_pressure vapour_partial_pressure 1e-3 1e-4
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1.000000.vtu solid_density solid_density 1e-3 1e-4
)

AddTest(
    NAME LARGE_TES_zeolite_discharge_Newton
    PATH Parabolic/TES/1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS tes-1D-zeolite-discharge-small-newton.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu pressure pressure 1.5e-3 1.5e-3
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu temperature temperature 1.5e-3 1.5e-3
    # tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu vapour_partial_pressure vapour_partial_pressure 1.5e-3 1.5e-3
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_newton_pcs_0_ts_32_t_0.100000.vtu solid_density solid_density 1.5e-3 1.5e-3
)


# SQUARE 1x1 TES TEST -- AXIALLY SYMMETRIC
# test results are compared to 3D simulation on a wedge-shaped domain
AddTest(
    NAME LARGE_TES_inert_axi
    PATH Parabolic/TES/2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS tes-inert-axi.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    # Note: Since the temperature and pressure only vary by a factor of ~ 1.e-6 in x-direction
    # the relative tolerance has to be much smaller than 1.e-6
    DIFF_DATA
    inert-wedge-extracted-surface-t-1s.vtu tes_inert_axi_pcs_0_ts_4_t_1.000000.vtu pressure pressure 1e-12 2e-9
    inert-wedge-extracted-surface-t-1s.vtu tes_inert_axi_pcs_0_ts_4_t_1.000000.vtu temperature temperature 1e-12 2e-9
    inert-wedge-extracted-surface-t-1s.vtu tes_inert_axi_pcs_0_ts_4_t_1.000000.vtu v_mass_frac v_mass_frac 1e-12 2e-9
)
# # WEDGE 1x1 TES TEST -- computes reference results for the above test
# AddTest(
#     NAME TES_inert_wedge
#     PATH Parabolic/TES/2D
#     EXECUTABLE ogs
#     EXECUTABLE_ARGS tes-inert-wedge.prj
# )

# For PETSc/MPI
AddTest(
    NAME TES_zeolite_discharge_small
    PATH Parabolic/TES/1D
    EXECUTABLE_ARGS tes-1D-zeolite-discharge-small.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0_100000_0.vtu pressure pressure 1e-7 5e-9
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0_100000_0.vtu temperature temperature 1e-7 5e-9
    tes_zeolite_discharge_small_ts_19_t_0_100000.vtu tes_zeolite_discharge_small_pcs_0_ts_19_t_0_100000_0.vtu v_mass_frac v_mass_frac 1e-7 5e-9
#        tes_zeolite_discharge_small_ts_19_t_0.100000.vtu solid_density solid_density 1e-7 5e-9
)

AddTest(
    NAME LARGE_TES_zeolite_discharge
    PATH Parabolic/TES/1D
    EXECUTABLE_ARGS tes-1D-zeolite-discharge-large.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000_0.vtu pressure pressure 1e-8 1e-8
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000_0.vtu temperature temperature 1e-8 1e-8
    tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000.vtu tes_zeolite_discharge_large_pcs_0_ts_28_t_1_000000_0.vtu v_mass_frac v_mass_frac 1e-8 1e-8
#        tes_zeolite_discharge_large_ts_28_t_1_0.vtu solid_density solid_density 1e-8 1e-8
)

