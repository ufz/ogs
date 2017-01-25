# Mechanics; Small deformations, linear (SDL)
AddTest(
    NAME Mechanics_SDL_square_1e0_displacementBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    square_1e0_expected_pcs_0_ts_4_t_1.000000.vtu square_1e0_pcs_0_ts_4_t_1.000000.vtu displacement displacement
)
AddTest(
    NAME Mechanics_SDL_square_1e2_tractionBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    square_1e2_expected_pcs_0_ts_4_t_1.000000.vtu square_1e2_pcs_0_ts_4_t_1.000000.vtu displacement displacement
)
AddTest(
    NAME LARGE_Mechanics_SDL_disc_with_hole
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    disc_with_hole_expected_pcs_0_ts_4_t_1.000000.vtu disc_with_hole_pcs_0_ts_4_t_1.000000.vtu displacement displacement
)
AddTest(
    NAME LARGE_Mechanics_SDL_square_1e5_tractionBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e5.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    square_1e5_expected_pcs_0_ts_4_t_1.000000.vtu square_1e5_pcs_0_ts_4_t_1.000000.vtu displacement displacement
)
AddTest(
    NAME Mechanics_SDL_square_1e2_quad8_traction_topBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_quad8_traction_top.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-14 RELTOL 1e-15
    DIFF_DATA
    expected_square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu displacement displacement
)

# Mechanics; Small deformations, Burgers (SDB)
AddTest(
    NAME Mechanics_SDB_cube_1e0_tractionBC
    PATH Mechanics/Burgers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    cube_1e0_expected_pcs_0_ts_1_t_0.000100.vtu cube_1e0_pcs_0_ts_1_t_0.000100.vtu displacement displacement
    cube_1e0_expected_pcs_0_ts_101_t_1.000000.vtu cube_1e0_pcs_0_ts_101_t_1.000000.vtu displacement displacement
)
AddTest(
    NAME LARGE_Mechanics_SDB_cube_1e3_tractionBC
    PATH Mechanics/Burgers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    cube_1e3_expected_pcs_0_ts_1_t_0.000100.vtu cube_1e3_pcs_0_ts_1_t_0.000100.vtu displacement displacement
    cube_1e3_expected_pcs_0_ts_101_t_1.000000.vtu cube_1e3_pcs_0_ts_101_t_1.000000.vtu displacement displacement
)

# Mechanics; Small deformations, Ehlers (SDE)
AddTest(
    NAME Mechanics_SDE_cube_1e0
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-13
    DIFF_DATA
    expected_cube_1e0_pcs_0_ts_101_t_2.550000.vtu cube_1e0_pcs_0_ts_101_t_2.550000.vtu displacement displacement
    expected_cube_1e0_pcs_0_ts_101_t_2.550000.vtu cube_1e0_pcs_0_ts_101_t_2.550000.vtu sigma_xx sigma_xx
    expected_cube_1e0_pcs_0_ts_203_t_5.100000.vtu cube_1e0_pcs_0_ts_203_t_5.100000.vtu displacement displacement
    expected_cube_1e0_pcs_0_ts_203_t_5.100000.vtu cube_1e0_pcs_0_ts_203_t_5.100000.vtu sigma_xx sigma_xx
)
AddTest(
    NAME Mechanics_SDE_cube_1e1
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e1.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-13
    DIFF_DATA
    expected_cube_1e1_pcs_0_ts_101_t_2.550000.vtu cube_1e1_pcs_0_ts_101_t_2.550000.vtu displacement displacement
    expected_cube_1e1_pcs_0_ts_101_t_2.550000.vtu cube_1e1_pcs_0_ts_101_t_2.550000.vtu sigma_xx sigma_xx
    expected_cube_1e1_pcs_0_ts_203_t_5.100000.vtu cube_1e1_pcs_0_ts_203_t_5.100000.vtu displacement displacement
    expected_cube_1e1_pcs_0_ts_203_t_5.100000.vtu cube_1e1_pcs_0_ts_203_t_5.100000.vtu sigma_xx sigma_xx
)
AddTest(
    NAME LARGE_Mechanics_SDE_cube_1e3
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 3e-13 RELTOL 1e-16
    DIFF_DATA
    expected_cube_1e3_pcs_0_ts_101_t_2.550000.vtu cube_1e3_pcs_0_ts_101_t_2.550000.vtu displacement displacement
    expected_cube_1e3_pcs_0_ts_101_t_2.550000.vtu cube_1e3_pcs_0_ts_101_t_2.550000.vtu sigma_xx sigma_xx
    expected_cube_1e3_pcs_0_ts_203_t_5.100000.vtu cube_1e3_pcs_0_ts_203_t_5.100000.vtu displacement displacement
    expected_cube_1e3_pcs_0_ts_203_t_5.100000.vtu cube_1e3_pcs_0_ts_203_t_5.100000.vtu sigma_xx sigma_xx
)

# Mechanics; Small deformations, Drucker-Prager (SDE-DP)
AddTest(
    NAME Mechanics_PlasticModel_SDE-DP_Ehlers_SpecialCase_DruckerPrager
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_dp.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 5e-15 RELTOL 1e-13
    DIFF_DATA
    expected_cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu displacement displacement
    expected_cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu sigma_xx sigma_xx
)


# SMALL DEFORMATION TEST -- AXIALLY SYMMETRIC
AddTest(
    NAME LARGE_SmallDeformation_ring_plane_strain_axi
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS ring_plane_strain.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 6e-4 RELTOL 1e-4
    DIFF_DATA
    ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu u displacement
    ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu sigma_rr sigma_rr
    ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu sigma_ff sigma_ff
)


