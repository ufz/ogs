# Mechanics; Small deformations, linear (SDL)
AddTest(
    NAME Mechanics_SDL_square_1e0_displacementBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1e0_expected_pcs_0_ts_4_t_1.000000.vtu square_1e0_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-16 1e-16
    square_1e0_expected_pcs_0_ts_4_t_1.000000.vtu square_1e0_pcs_0_ts_4_t_1.000000.vtu sigma sigma 1e-16 1e-16
)
AddTest(
    NAME Mechanics_SDL_square_1e2_tractionBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1e2_expected_pcs_0_ts_4_t_1.000000.vtu square_1e2_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-15 0
    square_1e2_expected_pcs_0_ts_4_t_1.000000.vtu square_1e2_pcs_0_ts_4_t_1.000000.vtu sigma sigma 1e-15 0
)
AddTest(
    NAME LARGE_Mechanics_SDL_disc_with_hole
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_expected_pcs_0_ts_4_t_1.000000.vtu disc_with_hole_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-16 1e-16
    VIS disc_with_hole_pcs_0_ts_4_t_1.000000.vtu
)
AddTest(
    NAME LARGE_Mechanics_SDL_square_1e5_tractionBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e5.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    square_1e5_expected_pcs_0_ts_4_t_1.000000.vtu square_1e5_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-16 1e-16
    square_1e5_expected_pcs_0_ts_4_t_1.000000.vtu square_1e5_pcs_0_ts_4_t_1.000000.vtu sigma sigma 1e-16 1e-16
)
AddTest(
    NAME Mechanics_SDL_square_1e2_quad8_traction_topBC
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2_quad8_traction_top.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu displacement displacement 2e-14 1e-15
    expected_square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu square_1e2_quad8_traction_topBC_pcs_0_ts_4_t_1.000000.vtu sigma sigma 2e-14 1e-15
)
# Mechanics; Small deformations, linear (SDL); Nodal forces
AddTest(
    NAME Mechanics_SDL_cube_1e0_nodal_forces
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e0_pcs_0_ts_4_t_1.000000.vtu cube_1e0_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-14 1e-15
    expected_cube_1e0_pcs_0_ts_4_t_1.000000.vtu cube_1e0_pcs_0_ts_4_t_1.000000.vtu NodalForces NodalForces 1e-14 1e-15
    expected_cube_1e0_pcs_0_ts_4_t_1.000000.vtu cube_1e0_pcs_0_ts_4_t_1.000000.vtu sigma sigma 1e-14 1e-15
    expected_cube_1e0_pcs_0_ts_4_t_1.000000.vtu cube_1e0_pcs_0_ts_4_t_1.000000.vtu epsilon epsilon 1e-14 1e-15
)
AddTest(
    NAME Mechanics_SDL_cube_1e0_simple_shear_nodal_forces
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_simple_shear.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-14 1e-15
    expected_cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu NodalForces NodalForces 1e-14 1e-15
    expected_cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu sigma sigma 1e-14 1e-15
    expected_cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu cube_1e0_simple_shear_pcs_0_ts_4_t_1.000000.vtu epsilon epsilon 1e-14 1e-15
)

# Mechanics; Small deformations, linear (SDL); Material forces
AddTest(
    NAME Mechanics_SDL_material_forces_bar
    PATH Mechanics/Linear/MaterialForces
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_bar_out_pcs_0_ts_2_t_1.000000.vtu bar_out_pcs_0_ts_2_t_1.000000.vtu displacement displacement 1e-14 1e-15
    expected_bar_out_pcs_0_ts_2_t_1.000000.vtu bar_out_pcs_0_ts_2_t_1.000000.vtu MaterialForces MaterialForces 1e-14 1e-15
    expected_bar_out_pcs_0_ts_2_t_1.000000.vtu bar_out_pcs_0_ts_2_t_1.000000.vtu free_energy_density free_energy_density 1e-14 1e-15
)

# Mechanics; Small deformations, linear (SDL); Material forces
AddTest(
    NAME Mechanics_SDL_material_forces_bar_3D
    PATH Mechanics/Linear/MaterialForces
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar_3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_bar_3D_out_pcs_0_ts_2_t_1.000000.vtu bar_3D_out_pcs_0_ts_2_t_1.000000.vtu displacement displacement 1e-14 1e-15
    expected_bar_3D_out_pcs_0_ts_2_t_1.000000.vtu bar_3D_out_pcs_0_ts_2_t_1.000000.vtu MaterialForces MaterialForces 1e-14 1e-15
    expected_bar_3D_out_pcs_0_ts_2_t_1.000000.vtu bar_3D_out_pcs_0_ts_2_t_1.000000.vtu free_energy_density free_energy_density 1e-14 1e-15
)

# Mechanics; Small deformations, Burgers (SDB)
AddTest(
    NAME Mechanics_SDB_cube_1e0_tractionBC
    PATH Mechanics/Burgers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e0_expected_pcs_0_ts_1_t_0.000100.vtu cube_1e0_pcs_0_ts_1_t_0.000100.vtu displacement displacement 1e-16 1e-16
    cube_1e0_expected_pcs_0_ts_101_t_1.000000.vtu cube_1e0_pcs_0_ts_101_t_1.000000.vtu displacement displacement 1e-16 1e-16
)
AddTest(
    NAME LARGE_Mechanics_SDB_cube_1e3_tractionBC
    PATH Mechanics/Burgers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    cube_1e3_expected_pcs_0_ts_1_t_0.000100.vtu cube_1e3_pcs_0_ts_1_t_0.000100.vtu displacement displacement 1e-16 1e-16
    cube_1e3_expected_pcs_0_ts_101_t_1.000000.vtu cube_1e3_pcs_0_ts_101_t_1.000000.vtu displacement displacement 1e-16 1e-16
)

# Mechanics; Small deformations, Ehlers (SDE)
AddTest(
    NAME Mechanics_SDE_cube_1e0
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e0_pcs_0_ts_101_t_2.550000.vtu cube_1e0_pcs_0_ts_101_t_2.550000.vtu displacement displacement 1e-17 0
    expected_cube_1e0_pcs_0_ts_101_t_2.550000.vtu cube_1e0_pcs_0_ts_101_t_2.550000.vtu sigma sigma 1e-15 0
    expected_cube_1e0_pcs_0_ts_203_t_5.100000.vtu cube_1e0_pcs_0_ts_203_t_5.100000.vtu displacement displacement 1e-17 0
    expected_cube_1e0_pcs_0_ts_203_t_5.100000.vtu cube_1e0_pcs_0_ts_203_t_5.100000.vtu sigma sigma 2e-15 0
)
AddTest(
    NAME Mechanics_SDE_cube_1e1
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e1.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e1_pcs_0_ts_101_t_2.550000.vtu cube_1e1_pcs_0_ts_101_t_2.550000.vtu displacement displacement 1e-11 0
    expected_cube_1e1_pcs_0_ts_101_t_2.550000.vtu cube_1e1_pcs_0_ts_101_t_2.550000.vtu sigma sigma 1e-11 0
    expected_cube_1e1_pcs_0_ts_203_t_5.100000.vtu cube_1e1_pcs_0_ts_203_t_5.100000.vtu displacement displacement 1e-11 0
    expected_cube_1e1_pcs_0_ts_203_t_5.100000.vtu cube_1e1_pcs_0_ts_203_t_5.100000.vtu sigma sigma 1e-11 0
)
AddTest(
    NAME LARGE_Mechanics_SDE_cube_1e3
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e3.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e3_pcs_0_ts_101_t_2.550000.vtu cube_1e3_pcs_0_ts_101_t_2.550000.vtu displacement displacement 3e-13 1e-16
    expected_cube_1e3_pcs_0_ts_101_t_2.550000.vtu cube_1e3_pcs_0_ts_101_t_2.550000.vtu sigma sigma 3e-13 1e-16
    expected_cube_1e3_pcs_0_ts_203_t_5.100000.vtu cube_1e3_pcs_0_ts_203_t_5.100000.vtu displacement displacement 3e-13 1e-16
    expected_cube_1e3_pcs_0_ts_203_t_5.100000.vtu cube_1e3_pcs_0_ts_203_t_5.100000.vtu sigma sigma 3e-13 1e-16
)

# Mechanics; Small deformations, Drucker-Prager (SDE-DP)
AddTest(
    NAME Mechanics_PlasticModel_SDE-DP_Ehlers_SpecialCase_DruckerPrager
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_dp.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu displacement displacement 1e-14 0
    expected_cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu sigma sigma 1e-14 0
    expected_cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu epsilon epsilon 1e-14 0
)


# SMALL DEFORMATION TEST -- AXIALLY SYMMETRIC
AddTest(
    NAME SmallDeformation_ring_plane_strain_axi
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS ring_plane_strain.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu displacement displacement 6e-4 0
    ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu sigma sigma 6e-4 0
)

#With PETSc
AddTest(
    NAME Parallel_Mechanics_SDL_disc_with_hole
    PATH Mechanics/Linear
    EXECUTABLE_ARGS disc_with_hole.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_pcs_0_ts_4_t_1_000000_0.vtu disc_with_hole_pcs_0_ts_4_t_1_000000_0.vtu displacement displacement 1e-16 1e-16
    VIS disc_with_hole_pcs_0_ts_0_t_0_000000.pvtu
)

# Pressure boundary condition
AddTest(
    NAME LARGE_SmallDeformation_PressureBC_hollow_sphere
    PATH Mechanics/Linear/PressureBC
    EXECUTABLE ogs
    EXECUTABLE_ARGS hollow_sphere.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_LIS
    DIFF_DATA
    expected_hollow_sphere_pcs_0_ts_1_t_1.000000.vtu hollow_sphere_pcs_0_ts_1_t_1.000000.vtu displacement displacement 1e-15 0
    expected_hollow_sphere_pcs_0_ts_1_t_1.000000.vtu hollow_sphere_pcs_0_ts_1_t_1.000000.vtu sigma sigma 1e-15 0
)

# Pressure boundary condition: elastic pipe plain strain
AddTest(
    NAME SmallDeformation_PressureBC_elastic_pipe_plain_strain
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS plain_strain_pipe.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_plain_strain_pipe_pcs_0_ts_1_t_1.000000.vtu plain_strain_pipe_pcs_0_ts_1_t_1.000000.vtu displacement displacement 2e-11 0
    ref_plain_strain_pipe_pcs_0_ts_1_t_1.000000.vtu plain_strain_pipe_pcs_0_ts_1_t_1.000000.vtu sigma sigma 2e-11 0
)

# Pressure boundary condition: elastic pipe axisymmetric
AddTest(
    NAME SmallDeformation_PressureBC_elastic_pipe_axisymmetric
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS axisymmetric_pipe.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ref_axisymmteric_pipe_pcs_0_ts_1_t_1.000000.vtu axisymmteric_pipe_pcs_0_ts_1_t_1.000000.vtu displacement displacement 1e-11 0
    ref_axisymmteric_pipe_pcs_0_ts_1_t_1.000000.vtu axisymmteric_pipe_pcs_0_ts_1_t_1.000000.vtu sigma sigma 1e-11 0
)

# Pressure boundary condition: elastic sphere axisymmetric
AddTest(
    NAME SmallDeformation_PressureBC_elastic_sphere_axisymmetric
    PATH Mechanics/Linear
    EXECUTABLE ogs
    EXECUTABLE_ARGS axisymmetric_sphere.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_LIS
    DIFF_DATA
    ref_axisymmteric_sphere_pcs_0_ts_1_t_1.000000.vtu axisymmteric_sphere_pcs_0_ts_1_t_1.000000.vtu displacement displacement 1e-15 0
    ref_axisymmteric_sphere_pcs_0_ts_1_t_1.000000.vtu axisymmteric_sphere_pcs_0_ts_1_t_1.000000.vtu sigma sigma 1e-15 0
)

# Pressure boundary condition: plastic sphere axisymmetric
AddTest(
    NAME SmallDeformation_PressureBC_plastic_sphere_axisymmetric
    PATH Mechanics/Ehlers
    EXECUTABLE ogs
    EXECUTABLE_ARGS axisymmetric_sphere_pl.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_LIS
    DIFF_DATA
    ref_axisymmetric_sphere_pcs_0_ts_100_t_1.000000.vtu axisymmetric_sphere_pcs_0_ts_100_t_1.000000.vtu displacement displacement 1e-15 0
    ref_axisymmetric_sphere_pcs_0_ts_100_t_1.000000.vtu axisymmetric_sphere_pcs_0_ts_100_t_1.000000.vtu sigma sigma 1e-15 0
)

# Two materials in gravity field
AddTest(
    NAME SmallDeformation_SDL_two_material_gravity
    PATH Mechanics/Linear
    EXECUTABLE_ARGS two_material_gravity.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_two_material_gravity.vtu two_material_gravity_pcs_0_ts_1_t_1.000000.vtu displacement displacement 5e-14 0
    expected_two_material_gravity.vtu two_material_gravity_pcs_0_ts_1_t_1.000000.vtu sigma sigma 5e-14 0
    expected_two_material_gravity.vtu two_material_gravity_pcs_0_ts_1_t_1.000000.vtu epsilon epsilon 5e-14 0
)
