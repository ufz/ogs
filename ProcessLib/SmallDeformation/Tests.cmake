if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e2.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/disc_with_hole.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e5.prj LARGE)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e2_quad8_traction_top.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/cube_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/cube_1e0_simple_shear.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MaterialForces/bar.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MaterialForces/bar_3D.prj LARGE)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/cube_1e0_orthotropic_xyz.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/cube_1e0_orthotropic_yzx.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/cube_1e0_orthotropic_zxy.prj)
    OgsTest(PROJECTFILE Mechanics/Burgers/cube_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Burgers/cube_1e3.prj LARGE)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e1.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e3.prj LARGE)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e0_dp.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/ring_plane_strain.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/plain_strain_pipe.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/two_material_gravity.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/two_material_gravity_Emodulus.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/PressureBC/axisymmetric_pipe.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/PressureBC/hollow_sphere.prj LARGE)
    OgsTest(PROJECTFILE Mechanics/Linear/PressureBC/axisymmetric_sphere.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_with_deactivated_hole.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/axisymmetric_sphere_pl.prj LARGE)
endif()

if (OGS_USE_PYTHON)
    OgsTest(PROJECTFILE Mechanics/Linear/PythonPiston/piston.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/PythonHertzContact/hertz_contact.prj RUNTIME 45)
endif()

if (OGS_USE_MPI)
    OgsTest(WRAPPER mpirun -np 4 PROJECTFILE Mechanics/Linear/disc_with_hole.prj)
endif()

# Basic test that MFront works at all.
# Linear elastic, no internal state variables, no external state variables.
AddTest(
    NAME Mechanics_SDL_disc_with_hole_mfront
    PATH Mechanics/Linear/MFront/disc_with_hole
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MFRONT AND NOT OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_expected_pcs_0_ts_4_t_1.000000.vtu disc_with_hole_pcs_0_ts_4_t_1.000000.vtu displacement displacement 1e-16 1e-16
    VIS disc_with_hole_pcs_0_ts_4_t_1.000000.vtu
)

# Tests that internal state variables work correcly.
AddTest(
    NAME Mechanics_DruckerPrager_mfront
    PATH Mechanics/Ehlers/MFront
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_dp.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MFRONT AND NOT OGS_USE_MPI
    # The reference solution has been computed by OGS's Ehlers model.
    # See also the prj file.
    DIFF_DATA
    cube_1e0_dp_ref_created_with_OGS_Ehlers.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu displacement displacement 1e-14 0
    cube_1e0_dp_ref_created_with_OGS_Ehlers.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu sigma sigma 1e-13 0
    cube_1e0_dp_ref_created_with_OGS_Ehlers.vtu cube_1e0_dp_pcs_0_ts_203_t_5.100000.vtu epsilon epsilon 1e-14 0
)

# Tests that axial symmetry works correctly.
# NB: Currently (2018-11-06) the plane strain hypothesis is used within MFront!
AddTest(
    NAME SmallDeformation_ring_plane_strain_axi_mfront
    PATH Mechanics/Linear/MFront/axisymm_ring
    EXECUTABLE ogs
    EXECUTABLE_ARGS ring_plane_strain.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MFRONT AND NOT OGS_USE_MPI
    DIFF_DATA
    ../../ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu displacement displacement 1e-16 0
    ../../ring_plane_strain_1e4_solution.vtu ring_plane_strain_pcs_0_ts_1_t_1.000000.vtu sigma sigma 1e-15 0
)
