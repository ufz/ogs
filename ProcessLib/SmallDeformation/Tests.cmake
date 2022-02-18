if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e2.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/disc_with_hole.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/ElementDeactivation3D/element_deactivation_M_3D.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e5.prj RUNTIME 200)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e2_quad8_traction_top.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/cube_1e0.prj)
    if(OGS_USE_MKL)
        OgsTest(PROJECTFILE Mechanics/Linear/cube_1e0_mkl.xml)
        set_tests_properties(ogs-Mechanics/Linear/cube_1e0_mkl PROPERTIES
            DEPENDS ogs-Mechanics/Linear/cube_1e0) # Prevent race condition
    endif()
    OgsTest(PROJECTFILE Mechanics/Linear/cube_1e0_simple_shear.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MaterialForces/bar.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MaterialForces/bar_3D.prj RUNTIME 40)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/cube_1e0_orthotropic_xyz.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/cube_1e0_orthotropic_yzx.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/cube_1e0_orthotropic_zxy.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/m3_3Dshearz.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/m3_3Dshearz_rot.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/m3_3Dtopload.prj RUNTIME 8)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/m3_3Dtoploadlc.prj RUNTIME 8)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/square_1e0_orthotropic_xyz.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/square_1e0_orthotropic_45xy_z.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/Orthotropy/square_1e0_orthotropic_y-xz.prj)
    OgsTest(PROJECTFILE Mechanics/Burgers/cube_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Burgers/cube_1e3.prj RUNTIME 85)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e0.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e1.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e3.prj RUNTIME 760)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e0_dp.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/ring_plane_strain.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/plain_strain_pipe.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/two_material_gravity.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/two_material_gravity_Emodulus.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/PressureBC/axisymmetric_pipe.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/PressureBC/hollow_sphere.prj RUNTIME 350)
    OgsTest(PROJECTFILE Mechanics/Linear/PressureBC/axisymmetric_sphere.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_with_deactivated_hole.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/axisymmetric_sphere_pl.prj RUNTIME 365)
    OgsTest(PROJECTFILE Mechanics/InitialStates/into_initial_state.prj)
    OgsTest(PROJECTFILE Mechanics/InitialStates/equilibrium_restart.prj)
    OgsTest(PROJECTFILE Mechanics/InitialStates/non_equilibrium_initial_state.prj)
    OgsTest(PROJECTFILE Mechanics/InitialStates/soil_column_nonequilibrium_sigma_elementwise.prj)
    OgsTest(PROJECTFILE Mechanics/Excavation/time_linear_excavation.prj RUNTIME 4)
endif()

if (OGS_USE_PYTHON)
    OgsTest(PROJECTFILE Mechanics/Linear/PythonPiston/piston.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/PythonHertzContact/hertz_contact.prj RUNTIME 70)
endif()

if (OGS_USE_MPI)
    # OgsTest(WRAPPER mpirun -np 4 PROJECTFILE Mechanics/Linear/disc_with_hole.prj)
    OgsTest(WRAPPER mpirun -np 2 PROJECTFILE Mechanics/InitialStates/soil_column_nonequilibrium_sigma_elementwise.prj)
    OgsTest(PROJECTFILE Mechanics/Ehlers/cube_1e0_SNES.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/square_1e0_SNES.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/cube_1e0_SNES.prj)
endif()

if (OGS_USE_MFRONT)
    OgsTest(PROJECTFILE Mechanics/MohrCoulombAbboSloan/slope.prj RUNTIME 50)
    if(TEST ogs-Mechanics/MohrCoulombAbboSloan/slope)
        set_tests_properties(ogs-Mechanics/MohrCoulombAbboSloan/slope PROPERTIES WILL_FAIL true)
    endif()
    OgsTest(PROJECTFILE Mechanics/MohrCoulombAbboSloan/load_test_mc.prj)
    OgsTest(PROJECTFILE Mechanics/MohrCoulombAbboSloanAnisotropic/triax_1e0_47.prj)
    OgsTest(PROJECTFILE Mechanics/MohrCoulombAbboSloanAnisotropic/triax_ortho_1e0_47.prj)
    OgsTest(PROJECTFILE Mechanics/MohrCoulombAbboSloanAnisotropic/triax_aniso_1e0_47.prj)
    #TODO (naumov) enable when output file format can be specified
    #OgsTest(PROJECTFILE Mechanics/MohrCoulombAbboSloan/oedometer.prj RUNTIME 80)
    OgsTest(PROJECTFILE Mechanics/Linear/MFront/cube_1e0_orthotropic_xyz.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MFront/cube_1e0_orthotropic_yzx.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MFront/cube_1e0_orthotropic_zxy.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MFront/square_1e0_orthotropic_xyz.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MFront/square_1e0_orthotropic_45xy_z.prj)
    OgsTest(PROJECTFILE Mechanics/Linear/MFront/square_1e0_orthotropic_y-xz.prj)
    OgsTest(PROJECTFILE Mechanics/Burgers/cube_1e0_mfront.prj)
    OgsTest(PROJECTFILE Mechanics/Burgers/cube_1e0_mfront_mod.prj)
    OgsTest(PROJECTFILE Mechanics/ModifiedCamClay/square_1e0_shear.prj)
    OgsTest(PROJECTFILE Mechanics/ModifiedCamClay/square_1e0_biax.prj)
    OgsTest(PROJECTFILE Mechanics/ModifiedCamClay/model_triaxtest.prj)
    OgsTest(PROJECTFILE Mechanics/GuentherSalzer/model_triaxtest.prj)

# Linear elastic, no internal state variables, no external state variables.
AddTest(
    NAME Mechanics_SDL_disc_with_hole_mfront
    PATH Mechanics/Linear/MFront/disc_with_hole
    EXECUTABLE ogs
    EXECUTABLE_ARGS disc_with_hole.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    disc_with_hole_expected_ts_4_t_1.000000.vtu disc_with_hole_ts_4_t_1.000000.vtu displacement displacement 1e-16 1e-16
    VIS disc_with_hole_ts_4_t_1.000000.vtu
)

# Tests that internal state variables work correctly.
AddTest(
    NAME Mechanics_DruckerPrager_mfront
    PATH Mechanics/Ehlers/MFront
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Mechanics/Ehlers/MFront
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0_dp.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    # The reference solution has been computed by OGS's Ehlers model.
    # See also the prj file.
    DIFF_DATA
    cube_1e0_dp_ref_created_with_OGS_Ehlers.vtu cube_1e0_dp_ts_203_t_5.100000.vtu displacement displacement 1e-14 0
    cube_1e0_dp_ref_created_with_OGS_Ehlers.vtu cube_1e0_dp_ts_203_t_5.100000.vtu sigma sigma 2e-13 0
    cube_1e0_dp_ref_created_with_OGS_Ehlers.vtu cube_1e0_dp_ts_203_t_5.100000.vtu epsilon epsilon 1e-14 0
)

# Tests that axial symmetry works correctly.
# NB: Currently (2018-11-06) the plane strain hypothesis is used within MFront!
AddTest(
    NAME SmallDeformation_ring_plane_strain_axi_mfront
    PATH Mechanics/Linear/MFront/axisymm_ring
    WORKING_DIRECTORY ${Data_SOURCE_DIR}/Mechanics/Linear/MFront/axisymm_ring
    EXECUTABLE ogs
    EXECUTABLE_ARGS ring_plane_strain.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    ../../ring_plane_strain_ts_1_t_1.000000.vtu ring_plane_strain_ts_1_t_1.000000.vtu displacement displacement 1e-16 0
    ../../ring_plane_strain_ts_1_t_1.000000.vtu ring_plane_strain_ts_1_t_1.000000.vtu sigma sigma 1e-15 0
)

endif()

AddTest(
    NAME Mechanics_m1_1Dload
    PATH Mechanics/m1_1Dload
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_1Dload.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    m1_1Dload_ts_1_t_1.000000.vtu m1_1Dload_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_1Dlozenge
    PATH Mechanics/m1_1Dlozenge
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_1Dlozenge.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    m1_1Dlozenge_ts_1_t_1.000000.vtu m1_1Dlozenge_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_2Dload
    PATH Mechanics/m1_2Dload
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_2Dload.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    m1_2Dload_ts_1_t_1.000000.vtu m1_2Dload_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_3Dbottom
    PATH Mechanics/m1_3Dbottom
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_3Dbottom.prj
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MFRONT AND NOT OGS_USE_MPI
    DIFF_DATA
    m1_3Dbottom_ts_1_t_1.000000.vtu m1_3Dbottom_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_3Dgravity
    PATH Mechanics/m1_3Dgravity
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_3Dgravity.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    m1_3Dgravity_ts_1_t_1.000000.vtu m1_3Dgravity_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_3Dload
    PATH Mechanics/m1_3Dload
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_3Dload.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    m1_3Dload_ts_1_t_1.000000.vtu m1_3Dload_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_3Dsquare
    PATH Mechanics/m1_3Dsquare
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_3Dsquare.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 27
    DIFF_DATA
    m1_3Dsquare_ts_1_t_1.000000.vtu m1_3Dsquare_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

AddTest(
    NAME Mechanics_m1_3Dtopload
    PATH Mechanics/m1_3Dtopload
    EXECUTABLE ogs
    EXECUTABLE_ARGS m1_3Dtopload.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 7
    DIFF_DATA
    m1_3Dtopload_ts_1_t_1.000000.vtu m1_3Dtopload_ts_1_t_1.000000.vtu displacement displacement 10e-12 0.0
)

# Tests for Principal Stress Output
AddTest(
    NAME Mechanics_hollow_sphere
    PATH Mechanics/Linear/PrincipalStress
    EXECUTABLE ogs
    EXECUTABLE_ARGS sphere.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    output_ts_1_t_1.000000.vtu output_ts_1_t_1.000000.vtu displacement displacement 0 1e-12
    output_ts_1_t_1.000000.vtu output_ts_1_t_1.000000.vtu principal_stress_values principal_stress_values 0 1e-10
)

AddTest(
    NAME Mechanics_creep_with_heterogeneous_reference_temperature
    PATH Mechanics/CreepWithHeterogeneousReferenceTemperature
    EXECUTABLE ogs
    EXECUTABLE_ARGS arehs-salt-M_gravity_only_element_refT.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 130
    DIFF_DATA
    arehs-salt-M_gravity_only_ts_10_t_31535999999.999996.vtu arehs-salt-M_gravity_only_ts_10_t_31535999999.999996.vtu displacement displacement 1.e-10 1.e-10
    arehs-salt-M_gravity_only_ts_10_t_31535999999.999996.vtu arehs-salt-M_gravity_only_ts_10_t_31535999999.999996.vtu sigma sigma 1.e-9 5.e-7
)

NotebookTest(NOTEBOOKFILE Notebooks/SimpleMechanics.ipynb RUNTIME 10)
