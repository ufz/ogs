if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE LIE/Mechanics/cohesive_zone_load_path.prj RUNTIME 2)
    OgsTest(PROJECTFILE LIE/Mechanics/coulomb_load_path.prj RUNTIME 2)
    OgsTest(PROJECTFILE LIE/Mechanics/elastic_push_pull_two_fractures.prj RUNTIME 1)
    OgsTest(PROJECTFILE LIE/Mechanics/mohr_coulomb_load_path_nu0p3.prj RUNTIME 15)
endif()

# LIE; Small deformation
AddTest(
    NAME LIE_M_single_joint_3D
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint_3D.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    single_joint_3D_ts_1_t_1.000000.vtu single_joint_3D_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    single_joint_3D_ts_1_t_1.000000.vtu single_joint_3D_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
)

AddTest(
    NAME LIE_M_single_joint
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    single_joint_ts_1_t_1.000000.vtu single_joint_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    single_joint_ts_1_t_1.000000.vtu single_joint_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
)

AddTest(
    NAME LIE_M_single_joint_inside
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint_inside.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    single_joint_inside_ts_1_t_1.000000.vtu single_joint_inside_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    single_joint_inside_ts_1_t_1.000000.vtu single_joint_inside_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
)

AddTest(
    NAME LIE_M_two_joints
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS two_joints.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    two_joints_ts_1_t_1.000000.vtu two_joints_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    two_joints_ts_1_t_1.000000.vtu two_joints_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
    two_joints_ts_1_t_1.000000.vtu two_joints_ts_1_t_1.000000.vtu displacement_jump2 displacement_jump2 1e-16 1e-16
)

AddTest(
    NAME LIE_M_single_joint_negative_aperture
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint_negative_aperture.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    single_joint_negative_aperture_ts_1_t_1.000000.vtu single_joint_negative_aperture_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    single_joint_negative_aperture_ts_1_t_1.000000.vtu single_joint_negative_aperture_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
    single_joint_negative_aperture_ts_1_t_1.000000.vtu single_joint_negative_aperture_ts_1_t_1.000000.vtu fracture_aperture fracture_aperture 1e-16 1e-16
    single_joint_negative_aperture_ts_1_t_1.000000.vtu single_joint_negative_aperture_ts_1_t_1.000000.vtu fracture_aperture_avg fracture_aperture_avg 1e-16 1e-16
)

AddTest(
    NAME LIE_M_single_joint_displacement_controlled
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint_displacement_controlled.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    single_joint_displacement_controlled_ts_10_t_1.000000.vtu single_joint_displacement_controlled_ts_10_t_1.000000.vtu displacement displacement 1e-16 0
    single_joint_displacement_controlled_ts_10_t_1.000000.vtu single_joint_displacement_controlled_ts_10_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 0
    single_joint_displacement_controlled_ts_10_t_1.000000.vtu single_joint_displacement_controlled_ts_10_t_1.000000.vtu fracture_aperture fracture_aperture 1e-16 1e-16
    single_joint_displacement_controlled_ts_10_t_1.000000.vtu single_joint_displacement_controlled_ts_10_t_1.000000.vtu fracture_aperture_avg fracture_aperture_avg 1e-16 1e-16
)

AddTest(
    NAME LIE_M_two_cracks_branch_pull
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS two_cracks_branch_pull.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu displacement_jump2 displacement_jump2 1e-16 1e-16
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu sigma sigma 1e-6 1e-6
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu sigma_avg sigma_avg 1e-6 1e-6
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu fracture_stress fracture_stress 1e-6 5e-6
    two_cracks_branch_pull_ts_1_t_1.000000.vtu two_cracks_branch_pull_ts_1_t_1.000000.vtu fracture_stress_avg fracture_stress_avg 1e-6 5e-6
)

AddTest(
    NAME LIE_M_two_cracks_junction_pull
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS two_cracks_junction_pull.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu displacement_jump2 displacement_jump2 1e-16 1e-16
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu displacement_jump3 displacement_jump3 1e-16 1e-16
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu sigma sigma 1e-6 1e-6
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu fracture_stress fracture_stress 1e-6 5e-6
    two_cracks_junction_pull_ts_1_t_1.000000.vtu two_cracks_junction_pull_ts_1_t_1.000000.vtu fracture_stress_avg fracture_stress_avg 1e-6 5e-6
)

AddTest(
    NAME LIE_M_sfrac_q
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS sfrac.q.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu displacement displacement 1e-16 1e-16
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-16 1e-16
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu displacement_jump2 displacement_jump2 1e-16 1e-16
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu displacement_jump3 displacement_jump3 1e-16 1e-16
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu displacement_jump4 displacement_jump4 1e-16 1e-16
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu displacement_jump4 displacement_jump4 1e-16 1e-16
    sfrac_ts_1_t_1.000000.vtu sfrac_ts_1_t_1.000000.vtu fracture_stress_avg fracture_stress_avg 1e-6 5e-6
)
