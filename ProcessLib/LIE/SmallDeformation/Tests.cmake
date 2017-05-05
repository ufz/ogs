# LIE; Small deformation
AddTest(
    NAME LIE_M_single_joint
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    single_joint_expected_pcs_0_ts_1_t_1.000000.vtu single_joint_pcs_0_ts_1_t_1.000000.vtu displacement displacement
    single_joint_expected_pcs_0_ts_1_t_1.000000.vtu single_joint_pcs_0_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1
)

AddTest(
    NAME LIE_M_single_joint_inside
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_joint_inside.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    single_joint_inside_expected_pcs_0_ts_1_t_1.000000.vtu single_joint_inside_pcs_0_ts_1_t_1.000000.vtu displacement displacement
    single_joint_inside_expected_pcs_0_ts_1_t_1.000000.vtu single_joint_inside_pcs_0_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1
)

AddTest(
    NAME LIE_M_two_joints
    PATH LIE/Mechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS two_joints.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-16 RELTOL 1e-16
    DIFF_DATA
    two_joints_expected_pcs_0_ts_1_t_1.000000.vtu two_joints_pcs_0_ts_1_t_1.000000.vtu displacement displacement
    two_joints_expected_pcs_0_ts_1_t_1.000000.vtu two_joints_pcs_0_ts_1_t_1.000000.vtu displacement_jump1 displacement_jump1
    two_joints_expected_pcs_0_ts_1_t_1.000000.vtu two_joints_pcs_0_ts_1_t_1.000000.vtu displacement_jump2 displacement_jump2
)
