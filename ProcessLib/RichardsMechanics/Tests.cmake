AddTest(
    NAME RichardsMechanics_square_1e2_gravity
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS gravity.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB gravity_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB gravity_pcs_0_ts_*.vtu sigma sigma 5e-10 0
    GLOB gravity_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB gravity_pcs_0_ts_*.vtu pressure pressure 0 2e-2
    GLOB gravity_pcs_0_ts_*.vtu velocity velocity 1e-7 1e-15
    GLOB gravity_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-5 0
    GLOB gravity_pcs_0_ts_*.vtu NodalForces NodalForces 1e-10 0
    GLOB gravity_pcs_0_ts_24_t_5.000000.vtu pressure pressure 1e-4 1e-15
)

AddTest(
    NAME RichardsMechanics_square_1e2_mechanics_linear
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS mechanics_linear.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB mechanics_linear_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB mechanics_linear_pcs_0_ts_*.vtu sigma sigma 1e-15 0
    GLOB mechanics_linear_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB mechanics_linear_pcs_0_ts_*.vtu pressure pressure 1e-15 1e-15
    GLOB mechanics_linear_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB mechanics_linear_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-15 0
    GLOB mechanics_linear_pcs_0_ts_*.vtu NodalForces NodalForces 1e-15 0
)

AddTest(
    NAME RichardsMechanics_square_1e2_confined_compression
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS confined_compression_fully_saturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu displacement displacement 1e-15 0
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu sigma sigma 1e-15 0
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu pressure pressure 1e-15 1e-15
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-15 0
    GLOB confined_compression_fully_saturated_pcs_0_ts_*.vtu NodalForces NodalForces 1e-15 0
)

AddTest(
    NAME RichardsMechanics_square_1e2_flow_fully_saturated
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_fully_saturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu displacement displacement 2e-14 0
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu sigma sigma 1e-14 0
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu pressure pressure 1e-15 1e-15
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-15 0
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu NodalForces NodalForces 1e-15 0
)

AddTest(
    NAME RichardsMechanics_RichardsFlow_2d_small
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu displacement displacement 2e-14 0
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu sigma sigma 1e-8 0
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu pressure pressure 1e-7 1e-15
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu saturation saturation 1e-11 1e-15
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu saturation_avg saturation_avg 1e-11 1e-15
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-13 0
    GLOB RichardsFlow_2d_small_pcs_0_ts_*.vtu NodalForces NodalForces 1e-9 0
)
AddTest(
    NAME RichardsMechanics_RichardsFlow_2d_quasinewton
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_quasinewton.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu displacement displacement 2e-14 0
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu sigma sigma 1e-8 0
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu epsilon epsilon 1e-15 0
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu pressure pressure 1e-10 1e-15
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu saturation saturation 1e-14 1e-15
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu saturation_avg saturation_avg 1e-14 1e-15
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu velocity velocity 1e-15 1e-15
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu HydraulicFlow HydraulicFlow 1e-15 0
    GLOB RichardsFlow_2d_quasinewton_pcs_0_ts_*.vtu NodalForces NodalForces 1e-15 0
)
AddTest(
    NAME RichardsMechanics_RichardsFlow_2d_richardsflow
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_richardsflow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB RichardsFlow_2d_richardsflow_pcs_0_ts_*.vtu pressure pressure 1e-11 1e-15
    GLOB RichardsFlow_2d_richardsflow_pcs_0_ts_*.vtu saturation saturation 1e-14 1e-15
)
