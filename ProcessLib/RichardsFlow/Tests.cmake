if (NOT (OGS_USE_MPI OR OGS_USE_LIS))
    # Comparison test for richards mechanics w/o deformations.
    if(NOT ENABLE_ASAN)
        OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_richardsflow.prj RUNTIME 2)
    endif()
    OgsTest(PROJECTFILE Parabolic/Richards/flow_fully_saturated.prj)
endif()

if (NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE Parabolic/Richards/richards-flow.py RUNTIME 5)
endif()

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_ogs5
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_compare_ogs5.prj
    TESTER vtkdiff
    DIFF_DATA
    h_us_quad_1000.vtu richards_ts_100_t_100.000000.vtu PRESSURE1 pressure 1e-1 1e-1
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small.prj
    TESTER vtkdiff
    RUNTIME 4
    DIFF_DATA
    ref_t_1600.000000.vtu richards_ts_1100_t_1600.000000.vtu pressure pressure 1e-8 1e-8
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
)

AddTest(
    NAME 2D_RichardsFlow_h_us_quad
    PATH Parabolic/Richards
    RUNTIME 18
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_large.prj
    TESTER vtkdiff
    DIFF_DATA
    ref_t_20000.000000.vtu richards_ts_18200_t_20000.000000.vtu pressure pressure 1e-8 1e-8
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
)

if(NOT OGS_USE_MPI)
    OgsTest(
        PROJECTFILE Parabolic/Richards/RichardsFlow_2d_small_PID_adaptive_dt.prj
    )
endif()

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small_iteration_adaptive_dt.prj
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    RUNTIME 8
    richards_pcs_PID_adaptive_dt_t_1600.vtu 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt_t_1600.000000.vtu 1e-8 1e-9
    richards_pcs_PID_adaptive_dt_t_1600.vtu 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt_t_1600.000000.vtu 1e-8 1e-9
)

# DEPENDS for preventing race condition writing to the same pvd-file
AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt_PiecewiseLinear
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS iteration_adaptive_dt_PiecewiseLinear.xml
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    RUNTIME 8
    DEPENDS ogs-2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt
    richards_pcs_PID_adaptive_dt_t_1600.vtu 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt_t_1600.000000.vtu 1e-8 1e-9
    richards_pcs_PID_adaptive_dt_t_1600.vtu 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt_t_1600.000000.vtu 1e-8 1e-9
)

#PETSc/MPI
#AddTest(
#    NAME 2D_RichardsFlow_h_us_quad_small_PID_adaptive_dt
#    PATH Parabolic/Richards
#    EXECUTABLE_ARGS RichardsFlow_2d_small_PID_adaptive_dt.prj
#    WRAPPER mpirun
#    WRAPPER_ARGS -np 1
#    TESTER vtkdiff
#    REQUIREMENTS OGS_USE_MPI
#    RUNTIME 220
#    DIFF_DATA
#    ref_t_1600.000000.vtu richards_ts_803_t_1600_000000_0.vtu pressure pressure 1e-8 1e-3
#)

# Comparison test for richards mechanics w/o deformations.
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_richardsflow_mpi.xml
            WRAPPER mpirun -np 1
    )
endif()
