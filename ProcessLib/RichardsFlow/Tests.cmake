if (NOT (OGS_USE_MPI OR OGS_USE_LIS))
    # Comparison test for richards mechanics w/o deformations.
    if(NOT ENABLE_ASAN)
        OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_richardsflow.prj RUNTIME 2)
    endif()
    OgsTest(PROJECTFILE Parabolic/Richards/flow_fully_saturated.prj)
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE Parabolic/Richards/richards-flow.py RUNTIME 5)
endif()
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/Richards/RichardsFlow_2d_compare_ogs5.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/Richards/RichardsFlow_2d_small.prj RUNTIME 4)
    OgsTest(PROJECTFILE Parabolic/Richards/RichardsFlow_2d_large.prj RUNTIME 18)
endif()

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/Richards/RichardsFlow_2d_small_PID_adaptive_dt.prj)
endif()
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE Parabolic/Richards/RichardsFlow_2d_small_iteration_adaptive_dt.prj
        RUNTIME 24
    )
    OgsTest(
        PROJECTFILE Parabolic/Richards/RichardsFlow_2d_small_iteration_adaptive_dt.prj
        PATCH_FILES iteration_adaptive_dt_PiecewiseLinear.xml
        NAME_SUFFIX PiecewiseLinear
        RUNTIME 24
    )
endif()

# Comparison test for richards mechanics w/o deformations.
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_richardsflow_mpi.xml
            WRAPPER mpirun -np 1
    )
endif()
