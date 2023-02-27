if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE LargeDeformation/RigidBody/square_1e0.prj RUNTIME 1)
endif()
