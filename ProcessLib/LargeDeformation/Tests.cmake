if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE LargeDeformation/RigidBody/square_1e0.prj RUNTIME 1)
    NotebookTest(NOTEBOOKFILE LargeDeformation/RigidBody/RigidBody.md RUNTIME 15)
endif()
