if(NOT OGS_USE_MPI)
    # TODO: re-enable this, see #3178. OgsTest(PROJECTFILE
    # Parabolic/RichardsComponentTransport/Padilla/Padilla_NaCl1/Padilla_NaCl1.prj
    # RUNTIME 7)
    OgsTest(
        PROJECTFILE
            Parabolic/RichardsComponentTransport/Padilla/Padilla_NaCl1/Padilla_NaCl1_quadratic.prj
        RUNTIME 7
    )
    OgsTest(
        PROJECTFILE
            Parabolic/RichardsComponentTransport/Padilla/Padilla_NaCl6/Padilla_NaCl6.prj
        RUNTIME 13
    )
endif()
