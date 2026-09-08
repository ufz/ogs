if(NOT OGS_USE_MFRONT)
    message(
        FATAL_ERROR
            "LargeDeformations process without MFront library is not available. Use enable "
            "MFront in the CMake settings via OGS_USE_MFRONT variable."
    )
endif()

ogs_add_process(LargeDeformation DIRECTORIES . ConstitutiveRelations)
