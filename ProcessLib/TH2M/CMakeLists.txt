if(WIN32)
    # Exceeds library limit of 65535 objects
    set(_static STATIC)
endif()

ogs_add_process(TH2M ${_static} DIRECTORIES . ConstitutiveRelations)
