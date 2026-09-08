if(WIN32)
    # Exceeds library limit of 65535 objects
    set(_static STATIC)
endif()

ogs_add_process(ComponentTransport ${_static} PUBLIC_LIBRARIES ChemistryLib)
