if(WIN32)
    # Exceeds library limit of 65535 objects
    set(_static STATIC)
endif()

ogs_add_process(
    LIE_HM
    ${_static}
    DIRECTORIES
    .
    LocalAssembler
    PUBLIC_LIBRARIES
    LIECommon
)
