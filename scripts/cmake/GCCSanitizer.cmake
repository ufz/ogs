option(OGS_ADDRESS_SANITIZER OFF "Use GCCs AddressSanitizer")
option(OGS_UNDEFINED_BEHAVIOR_SANITIZER OFF
       "Use GCCs UndefinedBehaviorSanitizer"
)

if(OGS_ADDRESS_SANITIZER)
    set(SANITIZE_FLAG_VALUE "address")
    set(ADDITIONAL_FLAGS "-fno-omit-frame-pointer")
endif()

if(OGS_UNDEFINED_BEHAVIOR_SANITIZER)
    set(SANITIZE_FLAG_VALUE
        "${SANITIZE_FLAG_VALUE},undefined,unreachable,integer-divide-by-zero,vla-bound,bounds,null"
    )
endif()

if(DEFINED SANITIZE_FLAG_VALUE)
    add_compile_options(-fsanitize=${SANITIZE_FLAG_VALUE} ${ADDITIONAL_FLAGS})
    link_libraries(-fsanitize=${SANITIZE_FLAG_VALUE} ${ADDITIONAL_FLAGS})
endif()
