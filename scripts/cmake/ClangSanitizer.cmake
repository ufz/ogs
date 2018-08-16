if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.6")
    return()
endif()

option(OGS_ADDRESS_SANITIZER OFF "Use Clangs AddressSanitizer")
option(OGS_UNDEFINED_BEHAVIOR_SANITIZER OFF "Use Clangs UndefinedBehaviorSanitizer")

if(OGS_ADDRESS_SANITIZER)
    set(SANITIZE_FLAG_VALUE "address")
    set(ADDITIONAL_FLAGS "-fno-omit-frame-pointer")
endif()

if(OGS_UNDEFINED_BEHAVIOR_SANITIZER)
    set(SANITIZE_FLAG_VALUE "${SANITIZE_FLAG_VALUE},undefined,integer;-fsanitize-blacklist=${CMAKE_CURRENT_SOURCE_DIR}/scripts/test/clang_sanitizer_blacklist.txt")
endif()

if(DEFINED SANITIZE_FLAG_VALUE)
    add_compile_options(-fsanitize=${SANITIZE_FLAG_VALUE} ${ADDITIONAL_FLAGS})
endif()
