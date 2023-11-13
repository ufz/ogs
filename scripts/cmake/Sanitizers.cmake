# Initial implementation from Professional CMake, 16th Edition, by Craig Scott
option(ENABLE_ASAN "Enable AddressSanitizer" NO)
if(MSVC)
    if(ENABLE_ASAN)
        string(REPLACE "/RTC1" "" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
        string(REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG
                       "${CMAKE_CXX_FLAGS_DEBUG}"
        )
        add_compile_options(
            /fsanitize=address /fsanitize-address-use-after-return
        )
    endif()
elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    option(ENABLE_LSAN "Enable LeakSanitizer" NO)
    option(ENABLE_TSAN "Enable ThreadSanitizer" NO)
    option(ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" NO)
    if(NOT APPLE)
        option(ENABLE_MSAN "Enable MemorySanitizer" NO)
    endif()
    if((ENABLE_ASAN AND (ENABLE_TSAN OR ENABLE_MSAN))
       OR (ENABLE_LSAN AND (ENABLE_TSAN OR ENABLE_MSAN)) OR (ENABLE_TSAN
                                                             AND ENABLE_MSAN)
    )
        message(
            FATAL_ERROR
                "Invalid sanitizer combination:\n"
                "  ENABLE_ASAN:  ${ENABLE_ASAN}\n"
                "  ENABLE_LSAN:  ${ENABLE_LSAN}\n"
                "  ENABLE_TSAN:  ${ENABLE_TSAN}\n"
                "  ENABLE_MSAN:  ${ENABLE_MSAN}"
        )
    endif()

    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        set(_ubsan_options
            "$<COMMA>unreachable$<COMMA>integer-divide-by-zero$<COMMA>vla-bound$<COMMA>bounds$<COMMA>null"
        )
    endif()
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(_ubsan_options "$<COMMA>integer")
    endif()
    add_compile_options(
        $<$<BOOL:${ENABLE_ASAN}>:-fsanitize=address>
        $<$<BOOL:${ENABLE_ASAN}>:-fno-omit-frame-pointer>
        $<$<BOOL:${ENABLE_LSAN}>:-fsanitize=leak>
        $<$<BOOL:${ENABLE_MSAN}>:-fsanitize=memory>
        $<$<BOOL:${ENABLE_TSAN}>:-fsanitize=thread>
        $<$<BOOL:${ENABLE_UBSAN}>:-fsanitize=undefined${_ubsan_options}>
        $<$<AND:$<BOOL:${ENABLE_UBSAN}>,$<C_COMPILER_ID:Clang>>:-fsanitize-blacklist=${CMAKE_CURRENT_SOURCE_DIR}/scripts/test/clang_sanitizer_blacklist.txt>
    )
    add_link_options(
        $<$<BOOL:${ENABLE_ASAN}>:-fsanitize=address>
        $<$<BOOL:${ENABLE_ASAN}>:-fno-omit-frame-pointer>
        $<$<BOOL:${ENABLE_LSAN}>:-fsanitize=leak>
        $<$<BOOL:${ENABLE_MSAN}>:-fsanitize=memory>
        $<$<BOOL:${ENABLE_TSAN}>:-fsanitize=thread>
        $<$<BOOL:${ENABLE_UBSAN}>:-fsanitize=undefined${_ubsan_options}>
        $<$<AND:$<BOOL:${ENABLE_UBSAN}>,$<C_COMPILER_ID:Clang>>:-fsanitize-blacklist=${CMAKE_CURRENT_SOURCE_DIR}/scripts/test/clang_sanitizer_blacklist.txt>
    )
endif()
