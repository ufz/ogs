if(WIN32 OR NOT CCACHE_TOOL_PATH OR OGS_DISABLE_CCACHE)
    return()
endif()

# Check ccache version
set(CCACHE_VERSION_REQUIRED 3.2.0)
execute_process(COMMAND ${CCACHE_TOOL_PATH} --version
    OUTPUT_VARIABLE CCACHE_VERSION
)
set(CCACHE_VERSION_GREATER_EQUAL_3_2_5 ON)
if("${CCACHE_VERSION}" MATCHES "ccache version ([0-9]\\.[0-9]\\.[0-9])")
    if(${CMAKE_MATCH_1} VERSION_LESS ${CCACHE_VERSION_REQUIRED})
        message(STATUS "CCache outdated. Installed: ${CMAKE_MATCH_1}, \
            required: ${CCACHE_VERSION_REQUIRED}. Caching disabled.")
        return()
    endif()
    if(${CMAKE_MATCH_1} VERSION_LESS 3.2.5)
        set(CCACHE_VERSION_GREATER_EQUAL_3_2_5 OFF)
    endif()
endif()

# Set ccache as the compiler launcher
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)

if(COMPILER_IS_CLANG)
    add_compile_options(-Qunused-arguments)
endif()

# Check ccache pre-compiled headers config
execute_process(COMMAND ${CCACHE_TOOL_PATH} -p
    OUTPUT_VARIABLE CCACHE_CONFIG
    ERROR_VARIABLE CCACHE_CONFIG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

string(REGEX MATCH ".*time_macros.*"
    COTIRE_CCACHE_CONFIG_TIME_MACROS ${CCACHE_CONFIG}
)

string(REGEX MATCH ".*pch_defines.*"
    COTIRE_CCACHE_CONFIG_PCH_DEFINES ${CCACHE_CONFIG}
)

# Regex should be "sloppiness.*time_macros.*pch_defines.*" but due to bug fixed
# in https://ccache.samba.org/releasenotes.html#_ccache_3_2_5 we have to leave
# out pch_defines if ccache version is older than 3.2.5.
# Ubuntu 16.04 comes with ccache 3.2.4 ...
if(NOT COTIRE_CCACHE_CONFIG_TIME_MACROS OR (CCACHE_VERSION_GREATER_EQUAL_3_2_5 AND NOT COTIRE_CCACHE_CONFIG_PCH_DEFINES))
    message(FATAL_ERROR "CCache configuration does not set sloppiness to pch_defines,time_macros. \
    Current options are: '${CCACHE_CONFIG}'. \
    See https://docs.opengeosys.org/docs/devguide/advanced/using-ccache")
endif()
