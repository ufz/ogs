if(NOT CCACHE_TOOL_PATH)
    return()
endif()

# Check ccache version
set(CCACHE_VERSION_REQUIRED 3.2.0)
execute_process(COMMAND ${CCACHE_TOOL_PATH} --version
    OUTPUT_VARIABLE CCACHE_VERSION
)
if("${CCACHE_VERSION}" MATCHES "ccache version ([0-9]\\.[0-9]\\.[0-9])")
    if(${CMAKE_MATCH_1} VERSION_LESS ${CCACHE_VERSION_REQUIRED})
        message(STATUS "CCache outdated. Installed: ${CMAKE_MATCH_1}, \
            required: ${CCACHE_VERSION_REQUIRED}. Caching disabled.")
        return()
    endif()
endif()

# Set ccache as the compiler launcher
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)

if(COMPILER_IS_CLANG)
    add_compile_options(-Qunused-arguments)
endif()

if($ENV{TRAVIS})
    return()
endif()

# Check ccache pre-compiled headers config
execute_process(COMMAND ${CCACHE_TOOL_PATH} --print-config
    OUTPUT_VARIABLE CCACHE_CONFIG
    ERROR_VARIABLE CCACHE_CONFIG
)

# Regex should be "sloppiness.*time_macros.*pch_defines.*" but due to bug fixed
# in https://ccache.samba.org/releasenotes.html#_ccache_3_2_5 we have to leave
# out pch_defines. Ubuntu 16.04 comes with ccache 3.2.4 ...
string(REGEX MATCH "sloppiness.*time_macros.*"
    COTIRE_CCACHE_CONFIG ${CCACHE_CONFIG}
)

if(NOT COTIRE_CCACHE_CONFIG)
    message(FATAL_ERROR "CCache not configured! You must set sloppiness to pch_defines,time_macros. See https://docs.opengeosys.org/docs/devguide/advanced/using-ccache")
endif()
