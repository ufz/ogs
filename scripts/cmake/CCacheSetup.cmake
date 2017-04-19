if(NOT CCACHE_TOOL_PATH)
    return()
endif()

# Set ccache as the compiler launcher
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)

if(COMPILER_IS_CLANG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
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
