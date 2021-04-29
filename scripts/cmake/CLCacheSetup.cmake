if(NOT WIN32 OR NOT CLCACHE_TOOL_PATH OR OGS_DISABLE_CLCACHE)
    return()
endif()

file(DOWNLOAD
    https://gitlab.com/timblechmann/clcache-launcher/uploads/c11842f963af3543b448cb400ed5fce9/clcache-launcher.exe
    ${PROJECT_BINARY_DIR}/clcache-launcher.exe
)

set(CMAKE_C_COMPILER_LAUNCHER ${PROJECT_BINARY_DIR}/clcache-launcher.exe)
set(CMAKE_CXX_COMPILER_LAUNCHER ${PROJECT_BINARY_DIR}/clcache-launcher.exe)

message(STATUS "clcache enabled! (${CLCACHE_TOOL_PATH})")
