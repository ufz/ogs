if(NOT OGS_USE_CONAN OR NOT OGS_USE_NETCDF)
    return()
endif()
find_program(CONAN_CMD conan)
if(NOT CONAN_CMD AND OGS_USE_NETCDF)
    message(WARNING "conan executable not found. Specify CMake option "
        "OGS_USE_CONAN=auto for automatic installation in the build directory "
        "OR install it system-wide (https://www.opengeosys.org/docs/devguide/"
        "getting-started/prerequisites/#step-install-conan-package-manager) "
        "OR disable this warning with OGS_USE_CONAN=OFF.")
    return()
endif()


if(CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "Multi-config generators are not yet supported when "
        "using Conan. Specify CMAKE_BUILD_TYPE, e.g. via cmd line: "
        "cmake . -DCMAKE_BUILD_TYPE=Release!")
endif()

# Treat Conan includes as system includes to suppress warnings
set(CONAN_SYSTEM_INCLUDES ON)

include(${PROJECT_SOURCE_DIR}/scripts/cmake/conan/conan.cmake)

if(OGS_USE_NETCDF)
    set(CONAN_REQUIRES ${CONAN_REQUIRES} netcdf-cxx/4.3.1-3@bilke/testing)
    set(CONAN_OPTIONS netcdf:dap=False)
endif()

if(NOT DEFINED CONAN_REQUIRES)
    return()
endif()

conan_check(VERSION ${ogs.minimum_version.conan})
conan_config_install(ITEM ${PROJECT_SOURCE_DIR}/scripts/cmake/conan/config)

set(CONAN_IMPORTS "")
if(APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.dylib* -> ./lib")
endif()
if(MSVC)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "bin, *.dll* -> ./bin")
endif()
if(UNIX AND NOT APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.so* -> ./lib")
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "plugins/platforms, *.so* -> ./bin/platforms")
endif()

file(TIMESTAMP ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt file_timestamp "%Y.%m.%d")
string(TIMESTAMP timestamp "%Y.%m.%d")

# Run conan install update only once a day
if("${file_timestamp}" VERSION_LESS ${timestamp} OR DEFINED ENV{CI})
    file(WRITE ${PROJECT_BINARY_DIR}/conan_install_timestamp.txt "${timestamp}\n")
    set(CONAN_UPDATE UPDATE)
    set(CONAN_COMMAND ${CONAN_CMD} CACHE INTERNAL "") # Speed up conan_add_remote
    conan_add_remote(NAME ogs INDEX 0
        URL https://ogs.jfrog.io/ogs/api/conan/conan)
    conan_add_remote(NAME bincrafters INDEX 1
        URL https://bincrafters.jfrog.io/artifactory/api/conan/public-conan)
    execute_process(COMMAND ${CONAN_COMMAND} config set general.revisions_enabled=1)
else()
    message(STATUS "Conan: Skipping update step.")
endif()

if(DEFINED OGS_CONAN_BUILD_TYPE)
    set(CONAN_BUILD_TYPE ${OGS_CONAN_BUILD_TYPE})
else()
    set(CONAN_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

# speed up conan_cmake_run
set(ARGUMENTS_CONAN_COMMAND ${CONAN_CMD} CACHE INTERNAL "")
conan_cmake_run(
    BASIC_SETUP
    ${CONAN_UPDATE}
    KEEP_RPATHS
    REQUIRES ${CONAN_REQUIRES}
    OPTIONS ${CONAN_OPTIONS}
    BUILD ${OGS_CONAN_BUILD}
    IMPORTS ${CONAN_IMPORTS}
    GENERATORS virtualrunenv cmake_find_package
    BUILD_TYPE ${CONAN_BUILD_TYPE}
)

if(NOT ${OGS_CONAN_BUILD} MATCHES "never|always|missing|outdated")
    message(STATUS "Warning: Resetting CMake variable OGS_CONAN_BUILD to its default value of 'missing'")
    set(OGS_CONAN_BUILD "missing" CACHE INTERNAL "")
endif()
