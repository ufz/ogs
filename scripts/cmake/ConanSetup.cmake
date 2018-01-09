if(NOT OGS_USE_CONAN)
    return()
endif()

if(CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "Multi-config generators are not yet supported when "
        "using Conan. Specify CMAKE_BUILD_TYPE!")
endif()

if(DEFINED OGS_LIB_Boost)
    if(${OGS_LIB_Boost} STREQUAL "Default")
        cmake_minimum_required(VERSION 3.4) # Conan Boost package requires this
    endif()
endif()

include(${PROJECT_SOURCE_DIR}/scripts/cmake/conan/conan.cmake)

set(CONAN_REQUIRES
    Boost/1.64.0@conan/stable
    Eigen3/3.2.9@bilke/stable
    VTK/7.1.0@bilke/stable
    CACHE INTERNAL ""
)

set(CONAN_OPTIONS
    Boost:header_only=True
    Qt:xmlpatterns=True
    CACHE INTERNAL ""
)

if(LINUX AND BUILD_SHARED_LIBS)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} VTK:fPIC=True)
endif()

if(OGS_USE_MPI)
    set(CONAN_OPTIONS ${CONAN_OPTIONS} VTK:mpi=True)
endif()

if(OGS_BUILD_GUI)
    set(CONAN_REQUIRES ${CONAN_REQUIRES}
        Shapelib/1.3.0@bilke/stable
        libgeotiff/1.4.2@bilke/stable
        Qt/5.9.2@osechet/stable
    )
endif()

# Find Conan and do version check
find_program(CONAN_CMD conan)
if(NOT CONAN_CMD)
    message(FATAL_ERROR "Conan executable not found!")
endif()
execute_process(COMMAND ${CONAN_CMD} --version
    OUTPUT_VARIABLE CONAN_VERSION_OUTPUT)
string(REGEX MATCH ".*Conan version ([0-9]+\.[0-9]+\.[0-9]+)" FOO "${CONAN_VERSION_OUTPUT}")
set(CONAN_VERSION_REQUIRED 0.26.0)
if(${CMAKE_MATCH_1} VERSION_LESS ${CONAN_VERSION_REQUIRED})
    message(FATAL_ERROR "Conan outdated. Installed: ${CONAN_VERSION}, \
        required: ${CONAN_VERSION_REQUIRED}. Consider updating via 'pip \
        install conan --upgrade'.")
endif()

execute_process(COMMAND ${CONAN_CMD} remote list OUTPUT_VARIABLE CONAN_REMOTES)

# Add ogs remote
if("${CONAN_REMOTES}" MATCHES ".*ogs:.*")
    # Make sure ogs repo is first
    execute_process(COMMAND ${CONAN_CMD} remote update -i 0 ogs
        https://ogs.jfrog.io/ogs/api/conan/conan)
else()
    # Add ogs repo as first
    message(STATUS "Conan adding ogs remote repositoy \
        (https://ogs.jfrog.io/ogs/api/conan/conan)")
    execute_process(COMMAND ${CONAN_CMD} remote add -i 0 ogs
        https://ogs.jfrog.io/ogs/api/conan/conan)
endif()

# Add conan-community remote
if("${CONAN_REMOTES}" MATCHES ".*conan-community:.*")
    execute_process(COMMAND ${CONAN_CMD} remote update -i 1 conan-community
        https://api.bintray.com/conan/conan-community/conan)
else()
    message(STATUS "Conan adding community remote repositoy \
        (https://api.bintray.com/conan/conan-community/conan)")
    execute_process(COMMAND ${CONAN_CMD} remote add -i 1 conan-community
        https://api.bintray.com/conan/conan-community/conan)
endif()

# Remove libraries from Conan which are set to "System"
message(STATUS "Third-party libraries:")
foreach(LIB ${OGS_LIBS})
    message("  - OGS_LIB_${LIB} = ${OGS_LIB_${LIB}}")
    if("${OGS_LIB_${LIB}}" STREQUAL System)
        list(FILTER CONAN_REQUIRES EXCLUDE REGEX ${LIB})
    endif()
endforeach()

set(CONAN_IMPORTS "")
if(APPLE)
    set(CONAN_IMPORTS ${CONAN_IMPORTS} "lib, *.dylib* -> ./bin")
endif()

conan_cmake_run(
    BASIC_SETUP
    UPDATE
    KEEP_RPATHS
    REQUIRES ${CONAN_REQUIRES}
    OPTIONS ${CONAN_OPTIONS}
    BUILD ${OGS_CONAN_BUILD}
    IMPORTS ${CONAN_IMPORTS}
)

if(NOT ${OGS_CONAN_BUILD} MATCHES "never|always|missing")
    message(STATUS "Warning: Resetting CMake variable OGS_CONAN_BUILD to its default value of 'missing'")
    set(OGS_CONAN_BUILD "missing" CACHE INTERNAL "")
endif()
